/**
 * @file AD_FluidEvaluator.h
 * @brief 自动微分统一流体状态评估器 (AD Fluid Evaluator) - 工业级增强版
 * @details
 * 将基于 double 的原生黑盒状态方程插值表 (Span-Wagner / IAPWS) 桥接到 ADVar 计算图中。
 * 覆盖全套 6 大物性参数：密度、黏度、定压比热、定容比热、比焓、导热系数。
 * 包含相界一致性检测 (Phase Mismatch Detection)、高精度微扰控制及显式内存安全组装。
 */

#pragma once

#include "ADVar.hpp"
#include "WaterPropertyTable.h"
#include "CO2PropertyTable.h"
#include "PropertiesSummary.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <exception>

namespace AD_Fluid
{
    // =========================================================
    // 1. 异常类与结构体定义
    // =========================================================

    /**
     * @brief 相态不一致异常 (用于拦截跨越相边界的中心差分)
     */
    class PhaseMismatchException : public std::exception {
    public:
        const char* what() const noexcept override {
            return "Fluid Phase Mismatch detected across numerical perturbation.";
        }
    };

    /**
     * @brief 差分模式枚举 (用于 V3 诊断补丁: 追踪物性评估是否逼近相边界或越界)
     */
    enum class DiffMode { Central, Forward, Backward, ZeroGrad };

    /**
     * @struct ADFluidProperties
     * @brief 携带雅可比梯度的全套流体物性集合
     * @tparam N 独立自变量的数量
     */
    template<int N>
    struct ADFluidProperties
    {
        ADVar<N> rho; ///< 密度 [kg/m^3]
        ADVar<N> mu;  ///< 动力黏度 [Pa*s]
        ADVar<N> cp;  ///< 定压比热 [J/(kg*K)]
        ADVar<N> cv;  ///< 定容比热 [J/(kg*K)]
        ADVar<N> h;   ///< 比焓 [J/kg]
        ADVar<N> k;   ///< 导热系数 [W/(m*K)]

        bool isFallback = false; ///< [新增] 标记此网格物性是否触发了兜底机制 (物理失效)

        DiffMode diff_mode_P = DiffMode::Central; ///< 压力求导所采用的差分模式
        DiffMode diff_mode_T = DiffMode::Central; ///< 温度求导所采用的差分模式
        bool near_bound = false;                  ///< 标记是否逼近或越界 (非中心差分均视为 near_bound)
    };

    // =========================================================
    // 2. 核心评估器类
    // =========================================================

    class Evaluator
    {
    public:
        /**
         * @brief 评估水相 (Water) 的 AD 物性
         * @tparam N 独立自变量的数量
         * @param P 水相压力 (ADVar, 单位: Pa)
         * @param T 温度 (ADVar, 单位: K)
         * @return ADFluidProperties<N>
         */
                template<int N>
        static ADFluidProperties<N> evaluateWater(const ADVar<N>& P, const ADVar<N>& T)
        {
            double p_val = P.val;
            double t_val = T.val;
            auto& wt = WaterPropertyTable::instance();
            const double p_clamped = wt.clampPressure(p_val);
            const double t_clamped = wt.clampTemperature(t_val);

            WaterProperties base_props;
            ADFluidProperties<N> res;
            res.isFallback = false;
            if (p_clamped != p_val || t_clamped != t_val) {
                res.near_bound = true;
            }

            // 1. 获取基准点物性
            try {
                base_props = wt.getProperties(p_clamped, t_clamped);
            }
            catch (...) {
                // 极端越界兜底参数
                base_props.rho = 1000.0;
                base_props.mu = 1e-3;
                base_props.cp = 4200.0;
                base_props.cv = 4182.0;
                base_props.h = 1.0e5;
                base_props.k = 0.6;
                res.isFallback = true; // 触发警告标记，供求解器诊断使用
            }

            // 2. 动态自适应差分步长 (采用 std::clamp 限制微扰范围 1e-6)
            // 压力扰动限制在 [10 Pa, 1000 Pa]；温度扰动限制在 [0.001 K, 0.1 K]
            double dP = std::max(10.0, std::min(std::abs(p_clamped) * 1e-6, 1000.0));
            double dT = std::max(0.001, std::min(std::abs(t_clamped) * 1e-6, 0.1));

            double dRho_dP = 0.0, dMu_dP = 0.0, dCp_dP = 0.0, dCv_dP = 0.0, dH_dP = 0.0, dK_dP = 0.0;
            double dRho_dT = 0.0, dMu_dT = 0.0, dCp_dT = 0.0, dCv_dT = 0.0, dH_dT = 0.0, dK_dT = 0.0;

            // 3. 计算对压力的偏导数 (固定 T)
            auto evalP = [&](double p_eval) { return wt.getProperties(p_eval, t_clamped); };
            res.diff_mode_P = computeRobustDerivatives(evalP, p_clamped, dP, base_props, dRho_dP, dMu_dP, dCp_dP, dCv_dP, dH_dP, dK_dP);

            // 4. 计算对温度的偏导数 (固定 P)
            auto evalT = [&](double t_eval) { return wt.getProperties(p_clamped, t_eval); };
            res.diff_mode_T = computeRobustDerivatives(evalT, t_clamped, dT, base_props, dRho_dT, dMu_dT, dCp_dT, dCv_dT, dH_dT, dK_dT);

            // 判定是否逼近相边界
            if (res.diff_mode_P != DiffMode::Central || res.diff_mode_T != DiffMode::Central) {
                res.near_bound = true;
            }

            // 5. 显式链式法则安全组装 (消除野指针或未初始化污染)
            AssembleADVar(res.rho, base_props.rho, dRho_dP, dRho_dT, P, T);
            AssembleADVar(res.mu, base_props.mu, dMu_dP, dMu_dT, P, T);
            AssembleADVar(res.cp, base_props.cp, dCp_dP, dCp_dT, P, T);
            AssembleADVar(res.cv, base_props.cv, dCv_dP, dCv_dT, P, T);
            AssembleADVar(res.h, base_props.h, dH_dP, dH_dT, P, T);
            AssembleADVar(res.k, base_props.k, dK_dP, dK_dT, P, T);

            return res;
        }

        /**
         * @brief 评估 CO2 相的 AD 物性
         */
                template<int N>
        static ADFluidProperties<N> evaluateCO2(const ADVar<N>& P, const ADVar<N>& T)
        {
            double p_val = P.val;
            double t_val = T.val;
            auto& gt = CO2PropertyTable::instance();
            const double p_clamped = gt.clampPressure(p_val);
            const double t_clamped = gt.clampTemperature(t_val);

            CO2Properties base_props;
            ADFluidProperties<N> res;
            res.isFallback = false;
            if (p_clamped != p_val || t_clamped != t_val) {
                res.near_bound = true;
            }

            try {
                base_props = gt.getProperties(p_clamped, t_clamped);
            }
            catch (...) {
                base_props.rho = 800.0;
                base_props.mu = 1.48e-5;
                base_props.cp = 1100.0;
                base_props.cv = 850.0;
                base_props.h = 3.0e5;
                base_props.k = 0.03;
                res.isFallback = true;
            }

            double dP = std::max(10.0, std::min(std::abs(p_clamped) * 1e-6, 1000.0));
            double dT = std::max(0.001, std::min(std::abs(t_clamped) * 1e-6, 0.1));

            double dRho_dP = 0.0, dMu_dP = 0.0, dCp_dP = 0.0, dCv_dP = 0.0, dH_dP = 0.0, dK_dP = 0.0;
            double dRho_dT = 0.0, dMu_dT = 0.0, dCp_dT = 0.0, dCv_dT = 0.0, dH_dT = 0.0, dK_dT = 0.0;

            auto evalP = [&](double p_eval) { return gt.getProperties(p_eval, t_clamped); };
            res.diff_mode_P = computeRobustDerivatives(evalP, p_clamped, dP, base_props, dRho_dP, dMu_dP, dCp_dP, dCv_dP, dH_dP, dK_dP);

            auto evalT = [&](double t_eval) { return gt.getProperties(p_clamped, t_eval); };
            res.diff_mode_T = computeRobustDerivatives(evalT, t_clamped, dT, base_props, dRho_dT, dMu_dT, dCp_dT, dCv_dT, dH_dT, dK_dT);

            // 判定是否逼近相边界
            if (res.diff_mode_P != DiffMode::Central || res.diff_mode_T != DiffMode::Central) {
                res.near_bound = true;
            }

            AssembleADVar(res.rho, base_props.rho, dRho_dP, dRho_dT, P, T);
            AssembleADVar(res.mu, base_props.mu, dMu_dP, dMu_dT, P, T);
            AssembleADVar(res.cp, base_props.cp, dCp_dP, dCp_dT, P, T);
            AssembleADVar(res.cv, base_props.cv, dCv_dP, dCv_dT, P, T);
            AssembleADVar(res.h, base_props.h, dH_dP, dH_dT, P, T);
            AssembleADVar(res.k, base_props.k, dK_dP, dK_dT, P, T);

            return res;
        }

    private:

        /**
         * @brief 显式安全的 AD 变量装配器 (消除未初始化污染风险)
         */
        template<int N>
        static void AssembleADVar(ADVar<N>& target, double val, double d_dP, double d_dT,
            const ADVar<N>& P, const ADVar<N>& T)
        {
            target = ADVar<N>(val);      // 构造函数已确保 grad 被彻底 setZero()
            target.grad += d_dP * P.grad; // Eigen 的 += 操作符是严格累加，无维度灾难
            target.grad += d_dT * T.grad;
        }

        /**
         * @brief 检查扰动点是否处于同一物理相态 (防止跨相界面的伪中心差分)
         */
        template<typename TProps>
        static bool isSamePhase(const TProps& p1, const TProps& p2)
        {
            // 启发式检测：如果极小微扰下的密度发生跳变 (>1%)，判定为跨越相界面（如液相切换至超临界）
            double diff = std::abs(p1.rho - p2.rho);
            double avg = 0.5 * (p1.rho + p2.rho);
            if (avg > 1e-6 && (diff / avg) > 0.01) {
                return false;
            }
            return true;
        }

        /**
         * @brief 高鲁棒性数值微分计算核 (四级降级防崩机制)
         * @return DiffMode 当前成功应用的数值微扰模式
         */
        template<typename Func, typename TProps>
        static DiffMode computeRobustDerivatives(Func evalFunc, double val, double delta, const TProps& base_props,
            double& drho, double& dmu, double& dcp, double& dcv, double& dh, double& dk)
        {
            try {
                // Tier 1: 中心差分 (Central Difference)
                auto prop_plus = evalFunc(val + delta);
                auto prop_minus = evalFunc(val - delta);

                // 【核心增强】物理一致性检测
                if (!isSamePhase(prop_plus, prop_minus)) {
                    throw PhaseMismatchException(); // 强制进入单侧差分捕获
                }

                double inv_2d = 1.0 / (2.0 * delta);
                drho = (prop_plus.rho - prop_minus.rho) * inv_2d;
                dmu = (prop_plus.mu - prop_minus.mu) * inv_2d;
                dcp = (prop_plus.cp - prop_minus.cp) * inv_2d;
                dcv = (prop_plus.cv - prop_minus.cv) * inv_2d;
                dh = (prop_plus.h - prop_minus.h) * inv_2d;
                dk = (prop_plus.k - prop_minus.k) * inv_2d;
                return DiffMode::Central;
            }
            catch (...) {
                try {
                    // Tier 2: 前向差分 (Forward Difference - 紧贴相边界)
                    auto prop_plus = evalFunc(val + delta);
                    double inv_d = 1.0 / delta;
                    drho = (prop_plus.rho - base_props.rho) * inv_d;
                    dmu = (prop_plus.mu - base_props.mu) * inv_d;
                    dcp = (prop_plus.cp - base_props.cp) * inv_d;
                    dcv = (prop_plus.cv - base_props.cv) * inv_d;
                    dh = (prop_plus.h - base_props.h) * inv_d;
                    dk = (prop_plus.k - base_props.k) * inv_d;
                    return DiffMode::Forward;
                }
                catch (...) {
                    try {
                        // Tier 3: 后向差分 (Backward Difference)
                        auto prop_minus = evalFunc(val - delta);
                        double inv_d = 1.0 / delta;
                        drho = (base_props.rho - prop_minus.rho) * inv_d;
                        dmu = (base_props.mu - prop_minus.mu) * inv_d;
                        dcp = (base_props.cp - prop_minus.cp) * inv_d;
                        dcv = (base_props.cv - prop_minus.cv) * inv_d;
                        dh = (base_props.h - prop_minus.h) * inv_d;
                        dk = (base_props.k - prop_minus.k) * inv_d;
                        return DiffMode::Backward;
                    }
                    catch (...) {
                        // Tier 4: 终极防御，赋予 0 梯度 (Zero Gradient)
                        drho = 0.0; dmu = 0.0; dcp = 0.0; dcv = 0.0; dh = 0.0; dk = 0.0;
                        return DiffMode::ZeroGrad;
                    }
                }
            }
        }
    };
}


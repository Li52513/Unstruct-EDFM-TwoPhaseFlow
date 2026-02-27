/**
 * @file CapRelPerm_HD_AD.h
 * @brief 自动微分扩展版：毛管力(Capillary Pressure)与相对渗透率(Relative Permeability)模型
 * @details
 * 无缝引入 ADVar<N> 模板类，实现了完全前向微分兼容的 VG 模型。
 * 采用 C++ 函数重载 (Function Overloading) 机制，与原生 double 接口完全同名且并存。
 * 核心逻辑与 CapRelPerm_HD.h 保持 100% 对应，完美继承了底层工业级平滑与截断防护。
 */

#pragma once

#include "CapRelPerm_HD.h"
#include "ADVar.hpp"
#include <algorithm>
 // 基础数学库已在 CapRelPerm_HD.h 中包含

namespace CapRelPerm {

    /**
     * @brief 计算有效饱和度 (Effective Saturation) - 自动微分重载版
     * @tparam N 独立自变量数量
     * @param Sw 当前水相饱和度 (ADVar<N>)
     * @param vg VG参数结构体
     * @return ADVar<N> 有效饱和度 Se (携带梯度信息)
     */
    template<int N>
    inline ADVar<N> calculate_Se(const ADVar<N>& Sw, const VGParams& vg)
    {
        double denom = 1.0 - vg.Swr - vg.Sgr;
        if (std::abs(denom) < kTiny) {
            // 防御性保护：处理分母趋于 0 的极端情况
            return (Sw.val >= vg.Swr) ? ADVar<N>(1.0) : ADVar<N>(0.0);
        }
        return (Sw - ADVar<N>(vg.Swr)) / denom;
    }

    /**
     * @brief 计算毛管力 Pc(Sw) - 自动微分重载版
     * @details
     * 复刻原生 double 版本的工业级保护逻辑（超限截断、kSe_Eps平滑线性延伸）。
     * 求解器在调用时，返回的 ADVar<N> pc 将通过操作符重载自动携带对主变量的解析雅可比偏导数。
     * @tparam N 独立自变量数量
     * @param Sw 当前水相饱和度 (ADVar<N>)
     * @param vg VG参数结构体
     * @return ADVar<N> 毛管力 Pc (单位: Pa)
     */
    template<int N>
    inline ADVar<N> pc_vG(const ADVar<N>& Sw, const VGParams& vg)
    {
        ADVar<N> Se = calculate_Se(Sw, vg);

        // 1. 物理上限截断 (完全润湿状态)
        if (Se.val >= 1.0) {
            return ADVar<N>(0.0);
        }

        // 2. 正常 VG 计算区间 (Se >= Threshold)
        if (Se.val >= kSe_Eps) {
            double m = vg.m();

            // val = Se^(-1/m) - 1.0; 
            // 注意: 此处无 std:: 前缀，依靠 C++ 的 ADL(参数依赖查找) 自动匹配 ADVar.hpp 中的 pow 重载
            ADVar<N> val = pow(Se, -1.0 / m) - ADVar<N>(1.0);

            // 防止 val 因底层浮点精度误差微小于0而导致后续分数次幂报错
            if (val.val < 0.0) {
                val = ADVar<N>(0.0);
            }

            ADVar<N> pc = ADVar<N>(1.0 / vg.alpha) * pow(val, 1.0 / vg.n);

            if (pc.val > kPcMax) {
                return ADVar<N>(kPcMax);
            }
            return pc;
        }
        // 3. 线性化扩展区间 (Smoothing Region) -> 抑制导数爆炸
        else {
            double Se_star = kSe_Eps;
            double m = vg.m();

            // 提前计算常数项斜率和截距，这部分是纯标量 (double) 的解析计算，不浪费 AD 的堆栈资源
            double val_star = std::pow(Se_star, -1.0 / m) - 1.0;
            double pc_star = (1.0 / vg.alpha) * std::pow(val_star, 1.0 / vg.n);

            double term1 = (1.0 / (vg.alpha * vg.n));
            double term2 = std::pow(val_star, (1.0 / vg.n) - 1.0);
            double term3 = (-1.0 / m) * std::pow(Se_star, (-1.0 / m) - 1.0);
            double slope = term1 * term2 * term3;

            // 线性延伸并利用 ADVar 的链式法则自动带出常数斜率下的梯度映射
            ADVar<N> pc_linear = ADVar<N>(pc_star) + (Se - ADVar<N>(Se_star)) * slope;

            // 安全截断
            if (pc_linear.val < 0.0) return ADVar<N>(0.0);
            if (pc_linear.val > kPcMax) return ADVar<N>(kPcMax);

            return pc_linear;
        }
    }

    /**
     * @brief 计算相对渗透率 (Mualem-vG Model) - 自动微分重载版
     * @details 输出的 krw 和 krg 会自动携带对网格主变量(如 P, Sw, T) 的精确偏导，消灭手动求导错误。
     * @tparam N 独立自变量数量
     * @param Sw 水相饱和度 (ADVar<N>)
     * @param vg VG参数结构体
     * @param rp 相渗参数(包含 L 等参数)
     * @param[out] krw 水相相对渗透率 (ADVar<N>)
     * @param[out] krg 气相相对渗透率 (ADVar<N>)
     */
    template<int N>
    inline void kr_Mualem_vG(const ADVar<N>& Sw, const VGParams& vg, const RelPermParams& rp, ADVar<N>& krw, ADVar<N>& krg)
    {
        ADVar<N> Se = calculate_Se(Sw, vg);

        // 边界截断与平滑防御
        if (Se.val >= 1.0) {
            krw = ADVar<N>(1.0);
            krg = ADVar<N>(0.0);
            return;
        }
        if (Se.val <= 0.0) {
            krw = ADVar<N>(0.0);
            krg = ADVar<N>(1.0);
            return;
        }

        double m = vg.m();

        // -------------------------
        // 1. Krw (水相相对渗透率) 计算
        // -------------------------
        ADVar<N> Se_inv_m = pow(Se, 1.0 / m);
        ADVar<N> one_minus_Se_inv_m = ADVar<N>(1.0) - Se_inv_m;

        // 保护幂运算底数非负
        if (one_minus_Se_inv_m.val < 0.0) {
            one_minus_Se_inv_m = ADVar<N>(0.0);
        }

        ADVar<N> term_m = pow(one_minus_Se_inv_m, m);
        ADVar<N> bracket_w = ADVar<N>(1.0) - term_m;
        krw = pow(Se, rp.L) * (bracket_w * bracket_w);

        // -------------------------
        // 2. Krg (气相相对渗透率) 计算
        // -------------------------
        ADVar<N> bracket_g = ADVar<N>(1.0) - Se_inv_m;
        if (bracket_g.val < 0.0) {
            bracket_g = ADVar<N>(0.0);
        }

        krg = pow(ADVar<N>(1.0) - Se, rp.L) * pow(bracket_g, 2.0 * m);
    }

} // namespace CapRelPerm
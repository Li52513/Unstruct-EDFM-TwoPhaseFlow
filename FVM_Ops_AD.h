/**
 * @file FVM_Ops_AD.h
 * @brief 全隐式动态离散算子库 (Fully Implicit Dynamic Discrete Operators)
 * @details 专为 EGS (增强型地热系统) 非结构化 EDFM 开发。
 * 完美适配项目原生的 Variable3D 向量类，直接复用其 operator* (点乘) 与 operator-
 */

#ifndef FVM_OPS_AD_H
#define FVM_OPS_AD_H
#include <array>
#include "UserDefineVarType.h" // Vector = Variable3D<double>
#include <cmath>

inline double GravityProj(const Vector& g, const Vector& xi, const Vector& xj) {
    return g * (xj - xi);
}

inline double GravityProj(const std::array<double, 3>& g,
    const std::array<double, 3>& xi,
    const std::array<double, 3>& xj) {
    return g[0] * (xj[0] - xi[0]) + g[1] * (xj[1] - xi[1]) + g[2] * (xj[2] - xi[2]);
}

namespace FVM_Ops {

    /**
     * @brief 计算单相或基准相势能差 (原生 Vector 适配版)
     * @tparam VecType 项目原生的向量类 (如 Variable3D<double>)，必须支持 operator- 和 operator* (点乘)
     */
    template <int N, typename ADVarType, typename VecType>
    inline ADVarType Compute_Potential_Diff(const ADVarType& P_i,
        const ADVarType& P_j,
        const ADVarType& rho_avg,
        const VecType& x_i,
        const VecType& x_j,
        const VecType& g_vec)
    {
        ADVarType delta_P = P_j - P_i;

        // 完美复用项目原生的 Variable3D 类的运算符重载：
        // (x_j - x_i) 触发向量减法，生成新向位移矢量
        // g_vec * (...) 触发向量点乘，返回 double
        double g_dot_dx = g_vec * (x_j - x_i);

        ADVarType gravity_term = rho_avg * g_dot_dx;

        return delta_P - gravity_term;
    }

    /**
     * @brief 计算两相流中携带毛管力 (Pc) 的相态势能差 (原生 Vector 适配版)
     */
    template <int N, typename ADVarType, typename VecType>
    inline ADVarType Compute_Potential_Diff(const ADVarType& P_i,
        const ADVarType& P_j,
        const ADVarType& Pc_i,
        const ADVarType& Pc_j,
        const ADVarType& rho_avg,
        const VecType& x_i,
        const VecType& x_j,
        const VecType& g_vec)
    {
        ADVarType PhaseP_i = P_i + Pc_i;
        ADVarType PhaseP_j = P_j + Pc_j;
        ADVarType delta_PhaseP = PhaseP_j - PhaseP_i;

        double g_dot_dx = GravityProj(g_vec, x_i, x_j);
        ADVarType gravity_term = rho_avg * g_dot_dx;

        return delta_PhaseP - gravity_term;
    }

    /**
     * @brief 迎风算子 (Upwind Operator) - 绝对锁定偏导数
     */
    template <int N, typename ADVarType>
    inline ADVarType Op_Upwind_AD(const ADVarType& delta_Phi,
        const ADVarType& var_i,
        const ADVarType& var_j)
    {
        if (delta_Phi.val < 0.0) {
            return var_i;
        }
        else {
            return var_j;
        }
    }

    /**
     * @brief 质量通量算子 (Mass Flux Operator)
     */
    template <int N, typename ADVarType>
    inline ADVarType Compute_Mass_Flux(double T_flow,
        const ADVarType& upwind_mob,
        const ADVarType& delta_Phi)
    {
        return upwind_mob * delta_Phi * T_flow;
    }

    /**
     * @brief 热通量算子 (Heat Flux Operator) - 单相流重载
     */
    template <int N, typename ADVarType>
    inline ADVarType Compute_Heat_Flux(double T_heat,
        const ADVarType& Temp_i,
        const ADVarType& Temp_j,
        const ADVarType& mass_flux,
        const ADVarType& upwind_enthalpy)
    {
        ADVarType q_cond = (Temp_j - Temp_i) * T_heat;
        ADVarType q_conv = mass_flux * upwind_enthalpy;
        return q_cond + q_conv;
    }

    /**
     * @brief 热通量算子 (Heat Flux Operator) - 两相流重载 (水 + CO2)
     */
    template <int N, typename ADVarType>
    inline ADVarType Compute_Heat_Flux(double T_heat,
        const ADVarType& Temp_i,
        const ADVarType& Temp_j,
        const ADVarType& mass_flux_w,
        const ADVarType& mass_flux_co2,
        const ADVarType& upwind_enthalpy_w,
        const ADVarType& upwind_enthalpy_co2)
    {
        ADVarType q_cond = (Temp_j - Temp_i) * T_heat;
        ADVarType q_conv = (mass_flux_w * upwind_enthalpy_w) +
            (mass_flux_co2 * upwind_enthalpy_co2);
        return q_cond + q_conv;
    }

} // namespace FVM_Ops

#endif // FVM_OPS_AD_H
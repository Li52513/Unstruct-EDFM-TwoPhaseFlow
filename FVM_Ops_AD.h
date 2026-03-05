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

    template <int N, typename ADVarType>
    inline ADVarType Op_Boundary_Dirichlet_AD(double T_face, const ADVarType& Phi_cell, double Phi_bc) {
        return T_face * (Phi_cell - Phi_bc);
    }

    template <int N, typename ADVarType>
    inline ADVarType Op_Boundary_Neumann_AD(double area, double Flux_bc) {
        ADVarType q_N; q_N.val = area * Flux_bc;
        for (int i = 0; i < N; ++i) q_N.grad(i) = 0.0;
        return q_N;
    }

    template <int N, typename ADVarType>
    inline ADVarType Op_Leakoff_Source_AD(bool enable, double C_L, const ADVarType& P_cell, double P_farfield) {
        if (!enable) {
            ADVarType zero; zero.val = 0.0;
            for (int i = 0; i < N; ++i) zero.grad(i) = 0.0;
            return zero;
        }
        return C_L * (P_cell - P_farfield);
    }

    template <int N, typename ADVarType>
    inline void Op_Leakoff_TwoPhase_AD(bool enable, double C_L_w, double C_L_g,
        const ADVarType& P_cell_w, const ADVarType& P_cell_g,
        double P_farfield, ADVarType& q_w_out, ADVarType& q_g_out) {
        if (!enable) {
            q_w_out.val = 0.0; q_g_out.val = 0.0;
            for (int i = 0; i < N; ++i) { q_w_out.grad(i) = 0.0; q_g_out.grad(i) = 0.0; }
            return;
        }
        q_w_out = C_L_w * (P_cell_w - P_farfield);
        q_g_out = C_L_g * (P_cell_g - P_farfield);
    }

    /**
     * @brief 定井底流压 (BHP) 井源项算子 (Outflow positive)
     * @param well_conductance 井指数与流度密度的乘积 (WI * rho * k_r / mu)
     * @param P_cell 网格中心压力 (ADVar)
     * @param P_bhp 目标井底流压 (Constant)
     */
    template <int N, typename ADVarType>
    inline ADVarType Op_Well_BHP_Source_AD(double well_conductance, const ADVarType& P_cell, double P_bhp) {
        return well_conductance * (P_cell - P_bhp);
    }

    /**
         * @brief 定流量 (Rate) 井源项算子 (Outflow positive)
         * @param q_target 目标流量 (采出为正，注入为负)
         * @note 导数显式置为 0，因为恒定流量对中心网格压力的偏导为零
         */
    template <int N, typename ADVarType>
    inline ADVarType Op_Well_Rate_Source_AD(double q_target) {
        // [Patch 8] 显式初始化后再赋值
        ADVarType q_mass;
        q_mass.val = q_target;
        for (int i = 0; i < N; ++i) {
            q_mass.grad(i) = 0.0;
        }
        return q_mass;
    }

    /**
     * @brief 能量方程井源项算子
     */
    template <int N, typename ADVarType>
    inline ADVarType Op_Well_Energy_Source_AD(const ADVarType& q_mass_w, double h_w, const ADVarType& q_mass_g, double h_g) {
        return (q_mass_w * h_w) + (q_mass_g * h_g);
    }

} // namespace FVM_Ops

#endif // FVM_OPS_AD_H
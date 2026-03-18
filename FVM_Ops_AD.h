/**
 * @file FVM_Ops_AD.h
 * @brief ȫ��ʽ��̬��ɢ���ӿ� (Fully Implicit Dynamic Discrete Operators)
 * @details רΪ EGS (��ǿ�͵���ϵͳ) �ǽṹ�� EDFM ������
 * ����������Ŀԭ���� Variable3D �����ֱ࣬�Ӹ����� operator* (���) �� operator-
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
     * @brief ���㵥����׼�����ܲ� (ԭ�� Vector �����)
     * @tparam VecType ��Ŀԭ���������� (�� Variable3D<double>)������֧�� operator- �� operator* (���)
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

        // ����������Ŀԭ���� Variable3D �����������أ�
        // (x_j - x_i) ����������������������λ��ʸ��
        // g_vec * (...) ����������ˣ����� double
        double g_dot_dx = g_vec * (x_j - x_i);

        ADVarType gravity_term = rho_avg * g_dot_dx;

        return delta_P - gravity_term;
    }

    /**
     * @brief ������������Я��ë���� (Pc) ����̬���ܲ� (ԭ�� Vector �����)
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
     * @brief ӭ������ (Upwind Operator) - ��������ƫ����
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
     * @brief ����ͨ������ (Mass Flux Operator)
     */
    template <int N, typename ADVarType>
    inline ADVarType Compute_Mass_Flux(double T_flow,
        const ADVarType& upwind_mob,
        const ADVarType& delta_Phi)
    {
        return upwind_mob * delta_Phi * T_flow;
    }

    /**
     * @brief ��ͨ������ (Heat Flux Operator) - ����������
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
     * @brief ��ͨ������ (Heat Flux Operator) - ���������� (ˮ + CO2)
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
    inline ADVarType Op_Boundary_ScaleFlux_AD(const ADVarType& q_raw, double coeff) {
        return q_raw * coeff;
    }

    template <int N, typename ADVarType>
    inline ADVarType Op_Boundary_ScaleFlux_AD(const ADVarType& q_raw, const ADVarType& coeff) {
        return q_raw * coeff;
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
     * @brief ��������ѹ (BHP) ��Դ������ (Outflow positive)
     * @param well_conductance ��ָ���������ܶȵĳ˻� (WI * rho * k_r / mu)
     * @param P_cell ��������ѹ�� (ADVar)
     * @param P_bhp Ŀ�꾮����ѹ (Constant)
     */
    template <int N, typename ADVarType>
    inline ADVarType Op_Well_BHP_Source_AD(double well_conductance, const ADVarType& P_cell, double P_bhp) {
        return well_conductance * (P_cell - P_bhp);
    }

    /**
         * @brief ������ (Rate) ��Դ������ (Outflow positive)
         * @param q_target Ŀ������ (�ɳ�Ϊ����ע��Ϊ��)
         * @note ������ʽ��Ϊ 0����Ϊ�㶨��������������ѹ����ƫ��Ϊ��
         */
    template <int N, typename ADVarType>
    inline ADVarType Op_Well_Rate_Source_AD(double q_target) {
        // [Patch 8] ��ʽ��ʼ�����ٸ�ֵ
        ADVarType q_mass;
        q_mass.val = q_target;
        for (int i = 0; i < N; ++i) {
            q_mass.grad(i) = 0.0;
        }
        return q_mass;
    }

    /**
     * @brief �������̾�Դ������
     */
    template <int N, typename ADVarType>
    inline ADVarType Op_Well_Energy_Source_AD(const ADVarType& q_mass_w, double h_w, const ADVarType& q_mass_g, double h_g) {
        return (q_mass_w * h_w) + (q_mass_g * h_g);
    }

    /**
     * @brief Non-orthogonal correction operator (deferred correction, scalar helper)
     * @details
     * Computes the tangential pressure correction:
     *   corr = grad_phi_f · vectorT
     * where grad_phi_f is the face-interpolated gradient and vectorT is the
     * non-orthogonal tangential vector (from Face::computeFaceVectors).
     *
     * This is an EXPLICIT scalar correction (constant wrt AD variables) — it
     * only contributes to the residual, not to the Jacobian.  The full flux
     * correction for the mass equation is:
     *   q_corr = K_eff_corr * mob_upwind * Op_NonOrthogonal_PressureCorr(...)
     * where K_eff_corr = T_Flow * aux_dist / aux_area.
     *
     * @param grad_phi_f  Face-interpolated cell-centre gradient [quantity/m]
     * @param vectorT     Non-orthogonal correction vector [m] (from Connection)
     * @return double     Dot product [quantity] (same units as pressure difference ΔP)
     */
    inline double Op_NonOrthogonal_PressureCorr(const Vector& grad_phi_f, const Vector& vectorT)
    {
        return grad_phi_f.m_x * vectorT.m_x
             + grad_phi_f.m_y * vectorT.m_y
             + grad_phi_f.m_z * vectorT.m_z;
    }

} // namespace FVM_Ops

#endif // FVM_OPS_AD_H
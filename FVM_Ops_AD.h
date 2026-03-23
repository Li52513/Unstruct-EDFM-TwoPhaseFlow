#ifndef FVM_OPS_AD_H
#define FVM_OPS_AD_H
#include <array>
#include "UserDefineVarType.h" // Vector = Variable3D<double>
#include <cmath>

/*
* brief: 计算重力势能差的工具函数，重载1：输入为Vector类型
* 输入参数：重力加速度向量g，单元中心坐标xi和xj
*/
inline double GravityProj(const Vector& g, const Vector& xi, const Vector& xj) {
    return g * (xj - xi);
}

/*
* brief: 计算重力势能差的工具函数，重载2：输入为std::array<double, 3>类型
* 输入参数：重力加速度向量g，单元中心坐标xi和xj
*/
inline double GravityProj(const std::array<double, 3>& g,
    const std::array<double, 3>& xi,
    const std::array<double, 3>& xj) {
    return g[0] * (xj[0] - xi[0]) + g[1] * (xj[1] - xi[1]) + g[2] * (xj[2] - xi[2]);
}

namespace FVM_Ops {

    /*
    * brief: 计算单相系统的势能差，考虑重力
	* 输入参数：owner单元压力P_i和P_j，平均密度rho_avg，单元中心坐标x_i和x_j，重力加速度向量g_vec
	* 输出：势能差 delta_Phi = phi_j - phi_i，其中 phi = P - rho * g * z
	* 公式： \Delta \Phi = (P_j - P_i) - \rho_{avg} * g \cdot (x_j - x_i)
	* 注： phi = P - rho * g * z (z为垂直坐标)，因此 delta_phi = (P_j - P_i) - rho_avg * g * (z_j - z_i)
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

        double g_dot_dx = g_vec * (x_j - x_i);

        ADVarType gravity_term = rho_avg * g_dot_dx;

        return delta_P - gravity_term;
    }

    /*
	*  brief: 计算两相系统的势能差，考虑毛细力和重力
	*  输入参数：owner单元压力P_i和P_j，毛细压力Pc_i和Pc_j，平均密度rho_avg，单元中心坐标x_i和x_j，重力加速度向量g_vec
	*  输出： 势能差 delta_Phi = phi_j - phi_i，其中 phi = P + Pc - rho * g * z
	*  公式：相态压力 = 绝对压力 + 毛细压力，因此势能差为 \Delta \Phi = (P_j + Pc_j) - (P_i + Pc_i) - \rho_{avg} * g \cdot (x_j - x_i)
	*  注： 在两相系统中，毛细压力的差异直接影响相态压力，从而改变迎风方向和通量计算。确保 AD 变量的导数正确传播对于 Newton-Raphson 收敛至关重要。
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

    /*
	*  brief: 迎风选择算子，基于势能差的符号判断流动方向，并返回对应的变量值
	*  输入参数：delta_Phi：势能差，var_i和var_j：对应单元i和j的变量值
	*  输出参数: 迎风变量值，根据 delta_Phi 的符号选择 var_i 或 var_j
	*  公式： 如果 delta_Phi < 0，流动方向为 i -> j，返回 var_i；否则流动方向为 j -> i，返回 var_j
	*  注：因为delta_Phi = phi_j - phi_i，delta_Phi < 0 意味着 phi_i > phi_j，流体从 i 流向 j，因此迎风变量是 var_i。反之亦然。确保 AD 变量的导数正确传播对于 Newton-Raphson 收敛至关重要。
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

    /*
	*  brief: 计算质量通量的工具函数，基于迎风流度和势能差
	*  输入参数：T_flow：流动传导系数，upwind_mob：迎风流度，delta_Phi：势能差
	*  输出参数：质量通量，根据 Darcy 定律计算，考虑迎风流度和势能差
	*  公式： Mass_Flux = T_flow * upwind_mob * delta_Phi
    */
    template <int N, typename ADVarType>
    inline ADVarType Compute_Mass_Flux(double T_flow,
        const ADVarType& upwind_mob,
        const ADVarType& delta_Phi)
    {
        return upwind_mob * delta_Phi * T_flow;
    }

    /*
	*  brief: 计算单相热通量的工具函数，基于温度差、传导系数、质量通量和迎风焓值
	*  输入参数: T_heat: 热传导系数，Temp_i和Temp_j: 单元i和j的温度，mass_flux: 质量通量，upwind_enthalpy: 迎风焓值
	*  输出参数：热通量，包含传导项和对流项两部分
	*  公式： 热通量 = 传导项 + 对流项 = T_heat * (Temp_j - Temp_i) + mass_flux * upwind_enthalpy
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

    /*
	*  brief: 计算两相热通量的工具函数，基于温度差、传导系数、各相质量通量和迎风焓值
	*  输入参数： T_heat: 热传导系数，Temp_i和Temp_j: 单元i和j的温度，mass_flux_w和mass_flux_co2: 水相和CO2相的质量通量，upwind_enthalpy_w和upwind_enthalpy_co2: 水相和CO2相的迎风焓值
	*  输出参数： 热通量，包含传导项和两相对流项
	*  公式： 热通量 = 传导项 + 对流项 = T_heat * (Temp_j - Temp_i) + (mass_flux_w * upwind_enthalpy_w) + (mass_flux_co2 * upwind_enthalpy_co2)
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

    /*
	*  brief: Dirichlet边界条件算子，基于单元内变量值和边界条件值计算边界通量
	*  输入参数： T_face: 边界传导系数，Phi_cell: 单元内变量值，Phi_bc: 边界条件值
	*  输出参数： 边界通量，根据 Dirichlet 边界条件计算，考虑边界传导系数和单元内变量与边界条件的差异
	*  公式： Boundary_Flux = T_face * (Phi_cell - Phi_bc)
    */
    template <int N, typename ADVarType>
    inline ADVarType Op_Boundary_Dirichlet_AD(double T_face, const ADVarType& Phi_cell, double Phi_bc) {
        return T_face * (Phi_cell - Phi_bc);
    }

    /*
	*  brief: Neumann边界条件算子，基于边界面积和通量值计算边界通量
	*  输入参数： area: 边界面积，Flux_bc: Neumann边界条件值（通量值）
	*  输出参数： 边界通量，根据 Neumann 边界条件计算，直接使用边界条件值乘以边界面积
	*  公式： Boundary_Flux = area * Flux_bc
    */
    template <int N, typename ADVarType>
    inline ADVarType Op_Boundary_Neumann_AD(double area, double Flux_bc) {
        ADVarType q_N; q_N.val = area * Flux_bc;
        for (int i = 0; i < N; ++i) q_N.grad(i) = 0.0;
        return q_N;
    }

    /*
	*  brief: 边界通量缩放算子，基于原始通量值和缩放系数计算新的边界通量
	*  输入参数： q_raw: 原始边界通量值，coeff: 缩放系数（可以是常数或AD变量）
	*  输出参数： 新的边界通量，根据缩放系数调整原始通量值
	*  公式： Scaled_Flux = q_raw * coeff
    */
    template <int N, typename ADVarType>
    inline ADVarType Op_Boundary_ScaleFlux_AD(const ADVarType& q_raw, double coeff) {
        
        return q_raw * coeff;
    }

    /*
	*  brief: 边界通量缩放算子重载，支持AD变量作为缩放系数
	*  输入参数： q_raw: 原始边界通量值，coeff: 缩放系数（AD变量）
	*  输出参数： 新的边界通量，根据缩放系数调整原始通量值
	*  公式： Scaled_Flux = q_raw * coeff
    */
    template <int N, typename ADVarType>
    inline ADVarType Op_Boundary_ScaleFlux_AD(const ADVarType& q_raw, const ADVarType& coeff) {
        return q_raw * coeff;
    }

    /*
	*  brief: Leakoff源项算子，基于单元内压力、远场压力和泄漏系数计算Leakoff通量
	*  输入参数： enable: 是否启用Leakoff，C_L: 泄漏系数，P_cell: 单元内压力（AD变量），P_farfield: 远场压力（常数）
	*  输出参数： Leakoff通量，根据Leakoff开关状态计算，如果enable=true，则根据单元内压力和远场压力的差异计算Leakoff通量；如果enable=false，则Leakoff通量为0
	*  公式： Leakoff_Flux = C_L * (P_cell - P_farfield)  (当enable=true时)，否则Leakoff_Flux=0
    */
    template <int N, typename ADVarType>
    inline ADVarType Op_Leakoff_Source_AD(bool enable, double C_L, const ADVarType& P_cell, double P_farfield) {
        if (!enable) {
            ADVarType zero; zero.val = 0.0;
            for (int i = 0; i < N; ++i) zero.grad(i) = 0.0;
            return zero;
        }
        return C_L * (P_cell - P_farfield);
    }

    /*
	*  brief: 两相Leakoff源项算子，基于单元内水相和气相压力、远场压力和泄漏系数计算水相和气相的Leakoff通量
	*  输入参数： enable: 是否启用Leakoff，C_L_w和C_L_g: 水相和气相的泄漏系数，P_cell_w和P_cell_g: 单元内水相和气相压力（AD变量），P_farfield: 远场压力（常数），q_w_out和q_g_out: 输出的水相和气相Leakoff通量（AD变量）
	*  输出参数： 水相和气相的Leakoff通量，根据Leakoff开关状态计算，如果enable=true，则根据单元内压力和远场压力的差异计算水相和气相的Leakoff通量；如果enable=false，则两相Leakoff通量均为0
	*  公式： 当enable=true时，水相Leakoff_Flux = C_L_w * (P_cell_w - P_farfield)，气相Leakoff_Flux = C_L_g * (P_cell_g - P_farfield)；否则两相Leakoff_Flux均为0
    */
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

    /*
	*  brief: 井算子（BHP源项），基于单元内压力、井底压力和井导纳计算井的质量通量
	*  输入参数： well_conductance: 井导纳（conductance），P_cell: 单元内压力（AD变量），P_bhp: 井底压力（常数）
	*  输出参数： 井的质量通量，根据 Darcy 定律计算，考虑井导纳和单元内压力与井底压力的差异 
	*  公式 ： Well_Flux = well_conductance * (P_cell - P_bhp) 正值表示生产井（流出），负值表示注入井（流入）。确保 AD 变量的导数正确传播对于 Newton-Raphson 收敛至关重要。
    */
    template <int N, typename ADVarType>
    inline ADVarType Op_Well_BHP_Source_AD(double well_conductance, const ADVarType& P_cell, double P_bhp) {
        return well_conductance * (P_cell - P_bhp);
    }

    /*
	*  brief: 井算子（Rate源项），基于目标注入/生产速率计算井的质量通量
	*  输入参数： q_target: 目标注入/生产速率（正值表示生产井，负值表示注入井）
	*  输出参数： 井的质量通量，直接使用目标速率值作为通量值，导数为0，因为速率是一个固定的边界条件，不随单元内变量变化而变化
	* 公式： Well_Flux = q_target，正值表示生产井（流出），负值表示注入井（流入），加到方程左边。确保 AD 变量的导数为零，以反映速率边界条件的固定性质。
    */
    template <int N, typename ADVarType>
    inline ADVarType Op_Well_Rate_Source_AD(double q_target) {
        ADVarType q_mass;
        q_mass.val = q_target;
        for (int i = 0; i < N; ++i) {
            q_mass.grad(i) = 0.0;
        }
        return q_mass;
    }

    /*
	*  brief: 井能量源项算子，基于水相和气相的质量通量以及对应的焓值计算井的热通量
    * 在调用 Op_Well_Energy_Source_AD 时，传入的 h_w 和 h_g 绝对不能在注入和生产时用同一个变量！必须根据 q_mass 的正负号进行严格的分支判断：当 q_mass > 0（生产井）：流体是从网格流向井筒的，因此携带的是网格内部的能量。你传入的 h 必须是根据当前网格的温度和压力计算出的原位热焓 $h_{cell}$。当 q_mass < 0（注入井）：流体是从地面沿着井筒打进网格的，携带的是外部的能量。
    * 你传入的 h 必须是根据地面注入温度（或者井底注入温度）计算出的外部热焓 $h_{inject}$。
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
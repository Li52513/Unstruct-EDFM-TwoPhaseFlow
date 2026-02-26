/**
 * @file Test_ADVar.h
 * @brief ADVar 自动微分内核生产级验证程序
 * @details 验证多项式、超越函数、以及综合物理关联式的解析梯度与有限差分梯度的一致性。
 */
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include "ADVar.hpp"

#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "ADVar.hpp"

 // =====================================================================
 // 场景一：微可压缩流体密度模型 (Exponential Compressibility Model)
 // 公式：rho(P, T) = rho_ref * exp( c_p * (P - P_ref) - c_T * (T - T_ref) )
 // =====================================================================

template<int N>
ADVar<N> CalcDensity_AD(const ADVar<N>& P, const ADVar<N>& T)
{
    double rho_ref = 1000.0;
    double P_ref = 1e5;
    double T_ref = 293.15;
    double c_p = 4.5e-10; // 液体压缩系数 (Pa^-1)
    double c_T = 2.0e-4;  // 热膨胀系数 (K^-1)

    return rho_ref * exp(c_p * (P - P_ref) - c_T * (T - T_ref));
}

double CalcDensity_Num(double P, double T)
{
    double rho_ref = 1000.0;
    double P_ref = 1e5;
    double T_ref = 293.15;
    double c_p = 4.5e-10;
    double c_T = 2.0e-4;

    return rho_ref * std::exp(c_p * (P - P_ref) - c_T * (T - T_ref));
}

// =====================================================================
// 场景二：Brooks-Corey 相对渗透率模型 (Relative Permeability)
// 公式：S_e = (S_w - S_wr) / (1 - S_wr - S_nr)
//       Krw = S_e ^ n_w
// =====================================================================

template<int N>
ADVar<N> CalcRelPerm_AD(const ADVar<N>& Sw)
{
    double S_wr = 0.2;  // 束缚水饱和度
    double S_nr = 0.05; // 残余非润湿相饱和度
    double n_w = 3.0;   // Corey 指数

    // 限制归一化饱和度在 [0, 1] 之间 (防越界)
    ADVar<N> Se = (Sw - S_wr) / (1.0 - S_wr - S_nr);
    Se = max(min(Se, 1.0), 0.0); // 嵌套使用 min/max 测试

    return pow(Se, n_w);
}

double CalcRelPerm_Num(double Sw)
{
    double S_wr = 0.2;
    double S_nr = 0.05;
    double n_w = 3.0;

    double Se = (Sw - S_wr) / (1.0 - S_wr - S_nr);
    Se = std::max(std::min(Se, 1.0), 0.0);

    return std::pow(Se, n_w);
}

// =====================================================================
// 场景三：迎风算子 / 惩罚函数 (Non-smooth Max/Min Operator)
// 公式：f(P) = max(P - P_threshold, 0.0)^2 * 1.5
// =====================================================================

template<int N>
ADVar<N> CalcUpwindPenalty_AD(const ADVar<N>& P)
{
    double P_threshold = 15e6; // 15 MPa
    ADVar<N> diff = P - P_threshold;
    ADVar<N> active_part = max(diff, 0.0);
    return pow(active_part, 2.0) * 1.5;
}

double CalcUpwindPenalty_Num(double P)
{
    double P_threshold = 15e6;
    double diff = P - P_threshold;
    double active_part = std::max(diff, 0.0);
    return std::pow(active_part, 2.0) * 1.5;
}

// =====================================================================
// 统一验证框架
// =====================================================================

void Run_ADVar_Comprehensive_Tests()
{
    std::cout << "\n=====================================================================" << std::endl;
    std::cout << "          STARTING ADVar COMPREHENSIVE VERIFICATION SUITE              " << std::endl;
    std::cout << "=====================================================================" << std::endl;

    constexpr int N_DOF = 3; // 假设系统变量为 P, S, T
    const int IDX_P = 0;
    const int IDX_S = 1;
    const int IDX_T = 2;

    std::cout << std::scientific << std::setprecision(8);

    // -----------------------------------------------------------------
    // 测试 A: 微可压缩密度模型 (P = 20 MPa, T = 350 K)
    // -----------------------------------------------------------------
    std::cout << "\n[Test A] Micro-Compressible Liquid Density EOS" << std::endl;
    double eval_P = 20e6;
    double eval_T = 350.0;

    ADVar<N_DOF> P_ad(eval_P, IDX_P);
    ADVar<N_DOF> T_ad(eval_T, IDX_T);

    ADVar<N_DOF> Rho_ad = CalcDensity_AD(P_ad, T_ad);

    double eps_P = eval_P * 1e-6;
    double eps_T = eval_T * 1e-6;
    double num_dRho_dP = (CalcDensity_Num(eval_P + eps_P, eval_T) - CalcDensity_Num(eval_P - eps_P, eval_T)) / (2.0 * eps_P);
    double num_dRho_dT = (CalcDensity_Num(eval_P, eval_T + eps_T) - CalcDensity_Num(eval_P, eval_T - eps_T)) / (2.0 * eps_T);

    double err_Rho_P = std::abs(Rho_ad.grad(IDX_P) - num_dRho_dP) / std::max(1e-10, std::abs(num_dRho_dP));
    double err_Rho_T = std::abs(Rho_ad.grad(IDX_T) - num_dRho_dT) / std::max(1e-10, std::abs(num_dRho_dT));

    std::cout << "  Value: rho(P, T) = " << Rho_ad.val << " kg/m^3" << std::endl;
    std::cout << "  dRho/dP -> Exact: " << Rho_ad.grad(IDX_P) << " | Num: " << num_dRho_dP << " | Err: " << err_Rho_P << std::endl;
    std::cout << "  dRho/dT -> Exact: " << Rho_ad.grad(IDX_T) << " | Num: " << num_dRho_dT << " | Err: " << err_Rho_T << std::endl;

    assert(err_Rho_P < 1e-6 && "Density Jacobian vs Pressure FAILED!");
    assert(err_Rho_T < 1e-6 && "Density Jacobian vs Temperature FAILED!");

    // -----------------------------------------------------------------
    // 测试 B: 相对渗透率 (Sw = 0.6)
    // -----------------------------------------------------------------
    std::cout << "\n[Test B] Brooks-Corey Relative Permeability" << std::endl;
    double eval_S = 0.6;
    ADVar<N_DOF> S_ad(eval_S, IDX_S);

    ADVar<N_DOF> Krw_ad = CalcRelPerm_AD(S_ad);

    double eps_S = 1e-6;
    double num_dKrw_dS = (CalcRelPerm_Num(eval_S + eps_S) - CalcRelPerm_Num(eval_S - eps_S)) / (2.0 * eps_S);

    double err_Krw_S = std::abs(Krw_ad.grad(IDX_S) - num_dKrw_dS) / std::max(1e-10, std::abs(num_dKrw_dS));

    std::cout << "  Value: Krw(Sw) = " << Krw_ad.val << std::endl;
    std::cout << "  dKrw/dSw -> Exact: " << Krw_ad.grad(IDX_S) << " | Num: " << num_dKrw_dS << " | Err: " << err_Krw_S << std::endl;

    assert(err_Krw_S < 1e-6 && "RelPerm Jacobian vs Saturation FAILED!");

    // --- 边界测试 (Sw = 0.1, 低于束缚水, 导数应严格为 0) ---
    ADVar<N_DOF> Krw_ad_bound = CalcRelPerm_AD(ADVar<N_DOF>(0.1, IDX_S));
    assert(Krw_ad_bound.val == 0.0 && "RelPerm value bounding FAILED!");
    assert(Krw_ad_bound.grad(IDX_S) == 0.0 && "RelPerm gradient bounding FAILED!");

    // -----------------------------------------------------------------
    // 测试 C: 非平滑迎风/惩罚函数 (P_high = 16 MPa, P_low = 14 MPa)
    // -----------------------------------------------------------------
    std::cout << "\n[Test C] Non-Smooth Max/Min Functions (Upwind Penalty)" << std::endl;

    // 情形 1：触发激活域 (P = 16 MPa > 15 MPa)
    double eval_P_high = 16e6;
    ADVar<N_DOF> P_high_ad(eval_P_high, IDX_P);
    ADVar<N_DOF> Pen_high_ad = CalcUpwindPenalty_AD(P_high_ad);

    double num_dPen_dP_high = (CalcUpwindPenalty_Num(eval_P_high + eps_P) - CalcUpwindPenalty_Num(eval_P_high - eps_P)) / (2.0 * eps_P);
    double err_Pen_high = std::abs(Pen_high_ad.grad(IDX_P) - num_dPen_dP_high) / std::max(1e-10, std::abs(num_dPen_dP_high));

    std::cout << "  Active (16MPa) dF/dP -> Exact: " << Pen_high_ad.grad(IDX_P) << " | Num: " << num_dPen_dP_high << " | Err: " << err_Pen_high << std::endl;
    assert(err_Pen_high < 1e-6 && "Penalty Jacobian (Active Region) FAILED!");

    // 情形 2：未触发激活域 (P = 14 MPa < 15 MPa)
    ADVar<N_DOF> P_low_ad(14e6, IDX_P);
    ADVar<N_DOF> Pen_low_ad = CalcUpwindPenalty_AD(P_low_ad);

    std::cout << "  Inactive(14MPa) dF/dP -> Exact: " << Pen_low_ad.grad(IDX_P) << " (Expected: 0.0)" << std::endl;
    assert(Pen_low_ad.val == 0.0 && "Penalty value (Inactive Region) FAILED!");
    assert(Pen_low_ad.grad(IDX_P) == 0.0 && "Penalty gradient (Inactive Region) FAILED!");

    std::cout << "\n=====================================================================" << std::endl;
    std::cout << "          [SUCCESS] ALL COMPREHENSIVE ADVar TESTS PASSED!            " << std::endl;
    std::cout << "=====================================================================" << std::endl;
}
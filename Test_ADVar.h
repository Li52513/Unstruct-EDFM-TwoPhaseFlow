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

 // =====================================================================
 // 目标测试函数: 模拟实际状态方程中复杂的非线性依赖
 // f(P, T) = P^2 / exp(T / 100) + sin(P/1e6) * T^1.5
 // =====================================================================

// 1. 基于 ADVar 的全隐式方程实现 (一套代码，同时算 Value 和 Jacobian)
template<int N>
ADVar<N> CalculateFluidProperty(const ADVar<N>& P, const ADVar<N>& T)
{
    // 利用 ADVar 的算子重载，代码外观与处理普通 double 完全一致
    ADVar<N> term1 = pow(P, 2.0) / exp(T / 100.0);
    ADVar<N> term2 = sin(P / 1e6) * pow(T, 1.5);
    return term1 + term2;
}

// 2. 基于纯 double 的相同方程 (用于有限差分)
double CalculateFluidProperty_Double(double P, double T)
{
    double term1 = std::pow(P, 2.0) / std::exp(T / 100.0);
    double term2 = std::sin(P / 1e6) * std::pow(T, 1.5);
    return term1 + term2;
}

// =====================================================================
// 验证模块
// =====================================================================
void Verify_AutomaticDifferentiation()
{
    std::cout << "\n========== [Test 3] ADVar Jacobian Precision Verification ==========" << std::endl;

    // 假设系统为主变量为 P, T，所以 N_dof = 2
    constexpr int N_DOF = 2;
    const int IDX_P = 0;
    const int IDX_T = 1;

    // 物理评估点 (例如 P = 10 MPa, T = 350 K)
    double eval_P = 10e6;
    double eval_T = 350.0;

    // A. 自动微分解析计算 (Analytical Jacobian via Forward-Mode AD)
    // -----------------------------------------------------------------
    // 构造自变量并注入种子梯度: P为第0个变量，T为第1个变量
    ADVar<N_DOF> P_ad(eval_P, IDX_P);
    ADVar<N_DOF> T_ad(eval_T, IDX_T);

    // 只需执行一次正向计算，不仅得到值，连梯度 (dF/dP, dF/dT) 也同时组装完毕！
    ADVar<N_DOF> F_ad = CalculateFluidProperty(P_ad, T_ad);

    // B. 传统有限差分计算 (Numerical Jacobian via Central Difference)
    // -----------------------------------------------------------------
    double eps_P = eval_P * 1e-6; // 压力扰动量
    double eps_T = eval_T * 1e-6; // 温度扰动量

    // dF/dP
    double F_P_plus = CalculateFluidProperty_Double(eval_P + eps_P, eval_T);
    double F_P_minus = CalculateFluidProperty_Double(eval_P - eps_P, eval_T);
    double num_dF_dP = (F_P_plus - F_P_minus) / (2.0 * eps_P);

    // dF/dT
    double F_T_plus = CalculateFluidProperty_Double(eval_P, eval_T + eps_T);
    double F_T_minus = CalculateFluidProperty_Double(eval_P, eval_T - eps_T);
    double num_dF_dT = (F_T_plus - F_T_minus) / (2.0 * eps_T);


    // C. 精度对齐验证
    // -----------------------------------------------------------------
    double err_dP = std::abs(F_ad.grad(IDX_P) - num_dF_dP) / std::max(1.0, std::abs(num_dF_dP));
    double err_dT = std::abs(F_ad.grad(IDX_T) - num_dF_dT) / std::max(1.0, std::abs(num_dF_dT));

    std::cout << std::scientific << std::setprecision(8);
    std::cout << "Function Evaluation : F(P, T) = " << F_ad.val << std::endl;
    std::cout << "-> dF/dP (AD Exact) : " << F_ad.grad(IDX_P) << " | (Numerical): " << num_dF_dP << " | RelErr: " << err_dP << std::endl;
    std::cout << "-> dF/dT (AD Exact) : " << F_ad.grad(IDX_T) << " | (Numerical): " << num_dF_dT << " | RelErr: " << err_dT << std::endl;

    // 断言容差检查 (由于中心差分本身有截断误差，且物理量级在 1e12 级别，相对误差通常在 1e-8 级别)
    assert(err_dP < 1e-7 && "Jacobian against Pressure fails precision test!");
    assert(err_dT < 1e-7 && "Jacobian against Temperature fails precision test!");

    // D. 稳定性与 NaN 防护测试
    // -----------------------------------------------------------------
    ADVar<N_DOF> zeroVar(0.0, IDX_P);
    ADVar<N_DOF> powZero = pow(zeroVar, 2.0); // 合法，梯标应为0
    ADVar<N_DOF> sqrtZero = sqrt(zeroVar);    // 极限合法，防除了 NaN

    assert(!powZero.hasNaN() && "Math stability check failed (pow)");
    assert(!sqrtZero.hasNaN() && "Math stability check failed (sqrt)");

    std::cout << "[PASS] Automatic Differentiation Kernel verified successfully." << std::endl;
}
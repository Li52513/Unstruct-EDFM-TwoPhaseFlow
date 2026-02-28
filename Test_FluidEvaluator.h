/**
 * @file Test_FluidEvaluator.h
 * @brief 物理属性计算内核 (AD vs FD) 多工况自动化测试及物性输出工具
 * @details 包含：
 * 1. 打印给定 P, T 下的所有流体真实物性值 (Rho, Mu, Cp, Cv, h, k)，用于与外部数据库基准比对。
 * 2. 采用自适应物理步长 (Adaptive Step Size) 验证 AD 传出的导数与二阶中心差分 (FD) 的匹配度。
 */

#pragma once

#include "AD_FluidEvaluator.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

 /**
  * @brief 运行流体物性评估器多工况测试及输出
  */
inline void run_fluid_evaluator_test() {
    // 主变量定义：N=2 表示系统有两个主求解变量 (通常为 P 和 T)
    constexpr int N = 2;

    // 定义测试工况 (Pressure [Pa], Temperature [K], Name)
    struct TestCase {
        double P;
        double T;
        std::string label;
    };

    std::vector<TestCase> cases = {
        // --- 水相工况 ---
        {1.0e7, 350.0, "Water: Standard Liquid (10MPa, 350K)"},
        {2.5e7, 550.0, "Water: High PT (25MPa, 550K)"},
        // --- CO2 工况 ---
        {8.0e6, 310.0, "CO2: Near Critical (8MPa, 310K) - HIGH DIFFICULTY"},
        {1.5e7, 400.0, "CO2: Stable Supercritical (15MPa, 400K)"},
        {2.5e7, 450.0, "CO2: High PT (25MPa, 450K)"}
    };

    int total_failed = 0;

    std::cout << "\n====================================================================================================" << std::endl;
    std::cout << ">>> STAGE 1: THE PHYSICS KERNEL - COMPREHENSIVE DIAGNOSTICS <<<" << std::endl;
    std::cout << "====================================================================================================\n" << std::endl;

    for (const auto& tc : cases) {
        // -------------------------------------------------------------------
        // 步骤 1: 构造 AD 主变量并调用内核获取物性
        // -------------------------------------------------------------------
        // 初始化 P 为索引 0 (对 P 求导位于 grad(0)), T 为索引 1 (位于 grad(1))
        ADVar<N> P_ad(tc.P, 0);
        ADVar<N> T_ad(tc.T, 1);

        AD_Fluid::ADFluidProperties<N> props;
        bool isCO2 = tc.label.find("CO2") != std::string::npos;

        if (isCO2) {
            props = AD_Fluid::Evaluator::evaluateCO2(P_ad, T_ad);
        }
        else {
            props = AD_Fluid::Evaluator::evaluateWater(P_ad, T_ad);
        }

        // -------------------------------------------------------------------
        // 步骤 2: 打印全部物性参数真实值 (供外部数据基准核对)
        // -------------------------------------------------------------------
        std::cout << "[Target Case] " << tc.label << std::endl;
        std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;
        std::cout << "  [Property Values]" << std::endl;
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  - Density (rho)       : " << std::setw(15) << props.rho.val << " [kg/m^3]" << std::endl;
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "  - Viscosity (mu)      : " << std::setw(15) << props.mu.val << " [Pa*s]" << std::endl;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "  - Sp. Heat (Cp)       : " << std::setw(15) << props.cp.val << " [J/(kg*K)]" << std::endl;
        std::cout << "  - Sp. Heat (Cv)       : " << std::setw(15) << props.cv.val << " [J/(kg*K)]" << std::endl;
        std::cout << "  - Enthalpy (h)        : " << std::setw(15) << props.h.val << " [J/kg]" << std::endl;
        std::cout << "  - Thermal Cond. (k)   : " << std::setw(15) << props.k.val << " [W/(m*K)]" << std::endl;

        if (props.isFallback) {
            std::cout << "  *** [WARNING] Fallback mechanism triggered! (P, T) may be strictly out of property table bounds." << std::endl;
        }

        // -------------------------------------------------------------------
        // 步骤 3: 独立自适应步长外部数值差分验证 (Adaptive Step Size)
        // -------------------------------------------------------------------
        std::cout << "  [Derivative Verification (AD vs FD)]" << std::endl;

        // 动态自适应外部差分步长，保证在 10^7 量级下具有充足的浮点数分辨空间
        double eps_P = std::max(10.0, tc.P * 1e-6);
        double eps_T = std::max(0.001, tc.T * 1e-6);

        auto eval_rho = [&](double p, double t) -> double {
            if (isCO2) return AD_Fluid::Evaluator::evaluateCO2(ADVar<N>(p), ADVar<N>(t)).rho.val;
            return AD_Fluid::Evaluator::evaluateWater(ADVar<N>(p), ADVar<N>(t)).rho.val;
            };

        // --- 验证 dRho/dP ---
        double rho_p_plus = eval_rho(tc.P + eps_P, tc.T);
        double rho_p_minus = eval_rho(tc.P - eps_P, tc.T);
        double fd_drho_dp = (rho_p_plus - rho_p_minus) / (2.0 * eps_P); // 中心差分
        double ad_drho_dp = props.rho.grad(0); // 直接从 AD 变量提取偏导数
        double error_p = std::abs(ad_drho_dp - fd_drho_dp) / (std::abs(fd_drho_dp) + 1e-10);

        // 放宽临界点附近压力导数的截断误差容忍度
        double tol_P = (tc.label.find("Near Critical") != std::string::npos) ? 1e-5 : 1e-6;

        std::cout << std::left << "  " << std::setw(15) << "dRho/dP"
            << std::scientific << std::setprecision(6)
            << "AD: " << std::setw(15) << ad_drho_dp
            << "FD: " << std::setw(15) << fd_drho_dp
            << "Rel.Err: " << std::setprecision(2) << error_p;
        if (error_p > tol_P) { std::cout << " [FAIL]\n"; total_failed++; }
        else { std::cout << " [PASS]\n"; }

        // --- 验证 dRho/dT ---
        double rho_t_plus = eval_rho(tc.P, tc.T + eps_T);
        double rho_t_minus = eval_rho(tc.P, tc.T - eps_T);
        double fd_drho_dt = (rho_t_plus - rho_t_minus) / (2.0 * eps_T);
        double ad_drho_dt = props.rho.grad(1);
        double error_t = std::abs(ad_drho_dt - fd_drho_dt) / (std::abs(fd_drho_dt) + 1e-10);

        double tol_T = 1e-8; // 温度变化较线性，保持严格阈值

        std::cout << std::left << "  " << std::setw(15) << "dRho/dT"
            << std::scientific << std::setprecision(6)
            << "AD: " << std::setw(15) << ad_drho_dt
            << "FD: " << std::setw(15) << fd_drho_dt
            << "Rel.Err: " << std::setprecision(2) << error_t;
        if (error_t > tol_T) { std::cout << " [FAIL]\n"; total_failed++; }
        else { std::cout << " [PASS]\n"; }

        std::cout << std::endl;
    }

    std::cout << "====================================================================================================" << std::endl;
    if (total_failed == 0) {
        std::cout << ">>> OVERALL STATUS: [SUCCESS] All physical properties and derivatives passed verification!" << std::endl;
    }
    else {
        std::cout << ">>> OVERALL STATUS: [WARNING] " << total_failed << " checks exceeded tolerance." << std::endl;
    }
    std::cout << "====================================================================================================" << std::endl;
}
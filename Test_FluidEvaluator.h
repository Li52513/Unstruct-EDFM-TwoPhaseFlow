/**
 * @file Test_FluidEvaluator.h
 * @brief 物理属性计算内核 (AD vs FD) 多工况自动化测试工具
 */

#pragma once

#include "AD_FluidEvaluator.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

 /**
  * @brief 运行流体物性评估器多工况测试
  * @details 验证 AD 传出的导数与二阶中心差分 (Finite Difference) 的匹配度
  */
inline void run_fluid_evaluator_test() {
    constexpr int N = 2; // 主变量：0-压力(P), 1-温度(T)
    const double eps = 1e-4; // 外部验证差分步长 (需略大于 Evaluator 内部微扰)
    const double tolerance = 1e-8; // 批判性收敛阈值

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

    std::cout << std::string(100, '=') << std::endl;
    std::cout << std::left << std::setw(45) << "Test Case Condition"
        << std::setw(15) << "Property"
        << std::setw(15) << "AD Deriv"
        << std::setw(15) << "FD Deriv"
        << "Rel.Error" << std::endl;
    std::cout << std::string(100, '-') << std::endl;

    int total_failed = 0;

    for (const auto& tc : cases) {
        // 1. 构造 AD 主变量 (初始化 P 为索引 0, T 为索引 1)
        ADVar<N> P_ad(tc.P, 0);
        ADVar<N> T_ad(tc.T, 1);

        // 2. 调用内核获取物性 (以密度 Rho 为例进行验证)
        AD_Fluid::ADFluidProperties<N> props;
        bool isCO2 = tc.label.find("CO2") != std::string::npos;

        if (isCO2) {
            props = AD_Fluid::Evaluator::evaluateCO2(P_ad, T_ad);
        }
        else {
            props = AD_Fluid::Evaluator::evaluateWater(P_ad, T_ad);
        }

        // 3. 执行外部数值差分 (验证 dRho/dP)
        auto eval_rho = [&](double p, double t) -> double {
            if (isCO2) return AD_Fluid::Evaluator::evaluateCO2(ADVar<N>(p), ADVar<N>(t)).rho.val;
            return AD_Fluid::Evaluator::evaluateWater(ADVar<N>(p), ADVar<N>(t)).rho.val;
            };

        // 二阶中心差分: [f(x+e) - f(x-e)] / 2e
        double rho_p_plus = eval_rho(tc.P + eps, tc.T);
        double rho_p_minus = eval_rho(tc.P - eps, tc.T);
        double fd_drho_dp = (rho_p_plus - rho_p_minus) / (2.0 * eps);

        // 4. 对比误差 (针对密度对压力的导数)
        double ad_drho_dp = props.rho.grad(0);
        double error_p = std::abs(ad_drho_dp - fd_drho_dp) / (std::abs(fd_drho_dp) + 1e-10);

        // 5. 打印结果
        std::cout << std::left << std::setw(45) << tc.label
            << std::setw(15) << "dRho/dP"
            << std::fixed << std::setprecision(6)
            << std::setw(15) << ad_drho_dp
            << std::setw(15) << fd_drho_dp
            << std::scientific << std::setprecision(2) << error_p;

        if (error_p > tolerance) {
            std::cout << " [FAIL]";
            total_failed++;
        }
        else {
            std::cout << " [PASS]";
        }
        std::cout << std::endl;

        // 也可以验证 dRho/dT (索引 1)
        double rho_t_plus = eval_rho(tc.P, tc.T + eps);
        double rho_t_minus = eval_rho(tc.P, tc.T - eps);
        double fd_drho_dt = (rho_t_plus - rho_t_minus) / (2.0 * eps);
        double ad_drho_dt = props.rho.grad(1);
        double error_t = std::abs(ad_drho_dt - fd_drho_dt) / (std::abs(fd_drho_dt) + 1e-10);

        std::cout << std::left << std::setw(45) << ""
            << std::setw(15) << "dRho/dT"
            << std::fixed << std::setprecision(6)
            << std::setw(15) << ad_drho_dt
            << std::setw(15) << fd_drho_dt
            << std::scientific << std::setprecision(2) << error_t;

        if (error_t > tolerance) { std::cout << " [FAIL]"; total_failed++; }
        else { std::cout << " [PASS]"; }
        std::cout << "\n" << std::endl;
    }

    std::cout << std::string(100, '=') << std::endl;
    if (total_failed == 0) {
        std::cout << ">>> CRITICAL SUCCESS: All physical property derivatives passed 10^-8 check!" << std::endl;
    }
    else {
        std::cout << ">>> CRITICAL FAILURE: " << total_failed << " checks exceeded tolerance." << std::endl;
        std::cout << ">>> Suggestion: Upgrade PropertyTable interpolation to Bicubic Spline (C1 continuous)." << std::endl;
    }
}
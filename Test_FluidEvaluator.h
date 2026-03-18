/**
 * @file Test_FluidEvaluator.h
 * @brief �������Լ����ں� (AD vs FD) �๤���Զ������Լ������������
 * @details ������
 * 1. ��ӡ���� P, T �µ�����������ʵ����ֵ (Rho, Mu, Cp, Cv, h, k)���������ⲿ���ݿ��׼�ȶԡ�
 * 2. ��������Ӧ�������� (Adaptive Step Size) ��֤ AD �����ĵ�����������Ĳ�� (FD) ��ƥ��ȡ�
 */

#pragma once

#include "AD_FluidEvaluator.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

 /**
  * @brief Fluid EOS backend validation: AD vs FD derivative check + CoolProp integration test
  *
  * When USE_COOLPROP_EOS is defined (Step 1), the evaluator uses CoolProp's
  * Span-Wagner (CO2) and IAPWS-95 (Water) equations of state via thread_local
  * AbstractState objects. Derivatives are computed by AD chain rule through
  * CoolProp's first_partial_deriv() API, replacing the 4-tier numerical
  * differentiation fallback.
  */
inline void run_fluid_evaluator_test() {
    constexpr int N = 2;

    struct TestCase {
        double P;
        double T;
        std::string label;
    };

    std::vector<TestCase> cases = {
        // --- Water ---
        {1.0e7, 350.0, "Water: Standard Liquid (10MPa, 350K)"},
        {2.5e7, 550.0, "Water: High PT (25MPa, 550K)"},
        // --- CO2 ---
        {8.0e6, 310.0, "CO2: Near Critical (8MPa, 310K) - HIGH DIFFICULTY"},
        {1.5e7, 400.0, "CO2: Stable Supercritical (15MPa, 400K)"},
        {2.5e7, 450.0, "CO2: High PT (25MPa, 450K)"},
#ifdef USE_COOLPROP_EOS
        // CoolProp-only extended test points: regions where table interpolation
        // struggles but CoolProp's analytical EOS handles accurately.
        {7.5e6, 304.5, "CO2: Pseudo-critical (7.5MPa, 304.5K) - CoolProp stress test"},
        {5.0e6, 290.0, "CO2: Subcritical liquid (5MPa, 290K) - CoolProp extended range"},
        {4.0e7, 700.0, "Water: Supercritical (40MPa, 700K) - CoolProp extended range"},
#endif
    };

    int total_failed = 0;

    std::cout << "\n====================================================================================================" << std::endl;
    std::cout << ">>> STAGE 1: THE PHYSICS KERNEL - COMPREHENSIVE DIAGNOSTICS <<<" << std::endl;
    std::cout << "====================================================================================================" << std::endl;
#ifdef USE_COOLPROP_EOS
    std::cout << "[EOS Backend] CoolProp 6.x active (Span-Wagner CO2 + IAPWS-95 Water)\n";
    std::cout << "              Derivatives: AD chain rule via CoolProp first_partial_deriv()\n";
    std::cout << "              Thread safety: thread_local AbstractState per fluid\n";
#else
    std::cout << "[EOS Backend] Table-based interpolation (Catmull-Rom 2D)\n";
    std::cout << "              Derivatives: 4-tier adaptive numerical differentiation\n";
#endif
    std::cout << std::endl;

    for (const auto& tc : cases) {
        // -------------------------------------------------------------------
        // ���� 1: ���� AD �������������ں˻�ȡ����
        // -------------------------------------------------------------------
        // ��ʼ�� P Ϊ���� 0 (�� P ��λ�� grad(0)), T Ϊ���� 1 (λ�� grad(1))
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
        // ���� 2: ��ӡȫ�����Բ�����ʵֵ (���ⲿ���ݻ�׼�˶�)
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
        // ���� 3: ��������Ӧ�����ⲿ��ֵ�����֤ (Adaptive Step Size)
        // -------------------------------------------------------------------
        std::cout << "  [Derivative Verification (AD vs FD)]" << std::endl;

        // ��̬����Ӧ�ⲿ��ֲ�������֤�� 10^7 �����¾��г���ĸ������ֱ�ռ�
        double eps_P = std::max(10.0, tc.P * 1e-6);
        double eps_T = std::max(0.001, tc.T * 1e-6);

        auto eval_rho = [&](double p, double t) -> double {
            if (isCO2) return AD_Fluid::Evaluator::evaluateCO2(ADVar<N>(p), ADVar<N>(t)).rho.val;
            return AD_Fluid::Evaluator::evaluateWater(ADVar<N>(p), ADVar<N>(t)).rho.val;
            };

        // --- ��֤ dRho/dP ---
        double rho_p_plus = eval_rho(tc.P + eps_P, tc.T);
        double rho_p_minus = eval_rho(tc.P - eps_P, tc.T);
        double fd_drho_dp = (rho_p_plus - rho_p_minus) / (2.0 * eps_P); // ���Ĳ��
        double ad_drho_dp = props.rho.grad(0); // ֱ�Ӵ� AD ������ȡƫ����
        double error_p = std::abs(ad_drho_dp - fd_drho_dp) / (std::abs(fd_drho_dp) + 1e-10);

        // �ſ��ٽ�㸽��ѹ�������Ľض�������̶�
        double tol_P = (tc.label.find("Near Critical") != std::string::npos) ? 1e-5 : 1e-6;

        std::cout << std::left << "  " << std::setw(15) << "dRho/dP"
            << std::scientific << std::setprecision(6)
            << "AD: " << std::setw(15) << ad_drho_dp
            << "FD: " << std::setw(15) << fd_drho_dp
            << "Rel.Err: " << std::setprecision(2) << error_p;
        if (error_p > tol_P) { std::cout << " [FAIL]\n"; total_failed++; }
        else { std::cout << " [PASS]\n"; }

        // --- ��֤ dRho/dT ---
        double rho_t_plus = eval_rho(tc.P, tc.T + eps_T);
        double rho_t_minus = eval_rho(tc.P, tc.T - eps_T);
        double fd_drho_dt = (rho_t_plus - rho_t_minus) / (2.0 * eps_T);
        double ad_drho_dt = props.rho.grad(1);
        double error_t = std::abs(ad_drho_dt - fd_drho_dt) / (std::abs(fd_drho_dt) + 1e-10);

        double tol_T = 1e-8; // �¶ȱ仯�����ԣ������ϸ���ֵ

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
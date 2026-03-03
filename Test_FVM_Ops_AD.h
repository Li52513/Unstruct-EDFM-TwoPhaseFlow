/**
 * @file Test_FVM_Ops_AD.h
 * @brief Day 2 测试排雷器 (Test Suite for FVM AD Operators)
 * @details 包含静水平衡、迎风导数阻断、多相热通量闭环以及毛细力驱动的严苛测试。
 */

#ifndef TEST_FVM_OPS_AD_H
#define TEST_FVM_OPS_AD_H

#include <iostream>
#include <cassert>
#include <cmath>
#include "FVM_Ops_AD.h"
#include "Variable3D.hpp" 

namespace Test_FVM {

    using Vector = Variable3D<double>;

    /**
     * @brief 自动化测试入口：重力静水平衡 (基于原生 Vector 点乘)
     */
    template <int N, typename ADVarType>
    void Benchmark_Hydrostatic() {
        ADVarType P_i; P_i.val = 100000.0;
        ADVarType P_j; P_j.val = 109810.0;
        ADVarType Pc_i; Pc_i.val = 5000.0;
        ADVarType Pc_j; Pc_j.val = 5000.0;
        ADVarType rho; rho.val = 1000.0;

        Vector x_i(0.0, 0.0, 0.0);
        Vector x_j(0.0, 0.0, -1.0);
        Vector g_vec(0.0, 0.0, -9.81);

        ADVarType delta_Phi = FVM_Ops::Compute_Potential_Diff<N, ADVarType, Vector>(P_i, P_j, Pc_i, Pc_j, rho, x_i, x_j, g_vec);
        assert(std::abs(delta_Phi.val) < 1e-8 && "Benchmark_Hydrostatic 失败：势能差不为零！");

        double T_flow = 1.5e-13;
        ADVarType upwind_mob; upwind_mob.val = 5.0e-4;
        ADVarType mass_flux = FVM_Ops::Compute_Mass_Flux<N, ADVarType>(T_flow, upwind_mob, delta_Phi);
        assert(std::abs(mass_flux.val) < 1e-8 && "Benchmark_Hydrostatic 失败：静水平衡下产生非物理质量通量！");

        std::cout << "  [PASS] Benchmark_Hydrostatic (Native Variable3D & Mass Flux)" << std::endl;
    }

    /**
     * @brief 自动化测试入口：毛细力驱动通量及导数链式法则测试
     * @details 设定：P_i = P_j, g*dx = 0, Pc_j > Pc_i。
     * 断言：引发非零通量，迎风方向正确，且对 Pc_i, Pc_j 的导数符号和绝对值严格匹配物理方程。
     */
    template <int N, typename ADVarType>
    void Test_Capillary_Driven_Flow() {
        // 1. 状态初始化 (静压相等)
        ADVarType P_i; P_i.val = 100000.0;
        ADVarType P_j; P_j.val = 100000.0;

        // 初始化毛细力，设 Pc_j (5000) > Pc_i (2000)
        ADVarType Pc_i; Pc_i.val = 2000.0;
        ADVarType Pc_j; Pc_j.val = 5000.0;

        // 赋予 AD 偏导数占位 (假设在 N=3 系统中，Pc_i/Pc_j 分别在索引 0 和 1 上有导数 1.0)
        for (int k = 0; k < N; ++k) {
            Pc_i.grad(k) = (k == 0) ? 1.0 : 0.0;
            Pc_j.grad(k) = (k == 1) ? 1.0 : 0.0;
            P_i.grad(k) = 0.0; // 本测试不关心绝对压力的导数
            P_j.grad(k) = 0.0;
        }

        // 屏蔽重力影响 (水平网格，g_vec点乘为0)
        ADVarType rho; rho.val = 1000.0;
        Vector x_i(0.0, 0.0, 0.0);
        Vector x_j(1.0, 0.0, 0.0);
        Vector g_vec(0.0, 0.0, -9.81);

        // 2. 计算相态势能差 ( \Delta \Phi = (P_j + Pc_j) - (P_i + Pc_i) )
        ADVarType delta_Phi = FVM_Ops::Compute_Potential_Diff<N, ADVarType, Vector>(P_i, P_j, Pc_i, Pc_j, rho, x_i, x_j, g_vec);

        // 断言势能差的数值: \Delta \Phi = 5000 - 2000 = 3000
        assert(std::abs(delta_Phi.val - 3000.0) < 1e-12 && "毛细力驱动势能差数值计算错误！");
        // 断言势能差的导数: \partial \Delta \Phi / \partial Pc_i = -1.0
        assert(std::abs(delta_Phi.grad(0) - (-1.0)) < 1e-12 && "毛细力驱动势能差对 Pc_i 偏导数错误！");
        // 断言势能差的导数: \partial \Delta \Phi / \partial Pc_j = 1.0
        assert(std::abs(delta_Phi.grad(1) - 1.0) < 1e-12 && "毛细力驱动势能差对 Pc_j 偏导数错误！");

        // 3. 迎风流度提取
        ADVarType mob_i; mob_i.val = 1.0;
        ADVarType mob_j; mob_j.val = 2.0;
        ADVarType upwind_mob = FVM_Ops::Op_Upwind_AD<N, ADVarType>(delta_Phi, mob_i, mob_j);

        // \Delta \Phi = 3000 > 0，流动方向为 j -> i，上游是 j
        assert(std::abs(upwind_mob.val - 2.0) < 1e-12 && "毛细力驱动迎风方向判断错误 (应为 j -> i)！");

        // 4. 组装通量并验证链式法则
        double T_flow = 1.0e-13;
        ADVarType flux = FVM_Ops::Compute_Mass_Flux<N, ADVarType>(T_flow, upwind_mob, delta_Phi);

        // 断言通量数值: Flux = 1.0e-13 * 2.0 * 3000 = 6.0e-10
        assert(std::abs(flux.val - 6.0e-10) < 1e-18 && "毛细力驱动通量数值错误！");

        // 断言通量的导数 (此时流度没有参与求导，仅作为常数 mob_j = 2.0 参与链式法则):
        // \partial F / \partial Pc_i = T_flow * mob_j * (-1.0) = -2.0e-13
        // \partial F / \partial Pc_j = T_flow * mob_j * (1.0)  =  2.0e-13
        assert(std::abs(flux.grad(0) - (-2.0e-13)) < 1e-20 && "通量对 Pc_i 导数(Jacobian)组装错误！");
        assert(std::abs(flux.grad(1) - (2.0e-13)) < 1e-20 && "通量对 Pc_j 导数(Jacobian)组装错误！");

        std::cout << "  [PASS] Test_Capillary_Driven_Flow (Pc-Driven Flux & AD Jacobian)" << std::endl;
    }

    /**
     * @brief 自动化测试入口：迎风导数锁定
     */
    template <int N, typename ADVarType>
    void Test_Upwind_Derivative() {
        ADVarType var_i; var_i.val = 1.0;
        ADVarType var_j; var_j.val = 2.0;

        for (int k = 0; k < N; ++k) {
            var_i.grad(k) = (k == 0) ? 1.0 : 0.0;
            var_j.grad(k) = (k == 1) ? 1.0 : 0.0;
        }

        ADVarType delta_Phi_neg; delta_Phi_neg.val = -5000.0;
        ADVarType upwind_var_A = FVM_Ops::Op_Upwind_AD<N, ADVarType>(delta_Phi_neg, var_i, var_j);
        assert(std::abs(upwind_var_A.val - var_i.val) < 1e-12 && "迎风判定流向错误 (i->j)！");
        for (int k = 0; k < N; ++k) {
            assert(std::abs(upwind_var_A.grad(k) - var_i.grad(k)) < 1e-12 && "[致命错误] i->j 下游变量(j)导数泄露污染！");
        }

        ADVarType delta_Phi_pos; delta_Phi_pos.val = 5000.0;
        ADVarType upwind_var_B = FVM_Ops::Op_Upwind_AD<N, ADVarType>(delta_Phi_pos, var_i, var_j);
        assert(std::abs(upwind_var_B.val - var_j.val) < 1e-12 && "迎风判定流向错误 (j->i)！");
        for (int k = 0; k < N; ++k) {
            assert(std::abs(upwind_var_B.grad(k) - var_j.grad(k)) < 1e-12 && "[致命错误] j->i 下游变量(i)导数泄露污染！");
        }
        std::cout << "  [PASS] Test_Upwind_Derivative (Eigen API Double-Checked)" << std::endl;
    }

    /**
     * @brief 自动化测试入口：多物理场热通量组装
     */
    template <int N, typename ADVarType>
    void Test_Heat_Flux_Assembly() {
        double T_heat = 2.0;
        ADVarType Temp_i; Temp_i.val = 300.0;
        ADVarType Temp_j; Temp_j.val = 350.0;

        ADVarType m_flux_1p; m_flux_1p.val = -0.5;
        ADVarType h_upwind_1p; h_upwind_1p.val = 4000.0;
        ADVarType q_heat_1p = FVM_Ops::Compute_Heat_Flux<N, ADVarType>(T_heat, Temp_i, Temp_j, m_flux_1p, h_upwind_1p);
        assert(std::abs(q_heat_1p.val - (-1900.0)) < 1e-12 && "单相热通量闭环计算错误！");

        ADVarType m_flux_w; m_flux_w.val = -0.5;
        ADVarType m_flux_co2; m_flux_co2.val = -0.2;
        ADVarType h_w; h_w.val = 4000.0;
        ADVarType h_co2; h_co2.val = 2000.0;
        ADVarType q_heat_2p = FVM_Ops::Compute_Heat_Flux<N, ADVarType>(T_heat, Temp_i, Temp_j, m_flux_w, m_flux_co2, h_w, h_co2);
        assert(std::abs(q_heat_2p.val - (-2300.0)) < 1e-12 && "两相热通量闭环计算错误！");

        std::cout << "  [PASS] Test_Heat_Flux_Assembly (Single & Two-Phase Fully Covered)" << std::endl;
    }

    /**
     * @brief 总控执行器
     */
    template <int N, typename ADVarType>
    void Run_All_Day2_Tests() {
        std::cout << "[Test_FVM_Ops_AD] 启动严苛测试套件..." << std::endl;
        Benchmark_Hydrostatic<N, ADVarType>();
        Test_Capillary_Driven_Flow<N, ADVarType>(); // <--- 新增的毛细力专项测试
        Test_Upwind_Derivative<N, ADVarType>();
        Test_Heat_Flux_Assembly<N, ADVarType>();
        std::cout << "[Test_FVM_Ops_AD] 所有自动化 CI 节点测试通过！可以安全合入主分支。\n" << std::endl;
    }

} // namespace Test_FVM

#endif // TEST_FVM_OPS_AD_H
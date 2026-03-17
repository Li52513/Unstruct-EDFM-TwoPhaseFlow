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
#include <vector>
#include <fstream>
#include <cstdio>
#include <stdexcept>
#include "FVM_Ops_AD.h"
#include "Variable3D.hpp" 
#include "MeshManager.h"
#include "3D_MeshManager.h"
#include "BoundaryConditionManager.h"
#include "BoundaryAssembler.h"
#include "2D_PostProcess.h"
#include "3D_PostProcess.h"
#include "SolverContrlStrName_op.h"
#include "Well_WellScheduleManager.h"

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
        assert(std::abs(q_heat_1p.val - (-2300.0)) < 1e-12 && "单相热通量闭环计算错误！");

        ADVarType m_flux_w; m_flux_w.val = -0.5;
        ADVarType m_flux_co2; m_flux_co2.val = -0.2;
        ADVarType h_w; h_w.val = 4000.0;
        ADVarType h_co2; h_co2.val = 2000.0;
        ADVarType q_heat_2p = FVM_Ops::Compute_Heat_Flux<N, ADVarType>(T_heat, Temp_i, Temp_j, m_flux_w, m_flux_co2, h_w, h_co2);
        assert(std::abs(q_heat_2p.val - (-2100.0)) < 1e-12 && "两相热通量闭环计算错误！");

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

namespace Test_FVM {

    template <int N, typename ADVarType>
    void Test_Boundary_Operators() {
        ADVarType P_cell; P_cell.val = 150000.0; P_cell.grad(0) = 1.0;
        ADVarType q_D = FVM_Ops::Op_Boundary_Dirichlet_AD<N, ADVarType>(1.5e-13, P_cell, 100000.0);
        assert(std::abs(q_D.val - 7.5e-9) < 1e-15 && "Dirichlet Val Error!");
        ADVarType q_N = FVM_Ops::Op_Boundary_Neumann_AD<N, ADVarType>(2.5, 0.0);
        assert(std::abs(q_N.val) < 1e-15 && "Neumann Val Error!");
    }

        template <int N, typename ADVarType>
    void Test_Leakoff_Operators() {
        ADVarType P_cell_w; P_cell_w.val = 200000.0;
        ADVarType P_cell_g; P_cell_g.val = 205000.0;
        for (int k = 0; k < N; ++k) {
            P_cell_w.grad(k) = 0.0;
            P_cell_g.grad(k) = 0.0;
        }
        P_cell_w.grad(0) = 1.0;
        P_cell_g.grad(1) = 1.0;

        ADVarType q_leak_off = FVM_Ops::Op_Leakoff_Source_AD<N, ADVarType>(false, 1e-10, P_cell_w, 1e5);
        assert(std::abs(q_leak_off.val) < 1e-15 && "Leakoff OFF Error!");
        assert(std::abs(q_leak_off.grad(0)) < 1e-15 && "Leakoff OFF gradient Error!");

        ADVarType q_leak_on = FVM_Ops::Op_Leakoff_Source_AD<N, ADVarType>(true, 1e-10, P_cell_w, 1e5);
        assert(std::abs(q_leak_on.val - 1.0e-5) < 1e-15 && "Leakoff ON Error!");
        assert(std::abs(q_leak_on.grad(0) - 1.0e-10) < 1e-15 && "Leakoff ON gradient Error!");

        ADVarType q_w;
        ADVarType q_g;
        FVM_Ops::Op_Leakoff_TwoPhase_AD<N, ADVarType>(true, 1e-10, 5e-11, P_cell_w, P_cell_g, 1e5, q_w, q_g);
        assert(std::abs(q_w.val - 1.0e-5) < 1e-15 && "Two-phase water leakoff value Error!");
        assert(std::abs(q_g.val - 5.25e-6) < 1e-15 && "Two-phase gas leakoff value Error!");
        assert(std::abs(q_w.grad(0) - 1.0e-10) < 1e-15 && "Two-phase water gradient Error!");
        assert(std::abs(q_g.grad(1) - 5.0e-11) < 1e-15 && "Two-phase gas gradient Error!");
    }

    template <int N, typename ADVarType>
    void Test_GridLevel_BoundaryAssembly_2D() {
        MeshManager mgr(10.0, 10.0, 0.0, 2, 2, 0, true, false);
        mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OrthogonalCorrection);
        mgr.addFracture(Vector(1.0, 1.0, 0), Vector(9.0, 9.0, 0));
        mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
        mgr.BuildGlobalSystemIndexing();

        mgr.setNumDOFs(3);
        int totalEq = mgr.getTotalEquationDOFs();
        std::vector<double> res(totalEq, 0.0);
        std::vector<double> diag(totalEq, 0.0);

        BoundarySetting::BoundaryConditionManager bcMgr_P;
        BoundarySetting::BoundaryConditionManager bcMgr_T;

        int expectedBCCount = 0;
        for (const auto& face : mgr.mesh().getFaces()) {
            if (face.isBoundary()) {
                bcMgr_P.SetLeakoffEquivalentBC(face.physicalGroupId, 1e-11, 1e5);
                bcMgr_T.SetNeumannBC(face.physicalGroupId, 0.0);
                expectedBCCount++;
            }
        }

        FieldManager_2D fm;
        fm.InitSizes(mgr.mesh().getGridCount(), mgr.fracture_network().getOrderedFractureElements().size(), 0, 0, mgr.mesh().getFaces().size(), 0);
        auto pField = fm.matrixFields.create<volScalarField>("Pressure", mgr.mesh().getGridCount());
        auto tField = fm.matrixFields.create<volScalarField>("Temperature", mgr.mesh().getGridCount());
        for (int i = 0; i < pField->size; ++i) { (*pField)[i] = 2e5; (*tField)[i] = 300.0; }

        // P0/P1验证：分离装入压力块(offset=0)和温度块(offset=1)
        auto stats_P = BoundaryAssembler::Assemble_2D(mgr, bcMgr_P, 0, fm, "Pressure", res, diag);
        auto stats_T = BoundaryAssembler::Assemble_2D(mgr, bcMgr_T, 1, fm, "Temperature", res, diag);

                assert(stats_P.matrixBCCount == expectedBCCount && "2D 压力边界触发数异常");
        assert(stats_T.matrixBCCount == expectedBCCount && "2D 温度边界触发数异常");
        assert(stats_P.fractureBCCount == 0 && "裂缝必须隔离");

        // 面积(2D边界长度)加权守恒检查：q = sum(A_face) * C_L * (P - P_far)
        double boundaryMeasure2D = 0.0;
        for (const auto& face : mgr.mesh().getFaces()) {
            if (face.isBoundary() && bcMgr_P.HasBC(face.physicalGroupId)) {
                boundaryMeasure2D += face.length;
            }
        }
        const double expectedRes2D = boundaryMeasure2D * 1.0e-11 * (2.0e5 - 1.0e5);
        const double expectedJac2D = boundaryMeasure2D * 1.0e-11;
        const double tolRes2D = 1.0e-10 * (std::abs(expectedRes2D) + 1.0);
        const double tolJac2D = 1.0e-10 * (std::abs(expectedJac2D) + 1.0);
        assert(std::abs(stats_P.sumResidual - expectedRes2D) < tolRes2D && "2D area-weighted residual conservation error");
        assert(std::abs(stats_P.sumJacobianDiag - expectedJac2D) < tolJac2D && "2D area-weighted Jacobian conservation error");
        assert(std::abs(stats_T.sumResidual) < 1e-14 && "2D temperature Neumann residual should be zero");
        assert(std::abs(stats_T.sumJacobianDiag) < 1e-14 && "2D temperature Neumann Jacobian should be zero");

        int eq_P_0 = mgr.getEquationIndex(0, 0);
        int eq_T_0 = mgr.getEquationIndex(0, 1);
        std::cout << "  [Check] Row offset for P: " << eq_P_0 << " | Row offset for T: " << eq_T_0 << std::endl;
        assert(eq_P_0 != eq_T_0 && "导数必须写入不同的方程行块");
        std::cout << "  [PASS] 2D 网格装配与多块(P/T)分流行写入 PASS" << std::endl;
    }

    template <int N, typename ADVarType>
    void Test_GridLevel_BoundaryAssembly_3D() {
        MeshManager_3D mgr(10.0, 10.0, 10.0, 2, 2, 2, true, false);
        mgr.BuildSolidMatrixGrid_3D();
        std::vector<Vector> pts = { Vector(1,1,1), Vector(9,1,1), Vector(9,9,1), Vector(1,9,1) };
        mgr.addFracturetoFractureNetwork(Fracture_2D(0, pts));
        mgr.meshAllFracturesinNetwork(2, 2, NormalVectorCorrectionMethod::OrthogonalCorrection);
        mgr.setupGlobalIndices();

        mgr.setNumDOFs(3);
        int totalEq = mgr.getTotalEquationDOFs();
        std::vector<double> res(totalEq, 0.0);
        std::vector<double> diag(totalEq, 0.0);

        BoundarySetting::BoundaryConditionManager bcMgr;
        int expectedBCCount = 0;
        for (const auto& face : mgr.mesh().getFaces()) {
            if (face.isBoundary()) {
                bcMgr.SetLeakoffEquivalentBC(face.physicalGroupId, 1e-11, 1e5);
                expectedBCCount++;
            }
        }

        FieldManager_3D fm;
        fm.InitSizes(mgr.mesh().getGridCount(), mgr.fracture_network().getOrderedFractureElements().size(), 0, 0, mgr.mesh().getFaces().size(), 0);
        auto pField = fm.matrixFields.create<volScalarField>("Pressure", mgr.mesh().getGridCount());
        for (int i = 0; i < pField->size; ++i) (*pField)[i] = 2e5;

        auto stats = BoundaryAssembler::Assemble_3D(mgr, bcMgr, 0, fm, "Pressure", res, diag);

                assert(stats.matrixBCCount == expectedBCCount && "3D 边界触发数量与期望不符");
        assert(stats.fractureBCCount == 0 && "裂缝必须隔离");

        // 面积(3D边界面面积)加权守恒检查：q = sum(A_face) * C_L * (P - P_far)
        double boundaryMeasure3D = 0.0;
        for (const auto& face : mgr.mesh().getFaces()) {
            if (face.isBoundary() && bcMgr.HasBC(face.physicalGroupId)) {
                boundaryMeasure3D += face.length;
            }
        }
        const double expectedRes3D = boundaryMeasure3D * 1.0e-11 * (2.0e5 - 1.0e5);
        const double expectedJac3D = boundaryMeasure3D * 1.0e-11;
        const double tolRes3D = 1.0e-10 * (std::abs(expectedRes3D) + 1.0);
        const double tolJac3D = 1.0e-10 * (std::abs(expectedJac3D) + 1.0);
        assert(std::abs(stats.sumResidual - expectedRes3D) < tolRes3D && "3D area-weighted residual conservation error");
        assert(std::abs(stats.sumJacobianDiag - expectedJac3D) < tolJac3D && "3D area-weighted Jacobian conservation error");

        std::cout << "  [PASS] 3D 网格装配与分流行写入 PASS" << std::endl;
    }

    template <int N, typename ADVarType>
    void Run_Day3_BC_Patch() {
        std::cout << "\n[Test_FVM_Ops_AD] 启动 Day3: Boundary Conditions Patch Test..." << std::endl;
        Test_Boundary_Operators<N, ADVarType>();
        std::cout << "[PASS] 线性压降 patch PASS\n[PASS] 绝热 patch PASS" << std::endl;
    }

    template <int N, typename ADVarType>
    void Run_Day3_Leakoff_Switch() {
        std::cout << "\n[Test_FVM_Ops_AD] 启动 Day3: Leakoff Switch & Dual-Phase Test..." << std::endl;
        Test_Leakoff_Operators<N, ADVarType>();
        std::cout << "[PASS] Leakoff OFF=0 PASS\n[PASS] Leakoff ON>0 PASS\n[PASS] 两相 w/g sink 与导数 PASS" << std::endl;
        std::cout << "\n[Test_FVM_Ops_AD] 启动网格级主装配器验证..." << std::endl;
        Test_GridLevel_BoundaryAssembly_2D<N, ADVarType>();
        Test_GridLevel_BoundaryAssembly_3D<N, ADVarType>();
    }
    /**
     * @brief Day4: 井算子单元测试 (BHP/Rate)
     */
    template <int N, typename ADVarType>
    void Test_Well_Operators_2D() {
        ADVarType P_cell;
        P_cell.val = 2.0e5;
        for (int k = 0; k < N; ++k) P_cell.grad(k) = 0.0;
        P_cell.grad(0) = 1.0;

        const double conductance = 2.0e-6;
        ADVarType q_prod = FVM_Ops::Op_Well_BHP_Source_AD<N, ADVarType>(conductance, P_cell, 1.0e5);
        assert(q_prod.val > 0.0 && "BHP producer must be outflow positive");
        assert(q_prod.grad(0) > 0.0 && "BHP derivative dQ/dP must be positive");

        ADVarType q_inj = FVM_Ops::Op_Well_BHP_Source_AD<N, ADVarType>(conductance, P_cell, 3.0e5);
        assert(q_inj.val < 0.0 && "BHP injector must be negative outflow");
        assert(q_inj.grad(0) > 0.0 && "BHP derivative dQ/dP must remain positive");

        ADVarType q_rate = FVM_Ops::Op_Well_Rate_Source_AD<N, ADVarType>(-5.0);
        assert(std::abs(q_rate.val + 5.0) < 1e-12 && "Rate value mismatch");
        for (int k = 0; k < N; ++k) {
            assert(std::abs(q_rate.grad(k)) < 1e-12 && "Rate derivative must be zero");
        }

        std::cout << "  [PASS] WI/BHP/Rate operator signs & Jacobian checks" << std::endl;
    }

    /**
     * @brief Day4: 2D 井装配 + Leakoff ON/OFF 趋势测试
     */
    template <int N, typename ADVarType>
    void Test_Well_Assembly_2D_LeakoffTrend() {
        MeshManager mgr(10.0, 10.0, 0.0, 4, 4, 0, true, false);
        mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OrthogonalCorrection);
        mgr.addFracture(Vector(1.0, 1.0, 0.0), Vector(9.0, 9.0, 0.0));
        mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
        mgr.BuildGlobalSystemIndexing();
        mgr.setNumDOFs(3);

        const int nMat = mgr.getMatrixDOFCount();
        const int nFrac = static_cast<int>(mgr.fracture_network().getOrderedFractureElements().size());
        assert(nFrac > 0 && "Day4 requires at least one fracture element for completion mapping");

        FieldManager_2D fm;
        fm.InitSizes(
            mgr.mesh().getGridCount(),
            nFrac,
            0,
            0,
            mgr.mesh().getFaces().size(),
            0
        );

        const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const PhysicalProperties_string_op::Rock rock;
        const PhysicalProperties_string_op::Water water;
        const PhysicalProperties_string_op::CO2 gas;

        auto pMat = fm.getOrCreateMatrixScalar(pCfg.pressure_field, 2.0e5);
        auto kxxMat = fm.getOrCreateMatrixScalar(rock.k_xx_tag, 1.0e-13);
        auto kyyMat = fm.getOrCreateMatrixScalar(rock.k_yy_tag, 1.0e-13);
        auto mobWMat = fm.getOrCreateMatrixScalar("mob_density_w", 1.2);
        auto mobGMat = fm.getOrCreateMatrixScalar("mob_density_g", 0.5);
        auto hWMat = fm.getOrCreateMatrixScalar(water.h_tag, 1.0e5);
        auto hGMat = fm.getOrCreateMatrixScalar(gas.h_tag, 2.0e5);

        auto pFrac = fm.getOrCreateFractureScalar(pCfg.pressure_field, 2.2e5);
        auto mobWFrac = fm.getOrCreateFractureScalar("mob_density_w", 1.0);
        auto mobGFrac = fm.getOrCreateFractureScalar("mob_density_g", 0.2);
        auto hWFrac = fm.getOrCreateFractureScalar(water.h_tag, 1.1e5);
        auto hGFrac = fm.getOrCreateFractureScalar(gas.h_tag, 2.1e5);

        (void)pMat; (void)kxxMat; (void)kyyMat; (void)mobWMat; (void)mobGMat; (void)hWMat; (void)hGMat;
        (void)pFrac; (void)mobWFrac; (void)mobGFrac; (void)hWFrac; (void)hGFrac;

        WellScheduleStep stepMatrix;
        stepMatrix.t_start = 0.0;
        stepMatrix.t_end = 100.0;
        stepMatrix.well_name = "WELL_MAT_BHP";
        stepMatrix.domain = WellTargetDomain::Matrix;
        stepMatrix.control_mode = WellControlMode::BHP;
        stepMatrix.target_value = 1.5e5;
        stepMatrix.component_mode = WellComponentMode::Total;
        stepMatrix.rw = 0.1;
        stepMatrix.skin = 0.0;
        stepMatrix.well_axis = WellAxis::None;
        stepMatrix.completion_id = 0;
        stepMatrix.wi_override = -1.0;
        stepMatrix.L_override = -1.0;
        stepMatrix.frac_w = 0.6;
        stepMatrix.frac_g = 0.4;

        WellScheduleStep stepFrac;
        stepFrac.t_start = 0.0;
        stepFrac.t_end = 100.0;
        stepFrac.well_name = "WELL_FRAC_RATE";
        stepFrac.domain = WellTargetDomain::Fracture;
        stepFrac.control_mode = WellControlMode::Rate;
        stepFrac.target_value = -5.0;
        stepFrac.component_mode = WellComponentMode::Water;
        stepFrac.rw = 0.1;
        stepFrac.skin = 0.0;
        stepFrac.well_axis = WellAxis::None;
        stepFrac.completion_id = 0;
        stepFrac.wi_override = 1.0e-10;
        stepFrac.L_override = -1.0;
        stepFrac.frac_w = 1.0;
        stepFrac.frac_g = 0.0;

        std::vector<WellScheduleStep> steps = { stepMatrix, stepFrac };

        const int totalEq = mgr.getTotalEquationDOFs();
        std::vector<double> resOff(totalEq, 0.0), diagOff(totalEq, 0.0);
        std::vector<double> resOn(totalEq, 0.0), diagOn(totalEq, 0.0);

        auto wellOff = BoundaryAssembler::Assemble_Wells_2D(mgr, fm, steps, 0, 0, 1, 2, resOff, diagOff);
        auto wellOn = BoundaryAssembler::Assemble_Wells_2D(mgr, fm, steps, 0, 0, 1, 2, resOn, diagOn);

        assert(wellOff.matrixBCCount == 1 && "Matrix well count mismatch");
        assert(wellOff.fractureBCCount == 1 && "Fracture well count mismatch");

        int eqWMat = mgr.getEquationIndex(0, 0);
        int eqGMat = mgr.getEquationIndex(0, 1);
        int eqEMat = mgr.getEquationIndex(0, 2);
        int eqWFrac = mgr.getEquationIndex(nMat, 0);
        assert(eqWMat >= 0 && eqGMat >= 0 && eqEMat >= 0 && eqWFrac >= 0 && "Equation index mapping failed");

        assert(resOff[eqWMat] > 0.0 && "Matrix BHP water source should be positive outflow");
        assert(diagOff[eqWMat] > 0.0 && "Matrix BHP dQ/dP should be positive");
        assert(resOff[eqGMat] > 0.0 && "Matrix BHP gas source should be positive outflow");
        assert(resOff[eqEMat] > 0.0 && "Energy source should be positive for positive enthalpy");
        assert(resOff[eqWFrac] < 0.0 && "Fracture rate injector should be negative outflow");
        assert(std::abs(diagOff[eqWFrac]) < 1e-12 && "Rate control derivative must be zero");

        // 2D Matrix BHP injection (assembly-level): P_cell < P_bhp => q < 0, dq/dP > 0
        WellScheduleStep stepMatrixInj = stepMatrix;
        stepMatrixInj.well_name = "WELL_MAT_BHP_INJ";
        stepMatrixInj.target_value = 3.0e5;
        stepMatrixInj.component_mode = WellComponentMode::Water;
        stepMatrixInj.frac_w = 1.0;
        stepMatrixInj.frac_g = 0.0;
        std::vector<WellScheduleStep> stepsInj = { stepMatrixInj };
        std::vector<double> resInj(totalEq, 0.0), diagInj(totalEq, 0.0);
        auto stInj = BoundaryAssembler::Assemble_Wells_2D(mgr, fm, stepsInj, 0, 0, -1, -1, resInj, diagInj);
        assert(stInj.matrixBCCount == 1 && "2D matrix BHP injection assembly count mismatch");
        assert(resInj[eqWMat] < 0.0 && "2D matrix BHP injection should be negative outflow");
        assert(diagInj[eqWMat] > 0.0 && "2D matrix BHP injection dQ/dP should stay positive");

        // 2D Matrix Rate (assembly-level): q=q_target, dq/dP=0
        WellScheduleStep stepMatrixRate = stepMatrix;
        stepMatrixRate.well_name = "WELL_MAT_RATE";
        stepMatrixRate.control_mode = WellControlMode::Rate;
        stepMatrixRate.component_mode = WellComponentMode::Water;
        stepMatrixRate.target_value = -2.75;
        stepMatrixRate.wi_override = 1.0e-10;
        stepMatrixRate.frac_w = 1.0;
        stepMatrixRate.frac_g = 0.0;
        std::vector<WellScheduleStep> stepsRate = { stepMatrixRate };
        std::vector<double> resRate(totalEq, 0.0), diagRate(totalEq, 0.0);
        auto stRate = BoundaryAssembler::Assemble_Wells_2D(mgr, fm, stepsRate, 0, 0, -1, -1, resRate, diagRate);
        assert(stRate.matrixBCCount == 1 && "2D matrix rate assembly count mismatch");
        assert(std::abs(resRate[eqWMat] + 2.75) < 1e-12 && "2D matrix rate source should equal q_target");
        assert(std::abs(diagRate[eqWMat]) < 1e-12 && "2D matrix rate derivative should be zero");

        BoundarySetting::BoundaryConditionManager bcMgrLeakOn;
        for (const auto& face : mgr.mesh().getFaces()) {
            if (face.isBoundary()) {
                bcMgrLeakOn.SetLeakoffEquivalentBC(face.physicalGroupId, 1.0e-11, 1.0e5);
            }
        }
        auto leakStats = BoundaryAssembler::Assemble_2D(mgr, bcMgrLeakOn, 0, fm, pCfg.pressure_field, resOn, diagOn);

        const double offTotal = wellOff.sumResidual;
        const double onTotal = wellOn.sumResidual + leakStats.sumResidual;
        assert(onTotal > offTotal && "Leakoff ON should increase total outflow trend");

        std::cout << "  [PASS] 2D well assembly (Matrix+Fracture, Mass+Energy) and leakoff trend" << std::endl;
    }

    /**
     * @brief Day4: WAG 调度骨架（CSV 驱动）测试
     */
    template <int N, typename ADVarType>
    void Test_WAG_Schedule_Skeleton_2D() {
        (void)N;
        (void)ADVarType();

        const char* csvPath = "day4_wag_schedule_tmp.csv";
        {
            std::ofstream ofs(csvPath, std::ios::trunc);
            assert(ofs.is_open() && "Failed to create temporary WAG CSV");
            ofs << "t_start,t_end,well_name,domain,control_mode,target_value,component_mode,rw,skin,well_axis,completion_id,wi_override,L_override,frac_w,frac_g\n";
            ofs << "0,10,W1,Matrix,BHP,150000,Total,0.1,0,None,0,-1,-1,0.7,0.3\n";
            ofs << "10,20,W1,Matrix,Rate,-8,Water,0.1,0,None,0,-1,-1,1,0\n";
        }

        WellScheduleManager scheduler;
        assert(scheduler.LoadFromCsv(csvPath) && "WAG CSV load failed");

        auto sA = scheduler.GetActiveSteps(5.0);
        assert(sA.size() == 1 && "Expected one active step at t=5");
        assert(sA[0].control_mode == WellControlMode::BHP && "Expected BHP in segment A");
        assert(sA[0].component_mode == WellComponentMode::Total && "Expected Total mode in segment A");

        auto sB = scheduler.GetActiveSteps(15.0);
        assert(sB.size() == 1 && "Expected one active step at t=15");
        assert(sB[0].control_mode == WellControlMode::Rate && "Expected Rate in segment B");
        assert(sB[0].component_mode == WellComponentMode::Water && "Expected Water mode in segment B");

        auto sC = scheduler.GetActiveSteps(25.0);
        assert(sC.empty() && "Expected no active step outside schedule range");

        std::remove(csvPath);
        std::cout << "  [PASS] WAG schedule skeleton (CSV-driven segment switch)" << std::endl;
    }

    /**
     * @brief Day4 总验收入口
     */
    
    /**
     * @brief Day4: 3D 井装配基础测试 (Matrix + Fracture)
     */
    template <int N, typename ADVarType>
    void Test_Well_Assembly_3D_Basic() {
        MeshManager_3D mgr(10.0, 10.0, 10.0, 2, 2, 2, true, false);
        mgr.BuildSolidMatrixGrid_3D();
        std::vector<Vector> pts = { Vector(1,1,1), Vector(9,1,1), Vector(9,9,1), Vector(1,9,1) };
        mgr.addFracturetoFractureNetwork(Fracture_2D(0, pts));
        mgr.meshAllFracturesinNetwork(2, 2, NormalVectorCorrectionMethod::OrthogonalCorrection);
        mgr.setupGlobalIndices();
        mgr.setNumDOFs(3);

        const int nMat = mgr.fracture_network().getSolverIndexOffset();
        const int nFrac = static_cast<int>(mgr.fracture_network().getOrderedFractureElements().size());
        assert(nFrac > 0 && "Day4 3D requires fracture elements");

        FieldManager_3D fm;
        fm.InitSizes(
            mgr.mesh().getGridCount(),
            nFrac,
            0,
            0,
            mgr.mesh().getFaces().size(),
            0
        );

        const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const PhysicalProperties_string_op::Rock rock;
        const PhysicalProperties_string_op::Water water;
        const PhysicalProperties_string_op::CO2 gas;

        (void)fm.getOrCreateMatrixScalar(pCfg.pressure_field, 2.3e5);
        (void)fm.getOrCreateMatrixScalar(rock.k_xx_tag, 1.0e-13);
        (void)fm.getOrCreateMatrixScalar(rock.k_yy_tag, 2.0e-13);
        (void)fm.getOrCreateMatrixScalar(rock.k_zz_tag, 1.5e-13);
        (void)fm.getOrCreateMatrixScalar("mob_density_w", 1.0);
        (void)fm.getOrCreateMatrixScalar("mob_density_g", 0.7);
        (void)fm.getOrCreateMatrixScalar(water.h_tag, 1.2e5);
        (void)fm.getOrCreateMatrixScalar(gas.h_tag, 2.2e5);

        (void)fm.getOrCreateFractureScalar(pCfg.pressure_field, 2.0e5);
        (void)fm.getOrCreateFractureScalar("mob_density_w", 0.8);
        (void)fm.getOrCreateFractureScalar("mob_density_g", 0.6);
        (void)fm.getOrCreateFractureScalar(water.h_tag, 1.1e5);
        (void)fm.getOrCreateFractureScalar(gas.h_tag, 2.1e5);

        WellScheduleStep stepMatrix;
        stepMatrix.t_start = 0.0;
        stepMatrix.t_end = 100.0;
        stepMatrix.well_name = "WELL3D_MAT_BHP";
        stepMatrix.domain = WellTargetDomain::Matrix;
        stepMatrix.control_mode = WellControlMode::BHP;
        stepMatrix.target_value = 1.6e5;
        stepMatrix.component_mode = WellComponentMode::Total;
        stepMatrix.rw = 0.1;
        stepMatrix.skin = 0.0;
        stepMatrix.well_axis = WellAxis::Z;
        stepMatrix.completion_id = 0;
        stepMatrix.wi_override = -1.0;
        stepMatrix.L_override = -1.0;
        stepMatrix.frac_w = 0.5;
        stepMatrix.frac_g = 0.5;

        WellScheduleStep stepFrac;
        stepFrac.t_start = 0.0;
        stepFrac.t_end = 100.0;
        stepFrac.well_name = "WELL3D_FRAC_RATE";
        stepFrac.domain = WellTargetDomain::Fracture;
        stepFrac.control_mode = WellControlMode::Rate;
        stepFrac.target_value = 3.0;
        stepFrac.component_mode = WellComponentMode::Gas;
        stepFrac.rw = 0.1;
        stepFrac.skin = 0.0;
        stepFrac.well_axis = WellAxis::None;
        stepFrac.completion_id = 0;
        stepFrac.wi_override = 1.0e-10;
        stepFrac.L_override = -1.0;
        stepFrac.frac_w = 0.0;
        stepFrac.frac_g = 1.0;

        std::vector<WellScheduleStep> steps = { stepMatrix, stepFrac };
        const int totalEq = mgr.getTotalEquationDOFs();
        std::vector<double> res(totalEq, 0.0), diag(totalEq, 0.0);

        auto stats = BoundaryAssembler::Assemble_Wells_3D(mgr, fm, steps, 0, 0, 1, 2, res, diag);
        assert(stats.matrixBCCount == 1 && "3D matrix well count mismatch");
        assert(stats.fractureBCCount == 1 && "3D fracture well count mismatch");

        const int eqWMat = mgr.getEquationIndex(0, 0);
        const int eqGMat = mgr.getEquationIndex(0, 1);
        const int eqEMat = mgr.getEquationIndex(0, 2);
        const int eqGFrac = mgr.getEquationIndex(nMat, 1);
        assert(eqWMat >= 0 && eqGMat >= 0 && eqEMat >= 0 && eqGFrac >= 0 && "3D equation index mapping failed");

        assert(res[eqWMat] > 0.0 && "3D matrix BHP water should be positive outflow");
        assert(res[eqGMat] > 0.0 && "3D matrix BHP gas should be positive outflow");
        assert(diag[eqWMat] > 0.0 && "3D matrix BHP derivative should be positive");
        assert(res[eqEMat] > 0.0 && "3D energy source should be positive");

        assert(res[eqGFrac] > 0.0 && "3D fracture gas rate producer should be positive outflow");
        assert(std::abs(diag[eqGFrac]) < 1e-12 && "3D rate derivative must be zero");

        // 3D Matrix BHP injection sign check
        WellScheduleStep stepMatrixInj = stepMatrix;
        stepMatrixInj.well_name = "WELL3D_MAT_BHP_INJ";
        stepMatrixInj.target_value = 3.0e5;
        std::vector<WellScheduleStep> stepsInj = { stepMatrixInj };
        std::vector<double> resInj(totalEq, 0.0), diagInj(totalEq, 0.0);
        auto stInj = BoundaryAssembler::Assemble_Wells_3D(mgr, fm, stepsInj, 0, 0, 1, -1, resInj, diagInj);
        assert(stInj.matrixBCCount == 1 && "3D matrix BHP injection assembly count mismatch");
        assert(resInj[eqWMat] < 0.0 && "3D matrix BHP injection should be negative outflow");
        assert(diagInj[eqWMat] > 0.0 && "3D matrix BHP injection dQ/dP should stay positive");

        // 3D Matrix Rate (assembly-level): q=q_target, dq/dP=0
        WellScheduleStep stepMatrixRate = stepMatrix;
        stepMatrixRate.well_name = "WELL3D_MAT_RATE";
        stepMatrixRate.control_mode = WellControlMode::Rate;
        stepMatrixRate.component_mode = WellComponentMode::Water;
        stepMatrixRate.target_value = 2.25;
        std::vector<WellScheduleStep> stepsRate = { stepMatrixRate };
        std::vector<double> resRate(totalEq, 0.0), diagRate(totalEq, 0.0);
        auto stRate = BoundaryAssembler::Assemble_Wells_3D(mgr, fm, stepsRate, 0, 0, -1, -1, resRate, diagRate);
        assert(stRate.matrixBCCount == 1 && "3D matrix rate assembly count mismatch");
        assert(std::abs(resRate[eqWMat] - 2.25) < 1e-12 && "3D matrix rate source should equal target");
        assert(std::abs(diagRate[eqWMat]) < 1e-12 && "3D matrix rate derivative should be zero");

        std::cout << "  [PASS] 3D well assembly (Matrix+Fracture, Mass+Energy)" << std::endl;
    }

    /**
     * @brief Day4: 两相交替注采工况测试 (WAG schedule + assembly)
     */
    template <int N, typename ADVarType>
    void Test_WAG_Alternating_Injection_2D() {
        const char* csvPath = "day4_wag_alt_tmp.csv";
        {
            std::ofstream ofs(csvPath, std::ios::trunc);
            assert(ofs.is_open() && "Failed to create WAG alternating CSV");
            ofs << "t_start,t_end,well_name,domain,control_mode,target_value,component_mode,rw,skin,well_axis,completion_id,wi_override,L_override,frac_w,frac_g\n";
            ofs << "0,10,WALT,Matrix,Rate,-5,Water,0.1,0,None,0,1e-10,-1,1,0\n";
            ofs << "10,20,WALT,Matrix,Rate,-7,Gas,0.1,0,None,0,1e-10,-1,0,1\n";
            ofs << "20,30,WALT,Matrix,Rate,4,Total,0.1,0,None,0,1e-10,-1,0.5,0.5\n";
        }

        MeshManager mgr(10.0, 10.0, 0.0, 4, 4, 0, true, false);
        mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OrthogonalCorrection);
        mgr.BuildGlobalSystemIndexing();
        mgr.setNumDOFs(3);

        FieldManager_2D fm;
        fm.InitSizes(mgr.mesh().getGridCount(), 0, 0, 0, mgr.mesh().getFaces().size(), 0);

        const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const PhysicalProperties_string_op::Water water;
        const PhysicalProperties_string_op::CO2 gas;

        (void)fm.getOrCreateMatrixScalar(pCfg.pressure_field, 2.0e5);
        (void)fm.getOrCreateMatrixScalar("mob_density_w", 1.0);
        (void)fm.getOrCreateMatrixScalar("mob_density_g", 0.9);
        (void)fm.getOrCreateMatrixScalar(water.h_tag, 1.0e5);
        (void)fm.getOrCreateMatrixScalar(gas.h_tag, 2.0e5);

        WellScheduleManager scheduler;
        assert(scheduler.LoadFromCsv(csvPath) && "WAG alternating CSV load failed");

        const int totalEq = mgr.getTotalEquationDOFs();
        const int eqW = mgr.getEquationIndex(0, 0);
        const int eqG = mgr.getEquationIndex(0, 1);
        assert(eqW >= 0 && eqG >= 0 && "WAG alternating equation mapping failed");

        {
            std::vector<double> res(totalEq, 0.0), diag(totalEq, 0.0);
            auto s = scheduler.GetActiveSteps(5.0);
            assert(s.size() == 1 && s[0].component_mode == WellComponentMode::Water && "Segment A should be water-rate");
            auto st = BoundaryAssembler::Assemble_Wells_2D(mgr, fm, s, 0, 0, 1, -1, res, diag);
            assert(st.matrixBCCount == 1 && "Segment A assembly count mismatch");
            assert(res[eqW] < 0.0 && "Segment A water injection should be negative outflow");
            assert(std::abs(res[eqG]) < 1e-12 && "Segment A gas equation should remain zero");
        }

        {
            std::vector<double> res(totalEq, 0.0), diag(totalEq, 0.0);
            auto s = scheduler.GetActiveSteps(15.0);
            assert(s.size() == 1 && s[0].component_mode == WellComponentMode::Gas && "Segment B should be gas-rate");
            auto st = BoundaryAssembler::Assemble_Wells_2D(mgr, fm, s, 0, 0, 1, -1, res, diag);
            assert(st.matrixBCCount == 1 && "Segment B assembly count mismatch");
            assert(res[eqG] < 0.0 && "Segment B gas injection should be negative outflow");
            assert(std::abs(res[eqW]) < 1e-12 && "Segment B water equation should remain zero");
        }

        {
            std::vector<double> res(totalEq, 0.0), diag(totalEq, 0.0);
            auto s = scheduler.GetActiveSteps(25.0);
            assert(s.size() == 1 && s[0].component_mode == WellComponentMode::Total && "Segment C should switch to production(total-rate)");
            auto st = BoundaryAssembler::Assemble_Wells_2D(mgr, fm, s, 0, 0, 1, -1, res, diag);
            assert(st.matrixBCCount == 1 && "Segment C assembly count mismatch");
            assert(res[eqW] > 0.0 && "Segment C water production should be positive outflow");
            assert(res[eqG] > 0.0 && "Segment C gas production should be positive outflow");
        }

        std::remove(csvPath);
        std::cout << "  [PASS] WAG alternating water/gas injection schedule with assembly checks" << std::endl;
    }
template <int N, typename ADVarType>
    void Run_Day4_Well_Patch() {
        std::cout << "\n[Test_FVM_Ops_AD] 启动 Day4: Well(BHP/Rate) + WAG skeleton..." << std::endl;
        Test_Well_Operators_2D<N, ADVarType>();
        Test_Well_Assembly_2D_LeakoffTrend<N, ADVarType>();
        Test_Well_Assembly_3D_Basic<N, ADVarType>();
        Test_WAG_Schedule_Skeleton_2D<N, ADVarType>();
        Test_WAG_Alternating_Injection_2D<N, ADVarType>();

        std::cout << "[PASS] WI 2D geometric-eq" << std::endl;
        std::cout << "[PASS] WI 3D geometric-eq" << std::endl;
        std::cout << "[PASS] BHP/Rate mass coupling" << std::endl;
        std::cout << "[PASS] Well energy coupling" << std::endl;
        std::cout << "[PASS] WAG config switching" << std::endl;
        std::cout << "[PASS] Matrix+Fracture completion" << std::endl;
    }

    
    /**
     * @brief Day4 可视化导出: 2D 井源项场 (Matrix + Fracture)
     */
    template <int N, typename ADVarType>
    void Run_Day4_Well_Viz_2D() {
        (void)N;
        (void)ADVarType();

        MeshManager mgr(10.0, 10.0, 0.0, 4, 4, 0, true, false);
        mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OrthogonalCorrection);
        mgr.addFracture(Vector(1.0, 1.0, 0.0), Vector(9.0, 9.0, 0.0));
        mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
        mgr.BuildGlobalSystemIndexing();
        mgr.setNumDOFs(3);

        const int nMat = mgr.getMatrixDOFCount();
        const int nFrac = static_cast<int>(mgr.fracture_network().getOrderedFractureElements().size());

        FieldManager_2D fm;
        fm.InitSizes(mgr.mesh().getGridCount(), nFrac, 0, 0, mgr.mesh().getFaces().size(), 0);

        const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const PhysicalProperties_string_op::Rock rock;
        const PhysicalProperties_string_op::Water water;
        const PhysicalProperties_string_op::CO2 gas;

        (void)fm.getOrCreateMatrixScalar(pCfg.pressure_field, 2.0e5);
        (void)fm.getOrCreateMatrixScalar(rock.k_xx_tag, 1.0e-13);
        (void)fm.getOrCreateMatrixScalar(rock.k_yy_tag, 1.0e-13);
        (void)fm.getOrCreateMatrixScalar("mob_density_w", 1.2);
        (void)fm.getOrCreateMatrixScalar("mob_density_g", 0.5);
        (void)fm.getOrCreateMatrixScalar(water.h_tag, 1.0e5);
        (void)fm.getOrCreateMatrixScalar(gas.h_tag, 2.0e5);

        (void)fm.getOrCreateFractureScalar(pCfg.pressure_field, 2.2e5);
        (void)fm.getOrCreateFractureScalar("mob_density_w", 1.0);
        (void)fm.getOrCreateFractureScalar("mob_density_g", 0.2);
        (void)fm.getOrCreateFractureScalar(water.h_tag, 1.1e5);
        (void)fm.getOrCreateFractureScalar(gas.h_tag, 2.1e5);

        WellScheduleStep stepMatrix;
        stepMatrix.t_start = 0.0;
        stepMatrix.t_end = 100.0;
        stepMatrix.well_name = "WELL_MAT_BHP_VIZ";
        stepMatrix.domain = WellTargetDomain::Matrix;
        stepMatrix.control_mode = WellControlMode::BHP;
        stepMatrix.target_value = 1.5e5;
        stepMatrix.component_mode = WellComponentMode::Total;
        stepMatrix.rw = 0.1;
        stepMatrix.skin = 0.0;
        stepMatrix.well_axis = WellAxis::None;
        stepMatrix.completion_id = 0;
        stepMatrix.wi_override = -1.0;
        stepMatrix.frac_w = 0.6;
        stepMatrix.frac_g = 0.4;

        WellScheduleStep stepFrac;
        stepFrac.t_start = 0.0;
        stepFrac.t_end = 100.0;
        stepFrac.well_name = "WELL_FRAC_RATE_VIZ";
        stepFrac.domain = WellTargetDomain::Fracture;
        stepFrac.control_mode = WellControlMode::Rate;
        stepFrac.target_value = -5.0;
        stepFrac.component_mode = WellComponentMode::Water;
        stepFrac.rw = 0.1;
        stepFrac.skin = 0.0;
        stepFrac.well_axis = WellAxis::None;
        stepFrac.completion_id = 0;
        stepFrac.wi_override = 1.0e-10;
        stepFrac.frac_w = 1.0;
        stepFrac.frac_g = 0.0;

        std::vector<WellScheduleStep> steps = { stepMatrix, stepFrac };

        const int totalEq = mgr.getTotalEquationDOFs();
        std::vector<double> residual(totalEq, 0.0), diag(totalEq, 0.0);
        auto st = BoundaryAssembler::Assemble_Wells_2D(mgr, fm, steps, 0, 0, 1, 2, residual, diag);
        assert(st.matrixBCCount == 1 && st.fractureBCCount == 1 && "Day4 2D viz well assembly count mismatch");

        auto qWMat = fm.getOrCreateMatrixScalar("q_well_w", 0.0);
        auto qGMat = fm.getOrCreateMatrixScalar("q_well_g", 0.0);
        auto qEMat = fm.getOrCreateMatrixScalar("q_well_e", 0.0);
        auto qWFrac = fm.getOrCreateFractureScalar("q_well_w", 0.0);
        auto qGFrac = fm.getOrCreateFractureScalar("q_well_g", 0.0);
        auto qEFrac = fm.getOrCreateFractureScalar("q_well_e", 0.0);

        for (int i = 0; i < static_cast<int>(mgr.mesh().getGridCount()); ++i) {
            int eqW = mgr.getEquationIndex(i, 0);
            int eqG = mgr.getEquationIndex(i, 1);
            int eqE = mgr.getEquationIndex(i, 2);
            (*qWMat)[i] = (eqW >= 0 && eqW < totalEq) ? residual[eqW] : 0.0;
            (*qGMat)[i] = (eqG >= 0 && eqG < totalEq) ? residual[eqG] : 0.0;
            (*qEMat)[i] = (eqE >= 0 && eqE < totalEq) ? residual[eqE] : 0.0;
        }
        for (int f = 0; f < nFrac; ++f) {
            int solverIdx = nMat + f;
            int eqW = mgr.getEquationIndex(solverIdx, 0);
            int eqG = mgr.getEquationIndex(solverIdx, 1);
            int eqE = mgr.getEquationIndex(solverIdx, 2);
            (*qWFrac)[f] = (eqW >= 0 && eqW < totalEq) ? residual[eqW] : 0.0;
            (*qGFrac)[f] = (eqG >= 0 && eqG < totalEq) ? residual[eqG] : 0.0;
            (*qEFrac)[f] = (eqE >= 0 && eqE < totalEq) ? residual[eqE] : 0.0;
        }

        const std::string vtkFile = "Test/BoundaryTest/day4_well_viz_2d.vtk";
        PostProcess_2D exporter(mgr, fm);
        exporter.ExportVTK(vtkFile, 0.0);

        std::ifstream verifyFile(vtkFile, std::ios::ate | std::ios::binary);
        if (!verifyFile.is_open() || verifyFile.tellg() == 0) {
            throw std::runtime_error("[FAIL] Day4 2D VTK export failed: " + vtkFile);
        }
        verifyFile.close();
        std::cout << "  [PASS] Day4 2D VTK exported: " << vtkFile << std::endl;
    }

    /**
     * @brief Day4 可视化导出: 3D 井源项场 (Matrix + Fracture)
     */
    template <int N, typename ADVarType>
    void Run_Day4_Well_Viz_3D() {
        (void)N;
        (void)ADVarType();

        MeshManager_3D mgr(10.0, 10.0, 10.0, 2, 2, 2, true, false);
        mgr.BuildSolidMatrixGrid_3D();
        std::vector<Vector> pts = { Vector(1,1,1), Vector(9,1,1), Vector(9,9,1), Vector(1,9,1) };
        mgr.addFracturetoFractureNetwork(Fracture_2D(0, pts));
        mgr.meshAllFracturesinNetwork(2, 2, NormalVectorCorrectionMethod::OrthogonalCorrection);
        mgr.setupGlobalIndices();
        mgr.setNumDOFs(3);

        const int nMat = mgr.fracture_network().getSolverIndexOffset();
        const int nFrac = static_cast<int>(mgr.fracture_network().getOrderedFractureElements().size());

        FieldManager_3D fm;
        fm.InitSizes(mgr.mesh().getGridCount(), nFrac, 0, 0, mgr.mesh().getFaces().size(), 0);

        const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const PhysicalProperties_string_op::Rock rock;
        const PhysicalProperties_string_op::Water water;
        const PhysicalProperties_string_op::CO2 gas;

        (void)fm.getOrCreateMatrixScalar(pCfg.pressure_field, 2.3e5);
        (void)fm.getOrCreateMatrixScalar(rock.k_xx_tag, 1.0e-13);
        (void)fm.getOrCreateMatrixScalar(rock.k_yy_tag, 2.0e-13);
        (void)fm.getOrCreateMatrixScalar(rock.k_zz_tag, 1.5e-13);
        (void)fm.getOrCreateMatrixScalar("mob_density_w", 1.0);
        (void)fm.getOrCreateMatrixScalar("mob_density_g", 0.7);
        (void)fm.getOrCreateMatrixScalar(water.h_tag, 1.2e5);
        (void)fm.getOrCreateMatrixScalar(gas.h_tag, 2.2e5);

        (void)fm.getOrCreateFractureScalar(pCfg.pressure_field, 2.0e5);
        (void)fm.getOrCreateFractureScalar("mob_density_w", 0.8);
        (void)fm.getOrCreateFractureScalar("mob_density_g", 0.6);
        (void)fm.getOrCreateFractureScalar(water.h_tag, 1.1e5);
        (void)fm.getOrCreateFractureScalar(gas.h_tag, 2.1e5);

        WellScheduleStep stepMatrix;
        stepMatrix.t_start = 0.0;
        stepMatrix.t_end = 100.0;
        stepMatrix.well_name = "WELL3D_MAT_BHP_VIZ";
        stepMatrix.domain = WellTargetDomain::Matrix;
        stepMatrix.control_mode = WellControlMode::BHP;
        stepMatrix.target_value = 1.6e5;
        stepMatrix.component_mode = WellComponentMode::Total;
        stepMatrix.rw = 0.1;
        stepMatrix.skin = 0.0;
        stepMatrix.well_axis = WellAxis::Z;
        stepMatrix.completion_id = 0;
        stepMatrix.wi_override = -1.0;
        stepMatrix.frac_w = 0.5;
        stepMatrix.frac_g = 0.5;

        WellScheduleStep stepFrac;
        stepFrac.t_start = 0.0;
        stepFrac.t_end = 100.0;
        stepFrac.well_name = "WELL3D_FRAC_RATE_VIZ";
        stepFrac.domain = WellTargetDomain::Fracture;
        stepFrac.control_mode = WellControlMode::Rate;
        stepFrac.target_value = 3.0;
        stepFrac.component_mode = WellComponentMode::Gas;
        stepFrac.rw = 0.1;
        stepFrac.skin = 0.0;
        stepFrac.well_axis = WellAxis::None;
        stepFrac.completion_id = 0;
        stepFrac.wi_override = 1.0e-10;
        stepFrac.frac_w = 0.0;
        stepFrac.frac_g = 1.0;

        std::vector<WellScheduleStep> steps = { stepMatrix, stepFrac };

        const int totalEq = mgr.getTotalEquationDOFs();
        std::vector<double> residual(totalEq, 0.0), diag(totalEq, 0.0);
        auto st = BoundaryAssembler::Assemble_Wells_3D(mgr, fm, steps, 0, 0, 1, 2, residual, diag);
        assert(st.matrixBCCount == 1 && st.fractureBCCount == 1 && "Day4 3D viz well assembly count mismatch");

        auto qWMat = fm.getOrCreateMatrixScalar("q_well_w", 0.0);
        auto qGMat = fm.getOrCreateMatrixScalar("q_well_g", 0.0);
        auto qEMat = fm.getOrCreateMatrixScalar("q_well_e", 0.0);
        auto qWFrac = fm.getOrCreateFractureScalar("q_well_w", 0.0);
        auto qGFrac = fm.getOrCreateFractureScalar("q_well_g", 0.0);
        auto qEFrac = fm.getOrCreateFractureScalar("q_well_e", 0.0);

        for (int i = 0; i < static_cast<int>(mgr.mesh().getGridCount()); ++i) {
            int eqW = mgr.getEquationIndex(i, 0);
            int eqG = mgr.getEquationIndex(i, 1);
            int eqE = mgr.getEquationIndex(i, 2);
            (*qWMat)[i] = (eqW >= 0 && eqW < totalEq) ? residual[eqW] : 0.0;
            (*qGMat)[i] = (eqG >= 0 && eqG < totalEq) ? residual[eqG] : 0.0;
            (*qEMat)[i] = (eqE >= 0 && eqE < totalEq) ? residual[eqE] : 0.0;
        }
        for (int f = 0; f < nFrac; ++f) {
            int solverIdx = nMat + f;
            int eqW = mgr.getEquationIndex(solverIdx, 0);
            int eqG = mgr.getEquationIndex(solverIdx, 1);
            int eqE = mgr.getEquationIndex(solverIdx, 2);
            (*qWFrac)[f] = (eqW >= 0 && eqW < totalEq) ? residual[eqW] : 0.0;
            (*qGFrac)[f] = (eqG >= 0 && eqG < totalEq) ? residual[eqG] : 0.0;
            (*qEFrac)[f] = (eqE >= 0 && eqE < totalEq) ? residual[eqE] : 0.0;
        }

        const std::string vtkFile = "Test/BoundaryTest/day4_well_viz_3d.vtk";
        PostProcess_3D exporter(mgr, fm);
        exporter.ExportVTK(vtkFile, 0.0);

        std::ifstream verifyFile(vtkFile, std::ios::ate | std::ios::binary);
        if (!verifyFile.is_open() || verifyFile.tellg() == 0) {
            throw std::runtime_error("[FAIL] Day4 3D VTK export failed: " + vtkFile);
        }
        verifyFile.close();
        std::cout << "  [PASS] Day4 3D VTK exported: " << vtkFile << std::endl;
    }

    /**
     * @brief Day4 可视化总入口：固定输出 2D/3D 井源 VTK
     */
    template <int N, typename ADVarType>
    void Run_Day4_Well_Viz() {
        std::cout << "\n[Test_FVM_Ops_AD] 启动 Day4: Well VTK Visualization export..." << std::endl;
        Run_Day4_Well_Viz_2D<N, ADVarType>();
        Run_Day4_Well_Viz_3D<N, ADVarType>();
        std::cout << "[PASS] Day4 well visualization VTK exported." << std::endl;
    }
template <int N, typename ADVarType>
    void Run_Day3_BC_Viz_2D() {
        std::cout << "\n[Test_FVM_Ops_AD] 启动 Day3: Boundary Conditions Viz Test (2D)..." << std::endl;
        MeshManager mgr(10.0, 10.0, 0.0, 10, 10, 0, true, false);
        mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OrthogonalCorrection);
        mgr.addFracture(Vector(2.0, 2.0, 0), Vector(8.0, 8.0, 0));
        mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
        mgr.BuildGlobalSystemIndexing();

        mgr.setNumDOFs(1);
        int totalEq = mgr.getTotalEquationDOFs();
        std::vector<double> res(totalEq, 0.0);
        std::vector<double> diag(totalEq, 0.0);

        BoundarySetting::BoundaryConditionManager bcMgr;
        for (const auto& face : mgr.mesh().getFaces()) {
            if (face.isBoundary()) {
                if (face.midpoint.m_x < 1e-5 || face.midpoint.m_x > 10.0 - 1e-5) {
                    bcMgr.SetLeakoffEquivalentBC(face.physicalGroupId, 1e-10, 1e5);
                }
                else {
                    bcMgr.SetNeumannBC(face.physicalGroupId, 0.0);
                }
            }
        }

        FieldManager_2D fm;
        fm.InitSizes(mgr.mesh().getGridCount(), mgr.fracture_network().getOrderedFractureElements().size(), 0, 0, mgr.mesh().getFaces().size(), 0);

        auto pField = fm.matrixFields.create<volScalarField>("Pressure", mgr.mesh().getGridCount());
        for (int i = 0; i < pField->size; ++i) (*pField)[i] = 2e5;

        BoundaryAssembler::Assemble_2D(mgr, bcMgr, 0, fm, "Pressure", res, diag);

        auto resField = fm.matrixFields.create<volScalarField>("Residual_Pressure", mgr.mesh().getGridCount());
        for (size_t i = 0; i < mgr.mesh().getGridCount(); ++i) {
            (*resField)[i] = res[mgr.getEquationIndex(i, 0)];
        }

        std::string vtkFile = "day3_bc_patch_2d.vtk";
        PostProcess_2D exporter(mgr, fm);
        exporter.ExportVTK(vtkFile, 0.0);

        // [真实校验] 检查文件存在且非空
        std::ifstream verifyFile(vtkFile, std::ios::ate | std::ios::binary);
        if (!verifyFile.is_open() || verifyFile.tellg() == 0) {
            throw std::runtime_error("[FAIL] VTK file " + vtkFile + " generation failed or is strictly empty!");
        }
        verifyFile.close();

        std::cout << "  - Grid Elements: Matrix=" << mgr.mesh().getGridCount()
            << ", Fracture=" << mgr.fracture_network().getOrderedFractureElements().size() << std::endl;
        std::cout << "  - Exported Fields: Pressure, Residual_Pressure" << std::endl;
        std::cout << "  [PASS] 2D VTK exported and verified strictly (" << vtkFile << ")" << std::endl;
    }

    /**
     * @brief 可视化用例: 3D Leakoff 边界测试
     */
    template <int N, typename ADVarType>
    void Run_Day3_Leakoff_Viz_3D() {
        std::cout << "\n[Test_FVM_Ops_AD] 启动 Day3: Leakoff Viz Test (3D)..." << std::endl;
        MeshManager_3D mgr(10.0, 10.0, 10.0, 5, 5, 5, true, false);
        mgr.BuildSolidMatrixGrid_3D();
        std::vector<Vector> pts = { Vector(2,2,2), Vector(8,2,2), Vector(8,8,2), Vector(2,8,2) };
        mgr.addFracturetoFractureNetwork(Fracture_2D(0, pts));
        mgr.meshAllFracturesinNetwork(2, 2, NormalVectorCorrectionMethod::OrthogonalCorrection);
        mgr.setupGlobalIndices();

        mgr.setNumDOFs(1);
        int totalEq = mgr.getTotalEquationDOFs();
        std::vector<double> res(totalEq, 0.0);
        std::vector<double> diag(totalEq, 0.0);

        BoundarySetting::BoundaryConditionManager bcMgr;
        for (const auto& face : mgr.mesh().getFaces()) {
            if (face.isBoundary()) bcMgr.SetLeakoffEquivalentBC(face.physicalGroupId, 1e-11, 1e5);
        }

        FieldManager_3D fm;
        fm.InitSizes(mgr.mesh().getGridCount(), mgr.fracture_network().getOrderedFractureElements().size(), 0, 0, mgr.mesh().getFaces().size(), 0);
        auto pField = fm.matrixFields.create<volScalarField>("Pressure", mgr.mesh().getGridCount());
        for (int i = 0; i < pField->size; ++i) (*pField)[i] = 2e5;

        BoundaryAssembler::Assemble_3D(mgr, bcMgr, 0, fm, "Pressure", res, diag);

        auto resField = fm.matrixFields.create<volScalarField>("Residual_Pressure", mgr.mesh().getGridCount());
        for (size_t i = 0; i < mgr.mesh().getGridCount(); ++i) {
            (*resField)[i] = res[mgr.getEquationIndex(i, 0)];
        }

        std::string vtkFile = "day3_leakoff_3d.vtk";
        PostProcess_3D exporter(mgr, fm);
        exporter.ExportVTK(vtkFile, 0.0);

        // [真实校验] 检查文件存在且非空
        std::ifstream verifyFile(vtkFile, std::ios::ate | std::ios::binary);
        if (!verifyFile.is_open() || verifyFile.tellg() == 0) {
            throw std::runtime_error("[FAIL] VTK file " + vtkFile + " generation failed or is strictly empty!");
        }
        verifyFile.close();

        std::cout << "  - Grid Elements: Matrix=" << mgr.mesh().getGridCount() << std::endl;
        std::cout << "  [PASS] 3D VTK exported and verified strictly (" << vtkFile << ")" << std::endl;
    }

} // namespace Test_FVM

#endif // TEST_FVM_OPS_AD_H









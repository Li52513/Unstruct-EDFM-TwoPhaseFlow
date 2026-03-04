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
#include "FVM_Ops_AD.h"
#include "Variable3D.hpp" 
#include "MeshManager.h"
#include "3D_MeshManager.h"
#include "BoundaryConditionManager.h"
#include "BoundaryAssembler.h"
#include "2D_PostProcess.h"
#include "3D_PostProcess.h"

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


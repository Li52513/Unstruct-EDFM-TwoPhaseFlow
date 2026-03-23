/**
 * @file Test_Day6_TransientSolver.cpp
 * @brief Day6 transient scenarios (thin wrappers over FIM_TransientCaseKit)
 *
 * Integrated features exercised by Run_Day6_Transient_2D_SP_InjProd():
 *   - CoolProp EOS (Step 1): Span-Wagner CO2 + IAPWS-95 Water via USE_COOLPROP_EOS
 *   - Non-orthogonal correction (Step 2): deferred T-vector gradient correction
 *   - Independent well DOF (Step 3): WellDOFManager adds P_wbh block per well
 *   - VTU/PVD post-processing (Step 4): ParaView XML Unstructured Grid export
 */

#include "Test_Day6_TransientSolver.h"
#include "FIM_TransientCaseKit.hpp"
#include "2D_PostProcess.h"

namespace Test_Day6 {

    void Run_Day6_Transient_2D_SP_InjProd() {
        std::cout << "\n=== Run_Day6_Transient_2D_SP_InjProd ===\n";

        // -- Feature flags summary --
#ifdef USE_COOLPROP_EOS
        std::cout << "  [Step 1] EOS backend: CoolProp 6.x (Span-Wagner CO2 + IAPWS-95 Water)\n";
#else
        std::cout << "  [Step 1] EOS backend: Table-based (Catmull-Rom + 4-tier FD)\n";
#endif
        std::cout << "  [Step 2] Non-orthogonal correction: ENABLED\n";
        std::cout << "  [Step 3] Well independent DOF: ENABLED (WellDOFManager)\n";
        std::cout << "  [Step 4] Post-processing: VTU/PVD (ParaView XML)\n\n";

        // =====================================================================
        // 1. Mesh & Topology Setup
        // =====================================================================
        MeshManager mgr(100.0, 100.0, 0.0, 10, 10, 0, true, false);
        mgr.BuildSolidMatrixGrid_2D();
        mgr.addFracture(Vector(10.0, 10.0, 0.0), Vector(90.0, 90.0, 0.0));
        mgr.addFracture(Vector(10.0, 90.0, 0.0), Vector(90.0, 10.0, 0.0));
        mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
        mgr.BuildGlobalSystemIndexing();
        mgr.BuildFracturetoFractureTopology();
        mgr.setNumDOFs(2);

        // =====================================================================
        // 2. Field Manager Initialization
        // =====================================================================
        FieldManager_2D fm;
        FIM_CaseKit::InitFieldManager(mgr, fm);

        // =====================================================================
        // 3. Rock/Fracture Properties Setup
        // =====================================================================
        auto preset = FIM_CaseKit::MakeDefaultPropertyPreset2D();
		preset.enable_rock_region = false;  // 2D case with uniform rock properties for simplicity
        preset.rock_bg.phi_r = 0.12;        // 孔隙度
        preset.rock_bg.kxx = 1.0e-13; 
        preset.rock_bg.kyy = 1.0e-13;
		preset.rock_bg.kzz = 1.0e-13;       // 渗透率，单位 m^2，约等于 0.5 mD
		preset.rock_bg.compressibility = 5.0e-9;  // Pa^-1, 约等于 5e-6 MPa^-1，保持适度压缩性以观察压力传播和孔隙体积变化的影响
		preset.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::Water; // 使用 CO2 作为单相流体，具有较高的可压缩性和温度敏感性，能够更明显地展示压力和温度变化的耦合效应

        // =====================================================================
        // 4. Initial Conditions (IC)
        // =====================================================================
        FIM_Engine::InitialConditions ic;
        ic.P_init = 30.0e6;    // 30 MPa  (深储层，>Pc 裕量×4)
        ic.T_init = 380.0;     // 380 K   (≈107°C，深3km地热，>Tc 裕量+76K)
        ic.Sw_init = 1.0;

        // =====================================================================
        // 5. External BC (Dirichlet / Neumann / Robin)
        // =====================================================================
		BoundarySetting::BoundaryConditionManager bcP;  // 压力 BC 管理器
		BoundarySetting::BoundaryConditionManager bcT;  // 温度 BC 管理器
        bcP.Clear();
        bcT.Clear();

        // Pressure: closed (no-flow) on all sides
        bcP.SetNeumannBC(MeshTags::LEFT, 0.0);          // 压力左边界
		bcP.SetNeumannBC(MeshTags::RIGHT, 0.0);         // 压力右边界
		bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);        // 压力底边界
		bcP.SetNeumannBC(MeshTags::TOP, 0.0);           // 压力上边界

        // Temperature: adiabatic on all sides
        bcT.SetNeumannBC(MeshTags::LEFT, 0.0);          // 左边界    
		bcT.SetNeumannBC(MeshTags::RIGHT, 0.0);         // 右边界
		bcT.SetNeumannBC(MeshTags::BOTTOM, 0.0);        // 底边界
		bcT.SetNeumannBC(MeshTags::TOP, 0.0);           // 上边界

		auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, nullptr); // 创建模型

        // =====================================================================
        // 6. Wells Setup
        // =====================================================================
        const auto& cells = mgr.mesh().getCells();
        auto findNearestCell = [&](double x, double y) {
            int best = 0;
            double bestD2 = std::numeric_limits<double>::max();
            for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
                const double dx = cells[i].center.m_x - x;
                const double dy = cells[i].center.m_y - y;
                const double d2 = dx * dx + dy * dy;
                if (d2 < bestD2) {
                    bestD2 = d2;
                    best = i;
                }
            }
            return best;
            };

        const int injCell = findNearestCell(15.0, 15.0);   // 移至单元中心(i=1,j=1)，避开y=10单元边界，且位于对角裂缝路径上
        const int prodCell = findNearestCell(85.0, 95.0);  // 明确落于顶行单元中心

        WellScheduleStep w1, w2;
        w1.well_name = "INJ_CO2_2D";
        w1.domain = WellTargetDomain::Matrix;
        w1.control_mode = WellControlMode::BHP;
        w1.target_value = 31.0e6;  // 31 MPa
        w1.component_mode = WellComponentMode::Water;
        w1.completion_id = injCell;
        w1.frac_w = 1.0;
        w1.frac_g = 0.0;
        w1.injection_temperature = 340.0;   // 340 K，比储层低 40K，观察冷注降温效果
        w1.injection_is_co2 = false;


        w2.well_name = "PROD_CO2_2D";
        w2.domain = WellTargetDomain::Matrix;
        w2.control_mode = WellControlMode::BHP;
        w2.target_value = 29.0e6;  // 29 MPa
        w2.component_mode = WellComponentMode::Water;
        w2.completion_id = prodCell;
        w2.frac_w = 1.0;
        w2.frac_g = 0.0;
        

        // =====================================================================
        // 7. Solver Parameters Setup
        // =====================================================================
        const double horizon_years = 1.0;
        // Startup: 30天，让温度斜坡在startup内完整走完，避免切换长阶段时的28.8K硬跳变
        const double startup_days = 30.0;
        const double startup_vtk_interval_days = 0.5;   // 30天startup内每5天出一帧
        const double long_vtk_interval_days = 5.0;

        auto params = FIM_CaseKit::BuildLongRunTemplate(
            FIM_CaseKit::LongRunScenario::S2D_SP,
            horizon_years,
            startup_days);

        // Explicit local controls for long-run setup in this case.
        params.target_end_time_s = horizon_years * 365.0 * 86400.0;
        params.startup_end_time_s = startup_days * 86400.0;
        params.startup_vtk_output_interval_s = startup_vtk_interval_days * 86400.0;
        params.long_vtk_output_interval_s = long_vtk_interval_days * 86400.0;
        params.gravity_vector = Vector(0.0, 0.0, 0.0);

        // Keep startup short and allow dt growth when iter~14.
		params.startup_profile.max_newton_iter = 18; // Startup阶段允许更深入的非线性迭代，以更快地适应初始条件和边界条件，促进时间步长增长
		params.startup_profile.ptc_lambda_init = 1.0; // 初始PTC lambda较大，提供更强的正则化以稳定初始迭代，促进时间步长增长
        params.startup_profile.ptc_lambda_decay = 0.70;
        params.startup_profile.ptc_lambda_min = 0.02;
        params.startup_profile.dt_relres_iter_grow_hi = 12;
        params.startup_profile.dt_relres_iter_neutral_hi = 18;
        params.startup_profile.dt_relres_iter_soft_shrink_hi = 24;
        params.startup_profile.dt_relres_grow_factor = 1.10;
        params.startup_profile.dt_relres_neutral_factor = 1.02;
        params.startup_profile.dt_relres_soft_shrink_factor = 0.995;
        params.startup_profile.dt_relres_hard_shrink_factor = 0.92;
        // 100步完成斜坡（~4天 @ dt_max=3600s），远早于30天startup结束，消除切换时硬跳变
        params.startup_profile.control_ramp_steps = 10;
		params.startup_profile.rel_res_tol = 1e-3; // Startup阶段使用较严格的相对残差容忍度，以确保初始非线性收敛质量，促进时间步长增长


        params.long_profile.max_newton_iter = 36; // Allow deeper nonlinear progress in late stiff periods.
        params.long_profile.rel_update_tol = 1e-3;
        params.long_profile.ptc_lambda_init = 0.5;
        params.long_profile.ptc_lambda_decay = 0.7;
        params.long_profile.ptc_lambda_min = 5.0e-2; // 加强long stage对角正则化，防止能量方程在高刚性时失稳
        params.long_profile.dt_relres_iter_grow_hi = 8;
        params.long_profile.dt_relres_iter_neutral_hi = 14;
        params.long_profile.dt_relres_iter_soft_shrink_hi = 18;
        params.long_profile.dt_relres_grow_factor = 1.05;
        params.long_profile.dt_relres_neutral_factor = 1.02;
        params.long_profile.dt_relres_soft_shrink_factor = 0.90;
        params.long_profile.dt_relres_hard_shrink_factor = 0.90;
        params.rollback_shrink_factor = 0.85;
        params.clamp_state_to_eos_bounds = true; // Keep single-phase CO2 state inside EOS table bounds.
        params.long_profile.rel_res_tol = 1e-2;
        
        // 降低Newton溫度更新上限：预设2K/iter过大，在偽临界區域（∂h/∂T極大）會過衝
        params.max_dT = 5;   // 从 0.5 放宽到 5K/iter，远离拟临界不再需要极小步

        // Issue#12: 切换至 CPR-AMG 求解器，验证全流程稳定性
        params.lin_solver = FIM_Engine::LinearSolverType::AMGCL_CPR;

        // Step 2: Enable deferred non-orthogonal correction (T-vector gradient term).
        // On this 10x10 orthogonal grid |vectorT|≈0, so the correction has near-zero effect
        // but exercises the full code path. Use skewed meshes for accuracy improvement.
        params.enable_non_orthogonal_correction = true;

        // =====================================================================
        // 8. Execution
        // =====================================================================
        // Step 3: Wells w1, w2 are passed to RunGenericFIMTransient, which internally
        // invokes WellDOFManager to add independent P_wbh DOF blocks per well.
        // The sparse matrix is extended to (nMatrix + nFrac + nWells) blocks;
        // well equations (Peaceman model) are assembled each Newton iteration.
        FIM_Engine::RunGenericFIMTransient<2>(
            "day6_transient_2d_sp_injprod", mgr, fm, ic, { w1, w2 }, params, FIM_Engine::SolverRoute::FIM, modules);

        // =====================================================================
        // 9. VTU/PVD Post-Processing (Step 4)
        // =====================================================================
        // Export final state as ParaView XML Unstructured Grid (.vtu).
        // RunGenericFIMTransient already synced field data at the last VTK snapshot,
        // so fm contains the final-time solution.
        {
            PostProcess_2D pp(mgr, fm);
            const std::string baseDir = "Test/Transient/Day6/day6_transient_2d_sp_injprod/";
            const std::string vtu_file = baseDir + "final_state.vtu";
            pp.ExportVTU(vtu_file, params.target_end_time_s);
            std::cout << "    [VTU] Exported: " << vtu_file << "\n";

            // PVD time-series collection (single-snapshot demo).
            // In production, RunGenericFIMTransient would collect all intermediate
            // VTU snapshots and produce a complete .pvd for ParaView animation.
            const std::string pvd_file = baseDir + "timeseries.pvd";
            PostProcess_2D::ExportPVD(pvd_file,
                                      { vtu_file },
                                      { params.target_end_time_s });
            std::cout << "    [PVD] Exported: " << pvd_file << "\n";
        }
    }

void Run_Day6_Transient_2D_TP_InjProd() {
    // =====================================================================
    // 1. Mesh & Topology Setup
    // =====================================================================
    MeshManager mgr(10.0, 10.0, 0.0, 10, 10, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    mgr.addFracture(Vector(1, 1, 0), Vector(9, 9, 0));
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(3);

    // =====================================================================
    // 2. Field Manager Initialization
    // =====================================================================
    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    // =====================================================================
    // 3. Rock/Fracture Properties Setup
    // =====================================================================
    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset2D();

    // =====================================================================
    // 4. Initial Conditions (IC)
    // =====================================================================
    FIM_Engine::InitialConditions ic;
    ic.P_init = 25.0e6;
    ic.T_init = 360.0;
    ic.Sw_init = 0.25;

    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    BoundarySetting::BoundaryConditionManager bcP;
    BoundarySetting::BoundaryConditionManager bcT;
    bcP.Clear();
    bcT.Clear();

    // Pressure: closed (no-flow) on all sides
    bcP.SetNeumannBC(MeshTags::LEFT, 0.0);
    bcP.SetNeumannBC(MeshTags::RIGHT, 0.0);
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcP.SetNeumannBC(MeshTags::TOP, 0.0);

    // Temperature: adiabatic on all sides
    bcT.SetNeumannBC(MeshTags::LEFT, 0.0);
    bcT.SetNeumannBC(MeshTags::RIGHT, 0.0);
    bcT.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcT.SetNeumannBC(MeshTags::TOP, 0.0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, nullptr);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    const int nMat = FIM_CaseKit::MatrixBlockCount(mgr);
    const int injCell = nMat / 5;
    const int prodCell = (4 * nMat) / 5;

    WellScheduleStep w1, w2;
    w1.well_name = "INJ_2D_TP";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::BHP;
    w1.target_value = 25.15e6;
    w1.component_mode = WellComponentMode::Water;
    w1.injection_temperature = 330.0;
    w1.completion_id = injCell;
    w1.frac_w = 1.0;
    w1.frac_g = 0.0;

    w2.well_name = "PROD_2D_TP";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 24.90e6;
    w2.component_mode = WellComponentMode::Total;
    w2.completion_id = prodCell;
    // Keep frac_w/frac_g unset for Total producer:
    // BoundaryAssembler now auto-splits by local mobility.

    // =====================================================================
    // 7. Solver Parameters Setup
    // =====================================================================
    auto params = FIM_CaseKit::BuildSolverParams(true);

    // =====================================================================
    // 8. Execution
    // =====================================================================
    FIM_Engine::RunGenericFIMTransient<3>(
        "day6_transient_2d_tp_injprod", mgr, fm, ic, { w1, w2 }, params, FIM_Engine::SolverRoute::FIM, modules);
}

void Run_Day6_Transient_2D_TP_Multiwell() {
    // =====================================================================
    // 1. Mesh & Topology Setup
    // =====================================================================
    MeshManager mgr(10.0, 10.0, 0.0, 10, 10, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    mgr.addFracture(Vector(1, 1, 0), Vector(9, 9, 0));
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(3);

    // =====================================================================
    // 2. Field Manager Initialization
    // =====================================================================
    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    // =====================================================================
    // 3. Rock/Fracture Properties Setup
    // =====================================================================
    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset2D();

    // =====================================================================
    // 4. Initial Conditions (IC)
    // =====================================================================
    FIM_Engine::InitialConditions ic;
    ic.P_init = 25.0e6;
    ic.T_init = 360.0;
    ic.Sw_init = 0.25;

    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    // Closed boundary for Day6 robustness: no external BC injection into solver.
    auto modules = FIM_CaseKit::BuildModules2D(preset, nullptr, nullptr, nullptr);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    const int nMat = FIM_CaseKit::MatrixBlockCount(mgr);
    const int injCell = nMat / 5;
    const int prodCell = (4 * nMat) / 5;
    const int midCell = nMat / 2;

    WellScheduleStep w1, w2, w3;
    w1.well_name = "INJ_W";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::Rate;
    w1.target_value = -0.04;
    w1.component_mode = WellComponentMode::Water;
    w1.injection_temperature = 330.0;
    w1.completion_id = injCell;
    w1.frac_w = 1.0;
    w1.frac_g = 0.0;

    w2.well_name = "PROD_1";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 24.92e6;
    w2.component_mode = WellComponentMode::Total;
    w2.completion_id = prodCell;

    w3.well_name = "PROD_2";
    w3.domain = WellTargetDomain::Matrix;
    w3.control_mode = WellControlMode::BHP;
    w3.target_value = 24.92e6;
    w3.component_mode = WellComponentMode::Total;
    w3.completion_id = midCell;

    // =====================================================================
    // 7. Solver Parameters Setup
    // =====================================================================
    auto params = FIM_CaseKit::BuildSolverParams(true);

    // =====================================================================
    // 8. Execution
    // =====================================================================
    FIM_Engine::RunGenericFIMTransient<3>(
        "day6_transient_2d_tp_multiwell", mgr, fm, ic, { w1, w2, w3 }, params, FIM_Engine::SolverRoute::FIM, modules);
}

void Run_Day6_Transient_3D_SP_InjProd() {
    // =====================================================================
    // 1. Mesh & Topology Setup
    // =====================================================================
    MeshManager_3D mgr(10.0, 10.0, 10.0, 20, 20, 20, true, false);
    mgr.BuildSolidMatrixGrid_3D();
    std::vector<Vector> pts = { Vector(1, 1, 1), Vector(9, 1, 1), Vector(9, 9, 1), Vector(1, 9, 1) };
    mgr.addFracturetoFractureNetwork(Fracture_2D(0, pts));
    mgr.meshAllFracturesinNetwork(2, 2);

    mgr.setupGlobalIndices();
    mgr.DetectFractureFractureIntersectionsInNetwork();
    mgr.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);
    mgr.fracture_network().rebuildEdgeProperties();
    mgr.removeDuplicateInteractions();
    mgr.resolveCoplanarInteractions();
    mgr.buildTopologyMaps();
    mgr.setNumDOFs(2);

    // =====================================================================
    // 2. Field Manager Initialization
    // =====================================================================
    FieldManager_3D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    // =====================================================================
    // 3. Rock/Fracture Properties Setup
    // =====================================================================
    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset3D();
    preset.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::CO2;

    // =====================================================================
    // 4. Initial Conditions (IC)
    // =====================================================================
    FIM_Engine::InitialConditions ic;
    ic.P_init = 25.0e6;
    ic.T_init = 360.0;
    ic.Sw_init = 1.0;

    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    // Closed boundary for Day6 robustness: no external BC injection into solver.
    auto modules = FIM_CaseKit::BuildModules3D(preset, nullptr, nullptr, nullptr);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    const int nMat = FIM_CaseKit::MatrixBlockCount(mgr);
    const int injCell = nMat / 5;
    const int prodCell = (4 * nMat) / 5;

    WellScheduleStep w1, w2;
    w1.well_name = "INJ_3D";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::Rate;
    w1.target_value = -0.08;
    w1.component_mode = WellComponentMode::Gas;
    w1.injection_temperature = 330.0;
    w1.injection_is_co2 = true;
    w1.completion_id = injCell;
    w1.frac_w = 0.0;
    w1.frac_g = 1.0;
    w1.well_axis = WellAxis::Z;

    w2.well_name = "PROD_3D";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 24.90e6;
    w2.component_mode = WellComponentMode::Gas;
    w2.completion_id = prodCell;
    w2.frac_w = 0.0;
    w2.frac_g = 1.0;
    w2.well_axis = WellAxis::Z;

    // =====================================================================
    // 7. Solver Parameters Setup
    // =====================================================================
    auto params = FIM_CaseKit::BuildSolverParams(false);

    // =====================================================================
    // 8. Execution
    // =====================================================================
    FIM_Engine::RunGenericFIMTransient<2>(
        "day6_transient_3d_sp_injprod", mgr, fm, ic, { w1, w2 }, params, FIM_Engine::SolverRoute::FIM, modules);
}

void Run_Day6_Transient_3D_TP_InjProd() {
    // =====================================================================
    // 1. Mesh & Topology Setup
    // =====================================================================
    MeshManager_3D mgr(10.0, 10.0, 10.0, 10, 10, 10, true, false);
    mgr.BuildSolidMatrixGrid_3D();
    std::vector<Vector> pts = { Vector(1, 1, 1), Vector(9, 1, 1), Vector(9, 9, 1), Vector(1, 9, 1) };
    mgr.addFracturetoFractureNetwork(Fracture_2D(0, pts));
    mgr.meshAllFracturesinNetwork(2, 2);

    mgr.setupGlobalIndices();
    mgr.DetectFractureFractureIntersectionsInNetwork();
    mgr.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);
    mgr.fracture_network().rebuildEdgeProperties();
    mgr.removeDuplicateInteractions();
    mgr.resolveCoplanarInteractions();
    mgr.buildTopologyMaps();
    mgr.setNumDOFs(3);

    // =====================================================================
    // 2. Field Manager Initialization
    // =====================================================================
    FieldManager_3D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    // =====================================================================
    // 3. Rock/Fracture Properties Setup
    // =====================================================================
    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset3D();
    preset.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::CO2;

    // =====================================================================
    // 4. Initial Conditions (IC)
    // =====================================================================
    FIM_Engine::InitialConditions ic;
    ic.P_init = 25.0e6;
    ic.T_init = 360.0;
    ic.Sw_init = 0.25;


    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    // Closed boundary for Day6 robustness: no external BC injection into solver.
    auto modules = FIM_CaseKit::BuildModules3D(preset, nullptr, nullptr, nullptr);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    const int nMat = FIM_CaseKit::MatrixBlockCount(mgr);
    const int injCell = nMat / 5;
    const int prodCell = (4 * nMat) / 5;

    WellScheduleStep w1, w2;
    w1.well_name = "INJ_3D_TP";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::BHP;
    w1.target_value = 25.15e6;
    w1.component_mode = WellComponentMode::Water;
    w1.injection_temperature = 330.0;
    w1.completion_id = injCell;
    w1.frac_w = 1.0;
    w1.frac_g = 0.0;
    w1.well_axis = WellAxis::Z;

    w2.well_name = "PROD_3D_TP";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 24.90e6;
    w2.component_mode = WellComponentMode::Total;
    w2.completion_id = prodCell;
    w2.well_axis = WellAxis::Z;


    // =====================================================================
    // 7. Solver Parameters Setup
    // =====================================================================
    auto params = FIM_CaseKit::BuildSolverParams(true);

    // =====================================================================
    // 8. Execution
    // =====================================================================
    FIM_Engine::RunGenericFIMTransient<3>(
        "day6_transient_3d_tp_injprod", mgr, fm, ic, { w1, w2 }, params, FIM_Engine::SolverRoute::FIM, modules);
}

void Run_Day6_Transient_3D_TP_Multiwell() {
    // =====================================================================
    // 1. Mesh & Topology Setup
    // =====================================================================
    MeshManager_3D mgr(10.0, 10.0, 10.0, 10, 10, 10, true, false);
    mgr.BuildSolidMatrixGrid_3D();
    std::vector<Vector> pts = { Vector(1, 1, 1), Vector(9, 1, 1), Vector(9, 9, 1), Vector(1, 9, 1) };
    mgr.addFracturetoFractureNetwork(Fracture_2D(0, pts));
    mgr.meshAllFracturesinNetwork(2, 2);

    mgr.setupGlobalIndices();
    mgr.DetectFractureFractureIntersectionsInNetwork();
    mgr.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);
    mgr.fracture_network().rebuildEdgeProperties();
    mgr.removeDuplicateInteractions();
    mgr.resolveCoplanarInteractions();
    mgr.buildTopologyMaps();
    mgr.setNumDOFs(3);

    // =====================================================================
    // 2. Field Manager Initialization
    // =====================================================================
    FieldManager_3D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    // =====================================================================
    // 3. Rock/Fracture Properties Setup
    // =====================================================================
    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset3D();
    preset.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::CO2;

    // =====================================================================
    // 4. Initial Conditions (IC)
    // =====================================================================
    FIM_Engine::InitialConditions ic;
    ic.P_init = 25.0e6;
    ic.T_init = 360.0;
    ic.Sw_init = 0.25;

    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    // Closed boundary for Day6 robustness: no external BC injection into solver.
    auto modules = FIM_CaseKit::BuildModules3D(preset, nullptr, nullptr, nullptr);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    const int nMat = FIM_CaseKit::MatrixBlockCount(mgr);
    const int injCell = nMat / 5;
    const int prodCell = (4 * nMat) / 5;
    const int midCell = nMat / 2;

    WellScheduleStep w1, w2, w3;
    w1.well_name = "INJ_W_3D";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::Rate;
    w1.target_value = -0.08;
    w1.component_mode = WellComponentMode::Water;
    w1.injection_temperature = 330.0;
    w1.completion_id = injCell;
    w1.frac_w = 1.0;
    w1.frac_g = 0.0;
    w1.well_axis = WellAxis::Z;

    w2.well_name = "PROD_1_3D";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 24.92e6;
    w2.component_mode = WellComponentMode::Total;
    w2.completion_id = prodCell;
    w2.well_axis = WellAxis::Z;

    w3.well_name = "PROD_2_3D";
    w3.domain = WellTargetDomain::Matrix;
    w3.control_mode = WellControlMode::BHP;
    w3.target_value = 24.92e6;
    w3.component_mode = WellComponentMode::Total;
    w3.completion_id = midCell;
    w3.well_axis = WellAxis::Z;

    // =====================================================================
    // 7. Solver Parameters Setup
    // =====================================================================
    auto params = FIM_CaseKit::BuildSolverParams(true);

    // =====================================================================
    // 8. Execution
    // =====================================================================
    FIM_Engine::RunGenericFIMTransient<3>(
        "day6_transient_3d_tp_multiwell", mgr, fm, ic, { w1, w2, w3 }, params, FIM_Engine::SolverRoute::FIM, modules);
}

// ---------- 新增 矩阵审计 测试代码 ----------
void Run_Day6_MatrixAudit_2D_EDFM() {
    MeshManager mgr(100.0, 100.0, 0.0, 10, 10, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    mgr.addFracture(Vector(10.0, 10.0, 0.0), Vector(90.0, 90.0, 0.0));
    mgr.addFracture(Vector(10.0, 90.0, 0.0), Vector(90.0, 10.0, 0.0));
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset2D();
    FIM_Engine::InitialConditions ic;
    auto modules = FIM_CaseKit::BuildModules2D(preset, nullptr, nullptr, nullptr);

    auto params = FIM_CaseKit::BuildSolverParams(false, 1, 1.0);
    params.enable_matrix_audit = true;
    params.matrix_audit_strict = true;
    params.matrix_audit_step = 1;
    params.matrix_audit_iter = 1;
    params.max_steps = 1;
    params.max_newton_iter = 1;
    params.abs_res_tol = 1e30; // 仅过一次残差及装配，不管牛顿收敛
    params.gravity_vector = Vector(0.0, 0.0, 0.0);

    FIM_Engine::RunGenericFIMTransient<2>(
        "day6_matrix_audit_2d_edfm", mgr, fm, ic, {}, params, FIM_Engine::SolverRoute::FIM, modules);
}

void Run_Day6_MatrixAudit_3D_EDFM() {
    MeshManager_3D mgr(10.0, 10.0, 10.0, 10, 10, 10, true, false);
    mgr.BuildSolidMatrixGrid_3D();
    std::vector<Vector> pts1 = { Vector(1, 1, 1), Vector(9, 1, 1), Vector(9, 9, 9), Vector(1, 9, 9) };
    std::vector<Vector> pts2 = { Vector(1, 9, 1), Vector(9, 9, 1), Vector(9, 1, 9), Vector(1, 1, 9) };
    mgr.addFracturetoFractureNetwork(Fracture_2D(0, pts1));
    mgr.addFracturetoFractureNetwork(Fracture_2D(1, pts2));
    mgr.meshAllFracturesinNetwork(2, 2);

    mgr.setupGlobalIndices();
    mgr.DetectFractureFractureIntersectionsInNetwork();
    mgr.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);
    mgr.fracture_network().rebuildEdgeProperties();
    mgr.removeDuplicateInteractions();
    mgr.resolveCoplanarInteractions();
    mgr.buildTopologyMaps();
    mgr.setNumDOFs(2);

    FieldManager_3D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset3D();
    preset.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::CO2;
    FIM_Engine::InitialConditions ic;
    auto modules = FIM_CaseKit::BuildModules3D(preset, nullptr, nullptr, nullptr);

    auto params = FIM_CaseKit::BuildSolverParams(false, 1, 1.0);
    params.enable_matrix_audit = true;
    params.matrix_audit_strict = true;
    params.matrix_audit_step = 1;
    params.matrix_audit_iter = 1;
    params.max_steps = 1;
    params.max_newton_iter = 1;
    params.abs_res_tol = 1e30;

    FIM_Engine::RunGenericFIMTransient<2>(
        "day6_matrix_audit_3d_edfm", mgr, fm, ic, {}, params, FIM_Engine::SolverRoute::FIM, modules);
}

} // namespace Test_Day6



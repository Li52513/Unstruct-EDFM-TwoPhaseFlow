/**
 * @file Test_Day6_TransientSolver.cpp
 * @brief Day6 transient scenarios (thin wrappers over FIM_TransientCaseKit)
 */

#include "Test_Day6_TransientSolver.h"
#include "FIM_TransientCaseKit.hpp"

namespace Test_Day6 {

    void Run_Day6_Transient_2D_SP_InjProd() {
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
        preset.enable_rock_region = false;
        preset.rock_bg.phi_r = 0.12;
        preset.rock_bg.kxx = 5.0e-14;
        preset.rock_bg.kyy = 5.0e-14;
        preset.rock_bg.kzz = 5.0e-14;
        preset.rock_bg.compressibility = 5.0e-9;
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

        const int injCell = findNearestCell(25.0, 25.0);
        const int prodCell = findNearestCell(75.0, 75.0);

        WellScheduleStep w1, w2;
        w1.well_name = "INJ_CO2_2D";
        w1.domain = WellTargetDomain::Matrix;
        w1.control_mode = WellControlMode::Rate;
        w1.target_value = -0.03;   // kg/s, softened startup stiffness for Day6 baseline
        w1.component_mode = WellComponentMode::Gas;
        w1.completion_id = injCell;
        w1.frac_w = 0.0;
        w1.frac_g = 1.0;
        w1.injection_temperature = 320.0;
        w1.injection_is_co2 = true;

        w2.well_name = "PROD_CO2_2D";
        w2.domain = WellTargetDomain::Matrix;
        w2.control_mode = WellControlMode::BHP;
        w2.target_value = 24.95e6; // Pa, keep drawdown moderate in closed boundary domain
        w2.component_mode = WellComponentMode::Gas;
        w2.completion_id = prodCell;
        w2.frac_w = 0.0;
        w2.frac_g = 1.0;

        // =====================================================================
        // 7. Solver Parameters Setup
        // =====================================================================
        const double horizon_years = 1.0;
        const double startup_seconds = 0.2;
        const double startup_days = startup_seconds / 86400.0;
        const double startup_vtk_interval_days = 0.5;
        const double long_vtk_interval_days = 30.0;

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
        params.startup_profile.max_newton_iter = 18;
        params.startup_profile.ptc_lambda_init = 1.0;
        params.startup_profile.ptc_lambda_decay = 0.70;
        params.startup_profile.ptc_lambda_min = 0.02;
        params.startup_profile.dt_relres_iter_grow_hi = 12;
        params.startup_profile.dt_relres_iter_neutral_hi = 18;
        params.startup_profile.dt_relres_iter_soft_shrink_hi = 24;
        params.startup_profile.dt_relres_grow_factor = 1.10;
        params.startup_profile.dt_relres_neutral_factor = 1.02;
        params.startup_profile.dt_relres_soft_shrink_factor = 0.995;
        params.startup_profile.dt_relres_hard_shrink_factor = 0.92;

		params.long_profile.max_newton_iter = 24; // Allow moderate Newton iterations for long-term stability, but keep it in check to prevent excessive runtime in this test scenario.
        params.long_profile.rel_update_tol = 1e-3;
        params.long_profile.ptc_lambda_init = 0.5;
        params.long_profile.ptc_lambda_decay = 0.7;
        params.long_profile.ptc_lambda_min = 1.0e-3;
        params.long_profile.dt_relres_iter_grow_hi = 8;
        params.long_profile.dt_relres_iter_neutral_hi = 14;
        params.long_profile.dt_relres_iter_soft_shrink_hi = 18;
        params.long_profile.dt_relres_grow_factor = 1.05;
        params.long_profile.dt_relres_neutral_factor = 1.02;
        params.long_profile.dt_relres_soft_shrink_factor = 0.90;
        params.long_profile.dt_relres_hard_shrink_factor = 0.90;
        params.rollback_shrink_factor = 0.7;

        // Issue#12: 切换至 CPR-AMG 求解器，验证全流程稳定性
        params.lin_solver = FIM_Engine::LinearSolverType::AMGCL_CPR;

        // =====================================================================
        // 8. Execution
        // =====================================================================
        FIM_Engine::RunGenericFIMTransient<2>(
            "day6_transient_2d_sp_injprod", mgr, fm, ic, { w1, w2 }, params, FIM_Engine::SolverRoute::FIM, modules);
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



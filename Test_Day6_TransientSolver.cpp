/**
 * @file Test_Day6_TransientSolver.cpp
 * @brief Day6 transient scenarios (thin wrappers over FIM_TransientCaseKit)
 */

#include "Test_Day6_TransientSolver.h"
#include "FIM_TransientCaseKit.hpp"

namespace Test_Day6 {

void Run_Day6_Transient_2D_SP_InjProd() {
    using FIM_CaseKit::BoundaryConditionManager;

    // =====================================================================
    // 1. Mesh & Topology Setup
    // =====================================================================
    MeshManager mgr(10.0, 10.0, 0.0, 10, 10, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    mgr.addFracture(Vector(1, 1, 0), Vector(9, 9, 0));
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

    // =====================================================================
    // 4. Initial Conditions (IC)
    // =====================================================================
    FIM_Engine::InitialConditions ic;
    ic.P_init = 6.0e6;
    ic.T_init = 300.0;
    ic.Sw_init = 1.0;

    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    BoundaryConditionManager bcP;
    BoundaryConditionManager bcT;
    FIM_CaseKit::ConfigurePressureBC2D(bcP, 7.0e6, 5.5e6, 6.0e6);
    FIM_CaseKit::ConfigureTemperatureBC2D(bcT, 300.0);
    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, nullptr);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    WellScheduleStep w1, w2;
    w1.well_name = "INJ_2D";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::BHP;
    w1.target_value = 7.0e6;
    w1.component_mode = WellComponentMode::Water;
    w1.completion_id = 0;
    w1.frac_w = 1.0;
    w1.frac_g = 0.0;

    w2.well_name = "PROD_2D";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 5.5e6;
    w2.component_mode = WellComponentMode::Total;
    w2.completion_id = FIM_CaseKit::MatrixBlockCount(mgr) - 1;
    w2.frac_w = 1.0;
    w2.frac_g = 0.0;

    // =====================================================================
    // 7. Solver Parameters Setup
    // =====================================================================
    auto params = FIM_CaseKit::BuildSolverParams(false);

    // =====================================================================
    // 8. Execution
    // =====================================================================
    FIM_Engine::RunGenericFIMTransient<2>(
        "day6_transient_2d_sp_injprod", mgr, fm, ic, { w1, w2 }, params, FIM_Engine::SolverRoute::FIM, modules);
}

void Run_Day6_Transient_2D_TP_InjProd() {
    using FIM_CaseKit::BoundaryConditionManager;

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
    ic.P_init = 6.0e6;
    ic.T_init = 300.0;
    ic.Sw_init = 0.2;

    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    BoundaryConditionManager bcP;
    BoundaryConditionManager bcT;
    BoundaryConditionManager bcSw;
    FIM_CaseKit::ConfigurePressureBC2D(bcP, 7.0e6, 5.5e6, 6.0e6);
    FIM_CaseKit::ConfigureTemperatureBC2D(bcT, 300.0);
    FIM_CaseKit::ConfigureSaturationBC2D(bcSw, 0.8, 0.2);
    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, &bcSw);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    WellScheduleStep w1, w2;
    w1.well_name = "INJ_2D_TP";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::BHP;
    w1.target_value = 7.0e6;
    w1.component_mode = WellComponentMode::Water;
    w1.completion_id = 0;
    w1.frac_w = 1.0;
    w1.frac_g = 0.0;

    w2.well_name = "PROD_2D_TP";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 5.5e6;
    w2.component_mode = WellComponentMode::Total;
    w2.completion_id = FIM_CaseKit::MatrixBlockCount(mgr) - 1;
    w2.frac_w = 0.5;
    w2.frac_g = 0.5;

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
    using FIM_CaseKit::BoundaryConditionManager;

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
    ic.P_init = 6.0e6;
    ic.T_init = 300.0;
    ic.Sw_init = 0.2;

    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    BoundaryConditionManager bcP;
    BoundaryConditionManager bcT;
    BoundaryConditionManager bcSw;
    FIM_CaseKit::ConfigurePressureBC2D(bcP, 7.0e6, 5.5e6, 6.0e6);
    FIM_CaseKit::ConfigureTemperatureBC2D(bcT, 300.0);
    FIM_CaseKit::ConfigureSaturationBC2D(bcSw, 0.8, 0.2);
    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, &bcSw);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    WellScheduleStep w1, w2, w3;
    w1.well_name = "INJ_W";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::Rate;
    w1.target_value = -0.5;
    w1.component_mode = WellComponentMode::Water;
    w1.completion_id = 0;
    w1.frac_w = 1.0;
    w1.frac_g = 0.0;

    w2.well_name = "PROD_1";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 5.5e6;
    w2.component_mode = WellComponentMode::Total;
    w2.completion_id = FIM_CaseKit::MatrixBlockCount(mgr) - 1;
    w2.frac_w = 0.5;
    w2.frac_g = 0.5;

    w3.well_name = "PROD_2";
    w3.domain = WellTargetDomain::Matrix;
    w3.control_mode = WellControlMode::BHP;
    w3.target_value = 5.5e6;
    w3.component_mode = WellComponentMode::Total;
    w3.completion_id = FIM_CaseKit::MatrixBlockCount(mgr) / 2;
    w3.frac_w = 0.5;
    w3.frac_g = 0.5;

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
    using FIM_CaseKit::BoundaryConditionManager;

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

    // =====================================================================
    // 4. Initial Conditions (IC)
    // =====================================================================
    FIM_Engine::InitialConditions ic;
    ic.P_init = 6.0e6;
    ic.T_init = 300.0;
    ic.Sw_init = 1.0;

    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    BoundaryConditionManager bcP;
    BoundaryConditionManager bcT;
    FIM_CaseKit::ConfigurePressureBC3D(bcP, 7.0e6, 5.5e6, 6.0e6);
    FIM_CaseKit::ConfigureTemperatureBC3D(bcT, 300.0);
    auto modules = FIM_CaseKit::BuildModules3D(preset, &bcP, &bcT, nullptr);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    WellScheduleStep w1, w2;
    w1.well_name = "INJ_3D";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::BHP;
    w1.target_value = 7.0e6;
    w1.component_mode = WellComponentMode::Water;
    w1.completion_id = 0;
    w1.frac_w = 1.0;
    w1.frac_g = 0.0;
    w1.well_axis = WellAxis::Z;

    w2.well_name = "PROD_3D";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 5.5e6;
    w2.component_mode = WellComponentMode::Total;
    w2.completion_id = FIM_CaseKit::MatrixBlockCount(mgr) - 1;
    w2.frac_w = 1.0;
    w2.frac_g = 0.0;
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
    using FIM_CaseKit::BoundaryConditionManager;

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

    // =====================================================================
    // 4. Initial Conditions (IC)
    // =====================================================================
    FIM_Engine::InitialConditions ic;
    ic.P_init = 6.0e6;
    ic.T_init = 300.0;
    ic.Sw_init = 0.2;

    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    BoundaryConditionManager bcP;
    BoundaryConditionManager bcT;
    BoundaryConditionManager bcSw;
    FIM_CaseKit::ConfigurePressureBC3D(bcP, 7.0e6, 5.5e6, 6.0e6);
    FIM_CaseKit::ConfigureTemperatureBC3D(bcT, 300.0);
    FIM_CaseKit::ConfigureSaturationBC3D(bcSw, 0.8, 0.2);
    auto modules = FIM_CaseKit::BuildModules3D(preset, &bcP, &bcT, &bcSw);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    WellScheduleStep w1, w2;
    w1.well_name = "INJ_3D_TP";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::BHP;
    w1.target_value = 7.0e6;
    w1.component_mode = WellComponentMode::Water;
    w1.completion_id = 0;
    w1.frac_w = 1.0;
    w1.frac_g = 0.0;
    w1.well_axis = WellAxis::Z;

    w2.well_name = "PROD_3D_TP";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 5.5e6;
    w2.component_mode = WellComponentMode::Total;
    w2.completion_id = FIM_CaseKit::MatrixBlockCount(mgr) - 1;
    w2.frac_w = 0.5;
    w2.frac_g = 0.5;
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
    using FIM_CaseKit::BoundaryConditionManager;

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

    // =====================================================================
    // 4. Initial Conditions (IC)
    // =====================================================================
    FIM_Engine::InitialConditions ic;
    ic.P_init = 6.0e6;
    ic.T_init = 300.0;
    ic.Sw_init = 0.2;

    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    BoundaryConditionManager bcP;
    BoundaryConditionManager bcT;
    BoundaryConditionManager bcSw;
    FIM_CaseKit::ConfigurePressureBC3D(bcP, 7.0e6, 5.5e6, 6.0e6);
    FIM_CaseKit::ConfigureTemperatureBC3D(bcT, 300.0);
    FIM_CaseKit::ConfigureSaturationBC3D(bcSw, 0.8, 0.2);
    auto modules = FIM_CaseKit::BuildModules3D(preset, &bcP, &bcT, &bcSw);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    WellScheduleStep w1, w2, w3;
    w1.well_name = "INJ_W_3D";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::Rate;
    w1.target_value = -0.5;
    w1.component_mode = WellComponentMode::Water;
    w1.completion_id = 0;
    w1.frac_w = 1.0;
    w1.frac_g = 0.0;
    w1.well_axis = WellAxis::Z;

    w2.well_name = "PROD_1_3D";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 5.5e6;
    w2.component_mode = WellComponentMode::Total;
    w2.completion_id = FIM_CaseKit::MatrixBlockCount(mgr) - 1;
    w2.frac_w = 0.5;
    w2.frac_g = 0.5;
    w2.well_axis = WellAxis::Z;

    w3.well_name = "PROD_2_3D";
    w3.domain = WellTargetDomain::Matrix;
    w3.control_mode = WellControlMode::BHP;
    w3.target_value = 5.5e6;
    w3.component_mode = WellComponentMode::Total;
    w3.completion_id = FIM_CaseKit::MatrixBlockCount(mgr) / 2;
    w3.frac_w = 0.5;
    w3.frac_g = 0.5;
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

} // namespace Test_Day6

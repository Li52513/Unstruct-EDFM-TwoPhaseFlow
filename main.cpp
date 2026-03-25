#include "3D_Improved_EDFM_test.h"
#include "2D_EDFM_test.h"
#include "2D_EDFM_MeshTest_Benchmark.h"

#include "3D_TwistFracIntersectionTest.h"
#include "3D_EDFM_GeomTest.h"
#include "3D_EDFM_MeshTest_Benchmark.h"
#include "3D_EDFM_Transmissibility_test.h"
#include "3D_InitializerTest.h"
#include "3D_InitializerTest_multiFrac.h"
#include "3D_PropTest_HD.h"

#include "3D_BoundarySetup_Test.h"

#include "test_FVM_Grad_Benchmark.h"
#include "Test_DOFMapper.h"
#include "Test_ADVar.h"
#include "3D_PropandInit_test.h"
#include "Test_FluidEvaluator.h"
#include "3D_Benchmark_ComplexFractureNetwork.h"

#include "test_Transmissibility_2D.h"
#include "test_Transmissibility_3D.h"
#include "Test_FVM_Ops_AD.h"
#include "Test_Day5_GlobalAssembly_Jacobian.h"
#include "Test_Day6_TransientSolver.h"
#include "FullCaseTest.h"
#include "Test_Issue11_FrozenMatrix.h"
#include "Test_Issue12_LinearSolverMemory.h"
#include "ADVar.hpp"
#include "test_MartixAssemble.h"

#include <array>
#include <exception>
#include <functional>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {
struct CaseEntry {
    std::string name;
    std::string desc;
    std::function<int()> run;
};

int RunDay2FvmAd() {
    std::cout << "\n=======================================\n";
    std::cout << ">>> Running Day2: FVM-AD operator tests <<<\n";
    std::cout << "=======================================\n\n";
    Test_FVM::Run_All_Day2_Tests<3, ADVar<3>>();
    return 0;
}

int RunDay1ArchitectureFreeze() {
    std::cout << "\n===============================================\n";
    std::cout << ">>> Running Day1: Architecture/Connection checks <<<\n";
    std::cout << "=== includes: R4(trans_2d) + R5(trans_3d) ===\n\n";

    Benchmark2D::run_TransmissibilityBenchmark_2D();
    std::cout << "\n-------------------------------------------\n" << std::endl;
    Benchmark3D::run_TransmissibilityBenchmark_3D();

    std::cout << "\n[Day1] Completed R4 + R5 sequence.\n";
    return 0;
}

int RunDay1ArchitectureFreezeRepro() {
    std::cout << "\n=========================================================\n";
    std::cout << ">>> Running Day1 Repro: same-input repeated sequence <<<\n";
    std::cout << "=== trans_2d x2 + trans_3d x2 (for log reproducibility) ===\n\n";

    Benchmark2D::run_TransmissibilityBenchmark_2D();
    std::cout << "\n[Day1-Repro] Repeat trans_2d...\n";
    Benchmark2D::run_TransmissibilityBenchmark_2D();

    std::cout << "\n-------------------------------------------\n" << std::endl;

    Benchmark3D::run_TransmissibilityBenchmark_3D();
    std::cout << "\n[Day1-Repro] Repeat trans_3d...\n";
    Benchmark3D::run_TransmissibilityBenchmark_3D();

    std::cout << "\n[Day1-Repro] Completed repeated sequence.\n";
    return 0;
}

void PrintUsage(const char* exeName) {
    std::cout << "Usage:\n";
    std::cout << "  " << exeName << " --list\n";
    std::cout << "  " << exeName << " --case=<case_name>\n";
    std::cout << "  " << exeName << " --case <case_name>\n";
    std::cout << "  " << exeName << " --help\n";
}

void PrintCases(const std::vector<CaseEntry>& cases) {
    std::cout << "Available cases:\n";
    for (const auto& entry : cases) {
        std::cout << "  - " << entry.name << " : " << entry.desc << '\n';
    }
}
} // namespace

int main (int argc, char** argv) {
    std::vector<CaseEntry> cases = {
        {"day1_arch_conn", "Day1 explicit gate: run R4(trans_2d) + R5(trans_3d)", []() { return RunDay1ArchitectureFreeze(); }},
        {"day1_arch_conn_repro", "Day1 reproducibility gate: trans_2d x2 + trans_3d x2", []() { return RunDay1ArchitectureFreezeRepro(); }},
        {"day2_fvm_ad", "Day2 FVM-AD operator acceptance tests (default)", []() { return RunDay2FvmAd(); }},
        {"day3_bc_patch", "Day3 explicit gate: Boundary condition operators & AD verification", []() { Test_FVM::Run_Day3_BC_Patch<3, ADVar<3>>(); return 0; }},
        {"day3_leakoff_switch", "Day3 explicit gate: Leakoff switch, AD ops & Grid-Level Assembly", []() { Test_FVM::Run_Day3_Leakoff_Switch<3, ADVar<3>>(); return 0; }},
        {"day3_viz", "Day3 explicit gate: Output VTK files for Boundary/Leakoff Assembly", []() { Test_FVM::Run_Day3_BC_Viz_2D<3, ADVar<3>>(); Test_FVM::Run_Day3_Leakoff_Viz_3D<3, ADVar<3>>(); return 0; }},
        {"day4_well_patch", "Day4 explicit gate: Well(BHP/Rate) + WAG schedule skeleton + matrix/fracture completion path", []() { Test_FVM::Run_Day4_Well_Patch<3, ADVar<3>>(); return 0; }},
        {"day4_well_viz", "Day4 visualization gate: export fixed-path 2D/3D well-source VTK files", []() { Test_FVM::Run_Day4_Well_Viz<3, ADVar<3>>(); return 0; }},
        {"day5_block_matrix_robust", "Day5 infra gate: FIM block sparse matrix robustness test", []() { Test_Day5::Test_FIM_BlockSparseMatrix_Robustness(); return 0; }},
        {"issue11_frozen_matrix", "Issue#11: FrozenMatrix CSR cache structure, numerical equivalence, pattern stability & timing", []() { Test_Issue11::Run_All(); return 0; }},
        {"issue12_linear_solver", "Issue#12: Linear solver memory (cache A_work) & CPR-AMG solvability", []() { Test_Issue12::Run_All(); return 0; }},
        { "day5_global_jac_2d","Day5 core gate: 2D FD vs AD Jacobian Assembly Verification",[]() { Test_Day5::Run_Day5_GlobalAssembly_Jacobian_2D(); return 0; }},
        { "day5_global_jac_3d","Day5 core gate: 3D FD vs AD Jacobian Assembly Verification",[]() { Test_Day5::Run_Day5_GlobalAssembly_Jacobian_3D(); return 0; }},
        {"day6_transient_2d_sp_injprod", "Day6: 2D single-phase transient stability + VTK export", []() { Test_Day6::Run_Day6_Transient_2D_SP_InjProd(); return 0; }},
        {"day6_transient_2d_tp_injprod", "Day6: 2D two-phase transient stability + VTK export", []() { Test_Day6::Run_Day6_Transient_2D_TP_InjProd(); return 0; }},
        {"day6_transient_2d_tp_multiwell", "Day6: 2D multi-well interference transient + VTK export", []() { Test_Day6::Run_Day6_Transient_2D_TP_Multiwell(); return 0; }},
        {"day6_transient_3d_sp_injprod", "Day6: 3D single-phase transient stability + VTK export", []() { Test_Day6::Run_Day6_Transient_3D_SP_InjProd(); return 0; }},
        {"day6_transient_3d_tp_injprod", "Day6: 3D two-phase transient stability + VTK export", []() { Test_Day6::Run_Day6_Transient_3D_TP_InjProd(); return 0; }},
        {"day6_transient_3d_tp_multiwell", "Day6: 3D multi-well interference transient + VTK export", []() { Test_Day6::Run_Day6_Transient_3D_TP_Multiwell(); return 0; }},
        {"day6_matrix_audit_2d_edfm", "Day6 explicit gate: 2D EDFM matrix assembly audit (NNC/FF)", []() { Test_Day6::Run_Day6_MatrixAudit_2D_EDFM(); return 0; }},
        {"day6_matrix_audit_3d_edfm", "Day6 explicit gate: 3D EDFM matrix assembly audit (NNC/FF)", []() { Test_Day6::Run_Day6_MatrixAudit_3D_EDFM(); return 0; }},
        {"day6_campaign_2d_all", "Day6 campaign: 2D full T01..T16 with mandatory F0/F1/F2 gating", []() { Test_Day6::Run_Day6_Campaign_2D_All_TxxFy(); return 0; }},
        {"day6_campaign_3d_all", "Day6 campaign: 3D full T01..T16 with mandatory F0/F1/F2 gating", []() { Test_Day6::Run_Day6_Campaign_3D_All_TxxFy(); return 0; }},
        {"day6_t01_f0", "Day6 T01 immediate split gate: F0 no-fracture", []() { Test_Day6::Run_Day6_Campaign_2D_T01_F0(); return 0; }},
        {"day6_t01_f1", "Day6 T01 immediate split gate: F1 single-fracture", []() { Test_Day6::Run_Day6_Campaign_2D_T01_F1(); return 0; }},
        {"day6_t01_f2", "Day6 T01 immediate split gate: F2 crossing-fracture", []() { Test_Day6::Run_Day6_Campaign_2D_T01_F2(); return 0; }},
        {"day6_t1_2d_sp_nowell_analytical", "Day6 T1 baseline: 2D no-well single-phase diffusion with Fourier analytical validation", []() { Test_Day6::Run_Day6_T1_2D_SP_NoWell_Analytical(); return 0; }},
        {"day6l1_2d_sp_co2_const_nowell_analytical", "Day6L1 (AD N=1): 2D pressure-only constant-baseline no-well analytical baseline (forced non-orth + VTK)", []() { FullCaseTest::RunN1L1ConstNoWellAnalytical(); return 0; }},
        {"day6l2_2d_sp_co2_const_nowell_grid", "Day6L2 (AD N=1): 2D pressure-only constant-baseline no-well grid convergence (forced non-orth + VTK)", []() { FullCaseTest::RunN1L2ConstNoWellGrid(); return 0; }},
        {"day6l3_2d_sp_co2_const_nowell_solver", "Day6L3 (AD N=1): 2D pressure-only constant-baseline no-well solver robustness (forced non-orth + VTK)", []() { FullCaseTest::RunN1L3ConstNoWellSolver(); return 0; }},
        {"day6l4_2d_sp_co2_varprop_nowell", "Day6L4 (AD N=1): 2D pressure-only CO2-EOS no-well nonlinear stress (forced non-orth + VTK)", []() { FullCaseTest::RunN1L4VarPropNoWell(); return 0; }},
        {"day6l4_2d_sp_co2_varprop_nowell_singlefrac", "Day6L4 (AD N=1): 2D pressure-only CO2-EOS no-well single-fracture nonlinear stress (forced non-orth + VTK)", []() { FullCaseTest::RunN1L4VarPropNoWellSingleFrac(); return 0; }},
        {"day6ladder_2d_sp_co2_const_nowell_all", "Day6 ladder (AD N=1): run L1+L2+L3+L4 pressure-only debug chain (forced non-orth + VTK)", []() { FullCaseTest::RunN1LadderAll(); return 0; }},
        {"day6l1_2d_sp_co2_const_nowell_analytical_legacy", "Day6L1 legacy: custom SparseLU pressure-only analytical baseline", []() { FullCaseTest::RunN1L1ConstNoWellAnalyticalLegacy(); return 0; }},
        {"day6l2_2d_sp_co2_const_nowell_grid_legacy", "Day6L2 legacy: custom SparseLU pressure-only grid convergence", []() { FullCaseTest::RunN1L2ConstNoWellGridLegacy(); return 0; }},
        {"day6l3_2d_sp_co2_const_nowell_solver_legacy", "Day6L3 legacy: custom SparseLU pressure-only solver robustness", []() { FullCaseTest::RunN1L3ConstNoWellSolverLegacy(); return 0; }},
        {"day6l4_2d_sp_co2_varprop_nowell_legacy", "Day6L4 legacy: custom SparseLU pressure-only variable-property stress", []() { FullCaseTest::RunN1L4VarPropNoWellLegacy(); return 0; }},
        {"day6ladder_2d_sp_co2_const_nowell_all_legacy", "Day6 ladder legacy: run L1+L2+L3+L4 legacy chain", []() { FullCaseTest::RunN1LadderAllLegacy(); return 0; }},
        {"full_n1_template_const_nowell_nofrac", "FullCase N=1 template: 2D single-phase CO2 const-property no-well no-fracture", []() { FullCaseTest::RunN1TemplateConstNoWellNoFrac(); return 0; }},
        {"full_n1_template_const_nowell_singlefrac", "FullCase N=1 template: 2D single-phase CO2 const-property no-well single-fracture", []() { FullCaseTest::RunN1TemplateConstNoWellSingleFrac(); return 0; }},
        {"full_n1_template_const_nowell_crossfrac", "FullCase N=1 template: 2D single-phase CO2 const-property no-well cross-fracture", []() { FullCaseTest::RunN1TemplateConstNoWellCrossFrac(); return 0; }},
        {"2d_edfm_single", "2D EDFM single-fracture end-to-end test", []() { return EDFM_test_2D(); }},
        {"2d_edfm_dfn", "2D EDFM DFN end-to-end test", []() { return EDFM_DFN_test_2D(); }},
        {"2d_geom_benchmark_dfn", "2D EDFM geometry benchmark with fixed DFN seed", []() { return EDFM_DFN_Geomtest_2D(); }},
        {"3d_twist_intersection", "3D twisted-fracture intersection test", []() { return TwistFracIntersectionTest_3D(); }},
        {"3d_edfm_improved", "3D improved EDFM end-to-end test", []() { return Improved_EDFM_test_3D(); }},
        {"3d_distance_accuracy", "3D matrix-to-fracture distance accuracy benchmark", []() { return RunBenchmark_Step1_Distance_Accuracy(); }},
        {"3d_nnc_ff_static", "3D static transmissibility check for NNC/FF", []() { RunBenchmark_Transmissibility_Static(); return 0; }},
        {"3d_mesh_benchmark", "3D EDFM mesh/geometry benchmark suite", []() { return Improved_3D_EDFM_MeshTest(); }},
        {"init_single_frac", "3D initializer visualization test (single fracture)", []() { RunTest_Initialization_And_Viz(); return 0; }},
        {"init_dfn", "3D initializer visualization test (DFN)", []() { RunTest_DFN_Initialization_And_Viz(); return 0; }},
        {"property_sweep", "3D CO2/water property accuracy sweep", []() { RunTest_Property_Accuracy_Sweep(); return 0; }},
        {"boundary_export", "3D boundary-tag setup/export verification", []() { RunTest_BoundaryCondition_Export(); return 0; }},
        {"grad_2d", "2D gradient operator benchmark", []() { run_Benchmark_2D_EDFM_Grad(); return 0; }},
        {"grad_3d", "3D gradient operator benchmark", []() { run_Benchmark_3D_EDFM_Grad(); return 0; }},
        {"grad_all", "2D + 3D gradient operator benchmarks", []() { run_Benchmark_2D_EDFM_Grad(); run_Benchmark_3D_EDFM_Grad(); return 0; }},
        {"dof_mapper", "DOF mapper verification", []() { RunDOFMapperVerification(); return 0; }},
        {"advar", "ADVar comprehensive tests", []() { Run_ADVar_Comprehensive_Tests(); return 0; }},
        {"prop_init_3d", "3D property+initializer benchmark", []() { RunBenchmark_3D_PropTest(); return 0; }},
        {"fluid_eval", "AD fluid evaluator tests", []() { run_fluid_evaluator_test(); return 0; }},
        {"complex_frac_benchmark", "3D complex fracture network benchmark", []() { RunBenchmark_ComplexFractureNetwork(); return 0; }},
        {"trans_2d", "2D transmissibility benchmark (MM/FI/NNC/FF)", []() { Benchmark2D::run_TransmissibilityBenchmark_2D(); return 0; }},
        {"trans_3d", "3D transmissibility benchmark (MM/FI/NNC/FF)", []() { Benchmark3D::run_TransmissibilityBenchmark_3D(); return 0; }},
        {"trans_all", "2D + 3D transmissibility benchmarks", []() { Benchmark2D::run_TransmissibilityBenchmark_2D(); Benchmark3D::run_TransmissibilityBenchmark_3D(); return 0; }},
        // ---- ValidationSuite migrated into FullCaseTest dispatcher ----
        {"val_t1_const_nofrac",    "ValidationSuite T1: Const NoFrac (analytical)",      []() { FullCaseTest::RunValidationCaseByKey("val_t1_const_nofrac"); return 0; }},
        {"val_t1_const_singlefrac","ValidationSuite T1: Const SingleFrac (stability)",   []() { FullCaseTest::RunValidationCaseByKey("val_t1_const_singlefrac"); return 0; }},
        {"val_t1_const_crossfrac", "ValidationSuite T1: Const CrossFrac (stability)",    []() { FullCaseTest::RunValidationCaseByKey("val_t1_const_crossfrac"); return 0; }},
        {"val_t1_var_nofrac",      "ValidationSuite T1: Var NoFrac (EOS cross-check)",   []() { FullCaseTest::RunValidationCaseByKey("val_t1_var_nofrac"); return 0; }},
        {"val_t1_var_singlefrac",  "ValidationSuite T1: Var SingleFrac (EOS stability)", []() { FullCaseTest::RunValidationCaseByKey("val_t1_var_singlefrac"); return 0; }},
        {"val_t1_var_crossfrac",   "ValidationSuite T1: Var CrossFrac (EOS stability)",  []() { FullCaseTest::RunValidationCaseByKey("val_t1_var_crossfrac"); return 0; }},
        {"val_t2_const_nofrac",    "ValidationSuite T2: Const NoFrac (analytical)",      []() { FullCaseTest::RunValidationCaseByKey("val_t2_const_nofrac"); return 0; }},
        {"val_t2_const_singlefrac","ValidationSuite T2: Const SingleFrac (stability)",   []() { FullCaseTest::RunValidationCaseByKey("val_t2_const_singlefrac"); return 0; }},
        {"val_t2_const_crossfrac", "ValidationSuite T2: Const CrossFrac (stability)",    []() { FullCaseTest::RunValidationCaseByKey("val_t2_const_crossfrac"); return 0; }},
        {"val_t2_var_nofrac",      "ValidationSuite T2: Var NoFrac (EOS cross-check)",   []() { FullCaseTest::RunValidationCaseByKey("val_t2_var_nofrac"); return 0; }},
        {"val_t2_var_singlefrac",  "ValidationSuite T2: Var SingleFrac (EOS stability)", []() { FullCaseTest::RunValidationCaseByKey("val_t2_var_singlefrac"); return 0; }},
        {"val_t2_var_crossfrac",   "ValidationSuite T2: Var CrossFrac (EOS stability)",  []() { FullCaseTest::RunValidationCaseByKey("val_t2_var_crossfrac"); return 0; }},
        {"val_t3_const_nofrac",    "ValidationSuite T3: Const NoFrac (analytical)",      []() { FullCaseTest::RunValidationCaseByKey("val_t3_const_nofrac"); return 0; }},
        {"val_t3_const_singlefrac","ValidationSuite T3: Const SingleFrac (stability)",   []() { FullCaseTest::RunValidationCaseByKey("val_t3_const_singlefrac"); return 0; }},
        {"val_t3_const_crossfrac", "ValidationSuite T3: Const CrossFrac (stability)",    []() { FullCaseTest::RunValidationCaseByKey("val_t3_const_crossfrac"); return 0; }},
        {"val_t3_var_nofrac",      "ValidationSuite T3: Var NoFrac (EOS cross-check)",   []() { FullCaseTest::RunValidationCaseByKey("val_t3_var_nofrac"); return 0; }},
        {"val_t3_var_singlefrac",  "ValidationSuite T3: Var SingleFrac (EOS stability)", []() { FullCaseTest::RunValidationCaseByKey("val_t3_var_singlefrac"); return 0; }},
        {"val_t3_var_crossfrac",   "ValidationSuite T3: Var CrossFrac (EOS stability)",  []() { FullCaseTest::RunValidationCaseByKey("val_t3_var_crossfrac"); return 0; }},
        {"val_t4_const_nofrac",    "ValidationSuite T4: Const NoFrac (analytical)",      []() { FullCaseTest::RunValidationCaseByKey("val_t4_const_nofrac"); return 0; }},
        {"val_t4_const_singlefrac","ValidationSuite T4: Const SingleFrac (stability)",   []() { FullCaseTest::RunValidationCaseByKey("val_t4_const_singlefrac"); return 0; }},
        {"val_t4_const_crossfrac", "ValidationSuite T4: Const CrossFrac (stability)",    []() { FullCaseTest::RunValidationCaseByKey("val_t4_const_crossfrac"); return 0; }},
        {"val_t4_var_nofrac",      "ValidationSuite T4: Var NoFrac (EOS cross-check)",   []() { FullCaseTest::RunValidationCaseByKey("val_t4_var_nofrac"); return 0; }},
        {"val_t4_var_singlefrac",  "ValidationSuite T4: Var SingleFrac (EOS stability)", []() { FullCaseTest::RunValidationCaseByKey("val_t4_var_singlefrac"); return 0; }},
        {"val_t4_var_crossfrac",   "ValidationSuite T4: Var CrossFrac (EOS stability)",  []() { FullCaseTest::RunValidationCaseByKey("val_t4_var_crossfrac"); return 0; }},
        {"val_t5_const_nofrac",    "ValidationSuite T5: Const NoFrac (analytical)",      []() { FullCaseTest::RunValidationCaseByKey("val_t5_const_nofrac"); return 0; }},
        {"val_t5_const_singlefrac","ValidationSuite T5: Const SingleFrac (stability)",   []() { FullCaseTest::RunValidationCaseByKey("val_t5_const_singlefrac"); return 0; }},
        {"val_t5_const_crossfrac", "ValidationSuite T5: Const CrossFrac (stability)",    []() { FullCaseTest::RunValidationCaseByKey("val_t5_const_crossfrac"); return 0; }},
        {"val_t5_var_nofrac",      "ValidationSuite T5: Var NoFrac (EOS cross-check)",   []() { FullCaseTest::RunValidationCaseByKey("val_t5_var_nofrac"); return 0; }},
        {"val_t5_var_singlefrac",  "ValidationSuite T5: Var SingleFrac (EOS stability)", []() { FullCaseTest::RunValidationCaseByKey("val_t5_var_singlefrac"); return 0; }},
        {"val_t5_var_crossfrac",   "ValidationSuite T5: Var CrossFrac (EOS stability)",  []() { FullCaseTest::RunValidationCaseByKey("val_t5_var_crossfrac"); return 0; }},
        {"val_t6_const_nofrac",    "ValidationSuite T6: Const NoFrac (analytical)",      []() { FullCaseTest::RunValidationCaseByKey("val_t6_const_nofrac"); return 0; }},
        {"val_t6_const_singlefrac","ValidationSuite T6: Const SingleFrac (stability)",   []() { FullCaseTest::RunValidationCaseByKey("val_t6_const_singlefrac"); return 0; }},
        {"val_t6_const_crossfrac", "ValidationSuite T6: Const CrossFrac (stability)",    []() { FullCaseTest::RunValidationCaseByKey("val_t6_const_crossfrac"); return 0; }},
        {"val_t6_var_nofrac",      "ValidationSuite T6: Var NoFrac (EOS cross-check)",   []() { FullCaseTest::RunValidationCaseByKey("val_t6_var_nofrac"); return 0; }},
        {"val_t6_var_singlefrac",  "ValidationSuite T6: Var SingleFrac (EOS stability)", []() { FullCaseTest::RunValidationCaseByKey("val_t6_var_singlefrac"); return 0; }},
        {"val_t6_var_crossfrac",   "ValidationSuite T6: Var CrossFrac (EOS stability)",  []() { FullCaseTest::RunValidationCaseByKey("val_t6_var_crossfrac"); return 0; }},
        {"val_t7_const_nofrac",    "ValidationSuite T7: Const NoFrac (analytical)",      []() { FullCaseTest::RunValidationCaseByKey("val_t7_const_nofrac"); return 0; }},
        {"val_t7_const_singlefrac","ValidationSuite T7: Const SingleFrac (stability)",   []() { FullCaseTest::RunValidationCaseByKey("val_t7_const_singlefrac"); return 0; }},
        {"val_t7_const_crossfrac", "ValidationSuite T7: Const CrossFrac (stability)",    []() { FullCaseTest::RunValidationCaseByKey("val_t7_const_crossfrac"); return 0; }},
        {"val_t7_var_nofrac",      "ValidationSuite T7: Var NoFrac (EOS cross-check)",   []() { FullCaseTest::RunValidationCaseByKey("val_t7_var_nofrac"); return 0; }},
        {"val_t7_var_singlefrac",  "ValidationSuite T7: Var SingleFrac (EOS stability)", []() { FullCaseTest::RunValidationCaseByKey("val_t7_var_singlefrac"); return 0; }},
        {"val_t7_var_crossfrac",   "ValidationSuite T7: Var CrossFrac (EOS stability)",  []() { FullCaseTest::RunValidationCaseByKey("val_t7_var_crossfrac"); return 0; }},
        {"val_t8_const_nofrac",    "ValidationSuite T8: Const NoFrac (analytical)",      []() { FullCaseTest::RunValidationCaseByKey("val_t8_const_nofrac"); return 0; }},
        {"val_t8_const_singlefrac","ValidationSuite T8: Const SingleFrac (stability)",   []() { FullCaseTest::RunValidationCaseByKey("val_t8_const_singlefrac"); return 0; }},
        {"val_t8_const_crossfrac", "ValidationSuite T8: Const CrossFrac (stability)",    []() { FullCaseTest::RunValidationCaseByKey("val_t8_const_crossfrac"); return 0; }},
        {"val_t8_var_nofrac",      "ValidationSuite T8: Var NoFrac (EOS cross-check)",   []() { FullCaseTest::RunValidationCaseByKey("val_t8_var_nofrac"); return 0; }},
        {"val_t8_var_singlefrac",  "ValidationSuite T8: Var SingleFrac (EOS stability)", []() { FullCaseTest::RunValidationCaseByKey("val_t8_var_singlefrac"); return 0; }},
        {"val_t8_var_crossfrac",   "ValidationSuite T8: Var CrossFrac (EOS stability)",  []() { FullCaseTest::RunValidationCaseByKey("val_t8_var_crossfrac"); return 0; }},
        {"val_t1", "ValidationSuite T1: no-well single-phase pressure diffusion (6 variants)", []() { FullCaseTest::RunValidationGroup(1); return 0; }},
        {"val_t2", "ValidationSuite T2: no-well single-phase thermal diffusion (6 variants)", []() { FullCaseTest::RunValidationGroup(2); return 0; }},
        {"val_t3", "ValidationSuite T3: no-well two-phase Buckley-Leverett (6 variants)", []() { FullCaseTest::RunValidationGroup(3); return 0; }},
        {"val_t4", "ValidationSuite T4: no-well two-phase + heat (6 variants)", []() { FullCaseTest::RunValidationGroup(4); return 0; }},
        {"val_t5", "ValidationSuite T5: with-well single-phase Theis (6 variants)", []() { FullCaseTest::RunValidationGroup(5); return 0; }},
        {"val_t6", "ValidationSuite T6: with-well single-phase thermal front (6 variants)", []() { FullCaseTest::RunValidationGroup(6); return 0; }},
        {"val_t7", "ValidationSuite T7: with-well two-phase radial BL (6 variants)", []() { FullCaseTest::RunValidationGroup(7); return 0; }},
        {"val_t8", "ValidationSuite T8: with-well two-phase radial BL + heat (6 variants)", []() { FullCaseTest::RunValidationGroup(8); return 0; }},
        {"val_suite_all", "ValidationSuite: run all 48 validation scenarios T1..T8", []() { FullCaseTest::RunValidationAllMigrated(); return 0; }},
    };

    const std::array<std::pair<const char*, FullCaseTest::TopologyVariant>, 3> topologyDefs = {{
        {"nofrac", FullCaseTest::TopologyVariant::NoFrac},
        {"singlefrac", FullCaseTest::TopologyVariant::SingleFrac},
        {"crossfrac", FullCaseTest::TopologyVariant::CrossFrac}
    }};

    for (int scenarioId = 1; scenarioId <= 8; ++scenarioId) {
        const std::string tLabel = "t" + std::to_string(scenarioId);
        for (const auto& topologyDef : topologyDefs) {
            const std::string case2D = "2d_" + tLabel + "_" + topologyDef.first;
            const std::string desc2D = "FullCase 2D " + tLabel + " " + topologyDef.first + " (Const topology split)";
            cases.push_back({ case2D, desc2D, [scenarioId, topology = topologyDef.second]() {
                FullCaseTest::Run2DCase(scenarioId, topology);
                return 0;
            }});

            const std::string case3D = "3d_" + tLabel + "_" + topologyDef.first;
            const std::string desc3D = "FullCase 3D " + tLabel + " " + topologyDef.first + " (Day6 3D mirrored)";
            cases.push_back({ case3D, desc3D, [scenarioId, topology = topologyDef.second]() {
                FullCaseTest::Run3DCase(scenarioId, topology);
                return 0;
            }});
        }
    }

    cases.push_back({"2d_all", "FullCase: run all 2D cases (T1..T8 x NoFrac/SingleFrac/CrossFrac)", []() {
        FullCaseTest::Run2DAll();
        return 0;
    }});
    cases.push_back({"3d_all", "FullCase: run all 3D cases (T1..T8 x NoFrac/SingleFrac/CrossFrac)", []() {
        FullCaseTest::Run3DAll();
        return 0;
    }});
    cases.push_back({"all", "FullCase: run all 2D+3D cases (48 total)", []() {
        FullCaseTest::RunAll();
        return 0;
    }});

    bool showHelp = false;
    bool showList = false;
    std::string caseName;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            showHelp = true;
            continue;
        }
        if (arg == "--list") {
            showList = true;
            continue;
        }
        if (arg == "--case") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --case requires a value.\n";
                return 2;
            }
            caseName = argv[++i];
            continue;
        }
        if (arg.rfind("--case=", 0) == 0) {
            caseName = arg.substr(7);
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        PrintUsage(argv[0]);
        return 2;
    }

    if (showHelp) {
        PrintUsage(argv[0]);
        std::cout << '\n';
        PrintCases(cases);
        return 0;
    }

    if (showList) {
        PrintCases(cases);
        return 0;
    }

    if (caseName.empty()) {
        caseName = "day2_fvm_ad";
        std::cout << "[Info] no --case provided. Using default: " << caseName << '\n';
    }

    for (const auto& entry : cases) {
        if (caseName == entry.name) {
            std::cout << "[Case] " << entry.name << " - " << entry.desc << "\n";
            try {
                return entry.run();
            }
            catch (const std::exception& e) {
                std::cerr << "[Error] exception: " << e.what() << '\n';
                return 1;
            }
            catch (...) {
                std::cerr << "[Error] unknown exception.\n";
                return 1;
            }
        }
    }

    std::cerr << "Unknown case: " << caseName << "\n\n";
    PrintCases(cases);
    return 2;
}

//int main()
//{
//    MatrixAssemblerTest::Test_BlockSparseMatrix();
//    return 0;
//}

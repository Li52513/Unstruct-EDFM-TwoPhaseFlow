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
#include "Test_Day3_BoundaryFixes.h"
#include "Test_Day5_GlobalAssembly_Jacobian.h"
#if __has_include("Test_Day6_TransientSolver.h")
#include "Test_Day6_TransientSolver.h"
#define CODEX_HAS_DAY6_TRANSIENT_SOLVER 1
#else
#define CODEX_HAS_DAY6_TRANSIENT_SOLVER 0
#endif
#include "Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.h"
#include "Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_NoWell.h"
#include "Test_2D_EDFM_H_T_CO2_VaryPP_NoFrac_NoWell.h"
#include "Test_2D_EDFM_H_T_CO2_ConstPP_ComplexFrac_NoWell.h"
#include "Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.h"
#include "Test_2D_EDFM_H_T_CO2_VaryPP_SingleFrac_NoWell.h"
#include "Test_2D_EDFM_H_T_CO2_VaryPP_ComplexFrac_NoWell.h"
#include "Test_2D_EDFM_H_CO2_ConstPP_ComplexFrac_NoWell.h"
#include "Test_2D_EDFM_H_CO2_ConstPP_SingleFrac_NoWell.h"
#include "Test_2D_EDFM_H_CO2_VaryPP_NoFrac_NoWell.h"
#include "Test_2D_EDFM_H_CO2_VaryPP_ComplexFrac_NoWell.h"
#include "Test_2D_EDFM_H_CO2_VaryPP_SingleFrac_NoWell.h"
#include "Test_2D_EDFM_H_TP_CO2H2O_VaryPP_NoFrac_NoWell.h"
#include "Test_2D_EDFM_H_TP_CO2H2O_ConstPP_SingleFrac_NoWell.h"
#include "Test_2D_EDFM_H_TP_CO2H2O_VaryPP_SingleFrac_NoWell.h"
#include "Test_2D_EDFM_H_TP_CO2H2O_ConstPP_ComplexFrac_NoWell.h"
#include "Test_2D_EDFM_H_TP_CO2H2O_VaryPP_ComplexFrac_NoWell.h"
#include "Test_Issue11_FrozenMatrix.h"
#include "Test_Issue12_LinearSolverMemory.h"
#include "CaseCommon_Catalog.h"
#include "ADVar.hpp"
#include "test_MartixAssemble.h"

#include <exception>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace {
struct CaseEntry {
    std::string name;
    std::string desc;
    std::function<int()> run;
};

struct CatalogAliasEntry {
    std::string name;
    std::string desc;
    std::string catalog_case;
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
    std::cout << "  " << exeName << " --case <case_name> --stage <solve_only|prepare_reference|validate_only|full_workflow>\n";
    std::cout << "  " << exeName << " --help\n";
}

void PrintCatalogCases() {
    CaseCommon::ValidateCaseCatalogOrThrow();
    std::cout << "A1-F12 catalog cases:\n";
    for (const auto& entry : CaseCommon::GetCaseCatalog()) {
        std::cout << "  - " << entry.metadata.case_code
                  << " [" << entry.metadata.implementation_status << "]"
                  << " : " << entry.metadata.dispatcher_key
                  << " : " << entry.metadata.description
                  << " : ref=" << entry.metadata.reference_mode << '\n';
    }
}

void PrintCases(const std::vector<CaseEntry>& cases) {
    std::cout << "Legacy/auxiliary cases:\n";
    for (const auto& entry : cases) {
        std::cout << "  - " << entry.name << " : " << entry.desc << '\n';
    }
}

void PrintCatalogAliases(const std::vector<CatalogAliasEntry>& aliases) {
    std::cout << "Catalog compatibility aliases:\n";
    for (const auto& entry : aliases) {
        std::cout << "  - " << entry.name
                  << " -> " << entry.catalog_case
                  << " : " << entry.desc << '\n';
    }
}

const CatalogAliasEntry* ResolveCatalogCaseAlias(const std::vector<CatalogAliasEntry>& aliases,
                                                 const std::string& name) {
    for (const auto& entry : aliases) {
        if (entry.name == name) {
            return &entry;
        }
    }
    return nullptr;
}
} // namespace

int main (int argc, char** argv) {
    std::vector<CatalogAliasEntry> legacyCatalogAliases = {
        {"test_h_co2_constpp_nofrac_nowell", "Compatibility alias for the canonical A1 no-well full workflow.", "A1"},
        {"test_h_t_co2_constpp_nofrac_nowell", "Compatibility alias for the canonical B1 no-well full workflow.", "B1"},
        {"test_h_tp_co2h2o_constpp_nofrac_nowell", "Compatibility alias for the canonical C1 no-well full workflow.", "C1"},
    };

    std::vector<CaseEntry> cases = {
        {"day1_arch_conn", "Day1 explicit gate: run R4(trans_2d) + R5(trans_3d)", []() { return RunDay1ArchitectureFreeze(); }},
        {"day1_arch_conn_repro", "Day1 reproducibility gate: trans_2d x2 + trans_3d x2", []() { return RunDay1ArchitectureFreezeRepro(); }},
        {"day2_fvm_ad", "Day2 FVM-AD operator acceptance tests (default)", []() { return RunDay2FvmAd(); }},
        {"day3_bc_patch", "Day3 explicit gate: Boundary condition operators & AD verification", []() { Test_FVM::Run_Day3_BC_Patch<3, ADVar<3>>(); Test_Day3_BoundaryFixes::RunPatchChecks(); return 0; }},
        {"day3_leakoff_switch", "Day3 explicit gate: Leakoff switch, AD ops & Grid-Level Assembly", []() { Test_FVM::Run_Day3_Leakoff_Switch<3, ADVar<3>>(); Test_Day3_BoundaryFixes::RunLeakoffChecks(); return 0; }},
        {"day3_viz", "Day3 explicit gate: Output VTK files for Boundary/Leakoff Assembly", []() { Test_FVM::Run_Day3_BC_Viz_2D<3, ADVar<3>>(); Test_FVM::Run_Day3_Leakoff_Viz_3D<3, ADVar<3>>(); return 0; }},
        {"day4_well_patch", "Day4 explicit gate: Well(BHP/Rate) + WAG schedule skeleton + matrix/fracture completion path", []() { Test_FVM::Run_Day4_Well_Patch<3, ADVar<3>>(); return 0; }},
        {"day4_well_viz", "Day4 visualization gate: export fixed-path 2D/3D well-source VTK files", []() { Test_FVM::Run_Day4_Well_Viz<3, ADVar<3>>(); return 0; }},
        {"day5_block_matrix_robust", "Day5 infra gate: FIM block sparse matrix robustness test", []() { Test_Day5::Test_FIM_BlockSparseMatrix_Robustness(); return 0; }},
        {"issue11_frozen_matrix", "Issue#11: FrozenMatrix CSR cache structure, numerical equivalence, pattern stability & timing", []() { Test_Issue11::Run_All(); return 0; }},
        {"issue12_linear_solver", "Issue#12: Linear solver memory (cache A_work) & CPR-AMG solvability", []() { Test_Issue12::Run_All(); return 0; }},
        { "day5_global_jac_2d","Day5 core gate: 2D FD vs AD Jacobian Assembly Verification",[]() { Test_Day5::Run_Day5_GlobalAssembly_Jacobian_2D(); return 0; }},
        { "day5_global_jac_3d","Day5 core gate: 3D FD vs AD Jacobian Assembly Verification",[]() { Test_Day5::Run_Day5_GlobalAssembly_Jacobian_3D(); return 0; }},
#if CODEX_HAS_DAY6_TRANSIENT_SOLVER
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
#endif
        {"test_h_t_co2_constpp_nofrac_nowell_fast", "Standalone fast-screen test: 2D single-phase CO2 const-property P-T coupled no-fracture no-well with t_end=1.0e7 s", []() { Test_H_T_CO2_ConstPP_NoFrac::ExecutePlanByKey("h_t_co2_constpp_nofrac_nowell_fast"); return 0; }},
        {"test_h_t_co2_constpp_nofrac_nowell_fast_grid", "Standalone fast-screen grid study: 2D single-phase CO2 const-property P-T coupled no-fracture no-well with t_end=1.0e7 s", []() { Test_H_T_CO2_ConstPP_NoFrac::ExecutePlanByKey("h_t_co2_constpp_nofrac_nowell_fast_grid"); return 0; }},
        {"test_h_t_co2_constpp_nofrac_nowell_grid", "Standalone validation-chain grid study: 2D single-phase CO2 const-property P-T coupled no-fracture no-well", []() { Test_H_T_CO2_ConstPP_NoFrac::ExecutePlanByKey("h_t_co2_constpp_nofrac_nowell_grid"); return 0; }},
        {"test_h_t_co2_constpp_nofrac_nowell_dt", "Standalone validation-chain time-step study: 2D single-phase CO2 const-property P-T coupled no-fracture no-well", []() { Test_H_T_CO2_ConstPP_NoFrac::ExecutePlanByKey("h_t_co2_constpp_nofrac_nowell_dt"); return 0; }},
        {"test_h_t_co2_constpp_nofrac_nowell_all", "Standalone validation-chain full study: 2D single-phase CO2 const-property P-T coupled no-fracture no-well", []() { Test_H_T_CO2_ConstPP_NoFrac::ExecutePlanByKey("h_t_co2_constpp_nofrac_nowell_all"); return 0; }},
        {"test_h_t_co2_constpp_nofrac_nowell_debug_96x12", "Standalone debug case: 96x12 fine-grid detailed-output run for the no-fracture P-T coupled validation case", []() { Test_H_T_CO2_ConstPP_NoFrac::ExecutePlanByKey("h_t_co2_constpp_nofrac_nowell_debug_96x12"); return 0; }},
        {"test_h_t_co2_constpp_nofrac_nowell_probe_96x12_dt5000", "Standalone long-time probe: 96x12, dt_init=5000 s no-fracture P-T coupled validation case", []() { Test_H_T_CO2_ConstPP_NoFrac::ExecutePlanByKey("h_t_co2_constpp_nofrac_nowell_probe_96x12_dt5000"); return 0; }},
        {"test_h_t_co2_varypp_nofrac_nowell", "Standalone test: 2D single-phase CO2 variable-property (EOS) P-T coupled no-fracture no-well", []() { Test_H_T_CO2_VaryPP_NoFrac::RunFullWorkflow(); return 0; }},
        {"test_h_t_co2_constpp_singlefrac_nowell", "Standalone test: 2D single-phase CO2 const-property P-T coupled single-fracture no-well", []() { Test_H_T_CO2_ConstPP_SingleFrac::RunFullWorkflow(); return 0; }},
        {"test_h_t_co2_constpp_singlefrac_nowell_prepare_reference", "Standalone validation-chain prepare phase: export engineering reference inputs for the N=2 single-fracture CO2 P-T case", []() { Test_H_T_CO2_ConstPP_SingleFrac::ExecutePlanByKey("h_t_co2_constpp_singlefrac_nowell_prepare_reference"); return 0; }},
        {"test_h_t_co2_constpp_singlefrac_nowell_validate_reference", "Standalone validation-chain validate phase: compare the N=2 single-fracture CO2 P-T case against prepared reference data", []() { Test_H_T_CO2_ConstPP_SingleFrac::ExecutePlanByKey("h_t_co2_constpp_singlefrac_nowell_validate_reference"); return 0; }},
        {"test_h_t_co2_constpp_singlefrac_nowell_full_comsol", "Standalone validation-chain full COMSOL mode: auto-prepare and validate the N=2 single-fracture CO2 P-T case", []() { Test_H_T_CO2_ConstPP_SingleFrac::ExecutePlanByKey("h_t_co2_constpp_singlefrac_nowell_full_comsol"); return 0; }},
        {"test_h_t_co2_varypp_singlefrac_nowell", "Standalone test: 2D single-phase CO2 variable-property (EOS) P-T coupled single-fracture no-well", []() { Test_H_T_CO2_VaryPP_SingleFrac::RunFullWorkflow(); return 0; }},
        {"test_h_t_co2_constpp_complexfrac_nowell", "Standalone test: 2D single-phase CO2 const-property P-T coupled complex-fracture no-well", []() { Test_H_T_CO2_ConstPP_ComplexFrac::RunFullWorkflow(); return 0; }},
        {"test_h_t_co2_varypp_complexfrac_nowell", "Standalone test: 2D single-phase CO2 variable-property (EOS) P-T coupled complex-fracture no-well", []() { Test_H_T_CO2_VaryPP_ComplexFrac::RunFullWorkflow(); return 0; }},
        {"test_h_tp_co2h2o_varypp_nofrac_nowell", "Standalone test: 2D two-phase CO2/H2O variable-property (EOS) P-T coupled no-fracture no-well", []() { Test_H_TP_CO2H2O_VaryPP_NoFrac::RunFullWorkflow(); return 0; }},
        {"test_h_tp_co2h2o_constpp_singlefrac_nowell", "Standalone test: 2D two-phase CO2/H2O const-property P-T coupled single-fracture no-well", []() { Test_H_TP_CO2H2O_ConstPP_SingleFrac::RunFullWorkflow(); return 0; }},
        {"test_h_tp_co2h2o_varypp_singlefrac_nowell", "Standalone test: 2D two-phase CO2/H2O variable-property (EOS) P-T coupled single-fracture no-well", []() { Test_H_TP_CO2H2O_VaryPP_SingleFrac::RunFullWorkflow(); return 0; }},
        {"test_h_tp_co2h2o_constpp_complexfrac_nowell", "Standalone test: 2D two-phase CO2/H2O const-property P-T coupled complex-fracture no-well", []() { Test_H_TP_CO2H2O_ConstPP_ComplexFrac::RunFullWorkflow(); return 0; }},
        {"test_h_tp_co2h2o_varypp_complexfrac_nowell", "Standalone test: 2D two-phase CO2/H2O variable-property (EOS) P-T coupled complex-fracture no-well", []() { Test_H_TP_CO2H2O_VaryPP_ComplexFrac::RunFullWorkflow(); return 0; }},
        {"test_h_co2_constpp_singlefrac_nowell", "Standalone test: 2D single-phase CO2 const-property single-fracture no-well", []() { Test_H_CO2_ConstPP_SingleFrac::RunFullWorkflow(); return 0; }},
        {"test_h_co2_constpp_complexfrac_nowell", "Standalone test: 2D single-phase CO2 const-property complex-fracture no-well", []() { Test_H_CO2_ConstPP_ComplexFrac::RunFullWorkflow(); return 0; }},
        {"test_h_co2_varypp_nofrac_nowell", "Standalone test: 2D single-phase CO2 variable-property (EOS) no-fracture no-well", []() { Test_H_CO2_VaryPP::RunFullWorkflow(); return 0; }},
        {"test_h_co2_varypp_singlefrac_nowell", "Standalone test: 2D single-phase CO2 variable-property (EOS) single-fracture no-well", []() { Test_H_CO2_VaryPP_SingleFrac::RunFullWorkflow(); return 0; }},
        {"test_h_co2_varypp_complexfrac_nowell", "Standalone test: 2D single-phase CO2 variable-property (EOS) complex-fracture no-well", []() { Test_H_CO2_VaryPP_ComplexFrac::RunFullWorkflow(); return 0; }},
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
    };

    bool showHelp = false;
    bool showList = false;
    std::string caseName;
    CaseCommon::CaseStage caseStage = CaseCommon::CaseStage::FullWorkflow;
    bool stageWasSet = false;

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
        if (arg == "--stage") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --stage requires a value.\n";
                return 2;
            }
            if (!CaseCommon::ParseCaseStage(argv[++i], &caseStage)) {
                std::cerr << "Error: unsupported stage value.\n";
                return 2;
            }
            stageWasSet = true;
            continue;
        }
        if (arg.rfind("--stage=", 0) == 0) {
            if (!CaseCommon::ParseCaseStage(arg.substr(8), &caseStage)) {
                std::cerr << "Error: unsupported stage value.\n";
                return 2;
            }
            stageWasSet = true;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        PrintUsage(argv[0]);
        return 2;
    }

    if (showHelp) {
        PrintUsage(argv[0]);
        std::cout << '\n';
        PrintCatalogCases();
        std::cout << '\n';
        PrintCatalogAliases(legacyCatalogAliases);
        std::cout << '\n';
        PrintCases(cases);
        return 0;
    }

    if (showList) {
        PrintCatalogCases();
        std::cout << '\n';
        PrintCatalogAliases(legacyCatalogAliases);
        std::cout << '\n';
        PrintCases(cases);
        return 0;
    }

    if (caseName.empty()) {
        caseName = "day2_fvm_ad";
        std::cout << "[Info] no --case provided. Using default: " << caseName << '\n';
    }

    const CatalogAliasEntry* resolvedAlias = ResolveCatalogCaseAlias(legacyCatalogAliases, caseName);
    const std::string dispatchCaseName = resolvedAlias ? resolvedAlias->catalog_case : caseName;

    try {
        CaseCommon::ValidateCaseCatalogOrThrow();
        if (const CaseCommon::CaseCatalogEntry* catalogEntry = CaseCommon::FindCaseCatalogEntry(dispatchCaseName)) {
            if (resolvedAlias) {
                std::cout << "[Catalog Alias] " << caseName
                          << " -> " << resolvedAlias->catalog_case
                          << " - stage=" << CaseCommon::ToString(caseStage) << "\n";
            }
            std::cout << "[Catalog Case] " << catalogEntry->metadata.case_code
                      << " - " << catalogEntry->metadata.dispatcher_key
                      << " - stage=" << CaseCommon::ToString(caseStage) << "\n";
            return CaseCommon::RunCatalogCase(*catalogEntry, caseStage);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "[Catalog Error] " << e.what() << '\n';
        return 1;
    }

    if (stageWasSet) {
        std::cerr << "Error: --stage is only supported for A1-F12 catalog cases.\n";
        return 2;
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
    PrintCatalogCases();
    std::cout << '\n';
    PrintCases(cases);
    return 2;
}

//int main()
//{
//    MatrixAssemblerTest::Test_BlockSparseMatrix();
//    return 0;
//}

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
#include "ADVar.hpp"

#include <exception>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace {
struct CaseEntry {
    const char* name;
    const char* desc;
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

int main(int argc, char** argv) {
    const std::vector<CaseEntry> cases = {
        {"day1_arch_conn", "Day1 explicit gate: run R4(trans_2d) + R5(trans_3d)", []() { return RunDay1ArchitectureFreeze(); }},
        {"day1_arch_conn_repro", "Day1 reproducibility gate: trans_2d x2 + trans_3d x2", []() { return RunDay1ArchitectureFreezeRepro(); }},
        {"day2_fvm_ad", "Day2 FVM-AD operator acceptance tests (default)", []() { return RunDay2FvmAd(); }},
        {"day3_bc_patch", "Day3 explicit gate: Boundary condition operators & AD verification", []() { Test_FVM::Run_Day3_BC_Patch<3, ADVar<3>>(); return 0; }},
        {"day3_leakoff_switch", "Day3 explicit gate: Leakoff switch, AD ops & Grid-Level Assembly", []() { Test_FVM::Run_Day3_Leakoff_Switch<3, ADVar<3>>(); return 0; }},
        {"day3_viz", "Day3 explicit gate: Output VTK files for Boundary/Leakoff Assembly", []() { Test_FVM::Run_Day3_BC_Viz_2D<3, ADVar<3>>(); Test_FVM::Run_Day3_Leakoff_Viz_3D<3, ADVar<3>>(); return 0; }},
        {"day4_well_patch", "Day4 explicit gate: Well(BHP/Rate) + WAG schedule skeleton + matrix/fracture completion path", []() { Test_FVM::Run_Day4_Well_Patch<3, ADVar<3>>(); return 0; }},
        {"day4_well_viz", "Day4 visualization gate: export fixed-path 2D/3D well-source VTK files", []() { Test_FVM::Run_Day4_Well_Viz<3, ADVar<3>>(); return 0; }},
        {"day5_block_matrix_robust", "Day5 infra gate: FIM block sparse matrix robustness test", []() { Test_Day5::Test_FIM_BlockSparseMatrix_Robustness(); return 0; }},
        { "day5_global_jac_2d","Day5 core gate: 2D FD vs AD Jacobian Assembly Verification",[]() { Test_Day5::Run_Day5_GlobalAssembly_Jacobian_2D(); return 0; }},
        { "day5_global_jac_3d","Day5 core gate: 3D FD vs AD Jacobian Assembly Verification",[]() { Test_Day5::Run_Day5_GlobalAssembly_Jacobian_3D(); return 0; }},
        {"day6_transient_2d_sp_injprod", "Day6: 2D single-phase transient stability + VTK export", []() { Test_Day6::Run_Day6_Transient_2D_SP_InjProd(); return 0; }},
        {"day6_transient_2d_tp_injprod", "Day6: 2D two-phase transient stability + VTK export", []() { Test_Day6::Run_Day6_Transient_2D_TP_InjProd(); return 0; }},
        {"day6_transient_2d_tp_multiwell", "Day6: 2D multi-well interference transient + VTK export", []() { Test_Day6::Run_Day6_Transient_2D_TP_Multiwell(); return 0; }},
        {"day6_transient_3d_sp_injprod", "Day6: 3D single-phase transient stability + VTK export", []() { Test_Day6::Run_Day6_Transient_3D_SP_InjProd(); return 0; }},
        {"day6_transient_3d_tp_injprod", "Day6: 3D two-phase transient stability + VTK export", []() { Test_Day6::Run_Day6_Transient_3D_TP_InjProd(); return 0; }},
        {"day6_transient_3d_tp_multiwell", "Day6: 3D multi-well interference transient + VTK export", []() { Test_Day6::Run_Day6_Transient_3D_TP_Multiwell(); return 0; }},
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


/**
 * @file Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_NoWell.cpp
 * @brief Standalone validation chain: 2D single-phase CO2 constant-property P-T coupled, no-fracture, no-well.
 */

#include "Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_NoWell.h"

#include "2D_PostProcess.h"
#include "BoundaryConditionManager.h"
#include "Case2D_Matlab.h"
#include "Case2D_ReferenceIO.h"
#include "CaseCommon_Artifacts.h"
#include "CaseCommon_Catalog.h"
#include "FVM_Grad.h"
#include "FIM_TransientCaseKit.hpp"
#include "MeshDefinitions.h"
#include "MeshManager.h"
#include "SolverContrlStrName_op.h"
#include "Well_WellControlTypes.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#include <windows.h>
#define TEST_MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define TEST_MKDIR(path) mkdir(path, 0777)
#endif

namespace Test_H_T_CO2_ConstPP_NoFrac {
namespace {

constexpr double kPi = 3.1415926535897932384626433832795;
constexpr double kPendingTimeTolerance = 1.0e-9;
constexpr double kGridTimeMonotoneAbsTol = 5.0e-3;
constexpr const char* kFamilyMatrixHorizontal = "matrix_horizontal";
constexpr const char* kFamilyMatrixVerticalMidline = "matrix_vertical_midline";
constexpr const char* kLocationMatrix = "matrix";
constexpr const char* kComsolRepresentationNonisothermal = "matrix_builtin_nonisothermal_flow_in_porous_media";
constexpr const char* kComsolRepresentationDlHtManual = "matrix_builtin_dl_ht_manual_coupling";
constexpr const char* kComsolRepresentationPdeFallback = "matrix_only_coupled_pdes_fallback";

struct SpatialSamplePoint {
    int id = -1;
    std::string label;
    std::string family;
    std::string location;
    double target_axis_m = 0.0;
    double target_x = 0.0;
    double target_y = 0.0;
    int block_id = -1;
    double actual_x = 0.0;
    double actual_y = 0.0;
};

struct SnapshotState {
    std::string tag;
    double requested_fraction = 0.0;
    double target_time_s = 0.0;
    double actual_time_s = 0.0;
    bool captured = false;
    std::vector<double> p_blocks;
    std::vector<double> t_blocks;
};

struct MonitorScheduleState {
    int sample_id = -1;
    double target_time_s = 0.0;
    double actual_time_s = 0.0;
    bool captured = false;
    std::vector<double> p_blocks;
    std::vector<double> t_blocks;
};

struct FieldErrorMetrics {
    int sample_count = 0;
    double l1_abs = 0.0;
    double l2_abs = 0.0;
    double linf_abs = 0.0;
    double l1_norm = 0.0;
    double l2_norm = 0.0;
    double linf_norm = 0.0;
};

struct PressureCellMetrics {
    std::string tag;
    double requested_fraction = 0.0;
    double target_time_s = 0.0;
    double actual_time_s = 0.0;
    FieldErrorMetrics pressure;
    int max_err_cell = -1;
    std::string compare_csv_path;
};

struct ProfileCompareMetrics {
    std::string family;
    std::string tag;
    double target_time_s = 0.0;
    double actual_time_s = 0.0;
    FieldErrorMetrics pressure;
    FieldErrorMetrics temperature;
    double pressure_uniformity_norm_std = 0.0;
    double temperature_uniformity_norm_std = 0.0;
    std::string compare_csv_path;
};

struct ReconstructedFieldValue {
    double value = 0.0;
    int carrier_cell_id = -1;
    double carrier_x = 0.0;
    double carrier_y = 0.0;
    std::string sample_mode;
};

struct ReconstructedSpatialPointSample {
    int point_id = -1;
    int sample_id = -1;
    std::string time_tag;
    std::string label;
    std::string family;
    std::string location;
    double target_axis_m = 0.0;
    double target_x = 0.0;
    double target_y = 0.0;
    int block_id = -1;
    double actual_x = 0.0;
    double actual_y = 0.0;
    double target_time_s = 0.0;
    double actual_time_s = 0.0;
    ReconstructedFieldValue pressure;
    ReconstructedFieldValue temperature;
};

using ProfileReferenceRow = Case2DReferenceIO::ProfileReferenceRow;
using ProfileReferenceTable = Case2DReferenceIO::ProfileReferenceTable;
using MonitorReferenceRow = Case2DReferenceIO::MonitorReferenceRow;
using MonitorReferenceSeries = Case2DReferenceIO::MonitorReferenceSeries;

struct SweepStudyRow {
    std::string label;
    std::string case_dir;
    int nx = 0;
    int ny = 0;
    int steps = 0;
    double dt_init = 0.0;
    double h_char = 0.0;
    double t_end = 0.0;
    double pressure_cell_l2_norm = 0.0;
    double pressure_horizontal_l2_norm = 0.0;
    double temperature_horizontal_l2_norm = 0.0;
    double pressure_vertical_uniformity_norm_std = 0.0;
    double temperature_vertical_uniformity_norm_std = 0.0;
    double pressure_time_self_l2_norm = std::numeric_limits<double>::quiet_NaN();
    double temperature_time_self_l2_norm = std::numeric_limits<double>::quiet_NaN();
    double pressure_order = std::numeric_limits<double>::quiet_NaN();
    double temperature_order = std::numeric_limits<double>::quiet_NaN();
    std::vector<ReconstructedSpatialPointSample> profile_samples;
};

struct TestCaseSpec {
    std::string case_name = "h_t_co2_constpp_nofrac_nowell";
    std::string output_base_dir = "Test/Transient/FullCaseTest";
    std::string sub_dir = "H_T_CO2_ConstPP";
    std::string comsol_reference_case_name;
    bool use_variable_properties = false;
    double lx = 400.0;
    double ly = 40.0;
    int nx = 48;
    int ny = 6;
    double matrix_phi = 0.10;
    double matrix_perm = 1.0e-13;
    double matrix_ct = 5.0e-9;
    double matrix_rho_r = 2600.0;
    double matrix_cp_r = 800.0;
    double matrix_lambda_r = 3.0;
    double co2_rho_const = 700.0;
    double co2_mu_const = 6.0e-5;
    double co2_cp_const = 1100.0;
    double co2_cv_const = 850.0;
    double co2_k_const = 0.03;
    double p_init = 10.0e6;
    double p_left = 12.0e6;
    double p_right = 8.0e6;
    double t_init = 380.0;
    double t_left = 440.0;
    double t_right = 320.0;
    double dt_init = 2.0e4;
    double dt_min = 100.0;
    double dt_max = 2.0e5;
    double target_end_time_s = 2.0e7;
    int max_steps = 5000;
    int max_newton_iter = 14;
    FIM_Engine::LinearSolverType lin_solver = FIM_Engine::LinearSolverType::SparseLU;
    bool amgcl_use_fallback_sparselu = true;
    double amgcl_tol = 1.0e-6;
    int amgcl_maxiter = 500;
    bool enable_non_orthogonal_correction = true;
    bool enable_row_scaling = true;
    double abs_res_tol = 1.0e-8;
    double rel_res_tol = 1.0e-6;
    double rel_update_tol = 1.0e-8;
    double max_dP = 2.0e7;
    double max_dT = 25.0;
    bool enable_armijo_line_search = false;
    double rollback_shrink_factor = 0.7;
    double dt_relres_grow_factor = 1.15;
    Vector gravity_vector = Vector(0.0, 0.0, 0.0);
    FIM_Engine::DiagLevel diag_level = FIM_Engine::DiagLevel::Off;
    int analytical_terms = 200;
    std::vector<double> report_time_fractions = {0.1, 0.5, 1.0};
    int matrix_horizontal_station_count = 81;
    int matrix_vertical_station_count = 25;
    std::vector<std::pair<int, int> > grid_sweep_cases = {
        std::make_pair(24, 3),
        std::make_pair(48, 6),
        std::make_pair(96, 12)
    };
    std::vector<double> time_step_sweep = {8.0e4, 2.0e4, 5.0e3};
    double pressure_l2_threshold = 5.0e-2;
    double pressure_linf_threshold = 1.5e-1;
    double temperature_l2_threshold = 8.0e-2;
    double temperature_linf_threshold = 2.0e-1;
    double vertical_uniformity_threshold = 2.0e-2;
    double grid_pressure_order_threshold = 0.9;
    double grid_temperature_order_threshold = 0.7;
    double time_pressure_order_threshold = 0.7;
    double time_temperature_order_threshold = 0.5;
    bool enable_grid_convergence_study = false;
    bool enable_time_sensitivity_study = false;
    bool export_vtk = true;
    bool emit_detailed_outputs = true;
    bool allow_full_workflow_comsol_autorun = false;
    std::string comsol_wrapper_relpath =
        "tools/COMSOL/H_T_CO2_ConstPP_NoFrac_NoWell/run_comsol_reference.ps1";
};

struct TestCaseSummary {
    std::string case_dir;
    std::string studies_dir;
    std::string figures_dir;
    std::string engineering_dir;
    std::string reference_dir;
    std::string analytic_dir;
    std::string comsol_input_dir;
    std::string comsol_output_dir;
    std::string report_dir;
    std::string report_scripts_dir;
    std::string convergence_log_path;
    std::string run_log_path;
    std::string metrics_csv_path;
    std::string validation_summary_path;
    std::string validation_summary_csv_path;
    std::string reference_spec_path;
    std::string analytical_note_path;
    std::string comsol_reference_spec_path;
    std::string matlab_script_path;
    std::string property_table_path;
    std::string comsol_property_table_path;
    std::string profile_station_definitions_path;
    std::string monitor_point_definitions_path;
    std::string profile_schedule_path;
    std::string monitor_schedule_path;
    std::string eng_monitor_timeseries_path;
    std::string grid_convergence_csv_path;
    std::string time_sensitivity_csv_path;
    int nx = 0;
    int ny = 0;
    int n_cells = 0;
    double h_char = 0.0;
    int steps = 0;
    int total_rollbacks = 0;
    double avg_iters = 0.0;
    int max_iters = 0;
    double t_end = 0.0;
    std::string resolved_reference_mode = "pressure_analytical_plus_temperature_comsol";
    std::string comsol_representation = "";
    bool comsol_fine_check_skipped = false;
    bool reference_ready = false;
    bool validation_performed = false;
    bool validation_passed = false;
    bool grid_convergence_ok = true;
    bool time_sensitivity_ok = true;
    std::string validation_status = "not_run";
    double final_pressure_cell_l1_norm = 0.0;
    double final_pressure_cell_l2_norm = 0.0;
    double final_pressure_cell_linf_norm = 0.0;
    double final_pressure_horizontal_l2_norm = 0.0;
    double final_pressure_horizontal_linf_norm = 0.0;
    double final_temperature_horizontal_l2_norm = std::numeric_limits<double>::quiet_NaN();
    double final_temperature_horizontal_linf_norm = std::numeric_limits<double>::quiet_NaN();
    double final_monitor_pressure_l2_norm = 0.0;
    double final_monitor_pressure_linf_norm = 0.0;
    double final_monitor_temperature_l2_norm = std::numeric_limits<double>::quiet_NaN();
    double final_monitor_temperature_linf_norm = std::numeric_limits<double>::quiet_NaN();
    double final_pressure_vertical_uniformity_norm_std = 0.0;
    double final_temperature_vertical_uniformity_norm_std = std::numeric_limits<double>::quiet_NaN();
    double grid_pressure_order_min = std::numeric_limits<double>::quiet_NaN();
    double grid_temperature_order_min = std::numeric_limits<double>::quiet_NaN();
    double time_pressure_order_min = std::numeric_limits<double>::quiet_NaN();
    double time_temperature_order_min = std::numeric_limits<double>::quiet_NaN();
    std::vector<PressureCellMetrics> pressure_cell_metrics;
    std::vector<ProfileCompareMetrics> report_metrics;
    std::vector<std::string> missing_reference_files;
};

struct CaseRunArtifacts {
    TestCaseSummary summary;
    std::vector<SpatialSamplePoint> profile_stations;
    std::vector<SpatialSamplePoint> monitor_points;
    std::vector<SnapshotState> snapshots;
    std::vector<MonitorScheduleState> monitor_schedule;
    std::vector<ReconstructedSpatialPointSample> sampled_profile_points;
    std::vector<ReconstructedSpatialPointSample> sampled_monitor_points;
};

struct TestCasePlan {
    std::string plan_key;
    TestCaseSpec spec;
};

double PressureScale(const TestCaseSpec& cfg);
double TemperatureScale(const TestCaseSpec& cfg);
FieldErrorMetrics BuildErrorMetrics(double sumAbs, double sumSq, double maxAbs, int count, double scale);

std::string BoolString(bool value) { return value ? "true" : "false"; }

void EnsureDirRecursive(const std::string& rawPath) {
    if (rawPath.empty()) return;
    std::string path = rawPath;
    for (char& ch : path) {
        if (ch == '\\') ch = '/';
    }

    std::stringstream ss(path);
    std::string token;
    std::string current;
    while (std::getline(ss, token, '/')) {
        if (token.empty() || token == ".") continue;
        if (!current.empty()) current += "/";
        current += token;
        TEST_MKDIR(current.c_str());
    }
}

std::string GetNoFracTemplateCaseCode(const TestCaseSpec& cfg) {
    return cfg.use_variable_properties ? "B2" : "B1";
}

const CaseCommon::CaseCatalogEntry& GetNoFracCatalogEntryOrThrow(const TestCaseSpec& cfg) {
    const std::string caseCode = GetNoFracTemplateCaseCode(cfg);
    const CaseCommon::CaseCatalogEntry* entry = CaseCommon::FindCaseCatalogEntry(caseCode);
    if (!entry) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] missing " + caseCode + " catalog entry.");
    }
    return *entry;
}

void ConfigureSummaryPaths(TestCaseSummary& summary,
                           const TestCaseSpec& cfg,
                           const std::string& caseDir) {
    summary.case_dir = caseDir;
    summary.studies_dir = summary.case_dir + "/studies";
    summary.figures_dir = summary.case_dir + "/figures";
    summary.engineering_dir = summary.case_dir + "/engineering";
    summary.reference_dir = summary.case_dir + "/reference";
    summary.analytic_dir = summary.reference_dir + "/analytic";
    summary.comsol_input_dir = summary.reference_dir + "/comsol_input";
    summary.comsol_output_dir = summary.reference_dir + "/comsol";
    if (!cfg.comsol_reference_case_name.empty()) {
        summary.comsol_output_dir =
            cfg.output_base_dir + "/" + cfg.sub_dir + "/" + cfg.comsol_reference_case_name + "/reference/comsol";
    }
    summary.report_dir = summary.case_dir + "/report";
    summary.report_scripts_dir = summary.report_dir + "/scripts";
    summary.convergence_log_path = summary.case_dir + "/convergence.log";
    summary.run_log_path = summary.case_dir + "/run.log";
    summary.metrics_csv_path = summary.report_dir + "/metrics.csv";
    summary.validation_summary_path = summary.report_dir + "/validation_summary.md";
    summary.validation_summary_csv_path = summary.report_dir + "/validation_summary.csv";
    summary.reference_spec_path = summary.engineering_dir + "/reference_spec.md";
    summary.analytical_note_path = summary.analytic_dir + "/pressure_analytical_note.md";
    summary.comsol_reference_spec_path = summary.reference_dir + "/comsol_reference_spec.md";
    summary.matlab_script_path = summary.report_scripts_dir + "/plot_validation_results.m";
    summary.property_table_path = summary.engineering_dir + "/property_table.csv";
    summary.comsol_property_table_path = summary.comsol_input_dir + "/property_table.csv";
    summary.profile_station_definitions_path = summary.engineering_dir + "/profile_station_definitions.csv";
    summary.monitor_point_definitions_path = summary.engineering_dir + "/monitor_point_definitions.csv";
    summary.profile_schedule_path = summary.engineering_dir + "/profile_report_schedule.csv";
    summary.monitor_schedule_path = summary.engineering_dir + "/monitor_sample_schedule.csv";
    summary.eng_monitor_timeseries_path = summary.engineering_dir + "/eng_monitor_timeseries.csv";
    summary.grid_convergence_csv_path = summary.studies_dir + "/grid_convergence.csv";
    summary.time_sensitivity_csv_path = summary.studies_dir + "/time_sensitivity.csv";
    summary.nx = cfg.nx;
    summary.ny = cfg.ny;
}

void EnsureSummaryDirs(const TestCaseSummary& summary) {
    EnsureDirRecursive(summary.case_dir);
    EnsureDirRecursive(summary.studies_dir);
    EnsureDirRecursive(summary.figures_dir);
    EnsureDirRecursive(summary.engineering_dir);
    EnsureDirRecursive(summary.reference_dir);
    EnsureDirRecursive(summary.analytic_dir);
    EnsureDirRecursive(summary.comsol_input_dir);
    EnsureDirRecursive(summary.comsol_output_dir);
    EnsureDirRecursive(summary.report_dir);
    EnsureDirRecursive(summary.report_scripts_dir);
}

std::string JoinFractions(const std::vector<double>& values) {
    std::ostringstream oss;
    oss << std::setprecision(6);
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) oss << ", ";
        oss << values[i];
    }
    return oss.str();
}

CaseCommon::CaseArtifactPaths BuildB1ArtifactPaths(const TestCaseSpec& cfg) {
    const CaseCommon::CaseCatalogEntry& entry = GetNoFracCatalogEntryOrThrow(cfg);
    return CaseCommon::BuildArtifactPaths(
        entry.metadata.output_root,
        entry.metadata.case_code,
        entry.metadata.case_slug);
}

void EnsureB1ArtifactContractDirs(const CaseCommon::CaseArtifactPaths& artifacts) {
    EnsureDirRecursive(artifacts.root_dir);
    EnsureDirRecursive(artifacts.case_dir);
    EnsureDirRecursive(artifacts.studies_dir);
    EnsureDirRecursive(artifacts.figures_dir);
    EnsureDirRecursive(artifacts.engineering_dir);
    EnsureDirRecursive(artifacts.reference_dir);
    EnsureDirRecursive(artifacts.report_dir);
    EnsureDirRecursive(artifacts.report_scripts_dir);
}

void WriteB1StageManifest(const CaseCommon::CaseArtifactPaths& artifacts,
                          const std::string& caseCode,
                          CaseCommon::CaseStage stage,
                          const std::string& status,
                          const std::string& outputDir) {
    Case2DReferenceIO::WriteAsciiFile(
        artifacts.engineering_stage_manifest_path,
        ("[Test_H_T_CO2_ConstPP_NoFrac] failed to write " + caseCode + " stage manifest").c_str(),
        [&](std::ofstream& out) {
            out << "case_code=" << caseCode << "\n";
            out << "stage=" << CaseCommon::ToString(stage) << "\n";
            out << "status=" << status << "\n";
            out << "output_dir=" << outputDir << "\n";
            out << "case_dir=" << artifacts.case_dir << "\n";
            out << "studies_dir=" << artifacts.studies_dir << "\n";
            out << "engineering_dir=" << artifacts.engineering_dir << "\n";
            out << "reference_dir=" << artifacts.reference_dir << "\n";
            out << "report_dir=" << artifacts.report_dir << "\n";
            out << "report_scripts_dir=" << artifacts.report_scripts_dir << "\n";
        });
}

void WriteB1ReferenceContract(const TestCaseSpec& cfg,
                              const CaseCommon::CaseArtifactPaths& artifacts,
                              const std::string& caseCode,
                              CaseCommon::CaseStage stage) {
    Case2DReferenceIO::WriteAsciiFile(
        artifacts.reference_contract_path,
        ("[Test_H_T_CO2_ConstPP_NoFrac] failed to write " + caseCode + " reference contract").c_str(),
        [&](std::ofstream& out) {
            out << "# " << caseCode << " Reference Contract\n\n";
            out << "- Stage: `" << CaseCommon::ToString(stage) << "`\n";
            out << "- Case dir: `" << artifacts.case_dir << "`\n";
            if (cfg.use_variable_properties) {
                out << "- Reference mode: `COMSOL`\n";
                out << "- Fluid model: `single-phase CO2 EOS`\n";
            } else {
                out << "- Reference mode: `pressure analytical + temperature COMSOL`\n";
                out << "- Fluid model: `single-phase CO2 constant-property`\n";
            }
            out << "- Temperature COMSOL wrapper: `" << cfg.comsol_wrapper_relpath << "`\n\n";
            out << "## Required directories\n";
            out << "- `engineering/`\n";
            out << "- `reference/`\n";
            out << "- `studies/`\n";
            out << "- `report/`\n";
            out << "- `report/scripts/`\n\n";
            out << "## Notes\n";
            out << "- `prepare_reference` writes property tables and COMSOL/reference specs without running the transient solve.\n";
            out << "- Profile and monitor engineering CSVs still materialize during solve/validate until a mesh-only prep path exists.\n";
        });
}

void WriteB1StageStatus(const TestCaseSpec& cfg,
                        const CaseCommon::CaseArtifactPaths& artifacts,
                        const std::string& caseCode,
                        CaseCommon::CaseStage stage,
                        const std::string& status,
                        const std::string& outputDir) {
    Case2DReferenceIO::WriteAsciiFile(
        artifacts.report_status_markdown_path,
        ("[Test_H_T_CO2_ConstPP_NoFrac] failed to write " + caseCode + " stage status").c_str(),
        [&](std::ofstream& out) {
            out << "# " << caseCode << " Template Status\n\n";
            out << "- Stage: `" << CaseCommon::ToString(stage) << "`\n";
            out << "- Status: `" << status << "`\n";
            out << "- Output dir: `" << outputDir << "`\n";
            out << "- Case dir: `" << artifacts.case_dir << "`\n\n";
            out << "## Stage semantics\n";
            out << "- `solve_only`: run the coupled single-case solve without grid/time sweep studies.\n";
            out << "- `prepare_reference`: emit COMSOL/reference contract files and property tables only.\n";
            out << "- `validate_only`: rerun validation workflow until engineering snapshot persistence exists.\n";
            out << "- `full_workflow`: run the legacy B1 validation chain under the template-system case root.\n\n";
            out << "## Numerical settings\n";
            out << "- Grid: `" << cfg.nx << "x" << cfg.ny << "`\n";
            out << "- End time: `" << std::setprecision(12) << cfg.target_end_time_s << " s`\n";
            out << "- Report fractions: `" << JoinFractions(cfg.report_time_fractions) << "`\n";
            out << "- Property mode: `" << (cfg.use_variable_properties ? "varypp" : "constpp") << "`\n";
        });
    Case2DReferenceIO::WriteAsciiFile(
        artifacts.report_scripts_dir + "/README.m",
        ("[Test_H_T_CO2_ConstPP_NoFrac] failed to write " + caseCode + " Matlab stub").c_str(),
        [&](std::ofstream& out) {
            out << "% " << caseCode << " Matlab script placeholder\n";
            out << "% Stage: " << CaseCommon::ToString(stage) << "\n";
            out << "% Output root: " << artifacts.case_dir << "\n";
        });
}

TestCaseSpec BuildStageSpec(const TestCaseSpec& base, CaseCommon::CaseStage stage) {
    TestCaseSpec cfg = base;
    cfg.allow_full_workflow_comsol_autorun = false;
    switch (stage) {
    case CaseCommon::CaseStage::SolveOnly:
        cfg.enable_grid_convergence_study = false;
        cfg.enable_time_sensitivity_study = false;
        cfg.emit_detailed_outputs = false;
        break;
    case CaseCommon::CaseStage::ValidateOnly:
        cfg.enable_grid_convergence_study = false;
        cfg.enable_time_sensitivity_study = false;
        cfg.export_vtk = false;
        cfg.emit_detailed_outputs = true;
        break;
    case CaseCommon::CaseStage::PrepareReference:
        cfg.enable_grid_convergence_study = false;
        cfg.enable_time_sensitivity_study = false;
        cfg.export_vtk = false;
        cfg.emit_detailed_outputs = false;
        break;
    case CaseCommon::CaseStage::FullWorkflow:
        break;
    default:
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] unsupported stage in BuildStageSpec.");
    }
    return cfg;
}

void PrintRunSummary(const std::string& banner, const TestCaseSummary& summary) {
    std::cout << "\n============================================\n";
    std::cout << "[Test_H_T_CO2_ConstPP_NoFrac] " << banner << "\n";
    std::cout << "  output_dir: " << summary.case_dir << "\n";
    std::cout << "  engineering_dir: " << summary.engineering_dir << "\n";
    std::cout << "  grid: " << summary.nx << " x " << summary.ny
              << " (" << summary.n_cells << " cells)\n";
    std::cout << "  steps: " << summary.steps
              << "  rollbacks: " << summary.total_rollbacks << "\n";
    std::cout << "  Newton iters: avg=" << std::fixed << std::setprecision(2) << summary.avg_iters
              << "  max=" << summary.max_iters << "\n";
    std::cout << "  final_time: " << std::scientific << std::setprecision(4) << summary.t_end << " s\n";
    std::cout << "  validation_status: " << summary.validation_status << "\n";
    std::cout << "  validation_summary: " << summary.validation_summary_path << "\n";
    std::cout << "============================================\n";
}

void PrintPrepareReferenceSummary(const CaseCommon::CaseArtifactPaths& artifacts) {
    std::cout << "\n============================================\n";
    std::cout << "[Test_H_T_CO2_ConstPP_NoFrac] prepare_reference completed\n";
    std::cout << "  output_dir: " << artifacts.reference_dir << "\n";
    std::cout << "  stage_manifest: " << artifacts.engineering_stage_manifest_path << "\n";
    std::cout << "  reference_contract: " << artifacts.reference_contract_path << "\n";
    std::cout << "  status_markdown: " << artifacts.report_status_markdown_path << "\n";
    std::cout << "============================================\n";
}

std::string PathToGenericString(const std::filesystem::path& path) {
    return path.generic_string();
}

std::string BuildOutputOpenDiagnostics(const std::filesystem::path& targetPath) {
    std::error_code ec;
    std::ostringstream oss;
    oss << "target=" << PathToGenericString(targetPath);
    if (!targetPath.parent_path().empty()) {
        oss << ", parent=" << PathToGenericString(targetPath.parent_path())
            << ", parent_exists=" << BoolString(std::filesystem::exists(targetPath.parent_path(), ec));
        ec.clear();
    }
    oss << ", target_exists=" << BoolString(std::filesystem::exists(targetPath, ec));
    ec.clear();
    const std::filesystem::path cwd = std::filesystem::current_path(ec);
    if (!ec) oss << ", cwd=" << PathToGenericString(cwd);
    return oss.str();
}

std::filesystem::path BuildShortAsciiStagingPath(const std::string& targetPath) {
    std::error_code ec;
    std::filesystem::path stagingRoot = std::filesystem::temp_directory_path(ec);
    if (ec || stagingRoot.empty()) stagingRoot = std::filesystem::path(".");
    stagingRoot /= "codex_nofrac_io";
    std::filesystem::create_directories(stagingRoot, ec);

    const std::size_t hashValue = std::hash<std::string>{}(targetPath);
    std::ostringstream oss;
    oss << std::hex << hashValue << "_" << std::filesystem::path(targetPath).filename().generic_string() << ".tmp";
    return stagingRoot / oss.str();
}

std::ofstream OpenAsciiStagingStream(const std::string& targetPath,
                                     const char* context,
                                     std::filesystem::path& stagingPath) {
    stagingPath = BuildShortAsciiStagingPath(targetPath);
    std::ofstream out;
    for (int attempt = 0; attempt < 8; ++attempt) {
        std::error_code ec;
        std::filesystem::remove(stagingPath, ec);
        out = std::ofstream(stagingPath, std::ios::out | std::ios::trunc);
        if (out.good()) return out;
        out.clear();
        std::this_thread::sleep_for(std::chrono::milliseconds(40));
    }
    throw std::runtime_error(
        std::string(context) + ": failed to open ASCII staging file for " + targetPath +
        " | staging=" + PathToGenericString(stagingPath) +
        " | " + BuildOutputOpenDiagnostics(std::filesystem::path(targetPath)));
}

void CommitAsciiStagingFile(std::ofstream& out,
                            const std::filesystem::path& stagingPath,
                            const std::string& targetPath,
                            const char* context) {
    out.flush();
    if (!out.good()) {
        out.close();
        std::error_code ec;
        std::filesystem::remove(stagingPath, ec);
        throw std::runtime_error(
            std::string(context) + ": failed while flushing staging file for " + targetPath +
            " | staging=" + PathToGenericString(stagingPath));
    }
    out.close();

    std::error_code ec;
    std::filesystem::path target(targetPath);
    std::filesystem::path absoluteTarget = std::filesystem::absolute(target, ec);
    if (ec) {
        ec.clear();
        absoluteTarget = std::filesystem::current_path(ec) / target;
    }
    if (!absoluteTarget.parent_path().empty()) {
        std::filesystem::create_directories(absoluteTarget.parent_path(), ec);
        ec.clear();
    }

    bool copied = false;
    for (int attempt = 0; attempt < 8; ++attempt) {
        std::error_code removeEc;
        if (std::filesystem::exists(absoluteTarget, removeEc)) {
            std::filesystem::remove(absoluteTarget, removeEc);
        }
        std::error_code copyEc;
        std::filesystem::copy_file(stagingPath, absoluteTarget, std::filesystem::copy_options::none, copyEc);
        if (!copyEc) {
            copied = true;
            break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(40));
    }
#ifdef _WIN32
    if (!copied) {
        const std::wstring stagingW = stagingPath.native();
        const std::wstring targetW = absoluteTarget.native();
        for (int attempt = 0; attempt < 8; ++attempt) {
            DeleteFileW(targetW.c_str());
            if (CopyFileW(stagingW.c_str(), targetW.c_str(), FALSE) != 0) {
                copied = true;
                break;
            }
            Sleep(40);
        }
    }
#endif
    if (!copied) {
        std::error_code cleanupEc;
        std::filesystem::remove(stagingPath, cleanupEc);
        throw std::runtime_error(
            std::string(context) + ": failed to commit staging file to " + targetPath +
            " | staging=" + PathToGenericString(stagingPath) +
            " | abs_target=" + PathToGenericString(absoluteTarget) +
            " | " + BuildOutputOpenDiagnostics(absoluteTarget));
    }
    std::filesystem::remove(stagingPath, ec);
}

void ApplyUniformScalarField(const std::shared_ptr<volScalarField>& field, double value) {
    if (!field) return;
    for (double& v : field->data) v = value;
}

void SyncPTFieldsToFM(MeshManager& mgr,
                      FieldManager_2D& fm,
                      const std::vector<double>& pBlocks,
                      const std::vector<double>& tBlocks,
                      double pFallback,
                      double tFallback) {
    if (pBlocks.empty()) return;

    const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
    const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();

    auto fPw = fm.getOrCreateMatrixScalar(pCfg.pressure_field, pFallback);
    auto fPviz = fm.getOrCreateMatrixScalar("P", pFallback);
    auto fTw = fm.getOrCreateMatrixScalar(tCfg.temperatue_field, tFallback);
    auto fTviz = fm.getOrCreateMatrixScalar("T", tFallback);

    auto fracPw = fm.getOrCreateFractureScalar(pCfg.pressure_field, pFallback);
    auto fracPviz = fm.getOrCreateFractureScalar("P", pFallback);
    auto fracTw = fm.getOrCreateFractureScalar(tCfg.temperatue_field, tFallback);
    auto fracTviz = fm.getOrCreateFractureScalar("T", tFallback);

    const int nMat = mgr.getMatrixDOFCount();
    const int nUse = std::min(static_cast<int>(pBlocks.size()), mgr.getTotalDOFCount());
    for (int i = 0; i < nUse; ++i) {
        const double p = pBlocks[static_cast<std::size_t>(i)];
        const double t = (i < static_cast<int>(tBlocks.size()))
            ? tBlocks[static_cast<std::size_t>(i)]
            : tFallback;
        if (i < nMat) {
            if (fPw && i < static_cast<int>(fPw->data.size())) fPw->data[static_cast<std::size_t>(i)] = p;
            if (fPviz && i < static_cast<int>(fPviz->data.size())) fPviz->data[static_cast<std::size_t>(i)] = p;
            if (fTw && i < static_cast<int>(fTw->data.size())) fTw->data[static_cast<std::size_t>(i)] = t;
            if (fTviz && i < static_cast<int>(fTviz->data.size())) fTviz->data[static_cast<std::size_t>(i)] = t;
        } else {
            const int fi = i - nMat;
            if (fracPw && fi < static_cast<int>(fracPw->data.size())) fracPw->data[static_cast<std::size_t>(fi)] = p;
            if (fracPviz && fi < static_cast<int>(fracPviz->data.size())) fracPviz->data[static_cast<std::size_t>(fi)] = p;
            if (fracTw && fi < static_cast<int>(fracTw->data.size())) fracTw->data[static_cast<std::size_t>(fi)] = t;
            if (fracTviz && fi < static_cast<int>(fracTviz->data.size())) fracTviz->data[static_cast<std::size_t>(fi)] = t;
        }
    }
}

BoundarySetting::BoundaryConditionManager BuildPressureBoundaryManager(const TestCaseSpec& cfg) {
    BoundarySetting::BoundaryConditionManager bcP;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT, cfg.p_left);
    bcP.SetDirichletBC(MeshTags::RIGHT, cfg.p_right);
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcP.SetNeumannBC(MeshTags::TOP, 0.0);
    return bcP;
}

BoundarySetting::BoundaryConditionManager BuildTemperatureBoundaryManager(const TestCaseSpec& cfg) {
    BoundarySetting::BoundaryConditionManager bcT;
    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT, cfg.t_left);
    bcT.SetDirichletBC(MeshTags::RIGHT, cfg.t_right);
    bcT.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcT.SetNeumannBC(MeshTags::TOP, 0.0);
    return bcT;
}

std::vector<Vector> BuildSortedCellPolygon2D(const Cell& cell, const std::unordered_map<int, Node>& nodes) {
    std::vector<Vector> polygon;
    polygon.reserve(cell.CellNodeIDs.size());
    for (int nodeId : cell.CellNodeIDs) {
        const auto it = nodes.find(nodeId);
        if (it == nodes.end()) {
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] missing node while building cell polygon.");
        }
        polygon.push_back(it->second.coord);
    }
    std::sort(polygon.begin(), polygon.end(), [&](const Vector& a, const Vector& b) {
        const double angleA = std::atan2(a.m_y - cell.center.m_y, a.m_x - cell.center.m_x);
        const double angleB = std::atan2(b.m_y - cell.center.m_y, b.m_x - cell.center.m_x);
        return angleA < angleB;
    });
    return polygon;
}

bool PointOnSegment2D(double px, double py, const Vector& a, const Vector& b, double tol) {
    const double cross = (px - a.m_x) * (b.m_y - a.m_y) - (py - a.m_y) * (b.m_x - a.m_x);
    if (std::abs(cross) > tol) return false;
    const double dot = (px - a.m_x) * (b.m_x - a.m_x) + (py - a.m_y) * (b.m_y - a.m_y);
    if (dot < -tol) return false;
    const double lenSq = (b.m_x - a.m_x) * (b.m_x - a.m_x) + (b.m_y - a.m_y) * (b.m_y - a.m_y);
    return dot <= lenSq + tol;
}

bool PointInPolygon2D(double px, double py, const std::vector<Vector>& polygon, double tol) {
    if (polygon.size() < 3) return false;
    bool inside = false;
    for (std::size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++) {
        const Vector& vi = polygon[i];
        const Vector& vj = polygon[j];
        if (PointOnSegment2D(px, py, vj, vi, tol)) return true;
        const bool intersects = ((vi.m_y > py) != (vj.m_y > py)) &&
            (px <= (vj.m_x - vi.m_x) * (py - vi.m_y) / ((vj.m_y - vi.m_y) + 1.0e-30) + vi.m_x + tol);
        if (intersects) inside = !inside;
    }
    return inside;
}

int FindNearestMatrixCell(const MeshManager& mgr, double targetX, double targetY) {
    const auto& cells = mgr.mesh().getCells();
    if (cells.empty()) return -1;
    int bestCell = -1;
    double bestDistSq = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < cells.size(); ++i) {
        const double dx = cells[i].center.m_x - targetX;
        const double dy = cells[i].center.m_y - targetY;
        const double distSq = dx * dx + dy * dy;
        if (distSq < bestDistSq) {
            bestDistSq = distSq;
            bestCell = static_cast<int>(i);
        }
    }
    return bestCell;
}

int FindContainingMatrixCell(const MeshManager& mgr, double targetX, double targetY) {
    const auto& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& nodes = mesh.getNodesMap();
    const double tol = 1.0e-8;
    int bestCell = -1;
    double bestDistSq = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < cells.size(); ++i) {
        const Cell& cell = cells[i];
        if (targetX < cell.boundingBox.min.m_x - tol || targetX > cell.boundingBox.max.m_x + tol ||
            targetY < cell.boundingBox.min.m_y - tol || targetY > cell.boundingBox.max.m_y + tol) {
            continue;
        }
        if (!PointInPolygon2D(targetX, targetY, BuildSortedCellPolygon2D(cell, nodes), tol)) continue;
        const double dx = cell.center.m_x - targetX;
        const double dy = cell.center.m_y - targetY;
        const double distSq = dx * dx + dy * dy;
        if (distSq < bestDistSq) {
            bestDistSq = distSq;
            bestCell = static_cast<int>(i);
        }
    }
    return bestCell;
}

std::shared_ptr<volVectorField> BuildLeastSquaresGradient(const MeshManager& mgr,
                                                          const std::vector<double>& blocks,
                                                          const std::string& fieldName,
                                                          double fallbackValue,
                                                          const BoundarySetting::BoundaryConditionManager& bc) {
    const std::size_t nCells = mgr.mesh().getCells().size();
    auto field = std::make_shared<volScalarField>(fieldName, nCells, fallbackValue);
    const std::size_t nUse = std::min(nCells, blocks.size());
    for (std::size_t i = 0; i < nUse; ++i) {
        field->data[i] = blocks[i];
    }
    FVM_Grad grad(mgr.mesh(), nullptr, nullptr, &bc);
    grad.precomputeLS();
    return grad.compute(*field, FVM_Grad::Method::LeastSquares);
}

ReconstructedFieldValue SampleFieldAtTarget(const MeshManager& mgr,
                                            const SpatialSamplePoint& point,
                                            const std::vector<double>& blocks,
                                            const std::shared_ptr<volVectorField>& gradient) {
    const auto& cells = mgr.mesh().getCells();
    if (cells.empty()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] cannot sample empty mesh.");

    const int containingCell = FindContainingMatrixCell(mgr, point.target_x, point.target_y);
    int carrierCell = containingCell;
    bool usedContaining = containingCell >= 0;
    if (carrierCell < 0) {
        if (point.block_id >= 0 && point.block_id < static_cast<int>(cells.size())) {
            carrierCell = point.block_id;
        } else {
            carrierCell = FindNearestMatrixCell(mgr, point.target_x, point.target_y);
        }
    }
    if (carrierCell < 0 || carrierCell >= static_cast<int>(cells.size()) ||
        carrierCell >= static_cast<int>(blocks.size())) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] invalid carrier cell during fixed-target sampling.");
    }

    ReconstructedFieldValue out;
    out.carrier_cell_id = carrierCell;
    out.carrier_x = cells[static_cast<std::size_t>(carrierCell)].center.m_x;
    out.carrier_y = cells[static_cast<std::size_t>(carrierCell)].center.m_y;
    out.value = blocks[static_cast<std::size_t>(carrierCell)];
    out.sample_mode = usedContaining ? "containing_cell_center" : "nearest_cell_center";

    if (gradient && carrierCell < static_cast<int>(gradient->data.size())) {
        const Vector delta(
            point.target_x - out.carrier_x,
            point.target_y - out.carrier_y,
            0.0);
        out.value += gradient->data[static_cast<std::size_t>(carrierCell)] * delta;
        out.sample_mode = usedContaining ? "containing_ls" : "nearest_ls";
    }
    return out;
}

bool IsLeftDirichletTarget(const SpatialSamplePoint& point, const TestCaseSpec& cfg) {
    return std::abs(point.target_x) <= 1.0e-8 && point.target_y >= -1.0e-8 && point.target_y <= cfg.ly + 1.0e-8;
}

bool IsRightDirichletTarget(const SpatialSamplePoint& point, const TestCaseSpec& cfg) {
    return std::abs(point.target_x - cfg.lx) <= 1.0e-8 && point.target_y >= -1.0e-8 && point.target_y <= cfg.ly + 1.0e-8;
}

std::vector<ReconstructedSpatialPointSample> BuildReconstructedSamples(const MeshManager& mgr,
                                                                       const std::vector<SpatialSamplePoint>& points,
                                                                       const std::vector<double>& pBlocks,
                                                                       const std::vector<double>& tBlocks,
                                                                       const std::string& timeTag,
                                                                       int sampleId,
                                                                       double targetTimeS,
                                                                       double actualTimeS,
                                                                       const TestCaseSpec& cfg) {
    std::vector<ReconstructedSpatialPointSample> rows;
    rows.reserve(points.size());
    const auto pressureBC = BuildPressureBoundaryManager(cfg);
    const auto temperatureBC = BuildTemperatureBoundaryManager(cfg);
    const std::shared_ptr<volVectorField> pGrad =
        BuildLeastSquaresGradient(mgr, pBlocks, "p_fixed_target_sample", cfg.p_init, pressureBC);
    const std::shared_ptr<volVectorField> tGrad =
        BuildLeastSquaresGradient(mgr, tBlocks, "t_fixed_target_sample", cfg.t_init, temperatureBC);

    for (const auto& point : points) {
        ReconstructedSpatialPointSample row;
        row.point_id = point.id;
        row.sample_id = sampleId;
        row.time_tag = timeTag;
        row.label = point.label;
        row.family = point.family;
        row.location = point.location;
        row.target_axis_m = point.target_axis_m;
        row.target_x = point.target_x;
        row.target_y = point.target_y;
        row.block_id = point.block_id;
        row.actual_x = point.actual_x;
        row.actual_y = point.actual_y;
        row.target_time_s = targetTimeS;
        row.actual_time_s = actualTimeS;
        row.pressure = SampleFieldAtTarget(mgr, point, pBlocks, pGrad);
        row.temperature = SampleFieldAtTarget(mgr, point, tBlocks, tGrad);
        if (IsLeftDirichletTarget(point, cfg)) {
            row.pressure.value = cfg.p_left;
            row.pressure.sample_mode = "boundary_dirichlet_left";
            row.temperature.value = cfg.t_left;
            row.temperature.sample_mode = "boundary_dirichlet_left";
        } else if (IsRightDirichletTarget(point, cfg)) {
            row.pressure.value = cfg.p_right;
            row.pressure.sample_mode = "boundary_dirichlet_right";
            row.temperature.value = cfg.t_right;
            row.temperature.sample_mode = "boundary_dirichlet_right";
        }
        rows.push_back(row);
    }
    return rows;
}

std::vector<ReconstructedSpatialPointSample> FilterReconstructedSamples(const std::vector<ReconstructedSpatialPointSample>& rows,
                                                                        const std::string& family,
                                                                        const std::string& timeTag) {
    std::vector<ReconstructedSpatialPointSample> out;
    for (const auto& row : rows) {
        if (row.family == family && row.time_tag == timeTag) out.push_back(row);
    }
    return out;
}

bool IsInteriorHorizontalTarget(double targetX, const TestCaseSpec& cfg) {
    const double tol = 1.0e-8;
    const double profileDx =
        (cfg.matrix_horizontal_station_count > 1)
            ? (cfg.lx / static_cast<double>(cfg.matrix_horizontal_station_count - 1))
            : 0.0;
    const double margin = std::max(2.0 * profileDx, tol);
    return targetX > margin && targetX < (cfg.lx - margin);
}

std::vector<ReconstructedSpatialPointSample> FilterInteriorHorizontalSamples(
    const std::vector<ReconstructedSpatialPointSample>& rows,
    const TestCaseSpec& cfg) {
    std::vector<ReconstructedSpatialPointSample> out;
    out.reserve(rows.size());
    for (const auto& row : rows) {
        if (IsInteriorHorizontalTarget(row.target_x, cfg)) out.push_back(row);
    }
    return out;
}

ProfileReferenceTable FilterInteriorHorizontalReference(const ProfileReferenceTable& table, const TestCaseSpec& cfg) {
    ProfileReferenceTable out;
    for (const auto& entry : table.rows_by_station_id) {
        if (IsInteriorHorizontalTarget(entry.second.target_x, cfg)) {
            out.rows_by_station_id.emplace(entry.first, entry.second);
        }
    }
    return out;
}

std::string MakeReconstructedSampleKey(const ReconstructedSpatialPointSample& sample) {
    return sample.time_tag + "#" + std::to_string(sample.point_id);
}

std::pair<FieldErrorMetrics, FieldErrorMetrics> EvaluateSelfConvergenceMetrics(
    const std::vector<ReconstructedSpatialPointSample>& coarseSamples,
    const std::vector<ReconstructedSpatialPointSample>& fineSamples,
    const TestCaseSpec& cfg) {
    if (coarseSamples.empty() || fineSamples.empty()) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] empty profile samples during self-convergence comparison.");
    }
    std::unordered_map<std::string, ReconstructedSpatialPointSample> fineByKey;
    for (const auto& row : fineSamples) {
        fineByKey[MakeReconstructedSampleKey(row)] = row;
    }

    double pSumAbs = 0.0;
    double pSumSq = 0.0;
    double pMaxAbs = 0.0;
    double tSumAbs = 0.0;
    double tSumSq = 0.0;
    double tMaxAbs = 0.0;
    int count = 0;

    for (const auto& coarse : coarseSamples) {
        const auto it = fineByKey.find(MakeReconstructedSampleKey(coarse));
        if (it == fineByKey.end()) {
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] missing fine-grid/fine-dt sample during self-convergence comparison.");
        }
        const double pAbsErr = std::abs(coarse.pressure.value - it->second.pressure.value);
        const double tAbsErr = std::abs(coarse.temperature.value - it->second.temperature.value);
        pSumAbs += pAbsErr;
        pSumSq += pAbsErr * pAbsErr;
        pMaxAbs = std::max(pMaxAbs, pAbsErr);
        tSumAbs += tAbsErr;
        tSumSq += tAbsErr * tAbsErr;
        tMaxAbs = std::max(tMaxAbs, tAbsErr);
        ++count;
    }

    return std::make_pair(
        BuildErrorMetrics(pSumAbs, pSumSq, pMaxAbs, count, PressureScale(cfg)),
        BuildErrorMetrics(tSumAbs, tSumSq, tMaxAbs, count, TemperatureScale(cfg)));
}

double ComputeMeshCharLength(const MeshManager& mgr) {
    const auto& cells = mgr.mesh().getCells();
    if (cells.empty()) return 0.0;
    double totalV = 0.0;
    for (const auto& c : cells) totalV += std::max(c.volume, 0.0);
    if (totalV <= 0.0) return 0.0;
    return std::sqrt(totalV / static_cast<double>(cells.size()));
}

double ClampFraction(double value) { return std::max(0.0, std::min(1.0, value)); }

std::string MakeReportTag(double fraction) {
    std::ostringstream oss;
    oss << "t" << std::setw(3) << std::setfill('0')
        << static_cast<int>(std::round(ClampFraction(fraction) * 100.0)) << "pct";
    return oss.str();
}

std::vector<double> BuildSortedFractions(const std::vector<double>& raw) {
    std::vector<double> out;
    for (double value : raw) {
        const double clamped = ClampFraction(value);
        bool duplicate = false;
        for (double existing : out) {
            if (std::abs(existing - clamped) <= 1.0e-12) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) out.push_back(clamped);
    }
    std::sort(out.begin(), out.end());
    return out;
}

double PressureScale(const TestCaseSpec& cfg) {
    return std::max(std::abs(cfg.p_left - cfg.p_right), 1.0);
}

double TemperatureScale(const TestCaseSpec& cfg) {
    return std::max(std::abs(cfg.t_left - cfg.t_right), 1.0);
}

double ComputePressureDiffusivity(const TestCaseSpec& cfg) {
    const double denom = cfg.co2_mu_const * cfg.matrix_phi * cfg.matrix_ct;
    if (cfg.matrix_perm <= 0.0 || denom <= 0.0) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] invalid constant-property coefficients for analytical diffusivity.");
    }
    return cfg.matrix_perm / denom;
}

double EvaluateAnalyticalPressure(const TestCaseSpec& cfg, double x, double timeS) {
    if (cfg.lx <= 0.0) throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] invalid domain length.");
    if (timeS <= 0.0) return cfg.p_init;
    const double xClamped = std::max(0.0, std::min(cfg.lx, x));
    const double deltaP = cfg.p_left - cfg.p_right;
    const double pSteady = cfg.p_left - deltaP * (xClamped / cfg.lx);
    const double a0 = cfg.p_init - cfg.p_left;
    const double D = ComputePressureDiffusivity(cfg);
    double transient = 0.0;
    for (int n = 1; n <= cfg.analytical_terms; ++n) {
        const double sign_n = (n % 2 == 0) ? 1.0 : -1.0;
        const double bn = (2.0 / (static_cast<double>(n) * kPi))
            * (a0 * (1.0 - sign_n) - deltaP * sign_n);
        const double lambda = static_cast<double>(n) * kPi / cfg.lx;
        transient += bn * std::sin(lambda * xClamped) * std::exp(-D * lambda * lambda * timeS);
    }
    return pSteady + transient;
}

std::string MakeLabel(const std::string& prefix, int index, int width) {
    std::ostringstream oss;
    oss << prefix << std::setw(width) << std::setfill('0') << index;
    return oss.str();
}

SpatialSamplePoint MakeNearestMatrixSample(const MeshManager& mgr,
                                          int id,
                                          const std::string& label,
                                          const std::string& family,
                                          double targetAxisM,
                                          double targetX,
                                          double targetY) {
    SpatialSamplePoint sample;
    sample.id = id;
    sample.label = label;
    sample.family = family;
    sample.location = kLocationMatrix;
    sample.target_axis_m = targetAxisM;
    sample.target_x = targetX;
    sample.target_y = targetY;
    const auto& cells = mgr.mesh().getCells();
    if (cells.empty()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] cannot build sample on empty mesh.");
    sample.block_id = FindContainingMatrixCell(mgr, targetX, targetY);
    if (sample.block_id < 0) {
        sample.block_id = FindNearestMatrixCell(mgr, targetX, targetY);
    }
    if (sample.block_id >= 0) {
        sample.actual_x = cells[static_cast<std::size_t>(sample.block_id)].center.m_x;
        sample.actual_y = cells[static_cast<std::size_t>(sample.block_id)].center.m_y;
    }
    return sample;
}

std::vector<SpatialSamplePoint> BuildProfileStations(const MeshManager& mgr, const TestCaseSpec& cfg) {
    if (cfg.matrix_horizontal_station_count < 2 || cfg.matrix_vertical_station_count < 2) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] invalid profile station counts.");
    }
    std::vector<SpatialSamplePoint> stations;
    int nextId = 0;
    for (int i = 0; i < cfg.matrix_horizontal_station_count; ++i) {
        const double ratio = static_cast<double>(i) / static_cast<double>(cfg.matrix_horizontal_station_count - 1);
        const double x = ratio * cfg.lx;
        const double y = 0.5 * cfg.ly;
        stations.push_back(MakeNearestMatrixSample(
            mgr, nextId++, MakeLabel("mh", i, 3), kFamilyMatrixHorizontal, x, x, y));
    }
    for (int i = 0; i < cfg.matrix_vertical_station_count; ++i) {
        const double ratio = static_cast<double>(i) / static_cast<double>(cfg.matrix_vertical_station_count - 1);
        const double x = 0.5 * cfg.lx;
        const double y = ratio * cfg.ly;
        stations.push_back(MakeNearestMatrixSample(
            mgr, nextId++, MakeLabel("mv", i, 3), kFamilyMatrixVerticalMidline, y, x, y));
    }
    return stations;
}

std::vector<SpatialSamplePoint> BuildMonitorPoints(const MeshManager& mgr, const TestCaseSpec& cfg) {
    const double xFractions[] = {0.1, 0.25, 0.5, 0.75, 0.9};
    std::vector<SpatialSamplePoint> points;
    int nextId = 0;
    for (std::size_t i = 0; i < sizeof(xFractions) / sizeof(xFractions[0]); ++i) {
        const double x = xFractions[i] * cfg.lx;
        const double y = 0.5 * cfg.ly;
        points.push_back(MakeNearestMatrixSample(
            mgr, nextId++, MakeLabel("mx", static_cast<int>(i + 1), 2), "monitor", x, x, y));
    }
    return points;
}

std::vector<SnapshotState> BuildSnapshots(const TestCaseSpec& cfg) {
    std::vector<SnapshotState> snapshots;
    const std::vector<double> fractions = BuildSortedFractions(cfg.report_time_fractions);
    snapshots.reserve(fractions.size());
    for (double fraction : fractions) {
        SnapshotState state;
        state.tag = MakeReportTag(fraction);
        state.requested_fraction = fraction;
        state.target_time_s = fraction * cfg.target_end_time_s;
        snapshots.push_back(state);
    }
    return snapshots;
}

std::vector<MonitorScheduleState> BuildMonitorSchedule(const TestCaseSpec& cfg) {
    std::vector<MonitorScheduleState> schedule;
    int sampleId = 0;
    schedule.push_back(MonitorScheduleState());
    schedule.back().sample_id = sampleId++;
    schedule.back().target_time_s = 0.0;
    for (double fraction : BuildSortedFractions(cfg.report_time_fractions)) {
        MonitorScheduleState state;
        state.sample_id = sampleId++;
        state.target_time_s = fraction * cfg.target_end_time_s;
        schedule.push_back(state);
    }
    return schedule;
}

void CaptureSnapshotIfCrossed(double prevTime,
                              const std::vector<double>& prevPBlocks,
                              const std::vector<double>& prevTBlocks,
                              double currentTime,
                              const std::vector<double>& currentPBlocks,
                              const std::vector<double>& currentTBlocks,
                              SnapshotState& snapshot) {
    if (snapshot.captured || currentTime + kPendingTimeTolerance < snapshot.target_time_s) return;
    const bool usePrev = !prevPBlocks.empty() && !prevTBlocks.empty() && prevTime >= 0.0
        && std::abs(prevTime - snapshot.target_time_s) < std::abs(currentTime - snapshot.target_time_s);
    snapshot.captured = true;
    snapshot.actual_time_s = usePrev ? prevTime : currentTime;
    snapshot.p_blocks = usePrev ? prevPBlocks : currentPBlocks;
    snapshot.t_blocks = usePrev ? prevTBlocks : currentTBlocks;
}

void CaptureMonitorIfCrossed(double prevTime,
                             const std::vector<double>& prevPBlocks,
                             const std::vector<double>& prevTBlocks,
                             double currentTime,
                             const std::vector<double>& currentPBlocks,
                             const std::vector<double>& currentTBlocks,
                             MonitorScheduleState& sample) {
    if (sample.captured || currentTime + kPendingTimeTolerance < sample.target_time_s) return;
    const bool usePrev = !prevPBlocks.empty() && !prevTBlocks.empty() && prevTime >= 0.0
        && std::abs(prevTime - sample.target_time_s) < std::abs(currentTime - sample.target_time_s);
    sample.captured = true;
    sample.actual_time_s = usePrev ? prevTime : currentTime;
    sample.p_blocks = usePrev ? prevPBlocks : currentPBlocks;
    sample.t_blocks = usePrev ? prevTBlocks : currentTBlocks;
}

void FinalizeMissingSnapshots(double finalTime,
                              const std::vector<double>& finalPBlocks,
                              const std::vector<double>& finalTBlocks,
                              std::vector<SnapshotState>& snapshots) {
    for (auto& snapshot : snapshots) {
        if (!snapshot.captured) {
            snapshot.captured = true;
            snapshot.actual_time_s = finalTime;
            snapshot.p_blocks = finalPBlocks;
            snapshot.t_blocks = finalTBlocks;
        }
    }
}

void FinalizeMissingMonitorSamples(double finalTime,
                                   const std::vector<double>& finalPBlocks,
                                   const std::vector<double>& finalTBlocks,
                                   std::vector<MonitorScheduleState>& schedule) {
    for (auto& sample : schedule) {
        if (!sample.captured) {
            sample.captured = true;
            sample.actual_time_s = finalTime;
            sample.p_blocks = finalPBlocks;
            sample.t_blocks = finalTBlocks;
        }
    }
}

std::vector<std::string> OrderedProfileFamilies() {
    return {kFamilyMatrixHorizontal, kFamilyMatrixVerticalMidline};
}

std::pair<int, int> GetFinestGridCase(const TestCaseSpec& cfg) {
    if (cfg.grid_sweep_cases.empty()) {
        return std::make_pair(cfg.nx, cfg.ny);
    }
    return *std::max_element(
        cfg.grid_sweep_cases.begin(),
        cfg.grid_sweep_cases.end(),
        [](const std::pair<int, int>& lhs, const std::pair<int, int>& rhs) {
            const int lhsCells = lhs.first * lhs.second;
            const int rhsCells = rhs.first * rhs.second;
            if (lhsCells != rhsCells) return lhsCells < rhsCells;
            return lhs.first < rhs.first;
        });
}

std::vector<SpatialSamplePoint> FilterStationsByFamily(const std::vector<SpatialSamplePoint>& stations,
                                                       const std::string& family) {
    std::vector<SpatialSamplePoint> out;
    for (const auto& station : stations) {
        if (station.family == family) out.push_back(station);
    }
    return out;
}

std::string Trim(std::string value) {
    const char* whitespace = " \t\r\n";
    const std::size_t begin = value.find_first_not_of(whitespace);
    if (begin == std::string::npos) return "";
    const std::size_t end = value.find_last_not_of(whitespace);
    return value.substr(begin, end - begin + 1);
}

FieldErrorMetrics BuildErrorMetrics(double sumAbs, double sumSq, double maxAbs, int count, double scale) {
    FieldErrorMetrics metrics;
    metrics.sample_count = count;
    if (count <= 0) return metrics;
    const double safeScale = std::max(scale, 1.0);
    metrics.l1_abs = sumAbs / static_cast<double>(count);
    metrics.l2_abs = std::sqrt(sumSq / static_cast<double>(count));
    metrics.linf_abs = maxAbs;
    metrics.l1_norm = metrics.l1_abs / safeScale;
    metrics.l2_norm = metrics.l2_abs / safeScale;
    metrics.linf_norm = metrics.linf_abs / safeScale;
    return metrics;
}

double ComputeNormalizedStd(const std::vector<double>& values, double scale) {
    if (values.empty()) return 0.0;
    const double mean = std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
    double sumSq = 0.0;
    for (double value : values) {
        const double diff = value - mean;
        sumSq += diff * diff;
    }
    return std::sqrt(sumSq / static_cast<double>(values.size())) / std::max(scale, 1.0);
}

PressureCellMetrics EvaluatePressureCellMetrics(const MeshManager& mgr,
                                                const SnapshotState& snapshot,
                                                const TestCaseSpec& cfg,
                                                const std::string& csvPath,
                                                bool writeCsv) {
    const auto& cells = mgr.mesh().getCells();
    const std::size_t nUse = std::min(cells.size(), snapshot.p_blocks.size());
    if (nUse == 0) throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] empty state for analytical pressure comparison.");
    const double scale = PressureScale(cfg);
    double sumAbs = 0.0;
    double sumSq = 0.0;
    double maxAbs = 0.0;
    int maxErrCell = -1;
    std::ofstream out;
    std::filesystem::path stagingPath;
    if (writeCsv) {
        out = OpenAsciiStagingStream(
            csvPath,
            "[Test_H_T_CO2_ConstPP_NoFrac] failed to write analytical cell csv",
            stagingPath);
        out << "cell_id,x_m,y_m,target_time_s,actual_time_s,p_num_pa,p_ref_pa,p_abs_err_pa,p_abs_err_over_dp\n";
    }
    for (std::size_t i = 0; i < nUse; ++i) {
        const double pRef = EvaluateAnalyticalPressure(cfg, cells[i].center.m_x, snapshot.actual_time_s);
        const double pAbsErr = std::abs(snapshot.p_blocks[i] - pRef);
        sumAbs += pAbsErr;
        sumSq += pAbsErr * pAbsErr;
        if (pAbsErr >= maxAbs) {
            maxAbs = pAbsErr;
            maxErrCell = static_cast<int>(i);
        }
        if (writeCsv) {
            out << i << ","
                << std::setprecision(12) << cells[i].center.m_x << ","
                << cells[i].center.m_y << ","
                << snapshot.target_time_s << ","
                << snapshot.actual_time_s << ","
                << snapshot.p_blocks[i] << ","
                << pRef << ","
                << pAbsErr << ","
                << (pAbsErr / scale) << "\n";
        }
    }
    if (writeCsv) {
        CommitAsciiStagingFile(
            out,
            stagingPath,
            csvPath,
            "[Test_H_T_CO2_ConstPP_NoFrac] failed to write analytical cell csv");
    }
    PressureCellMetrics metrics;
    metrics.tag = snapshot.tag;
    metrics.requested_fraction = snapshot.requested_fraction;
    metrics.target_time_s = snapshot.target_time_s;
    metrics.actual_time_s = snapshot.actual_time_s;
    metrics.pressure = BuildErrorMetrics(sumAbs, sumSq, maxAbs, static_cast<int>(nUse), scale);
    metrics.max_err_cell = maxErrCell;
    metrics.compare_csv_path = writeCsv ? csvPath : "";
    return metrics;
}

void WriteAnalyticalProfileReferenceCsv(const std::vector<SpatialSamplePoint>& stations,
                                        const SnapshotState& snapshot,
                                        const TestCaseSpec& cfg,
                                        const std::string& family,
                                        const std::string& path) {
    Case2DReferenceIO::WriteAnalyticalProfileReferenceCsv(
        stations,
        snapshot,
        [&](double x, double timeS) { return EvaluateAnalyticalPressure(cfg, x, timeS); },
        family,
        path);
}

void WriteAnalyticalMonitorReferenceCsv(const std::vector<SpatialSamplePoint>& points,
                                        const std::vector<MonitorScheduleState>& schedule,
                                        const TestCaseSpec& cfg,
                                        const std::string& path) {
    Case2DReferenceIO::WriteAnalyticalMonitorReferenceCsv(
        points,
        schedule,
        [&](double x, double timeS) { return EvaluateAnalyticalPressure(cfg, x, timeS); },
        path);
}

void WriteProfileStationDefinitions(const std::vector<SpatialSamplePoint>& stations, const std::string& path) {
    Case2DReferenceIO::WriteProfileStationDefinitions(stations, path);
}

void WriteMonitorPointDefinitions(const std::vector<SpatialSamplePoint>& points, const std::string& path) {
    Case2DReferenceIO::WriteMonitorPointDefinitions(points, path);
}

void WriteProfileSchedule(const std::vector<SnapshotState>& snapshots, const std::string& path) {
    Case2DReferenceIO::WriteProfileSchedule(snapshots, path);
}

void WriteMonitorSchedule(const std::vector<MonitorScheduleState>& schedule, const std::string& path) {
    Case2DReferenceIO::WriteMonitorSchedule(schedule, path);
}

void WriteEngineeringProfileCsv(const std::vector<SpatialSamplePoint>& stations,
                                const std::string& family,
                                const SnapshotState& snapshot,
                                const std::string& path) {
    Case2DReferenceIO::WriteEngineeringProfileCsv(stations, family, snapshot, path);
}

void WriteEngineeringMonitorCsv(const std::vector<SpatialSamplePoint>& points,
                                const std::vector<MonitorScheduleState>& schedule,
                                const std::string& path) {
    Case2DReferenceIO::WriteEngineeringMonitorCsv(points, schedule, path);
}

void WritePropertyTables(const TestCaseSpec& cfg,
                         const std::string& engineeringPath,
                         const std::string& comsolInputPath) {
    Case2DReferenceIO::WritePropertyTables(cfg, engineeringPath, comsolInputPath);
}

void WriteReferenceSpec(const TestCaseSpec& cfg, const TestCaseSummary& summary) {
    Case2DReferenceIO::WriteReferenceSpec(cfg, summary);
}

void WriteAnalyticalNote(const TestCaseSpec& cfg, const std::string& path) {
    Case2DReferenceIO::WriteAnalyticalNote(cfg, ComputePressureDiffusivity(cfg), PressureScale(cfg), path);
}

void WriteComsolReferenceSpec(const TestCaseSpec& cfg, const std::string& caseDir, const std::string& path) {
    Case2DReferenceIO::WriteComsolReferenceSpec(
        cfg,
        caseDir,
        path,
        kComsolRepresentationNonisothermal,
        kComsolRepresentationDlHtManual,
        kComsolRepresentationPdeFallback);
}

ProfileReferenceTable LoadProfileReference(const std::string& path) {
    return Case2DReferenceIO::LoadProfileReference(path);
}

MonitorReferenceSeries LoadMonitorReference(const std::string& path, const std::vector<SpatialSamplePoint>& monitorPoints) {
    return Case2DReferenceIO::LoadMonitorReference(path, monitorPoints);
}

ProfileCompareMetrics EvaluateProfileAgainstReference(const std::vector<ReconstructedSpatialPointSample>& sampledRows,
                                                      const ProfileReferenceTable& reference,
                                                      const TestCaseSpec& cfg,
                                                      const std::string& csvPath,
                                                      bool writeCsv) {
    if (sampledRows.empty()) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] empty sampled profile rows during reference comparison.");
    }
    double pSumAbs = 0.0;
    double pSumSq = 0.0;
    double pMaxAbs = 0.0;
    double tSumAbs = 0.0;
    double tSumSq = 0.0;
    double tMaxAbs = 0.0;
    std::vector<double> pValues;
    std::vector<double> tValues;
    std::ofstream out;
    if (writeCsv) {
        out.open(csvPath.c_str(), std::ios::out | std::ios::trunc);
        if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] failed to write profile compare csv: " + csvPath);
        out << "station_id,label,family,location,target_axis_m,target_x_m,target_y_m,block_id,actual_x_m,actual_y_m,target_time_s,actual_time_s,"
               "p_sample_mode,p_carrier_cell_id,p_carrier_x_m,p_carrier_y_m,p_num_pa,p_ref_pa,p_abs_err_pa,p_abs_err_over_dp,"
               "t_sample_mode,t_carrier_cell_id,t_carrier_x_m,t_carrier_y_m,t_num_k,t_ref_k,t_abs_err_k,t_abs_err_over_dt\n";
    }
    for (const auto& sample : sampledRows) {
        const auto refIt = reference.rows_by_station_id.find(sample.point_id);
        if (refIt == reference.rows_by_station_id.end()) {
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] missing station id in profile reference.");
        }
        const ProfileReferenceRow& refRow = refIt->second;
        const double pNum = sample.pressure.value;
        const double tNum = sample.temperature.value;
        const double pRef = EvaluateAnalyticalPressure(cfg, sample.target_x, sample.actual_time_s);
        const double pAbsErr = std::abs(pNum - pRef);
        const double tAbsErr = std::abs(tNum - refRow.t_ref);
        pSumAbs += pAbsErr;
        pSumSq += pAbsErr * pAbsErr;
        pMaxAbs = std::max(pMaxAbs, pAbsErr);
        tSumAbs += tAbsErr;
        tSumSq += tAbsErr * tAbsErr;
        tMaxAbs = std::max(tMaxAbs, tAbsErr);
        pValues.push_back(pNum);
        tValues.push_back(tNum);
        if (writeCsv) {
            out << sample.point_id << "," << sample.label << "," << sample.family << "," << sample.location << ","
                << std::setprecision(12) << sample.target_axis_m << "," << sample.target_x << "," << sample.target_y << ","
                << sample.block_id << "," << sample.actual_x << "," << sample.actual_y << ","
                << sample.target_time_s << "," << sample.actual_time_s << ","
                << sample.pressure.sample_mode << "," << sample.pressure.carrier_cell_id << ","
                << sample.pressure.carrier_x << "," << sample.pressure.carrier_y << ","
                << pNum << "," << pRef << "," << pAbsErr << "," << (pAbsErr / PressureScale(cfg)) << ","
                << sample.temperature.sample_mode << "," << sample.temperature.carrier_cell_id << ","
                << sample.temperature.carrier_x << "," << sample.temperature.carrier_y << ","
                << tNum << "," << refRow.t_ref << "," << tAbsErr << "," << (tAbsErr / TemperatureScale(cfg)) << "\n";
        }
    }
    ProfileCompareMetrics metrics;
    metrics.family = sampledRows.empty() ? "" : sampledRows.front().family;
    metrics.tag = sampledRows.empty() ? "" : sampledRows.front().time_tag;
    metrics.target_time_s = sampledRows.empty() ? 0.0 : sampledRows.front().target_time_s;
    metrics.actual_time_s = sampledRows.empty() ? 0.0 : sampledRows.front().actual_time_s;
    metrics.pressure = BuildErrorMetrics(pSumAbs, pSumSq, pMaxAbs, static_cast<int>(sampledRows.size()), PressureScale(cfg));
    metrics.temperature = BuildErrorMetrics(tSumAbs, tSumSq, tMaxAbs, static_cast<int>(sampledRows.size()), TemperatureScale(cfg));
    metrics.pressure_uniformity_norm_std = ComputeNormalizedStd(pValues, PressureScale(cfg));
    metrics.temperature_uniformity_norm_std = ComputeNormalizedStd(tValues, TemperatureScale(cfg));
    metrics.compare_csv_path = writeCsv ? csvPath : "";
    return metrics;
}

std::pair<FieldErrorMetrics, FieldErrorMetrics> WriteMonitorCompareCsv(const std::vector<SpatialSamplePoint>& monitorPoints,
                                                                       const std::vector<ReconstructedSpatialPointSample>& sampledRows,
                                                                       const MonitorReferenceSeries& reference,
                                                                       const TestCaseSpec& cfg,
                                                                       const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] failed to write monitor compare csv: " + path);
    out << "sample_id,target_time_s,actual_time_s";
    for (const auto& point : monitorPoints) out << ",carrier_cell_" << point.label;
    for (const auto& point : monitorPoints) out << ",p_sample_mode_" << point.label;
    for (const auto& point : monitorPoints) out << ",t_sample_mode_" << point.label;
    for (const auto& point : monitorPoints) out << ",p_num_" << point.label;
    for (const auto& point : monitorPoints) out << ",p_ref_" << point.label;
    for (const auto& point : monitorPoints) out << ",t_num_" << point.label;
    for (const auto& point : monitorPoints) out << ",t_ref_" << point.label;
    out << "\n";

    double pSumAbs = 0.0;
    double pSumSq = 0.0;
    double pMaxAbs = 0.0;
    double tSumAbs = 0.0;
    double tSumSq = 0.0;
    double tMaxAbs = 0.0;
    int count = 0;

    std::unordered_map<int, std::unordered_map<std::string, ReconstructedSpatialPointSample> > samplesBySchedule;
    std::unordered_map<int, std::pair<double, double> > sampleTimes;
    for (const auto& row : sampledRows) {
        samplesBySchedule[row.sample_id][row.label] = row;
        sampleTimes[row.sample_id] = std::make_pair(row.target_time_s, row.actual_time_s);
    }

    for (const auto& sampleEntry : samplesBySchedule) {
        const int sampleId = sampleEntry.first;
        const auto refIt = reference.sample_id_to_row.find(sampleId);
        if (refIt == reference.sample_id_to_row.end()) {
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] missing monitor sample in reference.");
        }
        const MonitorReferenceRow& refRow = reference.rows[refIt->second];
        const auto timeIt = sampleTimes.find(sampleId);
        if (timeIt == sampleTimes.end()) {
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] missing monitor sample time metadata.");
        }
        out << sampleId << "," << std::setprecision(12) << timeIt->second.first << "," << timeIt->second.second;
        for (const auto& point : monitorPoints) {
            const auto sampleIt = sampleEntry.second.find(point.label);
            if (sampleIt == sampleEntry.second.end()) {
                throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] missing sampled monitor point.");
            }
            out << "," << sampleIt->second.pressure.carrier_cell_id;
        }
        for (const auto& point : monitorPoints) {
            const auto sampleIt = sampleEntry.second.find(point.label);
            out << "," << sampleIt->second.pressure.sample_mode;
        }
        for (const auto& point : monitorPoints) {
            const auto sampleIt = sampleEntry.second.find(point.label);
            out << "," << sampleIt->second.temperature.sample_mode;
        }
        for (const auto& point : monitorPoints) {
            const auto sampleIt = sampleEntry.second.find(point.label);
            out << "," << sampleIt->second.pressure.value;
        }
        for (const auto& point : monitorPoints) {
            const auto sampleIt = sampleEntry.second.find(point.label);
            const double pRef = EvaluateAnalyticalPressure(cfg, sampleIt->second.target_x, sampleIt->second.actual_time_s);
            const double pAbsErr = std::abs(sampleIt->second.pressure.value - pRef);
            pSumAbs += pAbsErr;
            pSumSq += pAbsErr * pAbsErr;
            pMaxAbs = std::max(pMaxAbs, pAbsErr);
            out << "," << pRef;
        }
        for (const auto& point : monitorPoints) {
            const auto sampleIt = sampleEntry.second.find(point.label);
            out << "," << sampleIt->second.temperature.value;
        }
        for (const auto& point : monitorPoints) {
            const auto tIt = refRow.t_ref_by_label.find(point.label);
            if (tIt == refRow.t_ref_by_label.end()) {
                throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] missing monitor point label in temperature reference.");
            }
            const auto sampleIt = sampleEntry.second.find(point.label);
            const double tAbsErr = std::abs(sampleIt->second.temperature.value - tIt->second);
            tSumAbs += tAbsErr;
            tSumSq += tAbsErr * tAbsErr;
            tMaxAbs = std::max(tMaxAbs, tAbsErr);
            out << "," << tIt->second;
            ++count;
        }
        out << "\n";
    }

    const FieldErrorMetrics pMetrics = BuildErrorMetrics(
        pSumAbs, pSumSq, pMaxAbs, count, PressureScale(cfg));
    const FieldErrorMetrics tMetrics = BuildErrorMetrics(
        tSumAbs, tSumSq, tMaxAbs, count, TemperatureScale(cfg));
    return std::make_pair(pMetrics, tMetrics);
}

double ComputeOrder(double coarseError, double fineError, double coarseScale, double fineScale) {
    if (coarseError <= 0.0 || fineError <= 0.0 || coarseScale <= 0.0 || fineScale <= 0.0 || coarseScale == fineScale) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::log(coarseError / fineError) / std::log(coarseScale / fineScale);
}

void AnnotateObservedOrders(std::vector<SweepStudyRow>& rows, bool gridStudy) {
    for (std::size_t i = 1; i < rows.size(); ++i) {
        const double coarseScale = gridStudy ? rows[i - 1].h_char : rows[i - 1].dt_init;
        const double fineScale = gridStudy ? rows[i].h_char : rows[i].dt_init;
        const double coarsePressureError = gridStudy
            ? rows[i - 1].pressure_cell_l2_norm
            : rows[i - 1].pressure_time_self_l2_norm;
        const double finePressureError = gridStudy
            ? rows[i].pressure_cell_l2_norm
            : rows[i].pressure_time_self_l2_norm;
        const double coarseTemperatureError = gridStudy
            ? rows[i - 1].temperature_horizontal_l2_norm
            : rows[i - 1].temperature_time_self_l2_norm;
        const double fineTemperatureError = gridStudy
            ? rows[i].temperature_horizontal_l2_norm
            : rows[i].temperature_time_self_l2_norm;
        rows[i].pressure_order = ComputeOrder(
            coarsePressureError,
            finePressureError,
            coarseScale,
            fineScale);
        rows[i].temperature_order = ComputeOrder(
            coarseTemperatureError,
            fineTemperatureError,
            coarseScale,
            fineScale);
    }
}

bool IsMonotonicNonIncreasingWithTol(const std::vector<double>& values, double absTol) {
    for (std::size_t i = 1; i < values.size(); ++i) {
        if (values[i] > values[i - 1] + absTol) return false;
    }
    return true;
}

double ComputeMinFiniteOrder(const std::vector<SweepStudyRow>& rows, bool pressureOrder) {
    double minOrder = std::numeric_limits<double>::infinity();
    bool hasFinite = false;
    for (const auto& row : rows) {
        const double value = pressureOrder ? row.pressure_order : row.temperature_order;
        if (std::isfinite(value)) {
            minOrder = std::min(minOrder, value);
            hasFinite = true;
        }
    }
    return hasFinite ? minOrder : std::numeric_limits<double>::quiet_NaN();
}

void WriteStudyCsv(const std::vector<SweepStudyRow>& rows, const std::string& path, const std::string& studyTag) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] failed to write study csv: " + path);
    out << "study,label,case_dir,nx,ny,dt_init,h_char,t_end,steps,pressure_cell_l2_norm,pressure_horizontal_l2_norm,temperature_horizontal_l2_norm,pressure_vertical_uniformity_norm_std,temperature_vertical_uniformity_norm_std,pressure_time_self_l2_norm,temperature_time_self_l2_norm,pressure_order,temperature_order\n";
    for (const auto& row : rows) {
        out << studyTag << "," << row.label << "," << row.case_dir << ","
            << row.nx << "," << row.ny << ","
            << std::setprecision(12) << row.dt_init << "," << row.h_char << "," << row.t_end << "," << row.steps << ","
            << row.pressure_cell_l2_norm << "," << row.pressure_horizontal_l2_norm << ","
            << row.temperature_horizontal_l2_norm << "," << row.pressure_vertical_uniformity_norm_std << ","
            << row.temperature_vertical_uniformity_norm_std << "," << row.pressure_time_self_l2_norm << ","
            << row.temperature_time_self_l2_norm << "," << row.pressure_order << "," << row.temperature_order << "\n";
    }
}

std::string ReadTextFile(const std::string& path) {
    std::ifstream in(path.c_str(), std::ios::in);
    if (!in.good()) return "";
    std::ostringstream oss;
    oss << in.rdbuf();
    return oss.str();
}

std::string ExtractComsolRepresentation(const std::string& runSummaryPath) {
    std::ifstream in(runSummaryPath.c_str(), std::ios::in);
    if (!in.good()) return "";
    std::string line;
    const std::string token = "Representation:";
    while (std::getline(in, line)) {
        const std::size_t pos = line.find(token);
        if (pos == std::string::npos) continue;
        std::string value = Trim(line.substr(pos + token.size()));
        if (!value.empty() && value.front() == '`') value.erase(value.begin());
        while (!value.empty() && (value.back() == '`' || value.back() == '.' || value.back() == '\r')) value.pop_back();
        return Trim(value);
    }
    return "";
}

bool ComsolFineCheckWasSkipped(const std::string& meshCheckPath) {
    const std::string text = ReadTextFile(meshCheckPath);
    return text.find("fine_check_skipped=true") != std::string::npos;
}

bool ReferenceFilesReady(TestCaseSummary& summary, const std::vector<SnapshotState>& snapshots) {
    return Case2DReferenceIO::ReferenceFilesReady(summary, snapshots, OrderedProfileFamilies());
}

std::string ResolveComsolCaseDir(const TestCaseSummary& summary) {
    const std::filesystem::path comsolDir(summary.comsol_output_dir);
    const std::filesystem::path caseDir = comsolDir.parent_path().parent_path();
    if (caseDir.empty()) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] invalid COMSOL case dir derived from: " + summary.comsol_output_dir);
    }
    return caseDir.generic_string();
}

void InvokeComsolReferenceGeneration(const TestCaseSpec& cfg, const TestCaseSummary& summary) {
    std::ostringstream cmd;
    cmd << "\"C:\\Windows\\System32\\WindowsPowerShell\\v1.0\\powershell.exe\" -ExecutionPolicy Bypass -File \"" << cfg.comsol_wrapper_relpath
        << "\" -Mode All -CaseDir \"" << ResolveComsolCaseDir(summary) << "\"";
    const int code = std::system(cmd.str().c_str());
    if (code != 0) {
        throw std::runtime_error(
            "[Test_H_T_CO2_ConstPP_NoFrac] COMSOL reference generation failed with exit code " + std::to_string(code));
    }
}

void WriteMetricsCsv(const TestCaseSpec& cfg, const TestCaseSummary& summary) {
    Case2DReferenceIO::WriteMetricsCsv(cfg, summary);
}

void WriteValidationSummaryCsv(const TestCaseSpec& cfg, const TestCaseSummary& summary) {
    Case2DReferenceIO::WriteValidationSummaryCsv(
        cfg,
        summary,
        "pressure_abs_cell__temperature_fixed_target_abs_profile_interior_only",
        "fixed_target_self_against_finest_dt");
}

bool FinalMetricsWithinThreshold(const TestCaseSpec& cfg, const TestCaseSummary& summary) {
    return summary.final_pressure_cell_l2_norm <= cfg.pressure_l2_threshold &&
           summary.final_pressure_cell_linf_norm <= cfg.pressure_linf_threshold &&
           summary.final_temperature_horizontal_l2_norm <= cfg.temperature_l2_threshold &&
           summary.final_temperature_horizontal_linf_norm <= cfg.temperature_linf_threshold &&
           summary.final_monitor_temperature_l2_norm <= cfg.temperature_l2_threshold &&
           summary.final_monitor_temperature_linf_norm <= cfg.temperature_linf_threshold &&
           summary.final_pressure_vertical_uniformity_norm_std <= cfg.vertical_uniformity_threshold &&
           summary.final_temperature_vertical_uniformity_norm_std <= cfg.vertical_uniformity_threshold;
}

FIM_Engine::TransientSolverParams BuildSolverParams(const TestCaseSpec& cfg) {
    FIM_Engine::TransientSolverParams p;
    p.max_steps = cfg.max_steps;
    p.dt_init = cfg.dt_init;
    p.dt_min = cfg.dt_min;
    p.dt_max = cfg.dt_max;
    p.target_end_time_s = cfg.target_end_time_s;
    p.max_newton_iter = cfg.max_newton_iter;
    p.abs_res_tol = cfg.abs_res_tol;
    p.rel_res_tol = cfg.rel_res_tol;
    p.rel_update_tol = cfg.rel_update_tol;
    p.lin_solver = cfg.lin_solver;
    p.amgcl_use_fallback_sparselu = cfg.amgcl_use_fallback_sparselu;
    p.amgcl_tol = cfg.amgcl_tol;
    p.amgcl_maxiter = cfg.amgcl_maxiter;
    p.enable_non_orthogonal_correction = cfg.enable_non_orthogonal_correction;
    p.enable_row_scaling = cfg.enable_row_scaling;
    p.max_dP = cfg.max_dP;
    p.max_dT = cfg.max_dT;
    p.max_dSw = 0.1;
    p.min_alpha = 1.0e-8;
    p.enable_armijo_line_search = cfg.enable_armijo_line_search;
    p.rollback_shrink_factor = cfg.rollback_shrink_factor;
    p.dt_relres_grow_factor = cfg.dt_relres_grow_factor;
    p.gravity_vector = cfg.gravity_vector;
    p.diag_level = cfg.diag_level;
    return p;
}

CaseRunArtifacts RunSingleCaseCore(const TestCaseSpec& cfg, const std::string& outputDirOverride) {
    CaseRunArtifacts artifacts;
    const std::string caseDir = outputDirOverride.empty()
        ? (cfg.output_base_dir + "/" + cfg.sub_dir + "/" + cfg.case_name)
        : outputDirOverride;
    ConfigureSummaryPaths(artifacts.summary, cfg, caseDir);
    EnsureSummaryDirs(artifacts.summary);

    std::ofstream convergenceLog(artifacts.summary.convergence_log_path.c_str(), std::ios::out | std::ios::trunc);
    if (!convergenceLog.good()) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] failed to open convergence log: " + artifacts.summary.convergence_log_path);
    }
    std::ofstream runLog(artifacts.summary.run_log_path.c_str(), std::ios::out | std::ios::trunc);
    if (!runLog.good()) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] failed to open run log: " + artifacts.summary.run_log_path);
    }
    runLog << "case=" << cfg.case_name << "\n";
    runLog << "reference_mode=" << artifacts.summary.resolved_reference_mode << "\n";
    runLog << "output_dir=" << artifacts.summary.case_dir << "\n";

    MeshManager mgr(cfg.lx, cfg.ly, 0.0, cfg.nx, cfg.ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OverRelaxed);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    const std::size_t nCells = mgr.mesh().getCells().size();
    const int totalBlocks = mgr.getTotalDOFCount();
    if (nCells == 0) throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] matrix cell count is zero");
    artifacts.summary.n_cells = static_cast<int>(nCells);

    artifacts.profile_stations = BuildProfileStations(mgr, cfg);
    artifacts.monitor_points = BuildMonitorPoints(mgr, cfg);
    artifacts.snapshots = BuildSnapshots(cfg);
    artifacts.monitor_schedule = BuildMonitorSchedule(cfg);

    FIM_Engine::InitialConditions ic;
    ic.P_init = cfg.p_init;
    ic.T_init = cfg.t_init;
    ic.Sw_init = 1.0;

    BoundarySetting::BoundaryConditionManager bcP;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT, cfg.p_left);
    bcP.SetDirichletBC(MeshTags::RIGHT, cfg.p_right);
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcP.SetNeumannBC(MeshTags::TOP, 0.0);

    BoundarySetting::BoundaryConditionManager bcT;
    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT, cfg.t_left);
    bcT.SetDirichletBC(MeshTags::RIGHT, cfg.t_right);
    bcT.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcT.SetNeumannBC(MeshTags::TOP, 0.0);

    const auto pEqCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
    const auto tEqCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
    VTKBoundaryVisualizationContext bcVizCtx;
    bcVizCtx.water_family_policy = VTKBCWaterFamilyDerivePolicy::FollowPrimaryFluid;
    bcVizCtx.primary_fluid_model = VTKBCPrimaryFluidModel::CO2;
    bcVizCtx.bindings.push_back(VTKBCVariableBinding{pEqCfg.pressure_field, &bcP, VTKBCTransportKind::Pressure});
    bcVizCtx.bindings.push_back(VTKBCVariableBinding{tEqCfg.temperatue_field, &bcT, VTKBCTransportKind::Temperature});

    std::vector<double> pBlocksLatest(static_cast<std::size_t>(std::max(totalBlocks, 0)), cfg.p_init);
    std::vector<double> tBlocksLatest(static_cast<std::size_t>(std::max(totalBlocks, 0)), cfg.t_init);
    SyncPTFieldsToFM(mgr, fm, pBlocksLatest, tBlocksLatest, cfg.p_init, cfg.t_init);
    if (cfg.export_vtk) {
        PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(artifacts.summary.case_dir + "/initial.vtk", 0.0);
    }

    for (auto& sample : artifacts.monitor_schedule) {
        if (sample.target_time_s <= kPendingTimeTolerance) {
            sample.captured = true;
            sample.actual_time_s = 0.0;
            sample.p_blocks = pBlocksLatest;
            sample.t_blocks = tBlocksLatest;
        }
    }

    bool midExported = false;
    int iterSum = 0;
    int iterCount = 0;
    int maxIters = 0;
    double prevAcceptedTime = 0.0;
    std::vector<double> prevAcceptedPBlocks = pBlocksLatest;
    std::vector<double> prevAcceptedTBlocks = tBlocksLatest;

    FIM_Engine::TransientOptionalModules<MeshManager, FieldManager_2D> modules;
    modules.pressure_bc = &bcP;
    modules.temperature_bc = &bcT;
    modules.disable_default_vtk_output = true;

    FluidConstantProperties co2_props;
    co2_props.rho = cfg.co2_rho_const;
    co2_props.mu = cfg.co2_mu_const;
    co2_props.cp = cfg.co2_cp_const;
    co2_props.cv = cfg.co2_cv_const;
    co2_props.k = cfg.co2_k_const;
    if (cfg.use_variable_properties) {
        modules.SetFluidModelConfig(FIM_Engine::UnifiedFluidModelConfig::MakeSinglePhaseCO2EOS());
    } else {
        modules.SetFluidModelConfig(FIM_Engine::UnifiedFluidModelConfig::MakeSinglePhaseCO2Constant(co2_props));
    }

    modules.property_initializer = [&cfg](MeshManager&, FieldManager_2D& fld) {
        const auto rock = PhysicalProperties_string_op::Rock();
        const auto frac = PhysicalProperties_string_op::Fracture_string();
        const auto gas = PhysicalProperties_string_op::CO2();
        const auto water = PhysicalProperties_string_op::Water();

        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.k_xx_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.k_yy_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.k_zz_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.phi_tag, cfg.matrix_phi), cfg.matrix_phi);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.c_r_tag, cfg.matrix_ct), cfg.matrix_ct);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.rho_tag, cfg.matrix_rho_r), cfg.matrix_rho_r);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.cp_tag, cfg.matrix_cp_r), cfg.matrix_cp_r);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.lambda_tag, cfg.matrix_lambda_r), cfg.matrix_lambda_r);

        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.k_t_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.k_n_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.phi_tag, cfg.matrix_phi), cfg.matrix_phi);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.c_r_tag, cfg.matrix_ct), cfg.matrix_ct);

        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(gas.k_tag, cfg.co2_k_const), cfg.co2_k_const);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(gas.k_tag, cfg.co2_k_const), cfg.co2_k_const);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(water.k_tag, cfg.co2_k_const), cfg.co2_k_const);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(water.k_tag, cfg.co2_k_const), cfg.co2_k_const);
    };

    modules.on_step_accepted =
        [&](int step, double timeS, double dtUsedS, int newtonIters, double residualInf,
            int totalRollbacks, const std::string& convergeMode,
            const std::vector<double>& pVec, const std::vector<double>& tVec,
            const std::vector<double>*) {
            if (step <= 0) return;
            if (!pVec.empty()) pBlocksLatest = pVec;
            if (!tVec.empty()) tBlocksLatest = tVec;
            SyncPTFieldsToFM(mgr, fm, pBlocksLatest, tBlocksLatest, cfg.p_init, cfg.t_init);

            artifacts.summary.steps = step;
            artifacts.summary.t_end = timeS;
            artifacts.summary.total_rollbacks = totalRollbacks;
            iterSum += newtonIters;
            iterCount += 1;
            maxIters = std::max(maxIters, newtonIters);

            convergenceLog << "[Step " << step << "] t=" << std::scientific << std::setprecision(8) << timeS
                           << " dt=" << dtUsedS
                           << " iters=" << newtonIters
                           << " residual_inf=" << residualInf
                           << " rollbacks=" << totalRollbacks
                           << " mode=" << convergeMode << "\n";

            for (auto& snapshot : artifacts.snapshots) {
                CaptureSnapshotIfCrossed(
                    prevAcceptedTime, prevAcceptedPBlocks, prevAcceptedTBlocks,
                    timeS, pBlocksLatest, tBlocksLatest, snapshot);
            }
            for (auto& sample : artifacts.monitor_schedule) {
                CaptureMonitorIfCrossed(
                    prevAcceptedTime, prevAcceptedPBlocks, prevAcceptedTBlocks,
                    timeS, pBlocksLatest, tBlocksLatest, sample);
            }
            prevAcceptedTime = timeS;
            prevAcceptedPBlocks = pBlocksLatest;
            prevAcceptedTBlocks = tBlocksLatest;

            if (cfg.export_vtk && !midExported && timeS >= 0.5 * cfg.target_end_time_s - 1.0e-12) {
                PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(artifacts.summary.case_dir + "/mid.vtk", timeS);
                midExported = true;
            }
        };

    const auto params = BuildSolverParams(cfg);
    FIM_Engine::RunGenericFIMTransient<2>(cfg.case_name, mgr, fm, ic, {}, params, FIM_Engine::SolverRoute::FIM, modules);

    SyncPTFieldsToFM(mgr, fm, pBlocksLatest, tBlocksLatest, cfg.p_init, cfg.t_init);
    if (cfg.export_vtk) {
        if (!midExported) PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(artifacts.summary.case_dir + "/mid.vtk", artifacts.summary.t_end);
        PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(artifacts.summary.case_dir + "/final.vtk", artifacts.summary.t_end);
    }

    FinalizeMissingSnapshots(artifacts.summary.t_end, pBlocksLatest, tBlocksLatest, artifacts.snapshots);
    FinalizeMissingMonitorSamples(artifacts.summary.t_end, pBlocksLatest, tBlocksLatest, artifacts.monitor_schedule);

    artifacts.summary.max_iters = maxIters;
    artifacts.summary.avg_iters = (iterCount > 0) ? (static_cast<double>(iterSum) / static_cast<double>(iterCount)) : 0.0;
    artifacts.summary.h_char = ComputeMeshCharLength(mgr);

    runLog << "n_cells=" << artifacts.summary.n_cells << "\n";
    runLog << "steps=" << artifacts.summary.steps << "\n";
    runLog << "t_end=" << artifacts.summary.t_end << "\n";

    for (const auto& snapshot : artifacts.snapshots) {
        artifacts.summary.pressure_cell_metrics.push_back(EvaluatePressureCellMetrics(
            mgr,
            snapshot,
            cfg,
            artifacts.summary.analytic_dir + "/pressure_cells_" + snapshot.tag + ".csv",
            cfg.emit_detailed_outputs));
    }

    artifacts.sampled_profile_points.clear();
    for (const auto& snapshot : artifacts.snapshots) {
        const std::vector<ReconstructedSpatialPointSample> sampled =
            BuildReconstructedSamples(
                mgr,
                artifacts.profile_stations,
                snapshot.p_blocks,
                snapshot.t_blocks,
                snapshot.tag,
                -1,
                snapshot.target_time_s,
                snapshot.actual_time_s,
                cfg);
        artifacts.sampled_profile_points.insert(
            artifacts.sampled_profile_points.end(),
            sampled.begin(),
            sampled.end());
    }

    artifacts.sampled_monitor_points.clear();
    for (const auto& sample : artifacts.monitor_schedule) {
        const std::vector<ReconstructedSpatialPointSample> sampled =
            BuildReconstructedSamples(
                mgr,
                artifacts.monitor_points,
                sample.p_blocks,
                sample.t_blocks,
                "monitor",
                sample.sample_id,
                sample.target_time_s,
                sample.actual_time_s,
                cfg);
        artifacts.sampled_monitor_points.insert(
            artifacts.sampled_monitor_points.end(),
            sampled.begin(),
            sampled.end());
    }

    if (cfg.emit_detailed_outputs) {
        WriteProfileStationDefinitions(artifacts.profile_stations, artifacts.summary.profile_station_definitions_path);
        WriteMonitorPointDefinitions(artifacts.monitor_points, artifacts.summary.monitor_point_definitions_path);
        WriteProfileSchedule(artifacts.snapshots, artifacts.summary.profile_schedule_path);
        WriteMonitorSchedule(artifacts.monitor_schedule, artifacts.summary.monitor_schedule_path);
        for (const auto& snapshot : artifacts.snapshots) {
            for (const std::string& family : OrderedProfileFamilies()) {
                WriteEngineeringProfileCsv(
                    artifacts.profile_stations,
                    family,
                    snapshot,
                    artifacts.summary.engineering_dir + "/eng_profile_" + family + "_" + snapshot.tag + ".csv");
                WriteAnalyticalProfileReferenceCsv(
                    artifacts.profile_stations,
                    snapshot,
                    cfg,
                    family,
                    artifacts.summary.analytic_dir + "/pressure_profile_" + family + "_" + snapshot.tag + ".csv");
            }
        }
        WriteEngineeringMonitorCsv(artifacts.monitor_points, artifacts.monitor_schedule, artifacts.summary.eng_monitor_timeseries_path);
        WriteAnalyticalMonitorReferenceCsv(
            artifacts.monitor_points,
            artifacts.monitor_schedule,
            cfg,
            artifacts.summary.analytic_dir + "/pressure_monitor_timeseries.csv");
        WritePropertyTables(cfg, artifacts.summary.property_table_path, artifacts.summary.comsol_property_table_path);
        WriteReferenceSpec(cfg, artifacts.summary);
        WriteAnalyticalNote(cfg, artifacts.summary.analytical_note_path);
        WriteComsolReferenceSpec(cfg, artifacts.summary.case_dir, artifacts.summary.comsol_reference_spec_path);
    }

    if (!artifacts.summary.pressure_cell_metrics.empty()) {
        const PressureCellMetrics& finalMetric = artifacts.summary.pressure_cell_metrics.back();
        artifacts.summary.final_pressure_cell_l1_norm = finalMetric.pressure.l1_norm;
        artifacts.summary.final_pressure_cell_l2_norm = finalMetric.pressure.l2_norm;
        artifacts.summary.final_pressure_cell_linf_norm = finalMetric.pressure.linf_norm;
    }

    if (!artifacts.snapshots.empty()) {
        const SnapshotState& finalSnapshot = artifacts.snapshots.back();
        const std::vector<SpatialSamplePoint> horizontalStations =
            FilterStationsByFamily(artifacts.profile_stations, kFamilyMatrixHorizontal);
        const std::vector<SpatialSamplePoint> verticalStations =
            FilterStationsByFamily(artifacts.profile_stations, kFamilyMatrixVerticalMidline);

        double pHAbsSum = 0.0;
        double pHSumSq = 0.0;
        double pHMaxAbs = 0.0;
        for (const auto& station : horizontalStations) {
            const double pNum = finalSnapshot.p_blocks[static_cast<std::size_t>(station.block_id)];
            const double pRef = EvaluateAnalyticalPressure(cfg, station.actual_x, finalSnapshot.actual_time_s);
            const double absErr = std::abs(pNum - pRef);
            pHAbsSum += absErr;
            pHSumSq += absErr * absErr;
            pHMaxAbs = std::max(pHMaxAbs, absErr);
        }
        const FieldErrorMetrics horizontalPMetrics = BuildErrorMetrics(
            pHAbsSum, pHSumSq, pHMaxAbs, static_cast<int>(horizontalStations.size()), PressureScale(cfg));
        artifacts.summary.final_pressure_horizontal_l2_norm = horizontalPMetrics.l2_norm;
        artifacts.summary.final_pressure_horizontal_linf_norm = horizontalPMetrics.linf_norm;

        std::vector<double> verticalPValues;
        for (const auto& station : verticalStations) {
            verticalPValues.push_back(finalSnapshot.p_blocks[static_cast<std::size_t>(station.block_id)]);
        }
        artifacts.summary.final_pressure_vertical_uniformity_norm_std =
            ComputeNormalizedStd(verticalPValues, PressureScale(cfg));
    }

    return artifacts;
}

TestCaseSummary RunCase(const TestCaseSpec& cfg, const std::string& outputDirOverride) {
    CaseRunArtifacts mainArtifacts = RunSingleCaseCore(cfg, outputDirOverride);
    TestCaseSummary& summary = mainArtifacts.summary;

    const bool readyBefore = ReferenceFilesReady(summary, mainArtifacts.snapshots);
    if (!readyBefore && cfg.allow_full_workflow_comsol_autorun) {
        try {
            InvokeComsolReferenceGeneration(cfg, summary);
            ReferenceFilesReady(summary, mainArtifacts.snapshots);
        } catch (const std::exception& ex) {
            std::ofstream runLog(summary.run_log_path.c_str(), std::ios::out | std::ios::app);
            if (runLog.good()) runLog << "comsol_autorun_error=" << ex.what() << "\n";
        }
    }

    if (!summary.reference_ready) {
        summary.validation_status = "missing_reference";
        summary.validation_performed = false;
        summary.validation_passed = false;
        Case2DReferenceIO::WriteValidationSummary(
            cfg,
            summary,
            "pressure uses absolute cell `L2`; temperature uses final-time horizontal fixed-target profile `L2` on interior stations only, excluding two profile intervals adjacent to each Dirichlet boundary.",
            "pressure and temperature use fixed-target self-convergence against the finest configured `dt_init`; absolute reference errors remain diagnostic only.");
        WriteValidationSummaryCsv(cfg, summary);
        WriteMetricsCsv(cfg, summary);
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] COMSOL temperature reference files are missing.");
    }

    summary.comsol_representation = ExtractComsolRepresentation(summary.comsol_output_dir + "/comsol_run_summary.md");
    if (summary.comsol_representation.empty()) summary.comsol_representation = kComsolRepresentationDlHtManual;
    summary.comsol_fine_check_skipped =
        ComsolFineCheckWasSkipped(summary.comsol_output_dir + "/comsol_reference_mesh_check.txt");

    const MonitorReferenceSeries monitorReference =
        LoadMonitorReference(summary.comsol_output_dir + "/comsol_monitor_timeseries.csv", mainArtifacts.monitor_points);

    const SnapshotState* finalSnapshot = mainArtifacts.snapshots.empty() ? nullptr : &mainArtifacts.snapshots.back();
    if (!finalSnapshot) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] no engineering snapshots available.");
    }

    summary.report_metrics.clear();
    for (const auto& snapshot : mainArtifacts.snapshots) {
        for (const std::string& family : OrderedProfileFamilies()) {
            const std::vector<ReconstructedSpatialPointSample> sampledProfileRows =
                FilterReconstructedSamples(mainArtifacts.sampled_profile_points, family, snapshot.tag);
            const ProfileReferenceTable reference = LoadProfileReference(
                summary.comsol_output_dir + "/comsol_profile_" + family + "_" + snapshot.tag + ".csv");
            const std::string comparePath = summary.studies_dir + "/compare_profile_" + family + "_" + snapshot.tag + ".csv";
            summary.report_metrics.push_back(EvaluateProfileAgainstReference(
                sampledProfileRows,
                reference,
                cfg,
                comparePath,
                cfg.emit_detailed_outputs));
        }
    }

    summary.final_temperature_horizontal_l2_norm = std::numeric_limits<double>::quiet_NaN();
    summary.final_temperature_horizontal_linf_norm = std::numeric_limits<double>::quiet_NaN();
    summary.final_temperature_vertical_uniformity_norm_std = std::numeric_limits<double>::quiet_NaN();
    for (const auto& metric : summary.report_metrics) {
        if (metric.tag != finalSnapshot->tag) continue;
        if (metric.family == kFamilyMatrixHorizontal) {
            summary.final_pressure_horizontal_l2_norm = metric.pressure.l2_norm;
            summary.final_pressure_horizontal_linf_norm = metric.pressure.linf_norm;
            summary.final_temperature_horizontal_l2_norm = metric.temperature.l2_norm;
            summary.final_temperature_horizontal_linf_norm = metric.temperature.linf_norm;
        } else if (metric.family == kFamilyMatrixVerticalMidline) {
            summary.final_pressure_vertical_uniformity_norm_std = metric.pressure_uniformity_norm_std;
            summary.final_temperature_vertical_uniformity_norm_std = metric.temperature_uniformity_norm_std;
        }
    }

    const std::pair<FieldErrorMetrics, FieldErrorMetrics> monitorMetrics = WriteMonitorCompareCsv(
        mainArtifacts.monitor_points,
        mainArtifacts.sampled_monitor_points,
        monitorReference,
        cfg,
        summary.studies_dir + "/compare_monitor_timeseries.csv");
    summary.final_monitor_pressure_l2_norm = monitorMetrics.first.l2_norm;
    summary.final_monitor_pressure_linf_norm = monitorMetrics.first.linf_norm;
    summary.final_monitor_temperature_l2_norm = monitorMetrics.second.l2_norm;
    summary.final_monitor_temperature_linf_norm = monitorMetrics.second.linf_norm;

    std::vector<SweepStudyRow> gridRows;
    std::vector<SweepStudyRow> timeRows;

    if (cfg.enable_grid_convergence_study) {
        for (const auto& gridCase : cfg.grid_sweep_cases) {
            TestCaseSpec sweepCfg = cfg;
            sweepCfg.nx = gridCase.first;
            sweepCfg.ny = gridCase.second;
            sweepCfg.case_name = cfg.case_name + "_grid_" + std::to_string(sweepCfg.nx) + "x" + std::to_string(sweepCfg.ny);
            sweepCfg.enable_grid_convergence_study = false;
            sweepCfg.enable_time_sensitivity_study = false;
            sweepCfg.export_vtk = false;
            sweepCfg.emit_detailed_outputs = false;

            CaseRunArtifacts sweepArtifacts = RunSingleCaseCore(
                sweepCfg,
                summary.studies_dir + "/grid_" + std::to_string(sweepCfg.nx) + "x" + std::to_string(sweepCfg.ny));

            SweepStudyRow row;
            row.label = "grid_" + std::to_string(sweepCfg.nx) + "x" + std::to_string(sweepCfg.ny);
            row.case_dir = sweepArtifacts.summary.case_dir;
            row.nx = sweepCfg.nx;
            row.ny = sweepCfg.ny;
            row.dt_init = sweepCfg.dt_init;
            row.h_char = sweepArtifacts.summary.h_char;
            row.t_end = sweepArtifacts.summary.t_end;
            row.steps = sweepArtifacts.summary.steps;
            row.pressure_cell_l2_norm = sweepArtifacts.summary.final_pressure_cell_l2_norm;
            const std::vector<ReconstructedSpatialPointSample> horizontalSamples =
                FilterReconstructedSamples(sweepArtifacts.sampled_profile_points, kFamilyMatrixHorizontal, finalSnapshot->tag);
            const ProfileReferenceTable horizontalReference =
                LoadProfileReference(summary.comsol_output_dir + "/comsol_profile_matrix_horizontal_" + finalSnapshot->tag + ".csv");
            const ProfileCompareMetrics hMetric = EvaluateProfileAgainstReference(
                FilterInteriorHorizontalSamples(horizontalSamples, sweepCfg),
                FilterInteriorHorizontalReference(horizontalReference, sweepCfg),
                sweepCfg,
                "",
                false);
            const ProfileCompareMetrics vMetric = EvaluateProfileAgainstReference(
                FilterReconstructedSamples(sweepArtifacts.sampled_profile_points, kFamilyMatrixVerticalMidline, finalSnapshot->tag),
                LoadProfileReference(summary.comsol_output_dir + "/comsol_profile_matrix_vertical_midline_" + finalSnapshot->tag + ".csv"),
                sweepCfg,
                "",
                false);
            row.pressure_horizontal_l2_norm = hMetric.pressure.l2_norm;
            row.pressure_vertical_uniformity_norm_std = vMetric.pressure_uniformity_norm_std;
            row.temperature_horizontal_l2_norm = hMetric.temperature.l2_norm;
            row.temperature_vertical_uniformity_norm_std = vMetric.temperature_uniformity_norm_std;
            gridRows.push_back(row);
        }
        AnnotateObservedOrders(gridRows, true);
        WriteStudyCsv(gridRows, summary.grid_convergence_csv_path, "grid");
        std::vector<double> pressureErrors;
        std::vector<double> temperatureErrors;
        for (const auto& row : gridRows) {
            pressureErrors.push_back(row.pressure_cell_l2_norm);
            temperatureErrors.push_back(row.temperature_horizontal_l2_norm);
        }
        summary.grid_pressure_order_min = ComputeMinFiniteOrder(gridRows, true);
        summary.grid_temperature_order_min = ComputeMinFiniteOrder(gridRows, false);
        summary.grid_convergence_ok =
            IsMonotonicNonIncreasingWithTol(pressureErrors, kGridTimeMonotoneAbsTol) &&
            IsMonotonicNonIncreasingWithTol(temperatureErrors, kGridTimeMonotoneAbsTol) &&
            std::isfinite(summary.grid_pressure_order_min) &&
            std::isfinite(summary.grid_temperature_order_min) &&
            summary.grid_pressure_order_min >= cfg.grid_pressure_order_threshold &&
            summary.grid_temperature_order_min >= cfg.grid_temperature_order_threshold;
    }

    if (cfg.enable_time_sensitivity_study) {
        const std::pair<int, int> finestGrid = GetFinestGridCase(cfg);
        for (double dtInit : cfg.time_step_sweep) {
            TestCaseSpec sweepCfg = cfg;
            sweepCfg.nx = finestGrid.first;
            sweepCfg.ny = finestGrid.second;
            sweepCfg.dt_init = dtInit;
            sweepCfg.case_name = cfg.case_name + "_dt_" + std::to_string(static_cast<int>(std::round(dtInit))) + "s";
            sweepCfg.enable_grid_convergence_study = false;
            sweepCfg.enable_time_sensitivity_study = false;
            sweepCfg.export_vtk = false;
            sweepCfg.emit_detailed_outputs = false;

            CaseRunArtifacts sweepArtifacts = RunSingleCaseCore(
                sweepCfg,
                summary.studies_dir + "/dt_" + std::to_string(static_cast<int>(std::round(dtInit))) + "s");

            SweepStudyRow row;
            row.label = "dt_" + std::to_string(static_cast<int>(std::round(dtInit))) + "s";
            row.case_dir = sweepArtifacts.summary.case_dir;
            row.nx = sweepCfg.nx;
            row.ny = sweepCfg.ny;
            row.dt_init = dtInit;
            row.h_char = sweepArtifacts.summary.h_char;
            row.t_end = sweepArtifacts.summary.t_end;
            row.steps = sweepArtifacts.summary.steps;
            row.pressure_cell_l2_norm = sweepArtifacts.summary.final_pressure_cell_l2_norm;
            const std::vector<ReconstructedSpatialPointSample> horizontalSamples =
                FilterReconstructedSamples(sweepArtifacts.sampled_profile_points, kFamilyMatrixHorizontal, finalSnapshot->tag);
            const ProfileReferenceTable horizontalReference =
                LoadProfileReference(summary.comsol_output_dir + "/comsol_profile_matrix_horizontal_" + finalSnapshot->tag + ".csv");
            const ProfileCompareMetrics hMetric = EvaluateProfileAgainstReference(
                FilterInteriorHorizontalSamples(horizontalSamples, sweepCfg),
                FilterInteriorHorizontalReference(horizontalReference, sweepCfg),
                sweepCfg,
                "",
                false);
            const ProfileCompareMetrics vMetric = EvaluateProfileAgainstReference(
                FilterReconstructedSamples(sweepArtifacts.sampled_profile_points, kFamilyMatrixVerticalMidline, finalSnapshot->tag),
                LoadProfileReference(summary.comsol_output_dir + "/comsol_profile_matrix_vertical_midline_" + finalSnapshot->tag + ".csv"),
                sweepCfg,
                "",
                false);
            row.pressure_horizontal_l2_norm = hMetric.pressure.l2_norm;
            row.pressure_vertical_uniformity_norm_std = vMetric.pressure_uniformity_norm_std;
            row.temperature_horizontal_l2_norm = hMetric.temperature.l2_norm;
            row.temperature_vertical_uniformity_norm_std = vMetric.temperature_uniformity_norm_std;
            row.profile_samples = sweepArtifacts.sampled_profile_points;
            timeRows.push_back(row);
        }
        auto finestIt = std::min_element(timeRows.begin(), timeRows.end(), [](const SweepStudyRow& a, const SweepStudyRow& b) {
            return a.dt_init < b.dt_init;
        });
        if (finestIt == timeRows.end()) {
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] time sensitivity study produced no rows.");
        }
        for (auto& row : timeRows) {
            if (std::abs(row.dt_init - finestIt->dt_init) <= 1.0e-12) {
                row.pressure_time_self_l2_norm = 0.0;
                row.temperature_time_self_l2_norm = 0.0;
                continue;
            }
            const std::pair<FieldErrorMetrics, FieldErrorMetrics> selfMetrics =
                EvaluateSelfConvergenceMetrics(row.profile_samples, finestIt->profile_samples, cfg);
            row.pressure_time_self_l2_norm = selfMetrics.first.l2_norm;
            row.temperature_time_self_l2_norm = selfMetrics.second.l2_norm;
        }
        AnnotateObservedOrders(timeRows, false);
        WriteStudyCsv(timeRows, summary.time_sensitivity_csv_path, "time");
        std::vector<double> pressureErrors;
        std::vector<double> temperatureErrors;
        for (const auto& row : timeRows) {
            pressureErrors.push_back(row.pressure_time_self_l2_norm);
            temperatureErrors.push_back(row.temperature_time_self_l2_norm);
        }
        summary.time_pressure_order_min = ComputeMinFiniteOrder(timeRows, true);
        summary.time_temperature_order_min = ComputeMinFiniteOrder(timeRows, false);
        summary.time_sensitivity_ok =
            IsMonotonicNonIncreasingWithTol(pressureErrors, kGridTimeMonotoneAbsTol) &&
            IsMonotonicNonIncreasingWithTol(temperatureErrors, kGridTimeMonotoneAbsTol) &&
            std::isfinite(summary.time_pressure_order_min) &&
            std::isfinite(summary.time_temperature_order_min) &&
            summary.time_pressure_order_min >= cfg.time_pressure_order_threshold &&
            summary.time_temperature_order_min >= cfg.time_temperature_order_threshold;
    }

    summary.validation_performed = true;
    summary.validation_passed =
        FinalMetricsWithinThreshold(cfg, summary) &&
        (!cfg.enable_grid_convergence_study || summary.grid_convergence_ok) &&
        (!cfg.enable_time_sensitivity_study || summary.time_sensitivity_ok);
    summary.validation_status = summary.validation_passed ? "passed" : "failed";

    Case2DReferenceIO::WriteValidationSummary(
        cfg,
        summary,
        "pressure uses absolute cell `L2`; temperature uses final-time horizontal fixed-target profile `L2` on interior stations only, excluding two profile intervals adjacent to each Dirichlet boundary.",
        "pressure and temperature use fixed-target self-convergence against the finest configured `dt_init`; absolute reference errors remain diagnostic only.");
    WriteValidationSummaryCsv(cfg, summary);
    WriteMetricsCsv(cfg, summary);
    {
        Case2DMatlab::ValidationPlotScriptSpec spec;
        spec.script_path = summary.matlab_script_path;
        spec.profile_families = OrderedProfileFamilies();
        spec.final_tag = "t100pct";
        Case2DMatlab::WriteNoFracPTValidationPlotScript(spec);
    }

    if (!summary.validation_passed) {
        std::ostringstream oss;
        oss << "[Test_H_T_CO2_ConstPP_NoFrac] validation failed: "
            << "p_cell_l2=" << summary.final_pressure_cell_l2_norm
            << ", p_cell_linf=" << summary.final_pressure_cell_linf_norm
            << ", t_profile_l2=" << summary.final_temperature_horizontal_l2_norm
            << ", t_profile_linf=" << summary.final_temperature_horizontal_linf_norm
            << ", t_monitor_l2=" << summary.final_monitor_temperature_l2_norm
            << ", t_monitor_linf=" << summary.final_monitor_temperature_linf_norm
            << ", p_vertical_std=" << summary.final_pressure_vertical_uniformity_norm_std
            << ", t_vertical_std=" << summary.final_temperature_vertical_uniformity_norm_std
            << ", grid_ok=" << BoolString(summary.grid_convergence_ok)
            << ", time_ok=" << BoolString(summary.time_sensitivity_ok);
        throw std::runtime_error(oss.str());
    }
    return summary;
}

TestCasePlan BuildDefaultPlan() {
    TestCasePlan plan;
    plan.plan_key = "h_t_co2_constpp_nofrac_nowell";
    return plan;
}

TestCasePlan BuildGridPlan() {
    TestCasePlan plan = BuildDefaultPlan();
    plan.plan_key = "h_t_co2_constpp_nofrac_nowell_grid";
    plan.spec.enable_grid_convergence_study = true;
    return plan;
}

TestCasePlan BuildVaryPPPlan() {
    TestCasePlan plan = BuildDefaultPlan();
    plan.plan_key = "h_t_co2_varypp_nofrac_nowell";
    plan.spec.case_name = "h_t_co2_varypp_nofrac_nowell";
    plan.spec.sub_dir = "H_T_CO2_VaryPP";
    plan.spec.use_variable_properties = true;
    plan.spec.comsol_wrapper_relpath = "tools/COMSOL/VaryPP_NoFrac_NoWell/run_comsol_reference.ps1";
    return plan;
}

TestCasePlan BuildVaryPPSmokePlan() {
    TestCasePlan plan = BuildVaryPPPlan();
    plan.plan_key = "h_t_co2_varypp_nofrac_nowell_smoke";
    plan.spec.case_name = "h_t_co2_varypp_nofrac_nowell_smoke";
    plan.spec.nx = 24;
    plan.spec.ny = 3;
    plan.spec.dt_init = 1.0e3;
    plan.spec.dt_max = 5.0e3;
    plan.spec.target_end_time_s = 5.0e4;
    plan.spec.max_steps = 400;
    return plan;
}

TestCasePlan BuildFastPlan() {
    TestCasePlan plan = BuildDefaultPlan();
    plan.plan_key = "h_t_co2_constpp_nofrac_nowell_fast";
    plan.spec.case_name = "ht_nofrac_fast";
    plan.spec.target_end_time_s = 1.0e7;
    return plan;
}

TestCasePlan BuildFastGridPlan() {
    TestCasePlan plan = BuildFastPlan();
    plan.plan_key = "h_t_co2_constpp_nofrac_nowell_fast_grid";
    plan.spec.enable_grid_convergence_study = true;
    return plan;
}

TestCasePlan BuildDtPlan() {
    TestCasePlan plan = BuildDefaultPlan();
    plan.plan_key = "h_t_co2_constpp_nofrac_nowell_dt";
    plan.spec.enable_time_sensitivity_study = true;
    return plan;
}

TestCasePlan BuildAllPlan() {
    TestCasePlan plan = BuildDefaultPlan();
    plan.plan_key = "h_t_co2_constpp_nofrac_nowell_all";
    plan.spec.enable_grid_convergence_study = true;
    plan.spec.enable_time_sensitivity_study = true;
    return plan;
}

TestCasePlan BuildDebug96x12Plan() {
    TestCasePlan plan = BuildDefaultPlan();
    plan.plan_key = "h_t_co2_constpp_nofrac_nowell_debug_96x12";
    plan.spec.case_name = "ht_nofrac_dbg_96x12";
    plan.spec.comsol_reference_case_name = "h_t_co2_constpp_nofrac_nowell";
    plan.spec.nx = 96;
    plan.spec.ny = 12;
    return plan;
}

TestCasePlan BuildProbe96x12Dt5000Plan() {
    TestCasePlan plan = BuildDefaultPlan();
    plan.plan_key = "h_t_co2_constpp_nofrac_nowell_probe_96x12_dt5000";
    plan.spec.case_name = "ht_nf_p96_d5000";
    plan.spec.comsol_reference_case_name = "h_t_co2_constpp_nofrac_nowell";
    plan.spec.nx = 96;
    plan.spec.ny = 12;
    plan.spec.dt_init = 5.0e3;
    plan.spec.export_vtk = true;
    return plan;
}

using BuilderFn = TestCasePlan(*)();

const std::unordered_map<std::string, BuilderFn>& GetRegistry() {
    static const std::unordered_map<std::string, BuilderFn> registry = {
        {"h_t_co2_constpp_nofrac_nowell", &BuildDefaultPlan},
        {"h_t_co2_varypp_nofrac_nowell", &BuildVaryPPPlan},
        {"h_t_co2_varypp_nofrac_nowell_smoke", &BuildVaryPPSmokePlan},
        {"h_t_co2_constpp_nofrac_nowell_fast", &BuildFastPlan},
        {"h_t_co2_constpp_nofrac_nowell_fast_grid", &BuildFastGridPlan},
        {"h_t_co2_constpp_nofrac_nowell_grid", &BuildGridPlan},
        {"h_t_co2_constpp_nofrac_nowell_dt", &BuildDtPlan},
        {"h_t_co2_constpp_nofrac_nowell_all", &BuildAllPlan},
        {"h_t_co2_constpp_nofrac_nowell_debug_96x12", &BuildDebug96x12Plan},
        {"h_t_co2_constpp_nofrac_nowell_probe_96x12_dt5000", &BuildProbe96x12Dt5000Plan}
    };
    return registry;
}

void MaterializeB1ReferenceInputs(const TestCaseSpec& cfg, const CaseCommon::CaseArtifactPaths& artifacts) {
    TestCaseSummary summary;
    ConfigureSummaryPaths(summary, cfg, artifacts.case_dir);
    EnsureSummaryDirs(summary);
    WritePropertyTables(cfg, summary.property_table_path, summary.comsol_property_table_path);
    WriteReferenceSpec(cfg, summary);
    WriteAnalyticalNote(cfg, summary.analytical_note_path);
    WriteComsolReferenceSpec(cfg, summary.case_dir, summary.comsol_reference_spec_path);
}

void RunStageByKeyImpl(const std::string& key, CaseCommon::CaseStage stage) {
    const auto& registry = GetRegistry();
    const auto it = registry.find(key);
    if (it == registry.end()) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] unknown registry key: " + key);
    }

    const TestCasePlan plan = it->second();
    const std::string caseCode = GetNoFracTemplateCaseCode(plan.spec);
    const CaseCommon::CaseArtifactPaths artifacts = BuildB1ArtifactPaths(plan.spec);
    EnsureB1ArtifactContractDirs(artifacts);
    WriteB1StageManifest(artifacts, caseCode, stage, "started", artifacts.case_dir);
    WriteB1ReferenceContract(plan.spec, artifacts, caseCode, stage);
    WriteB1StageStatus(plan.spec, artifacts, caseCode, stage, "started", artifacts.case_dir);
    MaterializeB1ReferenceInputs(plan.spec, artifacts);

    if (stage == CaseCommon::CaseStage::PrepareReference) {
        WriteB1StageManifest(artifacts, caseCode, stage, "prepared_reference_inputs", artifacts.reference_dir);
        WriteB1StageStatus(plan.spec, artifacts, caseCode, stage, "prepared_reference_inputs", artifacts.reference_dir);
        PrintPrepareReferenceSummary(artifacts);
        return;
    }

    const TestCaseSpec stageSpec = BuildStageSpec(plan.spec, stage);
    const std::string defaultOutputDir =
        (stage == CaseCommon::CaseStage::SolveOnly) ? artifacts.engineering_dir : artifacts.case_dir;
    try {
        TestCaseSummary summary;
        switch (stage) {
        case CaseCommon::CaseStage::SolveOnly:
            summary = RunSingleCaseCore(stageSpec, artifacts.case_dir).summary;
            WriteB1StageManifest(artifacts, caseCode, stage, "completed", artifacts.engineering_dir);
            WriteB1StageStatus(stageSpec, artifacts, caseCode, stage, "completed", artifacts.engineering_dir);
            PrintRunSummary("solve_only completed", summary);
            return;
        case CaseCommon::CaseStage::ValidateOnly:
            summary = RunCase(stageSpec, artifacts.case_dir);
            WriteB1StageManifest(artifacts, caseCode, stage, "completed", artifacts.case_dir);
            WriteB1StageStatus(stageSpec, artifacts, caseCode, stage, "completed", artifacts.case_dir);
            PrintRunSummary("validate_only completed", summary);
            return;
        case CaseCommon::CaseStage::FullWorkflow:
            summary = RunCase(stageSpec, artifacts.case_dir);
            WriteB1StageManifest(artifacts, caseCode, stage, "completed", artifacts.case_dir);
            WriteB1StageStatus(stageSpec, artifacts, caseCode, stage, "completed", artifacts.case_dir);
            PrintRunSummary("full_workflow completed", summary);
            return;
        default:
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_NoFrac] unsupported stage in RunStageByKeyImpl.");
        }
    } catch (const std::exception& ex) {
        const std::string message = ex.what();
        const bool missingReference =
            message.find("missing") != std::string::npos &&
            message.find("reference") != std::string::npos;
        const std::string failureStatus = missingReference ? "missing_reference" : "failed";
        WriteB1StageManifest(artifacts, caseCode, stage, failureStatus, defaultOutputDir);
        WriteB1StageStatus(stageSpec, artifacts, caseCode, stage, failureStatus, defaultOutputDir);
        throw;
    }
}

void ExecutePlanByKeyImpl(const std::string& key) {
    RunStageByKeyImpl(key, CaseCommon::CaseStage::FullWorkflow);
}

} // namespace

void RunTestCase() {
    ExecutePlanByKeyImpl("h_t_co2_constpp_nofrac_nowell");
}

void ExecutePlanByKey(const std::string& key) {
    ExecutePlanByKeyImpl(key);
}

void RunStageByKey(const std::string& key, CaseCommon::CaseStage stage) {
    RunStageByKeyImpl(key, stage);
}

void RunSolveOnly() {
    RunStageByKeyImpl("h_t_co2_constpp_nofrac_nowell", CaseCommon::CaseStage::SolveOnly);
}

void RunPrepareReference() {
    RunStageByKeyImpl("h_t_co2_constpp_nofrac_nowell", CaseCommon::CaseStage::PrepareReference);
}

void RunValidateOnly() {
    RunStageByKeyImpl("h_t_co2_constpp_nofrac_nowell", CaseCommon::CaseStage::ValidateOnly);
}

void RunFullWorkflow() {
    RunStageByKeyImpl("h_t_co2_constpp_nofrac_nowell", CaseCommon::CaseStage::FullWorkflow);
}

} // namespace Test_H_T_CO2_ConstPP_NoFrac

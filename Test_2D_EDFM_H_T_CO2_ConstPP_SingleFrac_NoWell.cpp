/**
 * @file Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.cpp
 * @brief Standalone validation chain: 2D single-phase CO2 constant-property P-T coupled, single-fracture, no-well.
 */

#include "Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.h"

#include "2D_PostProcess.h"
#include "BoundaryConditionManager.h"
#include "FIM_TransientCaseKit.hpp"
#include "Fracture.h"
#include "FractureElement.h"
#include "MeshDefinitions.h"
#include "MeshManager.h"
#include "SolverContrlStrName_op.h"
#include "Well_WellControlTypes.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#define TEST_MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define TEST_MKDIR(path) mkdir(path, 0777)
#endif

namespace Test_H_T_CO2_ConstPP_SingleFrac {
namespace {

constexpr double kPendingTimeTolerance = 1.0e-9;
constexpr double kMonitorSampleDtS = 1.0e5;
constexpr double kGridTimeMonotoneAbsTol = 5.0e-3;
constexpr double kComsolThinBandWidthM = 1.0e-3;
constexpr double kFractureApertureM = 1.0e-3;
constexpr double kFinalTemperatureSpanThresholdK = 10.0;
constexpr double kBoundaryAdjacentDeltaTThresholdK = 20.0;
constexpr double kBoundaryPrescribedTolerance = 1.0e-8;
constexpr const char* kFamilyMatrixHorizontal = "matrix_horizontal";
constexpr const char* kFamilyFractureTangent = "fracture_tangent";
constexpr const char* kFamilyCrossNormal = "cross_normal";
constexpr const char* kLocationMatrix = "matrix";
constexpr const char* kLocationFracture = "fracture";
constexpr const char* kComsolRepresentationExplicit = "explicit_lower_dimensional_fracture";
constexpr const char* kComsolRepresentationThinBand = "thin_band_fallback";

enum class ReferenceMode { Auto, Analytical, Comsol };
enum class WorkflowMode { PrepareReference, ValidateAgainstReference, Full };

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

struct ScalarStats {
    int count = 0;
    double min_value = std::numeric_limits<double>::quiet_NaN();
    double max_value = std::numeric_limits<double>::quiet_NaN();
    double avg_value = std::numeric_limits<double>::quiet_NaN();
};

struct BoundaryDiagnosticsRow {
    std::string snapshot_tag;
    double time_s = 0.0;
    ScalarStats left_prescribed_t;
    ScalarStats right_prescribed_t;
    ScalarStats left_adjacent_t;
    ScalarStats right_adjacent_t;
    ScalarStats domain_t;
};

struct ProfileCompareMetrics {
    std::string family;
    std::string tag;
    double target_time_s = 0.0;
    double actual_time_s = 0.0;
    int sample_count = 0;
    FieldErrorMetrics pressure;
    FieldErrorMetrics temperature;
    std::string compare_csv_path;
};

struct SweepStudyRow {
    std::string label;
    std::string case_dir;
    int nx = 0;
    int ny = 0;
    int steps = 0;
    double dt_init = 0.0;
    double h_char = 0.0;
    double t_end = 0.0;
    double final_profile_p_l2_norm = 0.0;
    double final_profile_t_l2_norm = 0.0;
    double final_monitor_p_l2_norm = 0.0;
    double final_monitor_t_l2_norm = 0.0;
};

struct ProfileReferenceRow {
    int station_id = -1;
    double target_axis_m = 0.0;
    double target_x = 0.0;
    double target_y = 0.0;
    double target_time_s = 0.0;
    double p_ref = 0.0;
    double t_ref = 0.0;
};

struct ProfileReferenceTable {
    std::unordered_map<int, ProfileReferenceRow> rows_by_station_id;
};

struct MonitorReferenceRow {
    int sample_id = -1;
    double target_time_s = 0.0;
    std::unordered_map<std::string, double> p_ref_by_label;
    std::unordered_map<std::string, double> t_ref_by_label;
};

struct MonitorReferenceSeries {
    std::vector<MonitorReferenceRow> rows;
    std::unordered_map<int, std::size_t> sample_id_to_row;
};

struct FractureGeometry {
    Vector start;
    Vector end;
    Vector tangent;
    Vector normal;
    Vector midpoint;
    double length = 0.0;
};

struct FractureBlockSample {
    int block_id = -1;
    double s = 0.0;
    double x = 0.0;
    double y = 0.0;
};

struct TestCaseSpec {
    std::string case_name = "h_t_co2_constpp_singlefrac_nowell";
    std::string output_base_dir = "Test/Transient/FullCaseTest";
    std::string sub_dir = "H_T_CO2_ConstPP";
    double lx = 400.0;
    double ly = 40.0;
    int nx = 48;
    int ny = 6;
    double frac_x0_ratio = 0.2;
    double frac_y0_ratio = 0.1;
    double frac_x1_ratio = 0.8;
    double frac_y1_ratio = 0.9;
    double matrix_phi = 0.10;
    double matrix_perm = 1.0e-13;
    double matrix_ct = 5.0e-9;
    double matrix_rho_r = 2600.0;
    double matrix_cp_r = 800.0;
    double matrix_lambda_r = 3.0;
    double fracture_phi = 0.15;
    double fracture_kt = 5.0e-12;
    double fracture_kn = 5.0e-13;
    double fracture_ct = 5.0e-9;
    double fracture_rho_r = 2600.0;
    double fracture_cp_r = 800.0;
    double fracture_lambda_r = 6.0;
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
    ReferenceMode reference_mode = ReferenceMode::Auto;
    WorkflowMode workflow_mode = WorkflowMode::Full;
    std::vector<double> report_time_fractions = {0.1, 0.5, 1.0};
    std::vector<std::pair<int, int> > grid_sweep_cases = {
        std::make_pair(24, 3),
        std::make_pair(48, 6),
        std::make_pair(96, 12)
    };
    std::vector<double> time_step_sweep = {8.0e4, 2.0e4, 5.0e3};
    int matrix_horizontal_station_count = 81;
    int fracture_tangent_station_count = 81;
    int cross_normal_station_count = 61;
    double cross_normal_half_span_m = 12.0;
    double monitor_dt_s = kMonitorSampleDtS;
    double profile_l2_threshold = 5.0e-2;
    double profile_linf_threshold = 1.0e-1;
    double monitor_l2_threshold = 5.0e-2;
    double monitor_linf_threshold = 1.0e-1;
    bool enable_grid_convergence_study = true;
    bool enable_time_sensitivity_study = true;
    bool export_vtk = true;
    bool emit_detailed_outputs = true;
    bool allow_full_workflow_comsol_autorun = true;
    std::string comsol_wrapper_relpath =
        "tools/COMSOL/H_T_CO2_ConstPP_SingleFrac_NoWell/run_comsol_reference.ps1";
};

struct TestCaseSummary {
    std::string case_dir;
    std::string engineering_dir;
    std::string reference_dir;
    std::string comsol_input_dir;
    std::string comsol_output_dir;
    std::string convergence_log_path;
    std::string run_log_path;
    std::string metrics_csv_path;
    std::string validation_summary_path;
    std::string validation_summary_csv_path;
    std::string analytical_feasibility_path;
    std::string reference_spec_path;
    std::string property_table_path;
    std::string comsol_property_table_path;
    std::string boundary_diagnostics_path;
    std::string profile_station_definitions_path;
    std::string monitor_point_definitions_path;
    std::string profile_schedule_path;
    std::string monitor_schedule_path;
    std::string eng_monitor_timeseries_path;
    std::string comsol_reference_spec_path;
    std::string matlab_script_path;
    std::string grid_convergence_csv_path;
    std::string time_sensitivity_csv_path;
    int nx = 0;
    int ny = 0;
    int n_cells = 0;
    int n_fracture_dofs = 0;
    double h_char = 0.0;
    int steps = 0;
    int total_rollbacks = 0;
    double avg_iters = 0.0;
    int max_iters = 0;
    double t_end = 0.0;
    std::string resolved_reference_mode = "Auto";
    std::string comsol_representation = kComsolRepresentationExplicit;
    bool reference_is_explicit_lower_dimensional = false;
    bool comsol_fine_check_skipped = false;
    bool analytical_closed_form_available = false;
    std::string validation_status = "not_run";
    bool reference_ready = false;
    bool validation_performed = false;
    bool validation_passed = false;
    bool grid_convergence_ok = true;
    bool time_sensitivity_ok = true;
    std::vector<std::string> missing_reference_files;
    double final_profile_p_l2_norm = 0.0;
    double final_profile_p_linf_norm = 0.0;
    double final_profile_t_l2_norm = 0.0;
    double final_profile_t_linf_norm = 0.0;
    double final_monitor_p_l2_norm = 0.0;
    double final_monitor_p_linf_norm = 0.0;
    double final_monitor_t_l2_norm = 0.0;
    double final_monitor_t_linf_norm = 0.0;
    double final_temperature_span_k = 0.0;
    double final_boundary_adjacent_delta_t_k = 0.0;
    bool final_boundary_prescribed_t_ok = false;
    std::vector<ProfileCompareMetrics> report_metrics;
    std::vector<BoundaryDiagnosticsRow> boundary_diagnostics_rows;
};

struct CaseRunArtifacts {
    TestCaseSummary summary;
    std::vector<SpatialSamplePoint> profile_stations;
    std::vector<SpatialSamplePoint> monitor_points;
    std::vector<SnapshotState> snapshots;
    std::vector<MonitorScheduleState> monitor_schedule;
};

struct TestCasePlan {
    std::string plan_key;
    TestCaseSpec spec;
};

void EnsureDirRecursive(const std::string& rawPath) {
    if (rawPath.empty()) return;
    std::string path = rawPath;
    for (char& ch : path) if (ch == '\\') ch = '/';
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

bool FileExists(const std::string& path) {
    std::ifstream in(path.c_str(), std::ios::in);
    return in.good();
}

std::string BoolString(bool value) { return value ? "true" : "false"; }

std::string ReferenceModeString(ReferenceMode mode) {
    switch (mode) {
    case ReferenceMode::Auto: return "Auto";
    case ReferenceMode::Analytical: return "Analytical";
    case ReferenceMode::Comsol: return "Comsol";
    default: return "Unknown";
    }
}

std::string WorkflowModeString(WorkflowMode mode) {
    switch (mode) {
    case WorkflowMode::PrepareReference: return "PrepareReference";
    case WorkflowMode::ValidateAgainstReference: return "ValidateAgainstReference";
    case WorkflowMode::Full: return "Full";
    default: return "Unknown";
    }
}

std::string Trim(std::string value) {
    while (!value.empty() && (value.back() == '\r' || value.back() == '\n' || value.back() == ' ' || value.back() == '\t')) {
        value.pop_back();
    }
    std::size_t first = 0;
    while (first < value.size() && (value[first] == ' ' || value[first] == '\t')) ++first;
    return value.substr(first);
}

std::vector<std::string> SplitCsvLine(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) out.push_back(Trim(token));
    if (!line.empty() && line.back() == ',') out.push_back("");
    return out;
}

struct CsvTable {
    std::vector<std::string> headers;
    std::unordered_map<std::string, std::size_t> header_index;
    std::vector<std::vector<std::string> > rows;
};

CsvTable ReadCsvTable(const std::string& path) {
    std::ifstream in(path.c_str(), std::ios::in);
    if (!in.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to open csv: " + path);
    CsvTable table;
    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] empty csv: " + path);
    table.headers = SplitCsvLine(line);
    for (std::size_t i = 0; i < table.headers.size(); ++i) table.header_index[table.headers[i]] = i;
    while (std::getline(in, line)) {
        if (!Trim(line).empty()) table.rows.push_back(SplitCsvLine(line));
    }
    return table;
}

double CsvGetDouble(const CsvTable& table, std::size_t row, const std::string& column) {
    const auto it = table.header_index.find(column);
    if (it == table.header_index.end()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] missing csv column: " + column);
    if (row >= table.rows.size()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] csv row out of range.");
    const std::size_t col = it->second;
    if (col >= table.rows[row].size()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] csv column out of range.");
    return std::stod(table.rows[row][col]);
}

double ClampFraction(double fraction) { return std::max(0.0, std::min(1.0, fraction)); }

std::string MakeReportTag(double fraction) {
    std::ostringstream oss;
    oss << "t" << std::setw(3) << std::setfill('0')
        << static_cast<int>(std::round(ClampFraction(fraction) * 100.0)) << "pct";
    return oss.str();
}

std::vector<double> BuildSortedFractions(const std::vector<double>& rawFractions) {
    std::vector<double> fractions;
    for (double value : rawFractions) {
        const double clamped = ClampFraction(value);
        bool duplicate = false;
        for (double existing : fractions) {
            if (std::abs(existing - clamped) <= 1.0e-12) { duplicate = true; break; }
        }
        if (!duplicate) fractions.push_back(clamped);
    }
    std::sort(fractions.begin(), fractions.end());
    return fractions;
}

struct RunningStatsAccumulator {
    int count = 0;
    double sum = 0.0;
    double min_value = std::numeric_limits<double>::infinity();
    double max_value = -std::numeric_limits<double>::infinity();

    void add(double value) {
        if (!std::isfinite(value)) return;
        ++count;
        sum += value;
        min_value = std::min(min_value, value);
        max_value = std::max(max_value, value);
    }

    ScalarStats finalize() const {
        ScalarStats stats;
        stats.count = count;
        if (count <= 0) return stats;
        stats.min_value = min_value;
        stats.max_value = max_value;
        stats.avg_value = sum / static_cast<double>(count);
        return stats;
    }
};

void ApplyUniformScalarField(const std::shared_ptr<volScalarField>& field, double value) {
    if (!field) return;
    for (double& v : field->data) v = value;
}

void SyncPTFieldsToFM(MeshManager& mgr, FieldManager_2D& fm,
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

double ComputeMeshCharLength(const MeshManager& mgr) {
    const auto& cells = mgr.mesh().getCells();
    if (cells.empty()) return 0.0;
    double totalV = 0.0;
    for (const auto& c : cells) totalV += std::max(c.volume, 0.0);
    return (totalV > 0.0) ? std::sqrt(totalV / static_cast<double>(cells.size())) : 0.0;
}

double PressureScale(const TestCaseSpec& cfg) {
    return std::max(std::abs(cfg.p_left - cfg.p_right), 1.0);
}

double TemperatureScale(const TestCaseSpec& cfg) {
    return std::max(std::abs(cfg.t_left - cfg.t_right), 1.0);
}

std::string MakeLabel(const std::string& prefix, int index, int width) {
    std::ostringstream oss;
    oss << prefix << std::setw(width) << std::setfill('0') << index;
    return oss.str();
}

FractureGeometry BuildFractureGeometry(const TestCaseSpec& cfg) {
    FractureGeometry geom;
    geom.start = Vector(cfg.frac_x0_ratio * cfg.lx, cfg.frac_y0_ratio * cfg.ly, 0.0);
    geom.end = Vector(cfg.frac_x1_ratio * cfg.lx, cfg.frac_y1_ratio * cfg.ly, 0.0);
    const Vector diff = geom.end - geom.start;
    geom.length = diff.Mag();
    if (geom.length <= 1.0e-20) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] fracture length is zero.");
    }
    geom.tangent = diff / geom.length;
    geom.normal = Vector(-geom.tangent.m_y, geom.tangent.m_x, 0.0);
    geom.midpoint = geom.start + 0.5 * diff;
    return geom;
}

std::vector<FractureBlockSample> BuildFractureBlockSamples(const MeshManager& mgr) {
    const auto& ordered = mgr.fracture_network().getOrderedFractureElements();
    const auto& fractures = mgr.fracture_network().fractures;
    std::vector<FractureBlockSample> out;
    out.reserve(ordered.size());
    for (const FractureElement* elem : ordered) {
        if (!elem) continue;
        if (elem->parentFractureID < 0 || elem->parentFractureID >= static_cast<int>(fractures.size())) continue;
        const Fracture& frac = fractures[static_cast<std::size_t>(elem->parentFractureID)];
        const double s = 0.5 * (elem->param0 + elem->param1);
        const Vector pos = frac.start + (frac.end - frac.start) * s;
        FractureBlockSample sample;
        sample.block_id = elem->solverIndex;
        sample.s = s;
        sample.x = pos.m_x;
        sample.y = pos.m_y;
        out.push_back(sample);
    }
    if (out.empty()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] no fracture block samples available.");
    return out;
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
    if (cells.empty()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] cannot build sample on empty mesh.");
    double bestDistSq = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < cells.size(); ++i) {
        const double dx = cells[i].center.m_x - targetX;
        const double dy = cells[i].center.m_y - targetY;
        const double distSq = dx * dx + dy * dy;
        if (distSq < bestDistSq) {
            bestDistSq = distSq;
            sample.block_id = static_cast<int>(i);
            sample.actual_x = cells[i].center.m_x;
            sample.actual_y = cells[i].center.m_y;
        }
    }
    return sample;
}

SpatialSamplePoint MakeNearestFractureSample(const std::vector<FractureBlockSample>& fracBlocks,
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
    sample.location = kLocationFracture;
    sample.target_axis_m = targetAxisM;
    sample.target_x = targetX;
    sample.target_y = targetY;
    double bestDistSq = std::numeric_limits<double>::max();
    for (const auto& frac : fracBlocks) {
        const double dx = frac.x - targetX;
        const double dy = frac.y - targetY;
        const double distSq = dx * dx + dy * dy;
        if (distSq < bestDistSq) {
            bestDistSq = distSq;
            sample.block_id = frac.block_id;
            sample.actual_x = frac.x;
            sample.actual_y = frac.y;
        }
    }
    if (sample.block_id < 0) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to assign fracture sample block.");
    }
    return sample;
}

std::vector<SpatialSamplePoint> BuildProfileStations(const MeshManager& mgr, const TestCaseSpec& cfg) {
    if (cfg.matrix_horizontal_station_count < 2 ||
        cfg.fracture_tangent_station_count < 2 ||
        cfg.cross_normal_station_count < 2) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] invalid profile station counts.");
    }

    const FractureGeometry geom = BuildFractureGeometry(cfg);
    const std::vector<FractureBlockSample> fracBlocks = BuildFractureBlockSamples(mgr);
    std::vector<SpatialSamplePoint> stations;
    int nextId = 0;

    for (int i = 0; i < cfg.matrix_horizontal_station_count; ++i) {
        const double ratio = static_cast<double>(i) / static_cast<double>(cfg.matrix_horizontal_station_count - 1);
        const double x = ratio * cfg.lx;
        const double y = 0.5 * cfg.ly;
        stations.push_back(MakeNearestMatrixSample(
            mgr, nextId++, MakeLabel("mh", i, 3), kFamilyMatrixHorizontal, x, x, y));
    }

    for (int i = 0; i < cfg.fracture_tangent_station_count; ++i) {
        const double ratio = static_cast<double>(i) / static_cast<double>(cfg.fracture_tangent_station_count - 1);
        const Vector pos = geom.start + ratio * (geom.end - geom.start);
        stations.push_back(MakeNearestFractureSample(
            fracBlocks, nextId++, MakeLabel("ft", i, 3), kFamilyFractureTangent,
            ratio * geom.length, pos.m_x, pos.m_y));
    }

    for (int i = 0; i < cfg.cross_normal_station_count; ++i) {
        const double ratio = static_cast<double>(i) / static_cast<double>(cfg.cross_normal_station_count - 1);
        const double offset = -cfg.cross_normal_half_span_m + 2.0 * cfg.cross_normal_half_span_m * ratio;
        const Vector pos = geom.midpoint + offset * geom.normal;
        if (std::abs(offset) <= 1.0e-12) {
            stations.push_back(MakeNearestFractureSample(
                fracBlocks, nextId++, MakeLabel("cn", i, 3), kFamilyCrossNormal,
                offset, pos.m_x, pos.m_y));
        } else {
            stations.push_back(MakeNearestMatrixSample(
                mgr, nextId++, MakeLabel("cn", i, 3), kFamilyCrossNormal,
                offset, pos.m_x, pos.m_y));
        }
    }
    return stations;
}

std::vector<SpatialSamplePoint> BuildMonitorPoints(const MeshManager& mgr, const TestCaseSpec& cfg) {
    const FractureGeometry geom = BuildFractureGeometry(cfg);
    const std::vector<FractureBlockSample> fracBlocks = BuildFractureBlockSamples(mgr);
    std::vector<SpatialSamplePoint> points;
    int nextId = 0;

    points.push_back(MakeNearestMatrixSample(
        mgr, nextId++, "mx01", "monitor", 0.1 * cfg.lx, 0.1 * cfg.lx, 0.5 * cfg.ly));
    points.push_back(MakeNearestMatrixSample(
        mgr, nextId++, "mx02", "monitor", 0.3 * cfg.lx, 0.3 * cfg.lx, 0.5 * cfg.ly));
    points.push_back(MakeNearestMatrixSample(
        mgr, nextId++, "mx03", "monitor", 0.9 * cfg.lx, 0.9 * cfg.lx, 0.5 * cfg.ly));

    const double fracFractions[3] = {0.2, 0.5, 0.8};
    const char* fracLabels[3] = {"fr01", "fr02", "fr03"};
    for (int i = 0; i < 3; ++i) {
        const Vector pos = geom.start + fracFractions[i] * (geom.end - geom.start);
        points.push_back(MakeNearestFractureSample(
            fracBlocks, nextId++, fracLabels[i], "monitor",
            fracFractions[i] * geom.length, pos.m_x, pos.m_y));
    }

    const Vector posPlus = geom.midpoint + 0.1 * cfg.ly * geom.normal;
    const Vector posMinus = geom.midpoint - 0.1 * cfg.ly * geom.normal;
    points.push_back(MakeNearestMatrixSample(
        mgr, nextId++, "nf_p", "monitor", +0.1 * cfg.ly, posPlus.m_x, posPlus.m_y));
    points.push_back(MakeNearestMatrixSample(
        mgr, nextId++, "nf_m", "monitor", -0.1 * cfg.ly, posMinus.m_x, posMinus.m_y));
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
    for (double timeS = 0.0; timeS <= cfg.target_end_time_s + kPendingTimeTolerance; timeS += cfg.monitor_dt_s) {
        MonitorScheduleState state;
        state.sample_id = sampleId++;
        state.target_time_s = std::min(timeS, cfg.target_end_time_s);
        schedule.push_back(state);
        if (state.target_time_s >= cfg.target_end_time_s - kPendingTimeTolerance) break;
    }
    if (schedule.empty() || schedule.back().target_time_s < cfg.target_end_time_s - kPendingTimeTolerance) {
        MonitorScheduleState finalState;
        finalState.sample_id = sampleId;
        finalState.target_time_s = cfg.target_end_time_s;
        schedule.push_back(finalState);
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

void WriteProfileStationDefinitions(const std::vector<SpatialSamplePoint>& stations, const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write profile station definitions: " + path);
    out << "station_id,label,family,location,target_axis_m,target_x_m,target_y_m,block_id,actual_x_m,actual_y_m\n";
    for (const auto& station : stations) {
        out << station.id << "," << station.label << "," << station.family << "," << station.location << ","
            << std::setprecision(12) << station.target_axis_m << "," << station.target_x << "," << station.target_y << ","
            << station.block_id << "," << station.actual_x << "," << station.actual_y << "\n";
    }
}

void WriteMonitorPointDefinitions(const std::vector<SpatialSamplePoint>& points, const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write monitor point definitions: " + path);
    out << "point_id,label,location,target_axis_m,target_x_m,target_y_m,block_id,actual_x_m,actual_y_m\n";
    for (const auto& point : points) {
        out << point.id << "," << point.label << "," << point.location << ","
            << std::setprecision(12) << point.target_axis_m << "," << point.target_x << "," << point.target_y << ","
            << point.block_id << "," << point.actual_x << "," << point.actual_y << "\n";
    }
}

void WriteProfileSchedule(const std::vector<SnapshotState>& snapshots, const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write profile schedule: " + path);
    out << "tag,requested_fraction,target_time_s,actual_time_s\n";
    for (const auto& snapshot : snapshots) {
        out << snapshot.tag << ","
            << std::setprecision(12) << snapshot.requested_fraction << ","
            << snapshot.target_time_s << "," << snapshot.actual_time_s << "\n";
    }
}

void WriteMonitorSchedule(const std::vector<MonitorScheduleState>& schedule, const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write monitor schedule: " + path);
    out << "sample_id,target_time_s,actual_time_s\n";
    for (const auto& sample : schedule) {
        out << sample.sample_id << ","
            << std::setprecision(12) << sample.target_time_s << "," << sample.actual_time_s << "\n";
    }
}

void WriteEngineeringProfileCsv(const std::vector<SpatialSamplePoint>& stations,
                                const std::string& family,
                                const SnapshotState& snapshot,
                                const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write engineering profile csv: " + path);
    out << "station_id,label,family,location,target_axis_m,target_x_m,target_y_m,block_id,actual_x_m,actual_y_m,"
           "target_time_s,actual_time_s,p_num_pa,t_num_k\n";
    for (const auto& station : stations) {
        if (station.family != family) continue;
        if (station.block_id < 0 ||
            station.block_id >= static_cast<int>(snapshot.p_blocks.size()) ||
            station.block_id >= static_cast<int>(snapshot.t_blocks.size())) {
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] invalid block id during engineering profile export.");
        }
        out << station.id << "," << station.label << "," << station.family << "," << station.location << ","
            << std::setprecision(12) << station.target_axis_m << "," << station.target_x << "," << station.target_y << ","
            << station.block_id << "," << station.actual_x << "," << station.actual_y << ","
            << snapshot.target_time_s << "," << snapshot.actual_time_s << ","
            << snapshot.p_blocks[static_cast<std::size_t>(station.block_id)] << ","
            << snapshot.t_blocks[static_cast<std::size_t>(station.block_id)] << "\n";
    }
}

void WriteEngineeringMonitorCsv(const std::vector<SpatialSamplePoint>& points,
                                const std::vector<MonitorScheduleState>& schedule,
                                const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write engineering monitor csv: " + path);
    out << "sample_id,target_time_s,actual_time_s";
    for (const auto& point : points) out << ",p_num_" << point.label;
    for (const auto& point : points) out << ",t_num_" << point.label;
    out << "\n";
    for (const auto& sample : schedule) {
        out << sample.sample_id << ","
            << std::setprecision(12) << sample.target_time_s << "," << sample.actual_time_s;
        for (const auto& point : points) {
            if (point.block_id < 0 || point.block_id >= static_cast<int>(sample.p_blocks.size())) {
                throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] invalid pressure block id during engineering monitor export.");
            }
            out << "," << sample.p_blocks[static_cast<std::size_t>(point.block_id)];
        }
        for (const auto& point : points) {
            if (point.block_id < 0 || point.block_id >= static_cast<int>(sample.t_blocks.size())) {
                throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] invalid temperature block id during engineering monitor export.");
            }
            out << "," << sample.t_blocks[static_cast<std::size_t>(point.block_id)];
        }
        out << "\n";
    }
}

void WriteAnalyticalFeasibility(const TestCaseSpec& cfg, const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write analytical feasibility note: " + path);
    out << "# Analytical Feasibility Note\n\n";
    out << "- Reference strategy: COMSOL\n";
    out << "- Closed-form analytical solution available: false\n";
    out << "- Reason: the full target case contains coupled pressure-temperature transport and an embedded single fracture.\n";
    out << "- This workflow does not accept pseudo-analytical references derived by embedding the same numerical property evaluation inside an alleged exact solution.\n";
    out << "- Gravity vector: `(" << cfg.gravity_vector.m_x << ", " << cfg.gravity_vector.m_y << ", " << cfg.gravity_vector.m_z << ")`\n";
    out << "- Pressure normalization: `" << PressureScale(cfg) << " Pa`\n";
    out << "- Temperature normalization: `" << TemperatureScale(cfg) << " K`\n";
}

void WritePropertyTables(const TestCaseSpec& cfg,
                         const std::string& engineeringPath,
                         const std::string& comsolInputPath) {
    const auto writeOne = [&](const std::string& path) {
        std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
        if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write property table: " + path);
        out << "key,value,unit\n";
        out << "lx," << std::setprecision(12) << cfg.lx << ",m\n";
        out << "ly," << cfg.ly << ",m\n";
        out << "frac_x0_ratio," << cfg.frac_x0_ratio << ",1\n";
        out << "frac_y0_ratio," << cfg.frac_y0_ratio << ",1\n";
        out << "frac_x1_ratio," << cfg.frac_x1_ratio << ",1\n";
        out << "frac_y1_ratio," << cfg.frac_y1_ratio << ",1\n";
        out << "p_init," << cfg.p_init << ",Pa\n";
        out << "p_left," << cfg.p_left << ",Pa\n";
        out << "p_right," << cfg.p_right << ",Pa\n";
        out << "t_init," << cfg.t_init << ",K\n";
        out << "t_left," << cfg.t_left << ",K\n";
        out << "t_right," << cfg.t_right << ",K\n";
        out << "dt_init," << cfg.dt_init << ",s\n";
        out << "target_end_time_s," << cfg.target_end_time_s << ",s\n";
        out << "matrix_phi," << cfg.matrix_phi << ",1\n";
        out << "matrix_perm," << cfg.matrix_perm << ",m^2\n";
        out << "matrix_ct," << cfg.matrix_ct << ",1/Pa\n";
        out << "matrix_rho_r," << cfg.matrix_rho_r << ",kg/m^3\n";
        out << "matrix_cp_r," << cfg.matrix_cp_r << ",J/(kg*K)\n";
        out << "matrix_lambda_r," << cfg.matrix_lambda_r << ",W/(m*K)\n";
        out << "fracture_phi," << cfg.fracture_phi << ",1\n";
        out << "fracture_kt," << cfg.fracture_kt << ",m^2\n";
        out << "fracture_kn," << cfg.fracture_kn << ",m^2\n";
        out << "fracture_ct," << cfg.fracture_ct << ",1/Pa\n";
        out << "fracture_rho_r," << cfg.fracture_rho_r << ",kg/m^3\n";
        out << "fracture_cp_r," << cfg.fracture_cp_r << ",J/(kg*K)\n";
        out << "fracture_lambda_r," << cfg.fracture_lambda_r << ",W/(m*K)\n";
        out << "co2_rho_const," << cfg.co2_rho_const << ",kg/m^3\n";
        out << "co2_mu_const," << cfg.co2_mu_const << ",Pa*s\n";
        out << "co2_cp_const," << cfg.co2_cp_const << ",J/(kg*K)\n";
        out << "co2_cv_const," << cfg.co2_cv_const << ",J/(kg*K)\n";
        out << "co2_k_const," << cfg.co2_k_const << ",W/(m*K)\n";
        out << "gravity_x," << cfg.gravity_vector.m_x << ",m/s^2\n";
        out << "gravity_y," << cfg.gravity_vector.m_y << ",m/s^2\n";
        out << "gravity_z," << cfg.gravity_vector.m_z << ",m/s^2\n";
        out << "fracture_aperture_m," << kFractureApertureM << ",m\n";
        out << "comsol_thin_band_width," << kComsolThinBandWidthM << ",m\n";
    };
    writeOne(engineeringPath);
    writeOne(comsolInputPath);
}

void WriteComsolReferenceSpec(const TestCaseSpec& cfg, const std::string& caseDir, const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write COMSOL reference spec: " + path);
    out << "# COMSOL Reference Specification\n\n";
    out << "## Case\n";
    out << "- Case directory: `" << caseDir << "`\n";
    out << "- Pressure window: `" << cfg.p_left / 1.0e6 << " MPa` -> `" << cfg.p_right / 1.0e6 << " MPa`\n";
    out << "- Temperature window: `" << cfg.t_left << " K` -> `" << cfg.t_right << " K`\n";
    out << "- Geometry: `" << cfg.lx << " m x " << cfg.ly << " m`\n";
    out << "- Gravity: disabled\n";
    out << "- Non-orthogonal correction: enabled on engineering side\n";
    out << "- Requested fracture route: explicit lower-dimensional fracture\n";
    out << "- Implemented COMSOL route: `" << kComsolRepresentationExplicit << "` with aperture `" << kFractureApertureM << " m`\n\n";
    out << "## Required Inputs\n";
    out << "- `engineering/profile_station_definitions.csv`\n";
    out << "- `engineering/profile_report_schedule.csv`\n";
    out << "- `engineering/monitor_point_definitions.csv`\n";
    out << "- `engineering/monitor_sample_schedule.csv`\n";
    out << "- `reference/comsol_input/property_table.csv`\n\n";
    out << "## Expected COMSOL Outputs\n";
    out << "- `reference/comsol/comsol_monitor_timeseries.csv`\n";
    out << "- `reference/comsol/comsol_model.mph`\n";
    out << "- `reference/comsol/comsol_batch.log`\n";
    out << "- `reference/comsol/comsol_progress.log`\n";
    out << "- `reference/comsol/comsol_java_stdout.log`\n";
    out << "- `reference/comsol/comsol_java_stderr.log`\n\n";
    out << "## Commands\n";
    out << "```powershell\n";
    out << "powershell -ExecutionPolicy Bypass -File " << cfg.comsol_wrapper_relpath << " -Mode Compile\n";
    out << "powershell -ExecutionPolicy Bypass -File " << cfg.comsol_wrapper_relpath << " -Mode Run\n";
    out << "```\n";
}

void WriteStudyCSV(const std::vector<SweepStudyRow>& rows, const std::string& path, const std::string& studyTag) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write study csv: " + path);
    out << "study,label,case_dir,nx,ny,dt_init,h_char,t_end,steps,"
           "final_profile_p_l2_norm,final_profile_t_l2_norm,final_monitor_p_l2_norm,final_monitor_t_l2_norm\n";
    for (const auto& row : rows) {
        out << studyTag << "," << row.label << "," << row.case_dir << ","
            << row.nx << "," << row.ny << ","
            << std::setprecision(12) << row.dt_init << "," << row.h_char << ","
            << row.t_end << "," << row.steps << ","
            << row.final_profile_p_l2_norm << "," << row.final_profile_t_l2_norm << ","
            << row.final_monitor_p_l2_norm << "," << row.final_monitor_t_l2_norm << "\n";
    }
}

bool IsMonotonicNonIncreasingWithTol(const std::vector<double>& values, double absTol) {
    for (std::size_t i = 1; i < values.size(); ++i) {
        if (values[i] > values[i - 1] + absTol) return false;
    }
    return true;
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

std::vector<std::string> OrderedProfileFamilies() {
    return {kFamilyMatrixHorizontal, kFamilyFractureTangent, kFamilyCrossNormal};
}

std::vector<SpatialSamplePoint> FilterStationsByFamily(const std::vector<SpatialSamplePoint>& stations, const std::string& family) {
    std::vector<SpatialSamplePoint> out;
    for (const auto& station : stations) {
        if (station.family == family) out.push_back(station);
    }
    return out;
}

ProfileReferenceTable LoadProfileReference(const std::string& path) {
    const CsvTable table = ReadCsvTable(path);
    ProfileReferenceTable ref;
    for (std::size_t row = 0; row < table.rows.size(); ++row) {
        ProfileReferenceRow item;
        item.station_id = static_cast<int>(std::llround(CsvGetDouble(table, row, "station_id")));
        item.target_axis_m = CsvGetDouble(table, row, "target_axis_m");
        item.target_x = CsvGetDouble(table, row, "target_x_m");
        item.target_y = CsvGetDouble(table, row, "target_y_m");
        item.target_time_s = CsvGetDouble(table, row, "target_time_s");
        item.p_ref = CsvGetDouble(table, row, "p_ref_pa");
        item.t_ref = CsvGetDouble(table, row, "t_ref_k");
        ref.rows_by_station_id[item.station_id] = item;
    }
    return ref;
}

MonitorReferenceSeries LoadMonitorReference(const std::string& path, const std::vector<SpatialSamplePoint>& monitorPoints) {
    const CsvTable table = ReadCsvTable(path);
    MonitorReferenceSeries ref;
    for (std::size_t row = 0; row < table.rows.size(); ++row) {
        MonitorReferenceRow item;
        item.sample_id = static_cast<int>(std::llround(CsvGetDouble(table, row, "sample_id")));
        item.target_time_s = CsvGetDouble(table, row, "target_time_s");
        for (const auto& point : monitorPoints) {
            item.p_ref_by_label[point.label] = CsvGetDouble(table, row, "p_ref_" + point.label);
            item.t_ref_by_label[point.label] = CsvGetDouble(table, row, "t_ref_" + point.label);
        }
        ref.sample_id_to_row[item.sample_id] = ref.rows.size();
        ref.rows.push_back(item);
    }
    return ref;
}

ProfileCompareMetrics EvaluateProfileAgainstReference(const std::vector<SpatialSamplePoint>& stations,
                                                      const SnapshotState& snapshot,
                                                      const ProfileReferenceTable& reference,
                                                      const TestCaseSpec& cfg,
                                                      const std::string& csvPath,
                                                      bool writeCsv) {
    double pSumAbs = 0.0;
    double pSumSq = 0.0;
    double pMaxAbs = 0.0;
    double tSumAbs = 0.0;
    double tSumSq = 0.0;
    double tMaxAbs = 0.0;
    std::ofstream out;
    if (writeCsv) {
        out.open(csvPath.c_str(), std::ios::out | std::ios::trunc);
        if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write profile compare csv: " + csvPath);
        out << "station_id,label,family,location,target_axis_m,target_x_m,target_y_m,block_id,actual_x_m,actual_y_m,"
               "target_time_s,actual_time_s,p_num_pa,p_ref_pa,p_abs_err_pa,p_abs_err_over_dp,"
               "t_num_k,t_ref_k,t_abs_err_k,t_abs_err_over_dt\n";
    }
    for (const auto& station : stations) {
        const auto refIt = reference.rows_by_station_id.find(station.id);
        if (refIt == reference.rows_by_station_id.end()) {
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] missing station id in profile reference.");
        }
        if (station.block_id < 0 ||
            station.block_id >= static_cast<int>(snapshot.p_blocks.size()) ||
            station.block_id >= static_cast<int>(snapshot.t_blocks.size())) {
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] invalid block id during profile comparison.");
        }
        const ProfileReferenceRow& refRow = refIt->second;
        const double pNum = snapshot.p_blocks[static_cast<std::size_t>(station.block_id)];
        const double tNum = snapshot.t_blocks[static_cast<std::size_t>(station.block_id)];
        const double pAbsErr = std::abs(pNum - refRow.p_ref);
        const double tAbsErr = std::abs(tNum - refRow.t_ref);
        pSumAbs += pAbsErr;
        pSumSq += pAbsErr * pAbsErr;
        pMaxAbs = std::max(pMaxAbs, pAbsErr);
        tSumAbs += tAbsErr;
        tSumSq += tAbsErr * tAbsErr;
        tMaxAbs = std::max(tMaxAbs, tAbsErr);
        if (writeCsv) {
            out << station.id << "," << station.label << "," << station.family << "," << station.location << ","
                << std::setprecision(12) << station.target_axis_m << "," << station.target_x << "," << station.target_y << ","
                << station.block_id << "," << station.actual_x << "," << station.actual_y << ","
                << snapshot.target_time_s << "," << snapshot.actual_time_s << ","
                << pNum << "," << refRow.p_ref << "," << pAbsErr << "," << (pAbsErr / PressureScale(cfg)) << ","
                << tNum << "," << refRow.t_ref << "," << tAbsErr << "," << (tAbsErr / TemperatureScale(cfg)) << "\n";
        }
    }
    ProfileCompareMetrics metrics;
    metrics.family = stations.empty() ? "" : stations.front().family;
    metrics.tag = snapshot.tag;
    metrics.target_time_s = snapshot.target_time_s;
    metrics.actual_time_s = snapshot.actual_time_s;
    metrics.sample_count = static_cast<int>(stations.size());
    metrics.pressure = BuildErrorMetrics(pSumAbs, pSumSq, pMaxAbs, metrics.sample_count, PressureScale(cfg));
    metrics.temperature = BuildErrorMetrics(tSumAbs, tSumSq, tMaxAbs, metrics.sample_count, TemperatureScale(cfg));
    metrics.compare_csv_path = writeCsv ? csvPath : "";
    return metrics;
}

FieldErrorMetrics EvaluateFinalMonitorField(const std::vector<SpatialSamplePoint>& monitorPoints,
                                            const MonitorScheduleState& sample,
                                            const MonitorReferenceRow& refRow,
                                            const std::vector<double>& values,
                                            bool usePressure,
                                            double scale) {
    double sumAbs = 0.0;
    double sumSq = 0.0;
    double maxAbs = 0.0;
    for (const auto& point : monitorPoints) {
        if (point.block_id < 0 || point.block_id >= static_cast<int>(values.size())) {
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] invalid block id during monitor evaluation.");
        }
        const auto& refMap = usePressure ? refRow.p_ref_by_label : refRow.t_ref_by_label;
        const auto refIt = refMap.find(point.label);
        if (refIt == refMap.end()) {
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] missing monitor label in reference.");
        }
        const double err = std::abs(values[static_cast<std::size_t>(point.block_id)] - refIt->second);
        sumAbs += err;
        sumSq += err * err;
        maxAbs = std::max(maxAbs, err);
    }
    return BuildErrorMetrics(sumAbs, sumSq, maxAbs, static_cast<int>(monitorPoints.size()), scale);
}

FieldErrorMetrics WriteMonitorCompareCsvRow(std::ofstream& out,
                                            const std::vector<SpatialSamplePoint>& monitorPoints,
                                            const MonitorScheduleState& sample,
                                            const MonitorReferenceRow& refRow,
                                            const std::vector<double>& values,
                                            bool usePressure,
                                            double scale) {
    double sumAbs = 0.0;
    double sumSq = 0.0;
    double maxAbs = 0.0;
    for (const auto& point : monitorPoints) {
        const auto& refMap = usePressure ? refRow.p_ref_by_label : refRow.t_ref_by_label;
        const auto refIt = refMap.find(point.label);
        if (refIt == refMap.end()) {
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] missing monitor label in compare row.");
        }
        if (point.block_id < 0 || point.block_id >= static_cast<int>(values.size())) {
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] invalid block id during monitor compare row.");
        }
        const double numValue = values[static_cast<std::size_t>(point.block_id)];
        const double absErr = std::abs(numValue - refIt->second);
        sumAbs += absErr;
        sumSq += absErr * absErr;
        maxAbs = std::max(maxAbs, absErr);
        out << "," << numValue << "," << refIt->second << "," << absErr << "," << (absErr / scale);
    }
    return BuildErrorMetrics(sumAbs, sumSq, maxAbs, static_cast<int>(monitorPoints.size()), scale);
}

FieldErrorMetrics WriteMonitorCompareCsv(const std::vector<SpatialSamplePoint>& monitorPoints,
                                         const std::vector<MonitorScheduleState>& schedule,
                                         const MonitorReferenceSeries& reference,
                                         const TestCaseSpec& cfg,
                                         const std::string& path,
                                         FieldErrorMetrics& tMetricsOut) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write monitor compare csv: " + path);
    out << "sample_id,target_time_s,actual_time_s";
    for (const auto& point : monitorPoints) {
        out << ",p_num_" << point.label << ",p_ref_" << point.label
            << ",p_abs_err_" << point.label << ",p_abs_err_over_dP_" << point.label;
    }
    for (const auto& point : monitorPoints) {
        out << ",t_num_" << point.label << ",t_ref_" << point.label
            << ",t_abs_err_" << point.label << ",t_abs_err_over_dT_" << point.label;
    }
    out << ",p_l1_norm,p_l2_norm,p_linf_norm,t_l1_norm,t_l2_norm,t_linf_norm\n";

    FieldErrorMetrics finalPMetrics;
    FieldErrorMetrics finalTMetrics;
    for (const auto& sample : schedule) {
        const auto refIt = reference.sample_id_to_row.find(sample.sample_id);
        if (refIt == reference.sample_id_to_row.end()) {
            throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] missing monitor sample in reference.");
        }
        const MonitorReferenceRow& refRow = reference.rows[refIt->second];
        out << sample.sample_id << "," << std::setprecision(12) << sample.target_time_s << "," << sample.actual_time_s;
        const FieldErrorMetrics pMetrics = WriteMonitorCompareCsvRow(
            out, monitorPoints, sample, refRow, sample.p_blocks, true, PressureScale(cfg));
        const FieldErrorMetrics tMetrics = WriteMonitorCompareCsvRow(
            out, monitorPoints, sample, refRow, sample.t_blocks, false, TemperatureScale(cfg));
        out << "," << pMetrics.l1_norm << "," << pMetrics.l2_norm << "," << pMetrics.linf_norm
            << "," << tMetrics.l1_norm << "," << tMetrics.l2_norm << "," << tMetrics.linf_norm << "\n";
        finalPMetrics = pMetrics;
        finalTMetrics = tMetrics;
    }
    tMetricsOut = finalTMetrics;
    return finalPMetrics;
}

std::vector<std::string> RequiredReferenceFiles(const TestCaseSummary& summary, const TestCaseSpec& cfg) {
    std::vector<std::string> files;
    for (const std::string& family : OrderedProfileFamilies()) {
        const std::vector<double> fractions = BuildSortedFractions(cfg.report_time_fractions);
        for (double fraction : fractions) {
            files.push_back(summary.comsol_output_dir + "/comsol_profile_" + family + "_" + MakeReportTag(fraction) + ".csv");
        }
    }
    files.push_back(summary.comsol_output_dir + "/comsol_monitor_timeseries.csv");
    files.push_back(summary.comsol_output_dir + "/comsol_run_summary.md");
    files.push_back(summary.comsol_output_dir + "/comsol_reference_mesh_check.txt");
    return files;
}

bool ReferenceFilesReady(TestCaseSummary& summary, const TestCaseSpec& cfg) {
    summary.missing_reference_files.clear();
    const std::vector<std::string> required = RequiredReferenceFiles(summary, cfg);
    for (const auto& path : required) {
        if (!FileExists(path)) summary.missing_reference_files.push_back(path);
    }
    summary.reference_ready = summary.missing_reference_files.empty();
    return summary.reference_ready;
}

void WriteMetricsCsv(const TestCaseSpec& cfg, const TestCaseSummary& summary) {
    std::ofstream metrics(summary.metrics_csv_path.c_str(), std::ios::out | std::ios::trunc);
    if (!metrics.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write metrics csv: " + summary.metrics_csv_path);
    metrics << "case_name,nx,ny,n_cells,h_char,t_end,steps,total_rollbacks,avg_nonlinear_iters,max_nonlinear_iters,"
               "n_fracture_dofs,dt_init,dt_min,dt_max,property_mode,solver_route,reference_mode_requested,"
               "reference_mode_resolved,workflow_mode,validation_status,"
               "reference_ready,validation_performed,validation_passed,grid_convergence_ok,time_sensitivity_ok,"
               "final_profile_p_l2_norm,final_profile_p_linf_norm,final_profile_t_l2_norm,final_profile_t_linf_norm,"
               "final_monitor_p_l2_norm,final_monitor_p_linf_norm,final_monitor_t_l2_norm,final_monitor_t_linf_norm,"
               "comsol_representation,reference_is_explicit_lower_dimensional,comsol_fine_check_skipped,"
               "final_temperature_span_k,final_boundary_adjacent_delta_t_k,final_boundary_prescribed_t_ok,"
               "profile_l2_threshold,profile_linf_threshold,monitor_l2_threshold,monitor_linf_threshold\n";
    metrics << cfg.case_name << "," << cfg.nx << "," << cfg.ny << "," << summary.n_cells << ","
            << std::setprecision(12) << summary.h_char << "," << summary.t_end << ","
            << summary.steps << "," << summary.total_rollbacks << "," << summary.avg_iters << "," << summary.max_iters << ","
            << summary.n_fracture_dofs << "," << cfg.dt_init << "," << cfg.dt_min << "," << cfg.dt_max << ","
            << "ConstantSinglePhaseCO2,RunGenericFIMTransient<2>,"
            << ReferenceModeString(cfg.reference_mode) << "," << summary.resolved_reference_mode << ","
            << WorkflowModeString(cfg.workflow_mode) << ","
            << summary.validation_status << "," << BoolString(summary.reference_ready) << ","
            << BoolString(summary.validation_performed) << "," << BoolString(summary.validation_passed) << ","
            << BoolString(summary.grid_convergence_ok) << "," << BoolString(summary.time_sensitivity_ok) << ","
            << summary.final_profile_p_l2_norm << "," << summary.final_profile_p_linf_norm << ","
            << summary.final_profile_t_l2_norm << "," << summary.final_profile_t_linf_norm << ","
            << summary.final_monitor_p_l2_norm << "," << summary.final_monitor_p_linf_norm << ","
            << summary.final_monitor_t_l2_norm << "," << summary.final_monitor_t_linf_norm << ","
            << summary.comsol_representation << "," << BoolString(summary.reference_is_explicit_lower_dimensional) << ","
            << BoolString(summary.comsol_fine_check_skipped) << ","
            << summary.final_temperature_span_k << "," << summary.final_boundary_adjacent_delta_t_k << ","
            << BoolString(summary.final_boundary_prescribed_t_ok) << ","
            << cfg.profile_l2_threshold << "," << cfg.profile_linf_threshold << ","
            << cfg.monitor_l2_threshold << "," << cfg.monitor_linf_threshold << "\n";
}

void WriteReferenceSpec(const TestCaseSpec& cfg, const TestCaseSummary& summary) {
    std::ofstream out(summary.reference_spec_path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write reference spec: " + summary.reference_spec_path);
    out << "# Reference Specification\n\n";
    out << "## Minimum Physical Alignment\n";
    out << "- Gravity: disabled (`0, 0, 0`).\n";
    out << "- Non-orthogonal correction: " << BoolString(cfg.enable_non_orthogonal_correction) << ".\n";
    out << "- Fluid model: single-phase CO2, constant properties.\n";
    out << "- Geometry: `" << cfg.lx << " m x " << cfg.ly << " m`.\n";
    out << "- Fracture: single diagonal segment from `(" << cfg.frac_x0_ratio * cfg.lx << ", " << cfg.frac_y0_ratio * cfg.ly
        << ")` to `(" << cfg.frac_x1_ratio * cfg.lx << ", " << cfg.frac_y1_ratio * cfg.ly << ")` m.\n";
    out << "- Wells / sources: none.\n";
    out << "- Pressure BC: left Dirichlet `" << cfg.p_left << " Pa`, right Dirichlet `" << cfg.p_right << " Pa`, top/bottom Neumann `0`.\n";
    out << "- Temperature BC: left Dirichlet `" << cfg.t_left << " K`, right Dirichlet `" << cfg.t_right << " K`, top/bottom Neumann `0`.\n";
    out << "- Initial state: `p=" << cfg.p_init << " Pa`, T=" << cfg.t_init << " K`.\n";
    out << "- Pressure normalization: `Delta p = " << PressureScale(cfg) << " Pa`.\n";
    out << "- Temperature normalization: `Delta T = " << TemperatureScale(cfg) << " K`.\n\n";
    out << "## Sampling\n";
    out << "- Profiles: `" << kFamilyMatrixHorizontal << "`, `" << kFamilyFractureTangent << "`, `" << kFamilyCrossNormal << "`.\n";
    out << "- Profile times: `t = 0.1, 0.5, 1.0 * T_end`.\n";
    out << "- Monitor sample period: `" << cfg.monitor_dt_s << " s`.\n";
    out << "- Engineering directory: `" << summary.engineering_dir << "`.\n";
    out << "- Reference directory: `" << summary.reference_dir << "`.\n";
    out << "- Boundary diagnostics: `" << summary.boundary_diagnostics_path << "`.\n";
}

void WriteValidationSummary(const TestCaseSpec& cfg, const TestCaseSummary& summary) {
    std::ofstream out(summary.validation_summary_path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write validation summary: " + summary.validation_summary_path);
    out << "# Validation Summary\n\n";
    out << "## Status\n";
    out << "- Case: `" << cfg.case_name << "`\n";
    out << "- Output directory: `" << summary.case_dir << "`\n";
    out << "- Requested reference mode: `" << ReferenceModeString(cfg.reference_mode) << "`\n";
    out << "- Resolved reference mode: `" << summary.resolved_reference_mode << "`\n";
    out << "- Workflow mode: `" << WorkflowModeString(cfg.workflow_mode) << "`\n";
    out << "- Validation status: `" << summary.validation_status << "`\n";
    out << "- Reference ready: `" << BoolString(summary.reference_ready) << "`\n";
    out << "- Validation performed: `" << BoolString(summary.validation_performed) << "`\n";
    out << "- Validation passed: `" << BoolString(summary.validation_passed) << "`\n";
    out << "- Grid convergence ok: `" << BoolString(summary.grid_convergence_ok) << "`\n";
    out << "- Time sensitivity ok: `" << BoolString(summary.time_sensitivity_ok) << "`\n\n";

    out << "## Physical Alignment\n";
    out << "- Gravity vector: `(" << cfg.gravity_vector.m_x << ", " << cfg.gravity_vector.m_y << ", " << cfg.gravity_vector.m_z << ")`\n";
    out << "- Non-orthogonal correction: `" << BoolString(cfg.enable_non_orthogonal_correction) << "`\n";
    out << "- COMSOL representation: `" << summary.comsol_representation << "`\n";
    out << "- Reference uses explicit lower-dimensional fracture: `" << BoolString(summary.reference_is_explicit_lower_dimensional) << "`\n";
    out << "- COMSOL fine mesh check skipped: `" << BoolString(summary.comsol_fine_check_skipped) << "`\n";
    out << "- Analytical closed form available: `" << BoolString(summary.analytical_closed_form_available) << "`\n\n";

    out << "## Final Acceptance Metrics\n";
    out << "- Final profile pressure `L2_norm`: `" << summary.final_profile_p_l2_norm << "`\n";
    out << "- Final profile pressure `Linf_norm`: `" << summary.final_profile_p_linf_norm << "`\n";
    out << "- Final profile temperature `L2_norm`: `" << summary.final_profile_t_l2_norm << "`\n";
    out << "- Final profile temperature `Linf_norm`: `" << summary.final_profile_t_linf_norm << "`\n";
    out << "- Final monitor pressure `L2_norm`: `" << summary.final_monitor_p_l2_norm << "`\n";
    out << "- Final monitor pressure `Linf_norm`: `" << summary.final_monitor_p_linf_norm << "`\n";
    out << "- Final monitor temperature `L2_norm`: `" << summary.final_monitor_t_l2_norm << "`\n";
    out << "- Final monitor temperature `Linf_norm`: `" << summary.final_monitor_t_linf_norm << "`\n";
    out << "- Final temperature span `Tmax-Tmin`: `" << summary.final_temperature_span_k << " K`\n";
    out << "- Final left/right adjacent-cell average temperature difference: `" << summary.final_boundary_adjacent_delta_t_k << " K`\n";
    out << "- Prescribed boundary temperature consistency: `" << BoolString(summary.final_boundary_prescribed_t_ok) << "`\n";
    out << "- Profile thresholds: `L2 <= " << cfg.profile_l2_threshold << "`, `Linf <= " << cfg.profile_linf_threshold << "`\n";
    out << "- Monitor thresholds: `L2 <= " << cfg.monitor_l2_threshold << "`, `Linf <= " << cfg.monitor_linf_threshold << "`\n\n";

    out << "## Outputs\n";
    out << "- Engineering dir: `" << summary.engineering_dir << "`\n";
    out << "- Reference dir: `" << summary.reference_dir << "`\n";
    out << "- Property table: `" << summary.property_table_path << "`\n";
    out << "- COMSOL property table: `" << summary.comsol_property_table_path << "`\n";
    out << "- Reference spec: `" << summary.reference_spec_path << "`\n";
    out << "- Analytical feasibility note: `" << summary.analytical_feasibility_path << "`\n";
    out << "- COMSOL reference spec: `" << summary.comsol_reference_spec_path << "`\n";
    out << "- MATLAB script: `" << summary.matlab_script_path << "`\n";
    out << "- Boundary diagnostics csv: `" << summary.boundary_diagnostics_path << "`\n";
    if (!summary.grid_convergence_csv_path.empty()) out << "- Grid convergence csv: `" << summary.grid_convergence_csv_path << "`\n";
    if (!summary.time_sensitivity_csv_path.empty()) out << "- Time sensitivity csv: `" << summary.time_sensitivity_csv_path << "`\n";

    if (!summary.report_metrics.empty()) {
        out << "\n## Profile Compare Records\n";
        for (const auto& metric : summary.report_metrics) {
            out << "- `" << metric.family << " / " << metric.tag
                << "`: pL2=`" << metric.pressure.l2_norm
                << "`, pLinf=`" << metric.pressure.linf_norm
                << "`, tL2=`" << metric.temperature.l2_norm
                << "`, tLinf=`" << metric.temperature.linf_norm
                << "`, csv=`" << metric.compare_csv_path << "`\n";
        }
    }

    if (!summary.missing_reference_files.empty()) {
        out << "\n## Missing Reference Files\n";
        for (const auto& path : summary.missing_reference_files) out << "- `" << path << "`\n";
    }
}

ScalarStats ComputeFieldStats(const std::shared_ptr<volScalarField>& field) {
    RunningStatsAccumulator acc;
    if (!field) return acc.finalize();
    for (double value : field->data) acc.add(value);
    return acc.finalize();
}

ScalarStats ComputeBoundaryPrescribedStats(const MeshManager& mgr,
                                           const BoundarySetting::BoundaryConditionManager& bcMgr,
                                           int boundaryTag) {
    RunningStatsAccumulator acc;
    const auto& faces = mgr.mesh().getFaces();
    for (const auto& face : faces) {
        if (!face.isBoundary() || face.physicalGroupId != boundaryTag) continue;
        if (!bcMgr.HasBC(boundaryTag)) continue;
        const BoundarySetting::BCCoefficients bc = bcMgr.GetBCCoefficients(boundaryTag, face.midpoint);
        if (bc.type != BoundarySetting::BoundaryType::Dirichlet) continue;
        acc.add(bc.c);
    }
    return acc.finalize();
}

ScalarStats ComputeBoundaryAdjacentCellStats(const MeshManager& mgr,
                                             const std::shared_ptr<volScalarField>& field,
                                             int boundaryTag) {
    RunningStatsAccumulator acc;
    if (!field) return acc.finalize();
    const auto& faces = mgr.mesh().getFaces();
    for (const auto& face : faces) {
        if (!face.isBoundary() || face.physicalGroupId != boundaryTag) continue;
        const int owner = face.ownerCell_index;
        if (owner < 0 || owner >= static_cast<int>(field->data.size())) continue;
        acc.add(field->data[static_cast<std::size_t>(owner)]);
    }
    return acc.finalize();
}

BoundaryDiagnosticsRow BuildBoundaryDiagnosticsRow(const MeshManager& mgr,
                                                   const FieldManager_2D& fm,
                                                   const BoundarySetting::BoundaryConditionManager& bcT,
                                                   const std::string& snapshotTag,
                                                   double timeS) {
    BoundaryDiagnosticsRow row;
    row.snapshot_tag = snapshotTag;
    row.time_s = timeS;
    const auto tField = fm.getMatrixScalar("T");
    row.left_prescribed_t = ComputeBoundaryPrescribedStats(mgr, bcT, MeshTags::LEFT);
    row.right_prescribed_t = ComputeBoundaryPrescribedStats(mgr, bcT, MeshTags::RIGHT);
    row.left_adjacent_t = ComputeBoundaryAdjacentCellStats(mgr, tField, MeshTags::LEFT);
    row.right_adjacent_t = ComputeBoundaryAdjacentCellStats(mgr, tField, MeshTags::RIGHT);
    row.domain_t = ComputeFieldStats(tField);
    return row;
}

void WriteBoundaryDiagnosticsCsv(const std::vector<BoundaryDiagnosticsRow>& rows, const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write boundary diagnostics csv: " + path);
    }
    out << "snapshot_tag,time_s,"
           "left_prescribed_t_count,left_prescribed_t_min_k,left_prescribed_t_max_k,left_prescribed_t_avg_k,"
           "right_prescribed_t_count,right_prescribed_t_min_k,right_prescribed_t_max_k,right_prescribed_t_avg_k,"
           "left_adjacent_t_count,left_adjacent_t_min_k,left_adjacent_t_max_k,left_adjacent_t_avg_k,"
           "right_adjacent_t_count,right_adjacent_t_min_k,right_adjacent_t_max_k,right_adjacent_t_avg_k,"
           "domain_t_count,domain_t_min_k,domain_t_max_k,domain_t_avg_k,domain_t_span_k\n";
    for (const auto& row : rows) {
        const double span =
            (row.domain_t.count > 0 && std::isfinite(row.domain_t.max_value) && std::isfinite(row.domain_t.min_value))
            ? (row.domain_t.max_value - row.domain_t.min_value)
            : std::numeric_limits<double>::quiet_NaN();
        out << row.snapshot_tag << "," << std::setprecision(12) << row.time_s << ","
            << row.left_prescribed_t.count << "," << row.left_prescribed_t.min_value << "," << row.left_prescribed_t.max_value
            << "," << row.left_prescribed_t.avg_value << ","
            << row.right_prescribed_t.count << "," << row.right_prescribed_t.min_value << "," << row.right_prescribed_t.max_value
            << "," << row.right_prescribed_t.avg_value << ","
            << row.left_adjacent_t.count << "," << row.left_adjacent_t.min_value << "," << row.left_adjacent_t.max_value
            << "," << row.left_adjacent_t.avg_value << ","
            << row.right_adjacent_t.count << "," << row.right_adjacent_t.min_value << "," << row.right_adjacent_t.max_value
            << "," << row.right_adjacent_t.avg_value << ","
            << row.domain_t.count << "," << row.domain_t.min_value << "," << row.domain_t.max_value << ","
            << row.domain_t.avg_value << "," << span << "\n";
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
    const std::string token = "Fracture representation:";
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

bool StatsWithinTolerance(const ScalarStats& stats, double expected, double tol) {
    return stats.count > 0 &&
           std::isfinite(stats.min_value) &&
           std::isfinite(stats.max_value) &&
           std::isfinite(stats.avg_value) &&
           std::abs(stats.min_value - expected) <= tol &&
           std::abs(stats.max_value - expected) <= tol &&
           std::abs(stats.avg_value - expected) <= tol;
}

void WriteValidationSummaryCsv(const TestCaseSpec& cfg, const TestCaseSummary& summary) {
    std::ofstream out(summary.validation_summary_csv_path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write validation summary csv: " + summary.validation_summary_csv_path);
    out << "case_name,reference_mode_requested,reference_mode_resolved,workflow_mode,validation_status,"
           "reference_ready,validation_performed,validation_passed,grid_convergence_ok,time_sensitivity_ok,"
           "final_profile_p_l2_norm,final_profile_p_linf_norm,final_profile_t_l2_norm,final_profile_t_linf_norm,"
           "final_monitor_p_l2_norm,final_monitor_p_linf_norm,final_monitor_t_l2_norm,final_monitor_t_linf_norm,"
           "comsol_representation,reference_is_explicit_lower_dimensional,comsol_fine_check_skipped,"
           "final_temperature_span_k,final_boundary_adjacent_delta_t_k,final_boundary_prescribed_t_ok\n";
    out << cfg.case_name << "," << ReferenceModeString(cfg.reference_mode) << "," << summary.resolved_reference_mode << ","
        << WorkflowModeString(cfg.workflow_mode) << "," << summary.validation_status << ","
        << BoolString(summary.reference_ready) << "," << BoolString(summary.validation_performed) << ","
        << BoolString(summary.validation_passed) << "," << BoolString(summary.grid_convergence_ok) << ","
        << BoolString(summary.time_sensitivity_ok) << ","
        << summary.final_profile_p_l2_norm << "," << summary.final_profile_p_linf_norm << ","
        << summary.final_profile_t_l2_norm << "," << summary.final_profile_t_linf_norm << ","
        << summary.final_monitor_p_l2_norm << "," << summary.final_monitor_p_linf_norm << ","
        << summary.final_monitor_t_l2_norm << "," << summary.final_monitor_t_linf_norm << ","
        << summary.comsol_representation << "," << BoolString(summary.reference_is_explicit_lower_dimensional) << ","
        << BoolString(summary.comsol_fine_check_skipped) << ","
        << summary.final_temperature_span_k << "," << summary.final_boundary_adjacent_delta_t_k << ","
        << BoolString(summary.final_boundary_prescribed_t_ok) << "\n";
}

void WriteMatlabPlotScript(const TestCaseSummary& summary) {
    std::ofstream out(summary.matlab_script_path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to write MATLAB plot script: " + summary.matlab_script_path);
    out <<
"rootDir = fileparts(mfilename('fullpath'));\n"
"if isempty(rootDir)\n"
"    rootDir = pwd;\n"
"end\n"
"figDir = fullfile(rootDir, 'figures');\n"
"if ~exist(figDir, 'dir')\n"
"    mkdir(figDir);\n"
"end\n"
"\n"
"requiredRoot = {\n"
"    fullfile(rootDir, 'validation_summary.csv')\n"
"    fullfile(rootDir, 'grid_convergence.csv')\n"
"    fullfile(rootDir, 'time_sensitivity.csv')\n"
"    fullfile(rootDir, 'compare_monitor_timeseries.csv')\n"
"    fullfile(rootDir, 'boundary_diagnostics.csv')\n"
"};\n"
"for i = 1:numel(requiredRoot)\n"
"    if ~isfile(requiredRoot{i})\n"
"        error('Missing required input file: %s', requiredRoot{i});\n"
"    end\n"
"end\n"
"\n"
"families = {'matrix_horizontal', 'fracture_tangent', 'cross_normal'};\n"
"tags = {'t010pct', 't050pct', 't100pct'};\n"
"for iFam = 1:numel(families)\n"
"    fam = families{iFam};\n"
"    for iTag = 1:numel(tags)\n"
"        cmpFile = fullfile(rootDir, ['compare_profile_' fam '_' tags{iTag} '.csv']);\n"
"        if ~isfile(cmpFile)\n"
"            error('Missing required profile compare file: %s', cmpFile);\n"
"        end\n"
"    end\n"
"end\n"
"\n"
"for iFam = 1:numel(families)\n"
"    fam = families{iFam};\n"
"    cmpFile = fullfile(rootDir, ['compare_profile_' fam '_t100pct.csv']);\n"
"    tbl = readtable(cmpFile);\n"
"    f = figure('Color', 'w', 'Position', [100 100 1200 420]);\n"
"    subplot(1,2,1);\n"
"    plot(tbl.target_axis_m, tbl.p_num_pa, 'LineWidth', 1.6); hold on;\n"
"    plot(tbl.target_axis_m, tbl.p_ref_pa, '--', 'LineWidth', 1.6);\n"
"    xlabel('Target Axis (m)'); ylabel('Pressure (Pa)');\n"
"    title(['Profile Compare: ' strrep(fam, '_', ' ') ' (Final)']);\n"
"    legend({'Engineering', 'Reference'}, 'Location', 'best'); grid on; box on;\n"
"    subplot(1,2,2);\n"
"    plot(tbl.target_axis_m, tbl.t_num_k, 'LineWidth', 1.6); hold on;\n"
"    plot(tbl.target_axis_m, tbl.t_ref_k, '--', 'LineWidth', 1.6);\n"
"    xlabel('Target Axis (m)'); ylabel('Temperature (K)');\n"
"    title(['Temperature Compare: ' strrep(fam, '_', ' ') ' (Final)']);\n"
"    legend({'Engineering', 'Reference'}, 'Location', 'best'); grid on; box on;\n"
"    exportgraphics(f, fullfile(figDir, ['profile_compare_' fam '.pdf']), 'ContentType', 'vector');\n"
"    exportgraphics(f, fullfile(figDir, ['profile_compare_' fam '.png']), 'Resolution', 300);\n"
"    close(f);\n"
"end\n"
"\n"
"mon = readtable(fullfile(rootDir, 'compare_monitor_timeseries.csv'));\n"
"labelsP = mon.Properties.VariableNames(startsWith(mon.Properties.VariableNames, 'p_num_'));\n"
"labelsT = mon.Properties.VariableNames(startsWith(mon.Properties.VariableNames, 't_num_'));\n"
"f = figure('Color', 'w', 'Position', [120 120 1200 460]);\n"
"subplot(1,2,1); hold on;\n"
"for i = 1:numel(labelsP)\n"
"    refName = strrep(labelsP{i}, 'p_num_', 'p_ref_');\n"
"    plot(mon.target_time_s, mon.(labelsP{i}), 'LineWidth', 1.2);\n"
"    plot(mon.target_time_s, mon.(refName), '--', 'LineWidth', 1.0);\n"
"end\n"
"xlabel('Time (s)'); ylabel('Pressure (Pa)'); title('Monitor Compare: Pressure'); grid on; box on;\n"
"subplot(1,2,2); hold on;\n"
"for i = 1:numel(labelsT)\n"
"    refName = strrep(labelsT{i}, 't_num_', 't_ref_');\n"
"    plot(mon.target_time_s, mon.(labelsT{i}), 'LineWidth', 1.2);\n"
"    plot(mon.target_time_s, mon.(refName), '--', 'LineWidth', 1.0);\n"
"end\n"
"xlabel('Time (s)'); ylabel('Temperature (K)'); title('Monitor Compare: Temperature'); grid on; box on;\n"
"exportgraphics(f, fullfile(figDir, 'monitor_compare.pdf'), 'ContentType', 'vector');\n"
"exportgraphics(f, fullfile(figDir, 'monitor_compare.png'), 'Resolution', 300);\n"
"close(f);\n"
"\n"
"gridTbl = readtable(fullfile(rootDir, 'grid_convergence.csv'));\n"
"f = figure('Color', 'w', 'Position', [140 140 1000 420]);\n"
"subplot(1,2,1);\n"
"plot(gridTbl.h_char, gridTbl.final_profile_p_l2_norm, '-o', 'LineWidth', 1.6); hold on;\n"
"plot(gridTbl.h_char, gridTbl.final_profile_t_l2_norm, '-s', 'LineWidth', 1.6);\n"
"set(gca, 'XDir', 'reverse'); xlabel('Characteristic Grid Size'); ylabel('Normalized L2');\n"
"title('Grid Convergence: Profile'); legend({'Pressure', 'Temperature'}, 'Location', 'best'); grid on; box on;\n"
"subplot(1,2,2);\n"
"plot(gridTbl.h_char, gridTbl.final_monitor_p_l2_norm, '-o', 'LineWidth', 1.6); hold on;\n"
"plot(gridTbl.h_char, gridTbl.final_monitor_t_l2_norm, '-s', 'LineWidth', 1.6);\n"
"set(gca, 'XDir', 'reverse'); xlabel('Characteristic Grid Size'); ylabel('Normalized L2');\n"
"title('Grid Convergence: Monitor'); legend({'Pressure', 'Temperature'}, 'Location', 'best'); grid on; box on;\n"
"exportgraphics(f, fullfile(figDir, 'grid_convergence.pdf'), 'ContentType', 'vector');\n"
"exportgraphics(f, fullfile(figDir, 'grid_convergence.png'), 'Resolution', 300);\n"
"close(f);\n"
"\n"
"timeTbl = readtable(fullfile(rootDir, 'time_sensitivity.csv'));\n"
"f = figure('Color', 'w', 'Position', [160 160 1000 420]);\n"
"subplot(1,2,1);\n"
"plot(timeTbl.dt_init, timeTbl.final_profile_p_l2_norm, '-o', 'LineWidth', 1.6); hold on;\n"
"plot(timeTbl.dt_init, timeTbl.final_profile_t_l2_norm, '-s', 'LineWidth', 1.6);\n"
"set(gca, 'XDir', 'reverse'); xlabel('Initial Time Step (s)'); ylabel('Normalized L2');\n"
"title('Time Sensitivity: Profile'); legend({'Pressure', 'Temperature'}, 'Location', 'best'); grid on; box on;\n"
"subplot(1,2,2);\n"
"plot(timeTbl.dt_init, timeTbl.final_monitor_p_l2_norm, '-o', 'LineWidth', 1.6); hold on;\n"
"plot(timeTbl.dt_init, timeTbl.final_monitor_t_l2_norm, '-s', 'LineWidth', 1.6);\n"
"set(gca, 'XDir', 'reverse'); xlabel('Initial Time Step (s)'); ylabel('Normalized L2');\n"
"title('Time Sensitivity: Monitor'); legend({'Pressure', 'Temperature'}, 'Location', 'best'); grid on; box on;\n"
"exportgraphics(f, fullfile(figDir, 'time_sensitivity.pdf'), 'ContentType', 'vector');\n"
"exportgraphics(f, fullfile(figDir, 'time_sensitivity.png'), 'Resolution', 300);\n"
"close(f);\n"
"\n"
"diagTbl = readtable(fullfile(rootDir, 'boundary_diagnostics.csv'));\n"
"f = figure('Color', 'w', 'Position', [170 170 1100 420]);\n"
"subplot(1,2,1);\n"
"plot(diagTbl.time_s, diagTbl.left_prescribed_t_avg_k, '-o', 'LineWidth', 1.4); hold on;\n"
"plot(diagTbl.time_s, diagTbl.right_prescribed_t_avg_k, '-s', 'LineWidth', 1.4);\n"
"plot(diagTbl.time_s, diagTbl.left_adjacent_t_avg_k, '--o', 'LineWidth', 1.2);\n"
"plot(diagTbl.time_s, diagTbl.right_adjacent_t_avg_k, '--s', 'LineWidth', 1.2);\n"
"xlabel('Time (s)'); ylabel('Temperature (K)'); title('Boundary Temperature Diagnostics');\n"
"legend({'Left prescribed','Right prescribed','Left adjacent cells','Right adjacent cells'}, 'Location', 'best'); grid on; box on;\n"
"subplot(1,2,2);\n"
"plot(diagTbl.time_s, diagTbl.domain_t_span_k, '-d', 'LineWidth', 1.6); hold on;\n"
"yline(10.0, 'k--', 'Threshold');\n"
"xlabel('Time (s)'); ylabel('Temperature Span (K)'); title('Domain Temperature Span'); grid on; box on;\n"
"exportgraphics(f, fullfile(figDir, 'boundary_temperature_diagnostics.pdf'), 'ContentType', 'vector');\n"
"exportgraphics(f, fullfile(figDir, 'boundary_temperature_diagnostics.png'), 'Resolution', 300);\n"
"close(f);\n"
"\n"
"finalProfiles = cellfun(@(fam) readtable(fullfile(rootDir, ['compare_profile_' fam '_t100pct.csv'])), families, 'UniformOutput', false);\n"
"pNum = []; pRef = []; tNum = []; tRef = []; pErr = []; tErr = [];\n"
"for i = 1:numel(finalProfiles)\n"
"    pNum = [pNum; finalProfiles{i}.p_num_pa]; %#ok<AGROW>\n"
"    pRef = [pRef; finalProfiles{i}.p_ref_pa]; %#ok<AGROW>\n"
"    tNum = [tNum; finalProfiles{i}.t_num_k]; %#ok<AGROW>\n"
"    tRef = [tRef; finalProfiles{i}.t_ref_k]; %#ok<AGROW>\n"
"    pErr = [pErr; finalProfiles{i}.p_abs_err_over_dp]; %#ok<AGROW>\n"
"    tErr = [tErr; finalProfiles{i}.t_abs_err_over_dt]; %#ok<AGROW>\n"
"end\n"
"f = figure('Color', 'w', 'Position', [180 180 1000 420]);\n"
"subplot(1,2,1);\n"
"scatter(pRef, pNum, 18, 'filled'); hold on;\n"
"plot([min(pRef) max(pRef)], [min(pRef) max(pRef)], 'k--', 'LineWidth', 1.2);\n"
"xlabel('Reference Pressure (Pa)'); ylabel('Engineering Pressure (Pa)'); title('Parity Plot: Pressure'); grid on; box on;\n"
"subplot(1,2,2);\n"
"scatter(tRef, tNum, 18, 'filled'); hold on;\n"
"plot([min(tRef) max(tRef)], [min(tRef) max(tRef)], 'k--', 'LineWidth', 1.2);\n"
"xlabel('Reference Temperature (K)'); ylabel('Engineering Temperature (K)'); title('Parity Plot: Temperature'); grid on; box on;\n"
"exportgraphics(f, fullfile(figDir, 'parity.pdf'), 'ContentType', 'vector');\n"
"exportgraphics(f, fullfile(figDir, 'parity.png'), 'Resolution', 300);\n"
"close(f);\n"
"\n"
"f = figure('Color', 'w', 'Position', [200 200 1000 420]);\n"
"subplot(1,2,1); histogram(pErr, 20); xlabel('Normalized Pressure Error'); ylabel('Count'); title('Error Histogram: Pressure'); grid on; box on;\n"
"subplot(1,2,2); histogram(tErr, 20); xlabel('Normalized Temperature Error'); ylabel('Count'); title('Error Histogram: Temperature'); grid on; box on;\n"
"exportgraphics(f, fullfile(figDir, 'error_histogram.pdf'), 'ContentType', 'vector');\n"
"exportgraphics(f, fullfile(figDir, 'error_histogram.png'), 'Resolution', 300);\n"
"close(f);\n";
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
    artifacts.summary.case_dir = outputDirOverride.empty()
        ? (cfg.output_base_dir + "/" + cfg.sub_dir + "/" + cfg.case_name)
        : outputDirOverride;
    artifacts.summary.engineering_dir = artifacts.summary.case_dir + "/engineering";
    artifacts.summary.reference_dir = artifacts.summary.case_dir + "/reference";
    artifacts.summary.comsol_input_dir = artifacts.summary.reference_dir + "/comsol_input";
    artifacts.summary.comsol_output_dir = artifacts.summary.reference_dir + "/comsol";
    EnsureDirRecursive(artifacts.summary.case_dir);
    EnsureDirRecursive(artifacts.summary.engineering_dir);
    EnsureDirRecursive(artifacts.summary.reference_dir);
    EnsureDirRecursive(artifacts.summary.comsol_input_dir);
    EnsureDirRecursive(artifacts.summary.comsol_output_dir);

    artifacts.summary.convergence_log_path = artifacts.summary.case_dir + "/convergence.log";
    artifacts.summary.run_log_path = artifacts.summary.case_dir + "/run.log";
    artifacts.summary.metrics_csv_path = artifacts.summary.case_dir + "/metrics.csv";
    artifacts.summary.validation_summary_path = artifacts.summary.case_dir + "/validation_summary.md";
    artifacts.summary.validation_summary_csv_path = artifacts.summary.case_dir + "/validation_summary.csv";
    artifacts.summary.analytical_feasibility_path = artifacts.summary.reference_dir + "/analytical_feasibility_note.md";
    artifacts.summary.reference_spec_path = artifacts.summary.engineering_dir + "/reference_spec.md";
    artifacts.summary.property_table_path = artifacts.summary.engineering_dir + "/property_table.csv";
    artifacts.summary.comsol_property_table_path = artifacts.summary.comsol_input_dir + "/property_table.csv";
    artifacts.summary.boundary_diagnostics_path = artifacts.summary.case_dir + "/boundary_diagnostics.csv";
    artifacts.summary.profile_station_definitions_path = artifacts.summary.engineering_dir + "/profile_station_definitions.csv";
    artifacts.summary.monitor_point_definitions_path = artifacts.summary.engineering_dir + "/monitor_point_definitions.csv";
    artifacts.summary.profile_schedule_path = artifacts.summary.engineering_dir + "/profile_report_schedule.csv";
    artifacts.summary.monitor_schedule_path = artifacts.summary.engineering_dir + "/monitor_sample_schedule.csv";
    artifacts.summary.eng_monitor_timeseries_path = artifacts.summary.engineering_dir + "/eng_monitor_timeseries.csv";
    artifacts.summary.comsol_reference_spec_path = artifacts.summary.reference_dir + "/comsol_reference_spec.md";
    artifacts.summary.matlab_script_path = artifacts.summary.case_dir + "/plot_validation_results.m";
    artifacts.summary.grid_convergence_csv_path = artifacts.summary.case_dir + "/grid_convergence.csv";
    artifacts.summary.time_sensitivity_csv_path = artifacts.summary.case_dir + "/time_sensitivity.csv";
    artifacts.summary.nx = cfg.nx;
    artifacts.summary.ny = cfg.ny;

    std::ofstream convergenceLog(artifacts.summary.convergence_log_path.c_str(), std::ios::out | std::ios::trunc);
    if (!convergenceLog.good()) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to open convergence log: " + artifacts.summary.convergence_log_path);
    }
    std::ofstream runLog(artifacts.summary.run_log_path.c_str(), std::ios::out | std::ios::trunc);
    if (!runLog.good()) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to open run log: " + artifacts.summary.run_log_path);
    }
    runLog << "case=" << cfg.case_name << "\n";
    runLog << "workflow=" << WorkflowModeString(cfg.workflow_mode) << "\n";
    runLog << "reference_mode_requested=" << ReferenceModeString(cfg.reference_mode) << "\n";
    runLog << "output_dir=" << artifacts.summary.case_dir << "\n";

    MeshManager mgr(cfg.lx, cfg.ly, 0.0, cfg.nx, cfg.ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OverRelaxed);

    const FractureGeometry geom = BuildFractureGeometry(cfg);
    mgr.addFracture(geom.start, geom.end);
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    const std::size_t nCells = mgr.mesh().getCells().size();
    const int totalBlocks = mgr.getTotalDOFCount();
    if (nCells == 0) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] matrix cell count is zero");
    artifacts.summary.n_cells = static_cast<int>(nCells);
    artifacts.summary.n_fracture_dofs = std::max(0, totalBlocks - mgr.getMatrixDOFCount());

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

    const auto exportSnapshotArtifacts = [&](const std::string& snapshotTag, const std::string& vtkName, double timeS) {
        artifacts.summary.boundary_diagnostics_rows.push_back(
            BuildBoundaryDiagnosticsRow(mgr, fm, bcT, snapshotTag, timeS));
        PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(artifacts.summary.case_dir + "/" + vtkName, timeS);
    };

    std::vector<double> pBlocksLatest(static_cast<std::size_t>(std::max(totalBlocks, 0)), cfg.p_init);
    std::vector<double> tBlocksLatest(static_cast<std::size_t>(std::max(totalBlocks, 0)), cfg.t_init);
    SyncPTFieldsToFM(mgr, fm, pBlocksLatest, tBlocksLatest, cfg.p_init, cfg.t_init);
    if (cfg.export_vtk) {
        exportSnapshotArtifacts("initial", "initial.vtk", 0.0);
    }

    bool midExported = false;
    int iterSum = 0;
    int iterCount = 0;
    int maxIters = 0;
    double prevAcceptedTime = 0.0;
    std::vector<double> prevAcceptedPBlocks = pBlocksLatest;
    std::vector<double> prevAcceptedTBlocks = tBlocksLatest;

    for (auto& sample : artifacts.monitor_schedule) {
        if (sample.target_time_s <= kPendingTimeTolerance) {
            sample.captured = true;
            sample.actual_time_s = 0.0;
            sample.p_blocks = pBlocksLatest;
            sample.t_blocks = tBlocksLatest;
        }
    }

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
    modules.SetFluidModelConfig(FIM_Engine::UnifiedFluidModelConfig::MakeSinglePhaseCO2Constant(co2_props));

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

        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.k_t_tag, cfg.fracture_kt), cfg.fracture_kt);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.k_n_tag, cfg.fracture_kn), cfg.fracture_kn);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.phi_tag, cfg.fracture_phi), cfg.fracture_phi);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.c_r_tag, cfg.fracture_ct), cfg.fracture_ct);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.rho_tag, cfg.fracture_rho_r), cfg.fracture_rho_r);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.cp_tag, cfg.fracture_cp_r), cfg.fracture_cp_r);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.lambda_tag, cfg.fracture_lambda_r), cfg.fracture_lambda_r);

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
                exportSnapshotArtifacts("mid", "mid.vtk", timeS);
                midExported = true;
            }
        };

    const auto params = BuildSolverParams(cfg);
    FIM_Engine::RunGenericFIMTransient<2>(cfg.case_name, mgr, fm, ic, {}, params, FIM_Engine::SolverRoute::FIM, modules);

    SyncPTFieldsToFM(mgr, fm, pBlocksLatest, tBlocksLatest, cfg.p_init, cfg.t_init);
    if (cfg.export_vtk) {
        if (!midExported) exportSnapshotArtifacts("mid", "mid.vtk", artifacts.summary.t_end);
        exportSnapshotArtifacts("final", "final.vtk", artifacts.summary.t_end);
    }

    FinalizeMissingSnapshots(artifacts.summary.t_end, pBlocksLatest, tBlocksLatest, artifacts.snapshots);
    FinalizeMissingMonitorSamples(artifacts.summary.t_end, pBlocksLatest, tBlocksLatest, artifacts.monitor_schedule);

    artifacts.summary.max_iters = maxIters;
    artifacts.summary.avg_iters = (iterCount > 0) ? (static_cast<double>(iterSum) / static_cast<double>(iterCount)) : 0.0;
    artifacts.summary.h_char = ComputeMeshCharLength(mgr);

    runLog << "n_cells=" << artifacts.summary.n_cells << "\n";
    runLog << "n_fracture_dofs=" << artifacts.summary.n_fracture_dofs << "\n";
    runLog << "steps=" << artifacts.summary.steps << "\n";
    runLog << "t_end=" << artifacts.summary.t_end << "\n";

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
            }
        }
        WriteEngineeringMonitorCsv(artifacts.monitor_points, artifacts.monitor_schedule, artifacts.summary.eng_monitor_timeseries_path);
        WritePropertyTables(cfg, artifacts.summary.property_table_path, artifacts.summary.comsol_property_table_path);
        WriteBoundaryDiagnosticsCsv(artifacts.summary.boundary_diagnostics_rows, artifacts.summary.boundary_diagnostics_path);
        WriteReferenceSpec(cfg, artifacts.summary);
        WriteAnalyticalFeasibility(cfg, artifacts.summary.analytical_feasibility_path);
        WriteComsolReferenceSpec(cfg, artifacts.summary.case_dir, artifacts.summary.comsol_reference_spec_path);
    }

    return artifacts;
}
ReferenceMode ResolveReferenceMode(const TestCaseSpec& cfg, TestCaseSummary& summary) {
    summary.analytical_closed_form_available = false;
    summary.comsol_representation = kComsolRepresentationExplicit;
    summary.reference_is_explicit_lower_dimensional = false;
    if (cfg.reference_mode == ReferenceMode::Analytical) {
        summary.resolved_reference_mode = "Analytical";
        throw std::runtime_error(
            "[Test_H_T_CO2_ConstPP_SingleFrac] strict closed-form analytical solution is not available for this full single-fracture N=2 case.");
    }
    summary.resolved_reference_mode = "Comsol";
    return ReferenceMode::Comsol;
}

void InvokeComsolReferenceGeneration(const TestCaseSpec& cfg, const TestCaseSummary& summary) {
    std::ostringstream cmd;
    cmd << "powershell -ExecutionPolicy Bypass -File \"" << cfg.comsol_wrapper_relpath
        << "\" -Mode All -CaseDir \"" << summary.case_dir << "\"";
    const int code = std::system(cmd.str().c_str());
    if (code != 0) {
        throw std::runtime_error(
            "[Test_H_T_CO2_ConstPP_SingleFrac] COMSOL reference generation failed with exit code " + std::to_string(code));
    }
}

bool FinalMetricsWithinThreshold(const TestCaseSpec& cfg, const TestCaseSummary& summary) {
    return summary.final_profile_p_l2_norm <= cfg.profile_l2_threshold &&
           summary.final_profile_p_linf_norm <= cfg.profile_linf_threshold &&
           summary.final_profile_t_l2_norm <= cfg.profile_l2_threshold &&
           summary.final_profile_t_linf_norm <= cfg.profile_linf_threshold &&
           summary.final_monitor_p_l2_norm <= cfg.monitor_l2_threshold &&
           summary.final_monitor_p_linf_norm <= cfg.monitor_linf_threshold &&
           summary.final_monitor_t_l2_norm <= cfg.monitor_l2_threshold &&
           summary.final_monitor_t_linf_norm <= cfg.monitor_linf_threshold &&
           summary.reference_is_explicit_lower_dimensional &&
           summary.final_boundary_prescribed_t_ok &&
           summary.final_temperature_span_k >= kFinalTemperatureSpanThresholdK &&
           summary.final_boundary_adjacent_delta_t_k >= kBoundaryAdjacentDeltaTThresholdK;
}

TestCaseSummary RunCase(const TestCaseSpec& cfg) {
    CaseRunArtifacts mainArtifacts = RunSingleCaseCore(cfg, "");
    TestCaseSummary& summary = mainArtifacts.summary;
    const ReferenceMode resolvedMode = ResolveReferenceMode(cfg, summary);

    {
        std::ofstream runLog(summary.run_log_path.c_str(), std::ios::out | std::ios::app);
        runLog << "reference_mode_resolved=" << summary.resolved_reference_mode << "\n";
    }

    if (resolvedMode == ReferenceMode::Comsol) {
        const bool readyBefore = ReferenceFilesReady(summary, cfg);
        if (!readyBefore && cfg.workflow_mode == WorkflowMode::Full && cfg.allow_full_workflow_comsol_autorun) {
            InvokeComsolReferenceGeneration(cfg, summary);
            ReferenceFilesReady(summary, cfg);
        }
    }

    if (cfg.workflow_mode == WorkflowMode::PrepareReference) {
        summary.validation_status = summary.reference_ready ? "reference_ready" : "prepared_reference_inputs";
        summary.validation_performed = false;
        summary.validation_passed = false;
        WriteMatlabPlotScript(summary);
        WriteValidationSummary(cfg, summary);
        WriteValidationSummaryCsv(cfg, summary);
        WriteMetricsCsv(cfg, summary);
        return summary;
    }

    if (!summary.reference_ready) {
        summary.validation_status = "missing_reference";
        summary.validation_performed = false;
        summary.validation_passed = false;
        WriteMatlabPlotScript(summary);
        WriteValidationSummary(cfg, summary);
        WriteValidationSummaryCsv(cfg, summary);
        WriteMetricsCsv(cfg, summary);
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] reference files are not ready for validation.");
    }

    summary.comsol_representation = ExtractComsolRepresentation(summary.comsol_output_dir + "/comsol_run_summary.md");
    summary.reference_is_explicit_lower_dimensional =
        (summary.comsol_representation == kComsolRepresentationExplicit);
    summary.comsol_fine_check_skipped =
        ComsolFineCheckWasSkipped(summary.comsol_output_dir + "/comsol_reference_mesh_check.txt");
    if (!summary.reference_is_explicit_lower_dimensional) {
        summary.validation_status = "reference_not_acceptable";
        summary.validation_performed = false;
        summary.validation_passed = false;
        WriteMatlabPlotScript(summary);
        WriteValidationSummary(cfg, summary);
        WriteValidationSummaryCsv(cfg, summary);
        WriteMetricsCsv(cfg, summary);
        throw std::runtime_error(
            "[Test_H_T_CO2_ConstPP_SingleFrac] COMSOL reference is not explicit lower-dimensional fracture.");
    }

    summary.report_metrics.clear();
    const MonitorReferenceSeries monitorReference =
        LoadMonitorReference(summary.comsol_output_dir + "/comsol_monitor_timeseries.csv", mainArtifacts.monitor_points);

    const SnapshotState* finalSnapshot = mainArtifacts.snapshots.empty() ? nullptr : &mainArtifacts.snapshots.back();
    if (!finalSnapshot) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] no engineering snapshots available.");
    }

    for (const auto& snapshot : mainArtifacts.snapshots) {
        for (const std::string& family : OrderedProfileFamilies()) {
            const std::vector<SpatialSamplePoint> familyStations = FilterStationsByFamily(mainArtifacts.profile_stations, family);
            const ProfileReferenceTable reference = LoadProfileReference(
                summary.comsol_output_dir + "/comsol_profile_" + family + "_" + snapshot.tag + ".csv");
            const std::string comparePath = summary.case_dir + "/compare_profile_" + family + "_" + snapshot.tag + ".csv";
            summary.report_metrics.push_back(EvaluateProfileAgainstReference(
                familyStations,
                snapshot,
                reference,
                cfg,
                comparePath,
                cfg.emit_detailed_outputs));
        }
    }

    summary.final_profile_p_l2_norm = 0.0;
    summary.final_profile_p_linf_norm = 0.0;
    summary.final_profile_t_l2_norm = 0.0;
    summary.final_profile_t_linf_norm = 0.0;
    for (const auto& metric : summary.report_metrics) {
        if (metric.tag == finalSnapshot->tag) {
            summary.final_profile_p_l2_norm = std::max(summary.final_profile_p_l2_norm, metric.pressure.l2_norm);
            summary.final_profile_p_linf_norm = std::max(summary.final_profile_p_linf_norm, metric.pressure.linf_norm);
            summary.final_profile_t_l2_norm = std::max(summary.final_profile_t_l2_norm, metric.temperature.l2_norm);
            summary.final_profile_t_linf_norm = std::max(summary.final_profile_t_linf_norm, metric.temperature.linf_norm);
        }
    }

    FieldErrorMetrics finalMonitorTMetrics;
    const FieldErrorMetrics finalMonitorPMetrics = WriteMonitorCompareCsv(
        mainArtifacts.monitor_points,
        mainArtifacts.monitor_schedule,
        monitorReference,
        cfg,
        summary.case_dir + "/compare_monitor_timeseries.csv",
        finalMonitorTMetrics);
    summary.final_monitor_p_l2_norm = finalMonitorPMetrics.l2_norm;
    summary.final_monitor_p_linf_norm = finalMonitorPMetrics.linf_norm;
    summary.final_monitor_t_l2_norm = finalMonitorTMetrics.l2_norm;
    summary.final_monitor_t_linf_norm = finalMonitorTMetrics.linf_norm;
    if (!summary.boundary_diagnostics_rows.empty()) {
        const BoundaryDiagnosticsRow& finalDiag = summary.boundary_diagnostics_rows.back();
        if (finalDiag.domain_t.count > 0 &&
            std::isfinite(finalDiag.domain_t.max_value) &&
            std::isfinite(finalDiag.domain_t.min_value)) {
            summary.final_temperature_span_k = finalDiag.domain_t.max_value - finalDiag.domain_t.min_value;
        }
        if (std::isfinite(finalDiag.left_adjacent_t.avg_value) && std::isfinite(finalDiag.right_adjacent_t.avg_value)) {
            summary.final_boundary_adjacent_delta_t_k =
                std::abs(finalDiag.left_adjacent_t.avg_value - finalDiag.right_adjacent_t.avg_value);
        }
        summary.final_boundary_prescribed_t_ok =
            StatsWithinTolerance(finalDiag.left_prescribed_t, cfg.t_left, kBoundaryPrescribedTolerance) &&
            StatsWithinTolerance(finalDiag.right_prescribed_t, cfg.t_right, kBoundaryPrescribedTolerance);
    }

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
                summary.case_dir + "/studies/grid_" + std::to_string(sweepCfg.nx) + "x" + std::to_string(sweepCfg.ny));

            SweepStudyRow row;
            row.label = "grid_" + std::to_string(sweepCfg.nx) + "x" + std::to_string(sweepCfg.ny);
            row.case_dir = sweepArtifacts.summary.case_dir;
            row.nx = sweepCfg.nx;
            row.ny = sweepCfg.ny;
            row.dt_init = sweepCfg.dt_init;
            row.h_char = sweepArtifacts.summary.h_char;
            row.t_end = sweepArtifacts.summary.t_end;
            row.steps = sweepArtifacts.summary.steps;

            const SnapshotState& sweepFinalSnapshot = sweepArtifacts.snapshots.back();
            for (const std::string& family : OrderedProfileFamilies()) {
                const ProfileReferenceTable reference = LoadProfileReference(
                    summary.comsol_output_dir + "/comsol_profile_" + family + "_" + finalSnapshot->tag + ".csv");
                const ProfileCompareMetrics metric = EvaluateProfileAgainstReference(
                    FilterStationsByFamily(sweepArtifacts.profile_stations, family),
                    sweepFinalSnapshot,
                    reference,
                    sweepCfg,
                    "",
                    false);
                row.final_profile_p_l2_norm = std::max(row.final_profile_p_l2_norm, metric.pressure.l2_norm);
                row.final_profile_t_l2_norm = std::max(row.final_profile_t_l2_norm, metric.temperature.l2_norm);
            }

            const MonitorScheduleState& sweepFinalMonitor = sweepArtifacts.monitor_schedule.back();
            const auto refIt = monitorReference.sample_id_to_row.find(sweepFinalMonitor.sample_id);
            if (refIt == monitorReference.sample_id_to_row.end()) {
                throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] missing final monitor sample in reference during grid sweep.");
            }
            const MonitorReferenceRow& finalMonitorRef = monitorReference.rows[refIt->second];
            const FieldErrorMetrics sweepPMetrics = EvaluateFinalMonitorField(
                sweepArtifacts.monitor_points,
                sweepFinalMonitor,
                finalMonitorRef,
                sweepFinalMonitor.p_blocks,
                true,
                PressureScale(cfg));
            const FieldErrorMetrics sweepTMetrics = EvaluateFinalMonitorField(
                sweepArtifacts.monitor_points,
                sweepFinalMonitor,
                finalMonitorRef,
                sweepFinalMonitor.t_blocks,
                false,
                TemperatureScale(cfg));
            row.final_monitor_p_l2_norm = sweepPMetrics.l2_norm;
            row.final_monitor_t_l2_norm = sweepTMetrics.l2_norm;
            gridRows.push_back(row);
        }
        WriteStudyCSV(gridRows, summary.grid_convergence_csv_path, "grid");
        std::vector<double> gProfileP;
        std::vector<double> gProfileT;
        std::vector<double> gMonitorP;
        std::vector<double> gMonitorT;
        for (const auto& row : gridRows) {
            gProfileP.push_back(row.final_profile_p_l2_norm);
            gProfileT.push_back(row.final_profile_t_l2_norm);
            gMonitorP.push_back(row.final_monitor_p_l2_norm);
            gMonitorT.push_back(row.final_monitor_t_l2_norm);
        }
        summary.grid_convergence_ok =
            IsMonotonicNonIncreasingWithTol(gProfileP, kGridTimeMonotoneAbsTol) &&
            IsMonotonicNonIncreasingWithTol(gProfileT, kGridTimeMonotoneAbsTol) &&
            IsMonotonicNonIncreasingWithTol(gMonitorP, kGridTimeMonotoneAbsTol) &&
            IsMonotonicNonIncreasingWithTol(gMonitorT, kGridTimeMonotoneAbsTol);
    }

    if (cfg.enable_time_sensitivity_study) {
        for (double dtInit : cfg.time_step_sweep) {
            TestCaseSpec sweepCfg = cfg;
            sweepCfg.dt_init = dtInit;
            sweepCfg.case_name = cfg.case_name + "_dt_" + std::to_string(static_cast<int>(std::round(dtInit))) + "s";
            sweepCfg.enable_grid_convergence_study = false;
            sweepCfg.enable_time_sensitivity_study = false;
            sweepCfg.export_vtk = false;
            sweepCfg.emit_detailed_outputs = false;
            CaseRunArtifacts sweepArtifacts = RunSingleCaseCore(
                sweepCfg,
                summary.case_dir + "/studies/dt_" + std::to_string(static_cast<int>(std::round(dtInit))) + "s");

            SweepStudyRow row;
            row.label = "dt_" + std::to_string(static_cast<int>(std::round(dtInit))) + "s";
            row.case_dir = sweepArtifacts.summary.case_dir;
            row.nx = sweepCfg.nx;
            row.ny = sweepCfg.ny;
            row.dt_init = dtInit;
            row.h_char = sweepArtifacts.summary.h_char;
            row.t_end = sweepArtifacts.summary.t_end;
            row.steps = sweepArtifacts.summary.steps;

            const SnapshotState& sweepFinalSnapshot = sweepArtifacts.snapshots.back();
            for (const std::string& family : OrderedProfileFamilies()) {
                const ProfileReferenceTable reference = LoadProfileReference(
                    summary.comsol_output_dir + "/comsol_profile_" + family + "_" + finalSnapshot->tag + ".csv");
                const ProfileCompareMetrics metric = EvaluateProfileAgainstReference(
                    FilterStationsByFamily(sweepArtifacts.profile_stations, family),
                    sweepFinalSnapshot,
                    reference,
                    sweepCfg,
                    "",
                    false);
                row.final_profile_p_l2_norm = std::max(row.final_profile_p_l2_norm, metric.pressure.l2_norm);
                row.final_profile_t_l2_norm = std::max(row.final_profile_t_l2_norm, metric.temperature.l2_norm);
            }

            const MonitorScheduleState& sweepFinalMonitor = sweepArtifacts.monitor_schedule.back();
            const auto refIt = monitorReference.sample_id_to_row.find(sweepFinalMonitor.sample_id);
            if (refIt == monitorReference.sample_id_to_row.end()) {
                throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] missing final monitor sample in reference during time sweep.");
            }
            const MonitorReferenceRow& finalMonitorRef = monitorReference.rows[refIt->second];
            const FieldErrorMetrics sweepPMetrics = EvaluateFinalMonitorField(
                sweepArtifacts.monitor_points,
                sweepFinalMonitor,
                finalMonitorRef,
                sweepFinalMonitor.p_blocks,
                true,
                PressureScale(cfg));
            const FieldErrorMetrics sweepTMetrics = EvaluateFinalMonitorField(
                sweepArtifacts.monitor_points,
                sweepFinalMonitor,
                finalMonitorRef,
                sweepFinalMonitor.t_blocks,
                false,
                TemperatureScale(cfg));
            row.final_monitor_p_l2_norm = sweepPMetrics.l2_norm;
            row.final_monitor_t_l2_norm = sweepTMetrics.l2_norm;
            timeRows.push_back(row);
        }
        WriteStudyCSV(timeRows, summary.time_sensitivity_csv_path, "time");
        std::vector<double> tProfileP;
        std::vector<double> tProfileT;
        std::vector<double> tMonitorP;
        std::vector<double> tMonitorT;
        for (const auto& row : timeRows) {
            tProfileP.push_back(row.final_profile_p_l2_norm);
            tProfileT.push_back(row.final_profile_t_l2_norm);
            tMonitorP.push_back(row.final_monitor_p_l2_norm);
            tMonitorT.push_back(row.final_monitor_t_l2_norm);
        }
        summary.time_sensitivity_ok =
            IsMonotonicNonIncreasingWithTol(tProfileP, kGridTimeMonotoneAbsTol) &&
            IsMonotonicNonIncreasingWithTol(tProfileT, kGridTimeMonotoneAbsTol) &&
            IsMonotonicNonIncreasingWithTol(tMonitorP, kGridTimeMonotoneAbsTol) &&
            IsMonotonicNonIncreasingWithTol(tMonitorT, kGridTimeMonotoneAbsTol);
    }

    summary.validation_performed = true;
    summary.validation_passed = FinalMetricsWithinThreshold(cfg, summary) &&
        (!cfg.enable_grid_convergence_study || summary.grid_convergence_ok) &&
        (!cfg.enable_time_sensitivity_study || summary.time_sensitivity_ok);
    summary.validation_status = summary.validation_passed ? "passed" : "failed";

    WriteMatlabPlotScript(summary);
    WriteValidationSummary(cfg, summary);
    WriteValidationSummaryCsv(cfg, summary);
    WriteMetricsCsv(cfg, summary);

    {
        std::ofstream runLog(summary.run_log_path.c_str(), std::ios::out | std::ios::app);
        runLog << "validation_status=" << summary.validation_status << "\n";
        runLog << "validation_passed=" << BoolString(summary.validation_passed) << "\n";
    }

    if (!summary.validation_passed) {
        std::ostringstream oss;
        oss << "[Test_H_T_CO2_ConstPP_SingleFrac] validation failed: "
            << "profile_p_l2=" << summary.final_profile_p_l2_norm
            << ", profile_t_l2=" << summary.final_profile_t_l2_norm
            << ", monitor_p_l2=" << summary.final_monitor_p_l2_norm
            << ", monitor_t_l2=" << summary.final_monitor_t_l2_norm
            << ", comsol_representation=" << summary.comsol_representation
            << ", t_span=" << summary.final_temperature_span_k
            << ", boundary_delta_t=" << summary.final_boundary_adjacent_delta_t_k
            << ", prescribed_t_ok=" << BoolString(summary.final_boundary_prescribed_t_ok)
            << ", grid_ok=" << BoolString(summary.grid_convergence_ok)
            << ", time_ok=" << BoolString(summary.time_sensitivity_ok);
        throw std::runtime_error(oss.str());
    }
    return summary;
}
TestCasePlan BuildDefaultPlan() {
    TestCasePlan plan;
    plan.plan_key = "h_t_co2_constpp_singlefrac_nowell";
    return plan;
}

TestCasePlan BuildPrepareReferencePlan() {
    TestCasePlan plan = BuildDefaultPlan();
    plan.plan_key = "h_t_co2_constpp_singlefrac_nowell_prepare_reference";
    plan.spec.workflow_mode = WorkflowMode::PrepareReference;
    return plan;
}

TestCasePlan BuildValidateReferencePlan() {
    TestCasePlan plan = BuildDefaultPlan();
    plan.plan_key = "h_t_co2_constpp_singlefrac_nowell_validate_reference";
    plan.spec.reference_mode = ReferenceMode::Comsol;
    plan.spec.workflow_mode = WorkflowMode::ValidateAgainstReference;
    return plan;
}

TestCasePlan BuildFullComsolPlan() {
    TestCasePlan plan = BuildDefaultPlan();
    plan.plan_key = "h_t_co2_constpp_singlefrac_nowell_full_comsol";
    plan.spec.reference_mode = ReferenceMode::Comsol;
    plan.spec.workflow_mode = WorkflowMode::Full;
    return plan;
}

using BuilderFn = TestCasePlan(*)();

const std::unordered_map<std::string, BuilderFn>& GetRegistry() {
    static const std::unordered_map<std::string, BuilderFn> registry = {
        {"h_t_co2_constpp_singlefrac_nowell", &BuildDefaultPlan},
        {"h_t_co2_constpp_singlefrac_nowell_prepare_reference", &BuildPrepareReferencePlan},
        {"h_t_co2_constpp_singlefrac_nowell_validate_reference", &BuildValidateReferencePlan},
        {"h_t_co2_constpp_singlefrac_nowell_full_comsol", &BuildFullComsolPlan}
    };
    return registry;
}

void ExecutePlanByKeyImpl(const std::string& key) {
    const auto& registry = GetRegistry();
    const auto it = registry.find(key);
    if (it == registry.end()) throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] unknown registry key: " + key);
    const TestCasePlan plan = it->second();
    const TestCaseSummary summary = RunCase(plan.spec);
    std::cout << "\n============================================\n";
    std::cout << "[Test_H_T_CO2_ConstPP_SingleFrac] run completed\n";
    std::cout << "  output_dir: " << summary.case_dir << "\n";
    std::cout << "  grid: " << summary.nx << " x " << summary.ny << " (" << summary.n_cells << " cells)\n";
    std::cout << "  fracture_dofs: " << summary.n_fracture_dofs << "\n";
    std::cout << "  steps: " << summary.steps << "  rollbacks: " << summary.total_rollbacks << "\n";
    std::cout << "  Newton iters: avg=" << std::fixed << std::setprecision(2) << summary.avg_iters
              << "  max=" << summary.max_iters << "\n";
    std::cout << "  final_time: " << std::scientific << std::setprecision(4) << summary.t_end << " s\n";
    std::cout << "  reference_mode: " << summary.resolved_reference_mode << "\n";
    std::cout << "  validation_status: " << summary.validation_status << "\n";
    if (summary.validation_performed) {
        std::cout << "  final profile p/t L2: "
                  << summary.final_profile_p_l2_norm << " / " << summary.final_profile_t_l2_norm << "\n";
        std::cout << "  final monitor p/t L2: "
                  << summary.final_monitor_p_l2_norm << " / " << summary.final_monitor_t_l2_norm << "\n";
    }
    if (!summary.missing_reference_files.empty()) {
        std::cout << "  missing reference files:\n";
        for (const auto& path : summary.missing_reference_files) std::cout << "    - " << path << "\n";
    }
    std::cout << "  validation_summary: " << summary.validation_summary_path << "\n";
    std::cout << "============================================\n";
}

} // namespace

void RunTestCase() {
    ExecutePlanByKeyImpl("h_t_co2_constpp_singlefrac_nowell");
}

void ExecutePlanByKey(const std::string& key) {
    ExecutePlanByKeyImpl(key);
}

} // namespace Test_H_T_CO2_ConstPP_SingleFrac

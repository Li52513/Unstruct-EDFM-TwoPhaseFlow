/**
 * @file Test_2D_EDFM_H_CO2_VaryPP_NoFrac_NoWell.cpp
 * @brief Standalone validation chain: 2D single-phase CO2 variable-property (EOS), no-fracture, no-well.
 */

#include "Test_2D_EDFM_H_CO2_VaryPP_NoFrac_NoWell.h"

#include "2D_PostProcess.h"
#include "AD_FluidEvaluator.h"
#include "BoundaryConditionManager.h"
#include "FIM_TransientCaseKit.hpp"
#include "MeshDefinitions.h"
#include "MeshManager.h"
#include "SolverContrlStrName_op.h"
#include "Well_WellControlTypes.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
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

namespace Test_H_CO2_VaryPP {
namespace {

constexpr double kPendingTimeTolerance = 1.0e-9;
constexpr double kPressureTableMinPa = 7.5e6;
constexpr double kPressureTableMaxPa = 12.5e6;
constexpr int kPressureTableSamples = 2001;
constexpr double kMonitorSampleDtS = 500.0;

enum class ReferenceMode { Auto, Analytical, Comsol };
enum class WorkflowMode { PrepareReference, ValidateAgainstReference, Full };

struct SpatialSamplePoint {
    int id = -1;
    std::string label;
    double target_x = 0.0;
    double target_y = 0.0;
    int cell_id = -1;
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
};

struct MonitorScheduleState {
    int sample_id = -1;
    double target_time_s = 0.0;
    double actual_time_s = 0.0;
    bool captured = false;
    std::vector<double> p_blocks;
};

struct ProfileCompareMetrics {
    std::string tag;
    double requested_fraction = 0.0;
    double target_time_s = 0.0;
    double actual_time_s = 0.0;
    int sample_count = 0;
    double l1_abs = 0.0;
    double l2_abs = 0.0;
    double linf_abs = 0.0;
    double l1_norm = 0.0;
    double l2_norm = 0.0;
    double linf_norm = 0.0;
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
    double l1_norm = 0.0;
    double l2_norm = 0.0;
    double linf_norm = 0.0;
};

struct DenseProfileReference {
    std::vector<double> x;
    std::vector<double> p_ref;
};

struct MonitorReferenceSeries {
    std::vector<double> time_s;
    std::vector<std::string> point_labels;
    std::vector<std::vector<double> > p_ref_series;
};

struct TestCaseSpec {
    std::string case_name = "h_co2_varypp_nofrac_nowell";
    std::string output_base_dir = "Test/Transient/FullCaseTest";
    std::string sub_dir = "H_CO2_VaryPP";
    double lx = 400.0;
    double ly = 40.0;
    int nx = 48;
    int ny = 6;
    double matrix_phi = 0.10;
    double matrix_perm = 1.0e-13;
    double matrix_ct = 5.0e-9;
    double p_init = 10.0e6;
    double p_left = 12.0e6;
    double p_right = 8.0e6;
    double t_init = 360.0;
    double dt_init = 480.0;
    double dt_min = 1.0;
    double dt_max = 10000.0;
    double target_end_time_s = 1.0e5;
    int max_steps = 12000;
    int max_newton_iter = 12;
    FIM_Engine::LinearSolverType lin_solver = FIM_Engine::LinearSolverType::AMGCL;
    bool amgcl_use_fallback_sparselu = true;
    double amgcl_tol = 1.0e-6;
    int amgcl_maxiter = 500;
    bool enable_non_orthogonal_correction = true;
    bool enable_row_scaling = true;
    double abs_res_tol = 1.0e-10;
    double rel_res_tol = 1.0e-6;
    double rel_update_tol = 1.0e-8;
    double max_dP = 2.0e7;
    bool enable_armijo_line_search = true;
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
    std::vector<double> time_step_sweep = {480.0, 120.0, 30.0};
    int profile_station_count = 81;
    double monitor_dt_s = kMonitorSampleDtS;
    std::vector<std::pair<double, double> > monitor_targets = {
        std::make_pair(40.0, 20.0),
        std::make_pair(120.0, 20.0),
        std::make_pair(200.0, 20.0),
        std::make_pair(280.0, 20.0),
        std::make_pair(360.0, 20.0)
    };
    double profile_l2_threshold = 3.0e-2;
    double profile_linf_threshold = 5.0e-2;
    bool enable_grid_convergence_study = true;
    bool enable_time_sensitivity_study = true;
    bool export_vtk = true;
    bool emit_detailed_outputs = true;
};

struct TestCaseSummary {
    std::string case_dir;
    std::string convergence_log_path;
    std::string metrics_csv_path;
    std::string validation_summary_path;
    std::string analytical_feasibility_path;
    std::string coolprop_property_table_path;
    std::string profile_station_definitions_path;
    std::string monitor_point_definitions_path;
    std::string eng_monitor_timeseries_path;
    std::string comsol_reference_spec_path;
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
    double rho_min = 0.0;
    double rho_max = 0.0;
    double mu_min = 0.0;
    double mu_max = 0.0;
    std::string validation_status = "not_run";
    bool reference_ready = false;
    bool validation_performed = false;
    bool validation_passed = false;
    bool grid_convergence_ok = true;
    bool time_sensitivity_ok = true;
    std::vector<std::string> missing_reference_files;
    double final_profile_l1_abs = 0.0;
    double final_profile_l2_abs = 0.0;
    double final_profile_linf_abs = 0.0;
    double final_profile_l1_norm = 0.0;
    double final_profile_l2_norm = 0.0;
    double final_profile_linf_norm = 0.0;
    double final_monitor_l2_norm = 0.0;
    std::vector<ProfileCompareMetrics> report_metrics;
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
    if (!in.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to open csv: " + path);
    CsvTable table;
    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("[Test_H_CO2_VaryPP] empty csv: " + path);
    table.headers = SplitCsvLine(line);
    for (std::size_t i = 0; i < table.headers.size(); ++i) table.header_index[table.headers[i]] = i;
    while (std::getline(in, line)) {
        if (!Trim(line).empty()) table.rows.push_back(SplitCsvLine(line));
    }
    return table;
}

double CsvGetDouble(const CsvTable& table, std::size_t row, const std::string& column) {
    const auto it = table.header_index.find(column);
    if (it == table.header_index.end()) throw std::runtime_error("[Test_H_CO2_VaryPP] missing csv column: " + column);
    if (row >= table.rows.size()) throw std::runtime_error("[Test_H_CO2_VaryPP] csv row out of range.");
    const std::size_t col = it->second;
    if (col >= table.rows[row].size()) throw std::runtime_error("[Test_H_CO2_VaryPP] csv column out of range.");
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

void ApplyUniformScalarField(const std::shared_ptr<volScalarField>& field, double value) {
    if (!field) return;
    for (double& v : field->data) v = value;
}

void SyncPressureFieldsToFM(MeshManager& mgr, FieldManager_2D& fm,
                            const std::vector<double>& pBlocks, double tConst) {
    if (pBlocks.empty()) return;
    const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
    const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
    auto fPw = fm.getOrCreateMatrixScalar(pCfg.pressure_field, 0.0);
    auto fPviz = fm.getOrCreateMatrixScalar("P", 0.0);
    auto fT = fm.getOrCreateMatrixScalar(tCfg.temperatue_field, tConst);
    auto fracPw = fm.getOrCreateFractureScalar(pCfg.pressure_field, 0.0);
    auto fracPviz = fm.getOrCreateFractureScalar("P", 0.0);
    auto fracT = fm.getOrCreateFractureScalar(tCfg.temperatue_field, tConst);
    const int nMat = mgr.getMatrixDOFCount();
    const int nUse = std::min(static_cast<int>(pBlocks.size()), mgr.getTotalDOFCount());
    for (int i = 0; i < nUse; ++i) {
        const double p = pBlocks[static_cast<std::size_t>(i)];
        if (i < nMat) {
            if (fPw && i < static_cast<int>(fPw->data.size())) fPw->data[static_cast<std::size_t>(i)] = p;
            if (fPviz && i < static_cast<int>(fPviz->data.size())) fPviz->data[static_cast<std::size_t>(i)] = p;
            if (fT && i < static_cast<int>(fT->data.size())) fT->data[static_cast<std::size_t>(i)] = tConst;
        } else {
            const int fi = i - nMat;
            if (fracPw && fi < static_cast<int>(fracPw->data.size())) fracPw->data[static_cast<std::size_t>(fi)] = p;
            if (fracPviz && fi < static_cast<int>(fracPviz->data.size())) fracPviz->data[static_cast<std::size_t>(fi)] = p;
            if (fracT && fi < static_cast<int>(fracT->data.size())) fracT->data[static_cast<std::size_t>(fi)] = tConst;
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

SpatialSamplePoint MakeNearestSample(const MeshManager& mgr, double targetX, double targetY, int id, const std::string& label) {
    SpatialSamplePoint sample;
    sample.id = id;
    sample.label = label;
    sample.target_x = targetX;
    sample.target_y = targetY;
    const auto& cells = mgr.mesh().getCells();
    if (cells.empty()) throw std::runtime_error("[Test_H_CO2_VaryPP] cannot build sample on empty mesh.");
    double bestDistSq = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < cells.size(); ++i) {
        const double dx = cells[i].center.m_x - targetX;
        const double dy = cells[i].center.m_y - targetY;
        const double distSq = dx * dx + dy * dy;
        if (distSq < bestDistSq) {
            bestDistSq = distSq;
            sample.cell_id = static_cast<int>(i);
            sample.actual_x = cells[i].center.m_x;
            sample.actual_y = cells[i].center.m_y;
        }
    }
    return sample;
}

std::vector<SpatialSamplePoint> BuildProfileStations(const MeshManager& mgr, const TestCaseSpec& cfg) {
    if (cfg.profile_station_count < 2) throw std::runtime_error("[Test_H_CO2_VaryPP] profile_station_count must be >= 2.");
    std::vector<SpatialSamplePoint> stations;
    stations.reserve(static_cast<std::size_t>(cfg.profile_station_count));
    const double yMid = 0.5 * cfg.ly;
    const double dx = cfg.lx / static_cast<double>(cfg.profile_station_count - 1);
    for (int i = 0; i < cfg.profile_station_count; ++i) {
        std::ostringstream label;
        label << "st" << std::setw(3) << std::setfill('0') << i;
        stations.push_back(MakeNearestSample(mgr, dx * static_cast<double>(i), yMid, i, label.str()));
    }
    return stations;
}

std::vector<SpatialSamplePoint> BuildMonitorPoints(const MeshManager& mgr, const TestCaseSpec& cfg) {
    std::vector<SpatialSamplePoint> points;
    points.reserve(cfg.monitor_targets.size());
    for (std::size_t i = 0; i < cfg.monitor_targets.size(); ++i) {
        std::ostringstream label;
        label << "pt" << std::setw(2) << std::setfill('0') << (i + 1);
        points.push_back(MakeNearestSample(mgr, cfg.monitor_targets[i].first, cfg.monitor_targets[i].second,
                                           static_cast<int>(i), label.str()));
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

void CaptureSnapshotIfCrossed(double prevTime, const std::vector<double>& prevPBlocks,
                              double currentTime, const std::vector<double>& currentPBlocks,
                              SnapshotState& snapshot) {
    if (snapshot.captured || currentTime + kPendingTimeTolerance < snapshot.target_time_s) return;
    const bool usePrev = !prevPBlocks.empty() && prevTime >= 0.0
        && std::abs(prevTime - snapshot.target_time_s) < std::abs(currentTime - snapshot.target_time_s);
    snapshot.captured = true;
    snapshot.actual_time_s = usePrev ? prevTime : currentTime;
    snapshot.p_blocks = usePrev ? prevPBlocks : currentPBlocks;
}

void CaptureMonitorIfCrossed(double prevTime, const std::vector<double>& prevPBlocks,
                             double currentTime, const std::vector<double>& currentPBlocks,
                             MonitorScheduleState& sample) {
    if (sample.captured || currentTime + kPendingTimeTolerance < sample.target_time_s) return;
    const bool usePrev = !prevPBlocks.empty() && prevTime >= 0.0
        && std::abs(prevTime - sample.target_time_s) < std::abs(currentTime - sample.target_time_s);
    sample.captured = true;
    sample.actual_time_s = usePrev ? prevTime : currentTime;
    sample.p_blocks = usePrev ? prevPBlocks : currentPBlocks;
}

void FinalizeMissingSnapshots(double finalTime, const std::vector<double>& finalPBlocks, std::vector<SnapshotState>& snapshots) {
    for (auto& snapshot : snapshots) {
        if (!snapshot.captured) {
            snapshot.captured = true;
            snapshot.actual_time_s = finalTime;
            snapshot.p_blocks = finalPBlocks;
        }
    }
}

void FinalizeMissingMonitorSamples(double finalTime, const std::vector<double>& finalPBlocks,
                                   std::vector<MonitorScheduleState>& schedule) {
    for (auto& sample : schedule) {
        if (!sample.captured) {
            sample.captured = true;
            sample.actual_time_s = finalTime;
            sample.p_blocks = finalPBlocks;
        }
    }
}

void WriteProfileStationDefinitions(const std::vector<SpatialSamplePoint>& stations, const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to write profile station definitions: " + path);
    out << "station_id,label,target_x_m,target_y_m,cell_id,actual_x_m,actual_y_m\n";
    for (const auto& station : stations) {
        out << station.id << "," << station.label << ","
            << std::setprecision(12) << station.target_x << "," << station.target_y << ","
            << station.cell_id << "," << station.actual_x << "," << station.actual_y << "\n";
    }
}

void WriteMonitorPointDefinitions(const std::vector<SpatialSamplePoint>& points, const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to write monitor point definitions: " + path);
    out << "point_id,label,target_x_m,target_y_m,cell_id,actual_x_m,actual_y_m\n";
    for (const auto& point : points) {
        out << point.id << "," << point.label << ","
            << std::setprecision(12) << point.target_x << "," << point.target_y << ","
            << point.cell_id << "," << point.actual_x << "," << point.actual_y << "\n";
    }
}

void WriteProfileSchedule(const std::vector<SnapshotState>& snapshots, const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to write profile schedule: " + path);
    out << "tag,requested_fraction,target_time_s,actual_time_s\n";
    for (const auto& snapshot : snapshots) {
        out << snapshot.tag << ","
            << std::setprecision(12) << snapshot.requested_fraction << ","
            << snapshot.target_time_s << "," << snapshot.actual_time_s << "\n";
    }
}

void WriteMonitorSchedule(const std::vector<MonitorScheduleState>& schedule, const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to write monitor schedule: " + path);
    out << "sample_id,target_time_s,actual_time_s\n";
    for (const auto& sample : schedule) {
        out << sample.sample_id << ","
            << std::setprecision(12) << sample.target_time_s << "," << sample.actual_time_s << "\n";
    }
}

void WriteEngineeringProfileCsv(const std::vector<SpatialSamplePoint>& stations,
                                const SnapshotState& snapshot,
                                const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to write engineering profile csv: " + path);
    out << "station_id,label,target_x_m,target_y_m,cell_id,actual_x_m,actual_y_m,target_time_s,actual_time_s,p_num_pa\n";
    for (const auto& station : stations) {
        if (station.cell_id < 0 || station.cell_id >= static_cast<int>(snapshot.p_blocks.size())) {
            throw std::runtime_error("[Test_H_CO2_VaryPP] invalid station cell id during engineering profile export.");
        }
        out << station.id << "," << station.label << ","
            << std::setprecision(12) << station.target_x << "," << station.target_y << ","
            << station.cell_id << "," << station.actual_x << "," << station.actual_y << ","
            << snapshot.target_time_s << "," << snapshot.actual_time_s << ","
            << snapshot.p_blocks[static_cast<std::size_t>(station.cell_id)] << "\n";
    }
}

void WriteEngineeringMonitorCsv(const std::vector<SpatialSamplePoint>& points,
                                const std::vector<MonitorScheduleState>& schedule,
                                const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to write engineering monitor csv: " + path);
    out << "sample_id,target_time_s,actual_time_s";
    for (const auto& point : points) out << ",p_num_" << point.label;
    out << "\n";
    for (const auto& sample : schedule) {
        out << sample.sample_id << ","
            << std::setprecision(12) << sample.target_time_s << "," << sample.actual_time_s;
        for (const auto& point : points) {
            if (point.cell_id < 0 || point.cell_id >= static_cast<int>(sample.p_blocks.size())) {
                throw std::runtime_error("[Test_H_CO2_VaryPP] invalid monitor cell id during engineering monitor export.");
            }
            out << "," << sample.p_blocks[static_cast<std::size_t>(point.cell_id)];
        }
        out << "\n";
    }
}

void WriteAnalyticalFeasibility(const TestCaseSpec& cfg, const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to write analytical feasibility note: " + path);
    out << "reference_strategy=COMSOL\n";
    out << "analytical_closed_form_available=false\n";
    out << "reason=Pressure-dependent rho(p), mu(p), and cf(p) from the real CO2 EOS make the diffusion problem nonlinear. "
           "This validation chain does not use a pseudo-analytical solution that re-calls the same EOS backend.\n";
    out << "fixed_temperature_k=" << std::setprecision(12) << cfg.t_init << "\n";
    out << "gravity_vector=(" << cfg.gravity_vector.m_x << "," << cfg.gravity_vector.m_y << "," << cfg.gravity_vector.m_z << ")\n";
#ifdef USE_COOLPROP_EOS
    out << "compiled_with_coolprop=true\n";
#else
    out << "compiled_with_coolprop=false\n";
#endif
}

void WriteCoolPropPropertyTables(const TestCaseSpec& cfg,
                                 const std::string& combinedCsvPath,
                                 const std::string& splitDir) {
    EnsureDirRecursive(splitDir);
    std::ofstream combined(combinedCsvPath.c_str(), std::ios::out | std::ios::trunc);
    if (!combined.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to open combined property table: " + combinedCsvPath);
    combined << "p_pa,rho,mu,cf,lambda,S\n";

    const std::array<std::string, 5> splitNames = {{
        splitDir + "/rho_fun.txt",
        splitDir + "/mu_fun.txt",
        splitDir + "/cf_fun.txt",
        splitDir + "/lambda_fun.txt",
        splitDir + "/S_fun.txt"
    }};
    std::array<std::ofstream, 5> splitFiles;
    for (std::size_t i = 0; i < splitNames.size(); ++i) {
        splitFiles[i].open(splitNames[i].c_str(), std::ios::out | std::ios::trunc);
        if (!splitFiles[i].good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to open split property table: " + splitNames[i]);
    }

    for (int i = 0; i < kPressureTableSamples; ++i) {
        const double ratio = (kPressureTableSamples <= 1) ? 0.0
            : static_cast<double>(i) / static_cast<double>(kPressureTableSamples - 1);
        const double pPa = kPressureTableMinPa + (kPressureTableMaxPa - kPressureTableMinPa) * ratio;
        ADVar<1> pVar(pPa, 0);
        ADVar<1> tVar(cfg.t_init);
        const AD_Fluid::ADFluidProperties<1> props = AD_Fluid::Evaluator::evaluateCO2<1>(pVar, tVar);
        const double rho = props.rho.val;
        const double mu = props.mu.val;
        const double cf = (std::abs(rho) > 1.0e-30) ? (props.rho.grad(0) / rho) : 0.0;
        const double lambda = (std::abs(mu) > 1.0e-30) ? (cfg.matrix_perm * rho / mu) : 0.0;
        const double storage = cfg.matrix_phi * rho * (cf + cfg.matrix_ct);

        combined << std::setprecision(12) << pPa << "," << rho << "," << mu << ","
                 << cf << "," << lambda << "," << storage << "\n";
        splitFiles[0] << std::setprecision(12) << pPa << " " << rho << "\n";
        splitFiles[1] << std::setprecision(12) << pPa << " " << mu << "\n";
        splitFiles[2] << std::setprecision(12) << pPa << " " << cf << "\n";
        splitFiles[3] << std::setprecision(12) << pPa << " " << lambda << "\n";
        splitFiles[4] << std::setprecision(12) << pPa << " " << storage << "\n";
    }
}

void WriteComsolReferenceSpec(const TestCaseSpec& cfg, const std::string& caseDir, const std::string& path) {
    const std::string javaSource = "Tools/COMSOL/VaryPP_NoFrac_NoWell/ComsolVaryPPNoFracNoWell.java";
    const std::string javaClass = "Tools/COMSOL/VaryPP_NoFrac_NoWell/ComsolVaryPPNoFracNoWell.class";
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to write COMSOL reference spec: " + path);
    out << "# COMSOL Reference Specification\n\n";
    out << "## Case\n";
    out << "- Case directory: `" << caseDir << "`\n";
    out << "- Temperature: `" << cfg.t_init << " K`\n";
    out << "- Pressure window: `" << cfg.p_left / 1.0e6 << " MPa` -> `" << cfg.p_right / 1.0e6 << " MPa`\n";
    out << "- Geometry: `" << cfg.lx << " m x " << cfg.ly << " m`\n";
    out << "- Gravity: disabled\n";
    out << "- Non-orthogonal correction: enabled on engineering side\n\n";
    out << "## Required Inputs\n";
    out << "- `coolprop_property_table.csv`\n";
    out << "- `reference/comsol_input/rho_fun.txt`\n";
    out << "- `reference/comsol_input/mu_fun.txt`\n";
    out << "- `reference/comsol_input/cf_fun.txt`\n";
    out << "- `reference/comsol_input/lambda_fun.txt`\n";
    out << "- `reference/comsol_input/S_fun.txt`\n";
    out << "- `profile_station_definitions.csv`\n";
    out << "- `profile_report_schedule.csv`\n";
    out << "- `monitor_point_definitions.csv`\n\n";
    out << "- `monitor_sample_schedule.csv`\n\n";
    out << "## Expected COMSOL Outputs\n";
    out << "- `reference/comsol/comsol_profile_t010pct.csv`\n";
    out << "- `reference/comsol/comsol_profile_t050pct.csv`\n";
    out << "- `reference/comsol/comsol_profile_t100pct.csv`\n";
    out << "- `reference/comsol/comsol_monitor_timeseries.csv`\n";
    out << "- `reference/comsol/comsol_model.mph`\n";
    out << "- `reference/comsol/comsol_batch.log`\n\n";
    out << "## Compile Command\n";
    out << "```powershell\n";
    out << "& \"D:\\1Comsol6.3\\Comosol6.3\\COMSOL63\\Multiphysics\\bin\\win64\\comsolcompile.exe\" \\\n";
    out << "  \"" << javaSource << "\"\n";
    out << "```\n\n";
    out << "## Batch Run Command\n";
    out << "```powershell\n";
    out << "& \"D:\\1Comsol6.3\\Comosol6.3\\COMSOL63\\Multiphysics\\bin\\win64\\comsolbatch.exe\" \\\n";
    out << "  -inputfile \"" << javaClass << "\" \\\n";
    out << "  -batchlog \"" << caseDir << "/reference/comsol/comsol_batch.log\"\n";
    out << "```\n\n";
    out << "## PDE\n";
    out << "- Solve `S(p) * dp/dt - div(lambda(p) * grad(p)) = 0`\n";
    out << "- Use `lambda(p) = k * rho(p) / mu(p)`\n";
    out << "- Use `S(p) = phi * rho(p) * (cf(p) + c_r)`\n";
    out << "- Use interpolation functions built from the exported split tables\n\n";
    out << "## Solver Settings\n";
    out << "- Time-dependent study, BDF, relative tolerance `1e-6`\n";
    out << "- Output times: `0:500:100000 s`\n";
    out << "- Dense profile exports at `t = 10000, 50000, 100000 s`\n";
    out << "- Centerline profile: `y = " << (0.5 * cfg.ly) << " m`, dense x-sampling over `[0, " << cfg.lx << "]`\n";
}

void WriteStudyCSV(const std::vector<SweepStudyRow>& rows, const std::string& path, const std::string& studyTag) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to write study csv: " + path);
    out << "study,label,case_dir,nx,ny,dt_init,h_char,t_end,steps,l1_norm,l2_norm,linf_norm\n";
    for (const auto& row : rows) {
        out << studyTag << "," << row.label << "," << row.case_dir << ","
            << row.nx << "," << row.ny << ","
            << std::setprecision(12) << row.dt_init << "," << row.h_char << ","
            << row.t_end << "," << row.steps << ","
            << row.l1_norm << "," << row.l2_norm << "," << row.linf_norm << "\n";
    }
}

bool IsMonotonicNonIncreasingL2(const std::vector<SweepStudyRow>& rows) {
    const double absTol = 1.0e-10;
    const double relTol = 1.0e-6;
    for (std::size_t i = 1; i < rows.size(); ++i) {
        const double prev = rows[i - 1].l2_norm;
        const double allowed = std::max(absTol, relTol * std::max(prev, 1.0));
        if (rows[i].l2_norm > prev + allowed) return false;
    }
    return true;
}

double InterpolateLinear(const std::vector<double>& x, const std::vector<double>& y, double xQuery) {
    if (x.empty() || y.empty() || x.size() != y.size()) {
        throw std::runtime_error("[Test_H_CO2_VaryPP] invalid interpolation data.");
    }
    if (x.size() == 1) return y.front();
    if (xQuery <= x.front()) return y.front();
    if (xQuery >= x.back()) return y.back();
    const auto upper = std::upper_bound(x.begin(), x.end(), xQuery);
    const std::size_t idxHi = static_cast<std::size_t>(upper - x.begin());
    const std::size_t idxLo = idxHi - 1;
    const double x0 = x[idxLo];
    const double x1 = x[idxHi];
    const double y0 = y[idxLo];
    const double y1 = y[idxHi];
    const double weight = (std::abs(x1 - x0) > 1.0e-30) ? ((xQuery - x0) / (x1 - x0)) : 0.0;
    return y0 + (y1 - y0) * weight;
}

DenseProfileReference LoadProfileReference(const std::string& path,
                                           const std::vector<SpatialSamplePoint>& referenceStations) {
    const CsvTable table = ReadCsvTable(path);
    DenseProfileReference ref;
    ref.x.reserve(table.rows.size());
    ref.p_ref.reserve(table.rows.size());
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        double xRef = CsvGetDouble(table, i, "x_m");
        const int stationId = static_cast<int>(std::llround(CsvGetDouble(table, i, "station_id")));
        if (stationId >= 0 && stationId < static_cast<int>(referenceStations.size())) {
            xRef = referenceStations[static_cast<std::size_t>(stationId)].target_x;
        }
        ref.x.push_back(xRef);
        ref.p_ref.push_back(CsvGetDouble(table, i, "p_ref_pa"));
    }
    std::vector<std::size_t> order(ref.x.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) { return ref.x[lhs] < ref.x[rhs]; });
    DenseProfileReference sorted;
    sorted.x.reserve(order.size());
    sorted.p_ref.reserve(order.size());
    for (std::size_t idx : order) {
        sorted.x.push_back(ref.x[idx]);
        sorted.p_ref.push_back(ref.p_ref[idx]);
    }
    return sorted;
}

MonitorReferenceSeries LoadMonitorReference(const std::string& path, const std::vector<SpatialSamplePoint>& monitorPoints) {
    const CsvTable table = ReadCsvTable(path);
    MonitorReferenceSeries ref;
    ref.time_s.reserve(table.rows.size());
    ref.point_labels.reserve(monitorPoints.size());
    ref.p_ref_series.resize(monitorPoints.size());
    for (const auto& point : monitorPoints) ref.point_labels.push_back(point.label);
    for (std::size_t row = 0; row < table.rows.size(); ++row) {
        ref.time_s.push_back(CsvGetDouble(table, row, "time_s"));
        for (std::size_t i = 0; i < monitorPoints.size(); ++i) {
            ref.p_ref_series[i].push_back(CsvGetDouble(table, row, "p_ref_" + monitorPoints[i].label));
        }
    }
    return ref;
}

ProfileCompareMetrics EvaluateProfileAgainstReference(const std::vector<SpatialSamplePoint>& stations,
                                                      const SnapshotState& snapshot,
                                                      const DenseProfileReference& reference,
                                                      const TestCaseSpec& cfg,
                                                      const std::string& csvPath,
                                                      bool writeCsv) {
    const double deltaP = PressureScale(cfg);
    double sumAbs = 0.0;
    double sumSq = 0.0;
    double maxAbs = 0.0;
    std::ofstream out;
    if (writeCsv) {
        out.open(csvPath.c_str(), std::ios::out | std::ios::trunc);
        if (!out.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to write profile compare csv: " + csvPath);
        out << "station_id,label,target_x_m,target_y_m,cell_id,actual_x_m,actual_y_m,ref_query_x_m,target_time_s,actual_time_s,"
               "p_num_pa,p_ref_pa,abs_err_pa,abs_err_over_dP\n";
    }
    for (const auto& station : stations) {
        if (station.cell_id < 0 || station.cell_id >= static_cast<int>(snapshot.p_blocks.size())) {
            throw std::runtime_error("[Test_H_CO2_VaryPP] invalid station cell id during profile comparison.");
        }
        const double pNum = snapshot.p_blocks[static_cast<std::size_t>(station.cell_id)];
        const double pRef = InterpolateLinear(reference.x, reference.p_ref, station.target_x);
        const double absErr = std::abs(pNum - pRef);
        sumAbs += absErr;
        sumSq += absErr * absErr;
        maxAbs = std::max(maxAbs, absErr);
        if (writeCsv) {
            out << station.id << "," << station.label << ","
                << std::setprecision(12) << station.target_x << "," << station.target_y << ","
                << station.cell_id << "," << station.actual_x << "," << station.actual_y << "," << station.target_x << ","
                << snapshot.target_time_s << "," << snapshot.actual_time_s << ","
                << pNum << "," << pRef << "," << absErr << "," << (absErr / deltaP) << "\n";
        }
    }
    ProfileCompareMetrics metrics;
    metrics.tag = snapshot.tag;
    metrics.requested_fraction = snapshot.requested_fraction;
    metrics.target_time_s = snapshot.target_time_s;
    metrics.actual_time_s = snapshot.actual_time_s;
    metrics.sample_count = static_cast<int>(stations.size());
    metrics.l1_abs = sumAbs / static_cast<double>(stations.size());
    metrics.l2_abs = std::sqrt(sumSq / static_cast<double>(stations.size()));
    metrics.linf_abs = maxAbs;
    metrics.l1_norm = metrics.l1_abs / deltaP;
    metrics.l2_norm = metrics.l2_abs / deltaP;
    metrics.linf_norm = metrics.linf_abs / deltaP;
    metrics.compare_csv_path = writeCsv ? csvPath : "";
    return metrics;
}

double InterpolateMonitorPressure(const MonitorReferenceSeries& reference, std::size_t pointIndex, double timeQuery) {
    if (pointIndex >= reference.p_ref_series.size()) throw std::runtime_error("[Test_H_CO2_VaryPP] invalid monitor reference point index.");
    return InterpolateLinear(reference.time_s, reference.p_ref_series[pointIndex], timeQuery);
}

double WriteMonitorCompareCsv(const std::vector<SpatialSamplePoint>& monitorPoints,
                              const std::vector<MonitorScheduleState>& schedule,
                              const MonitorReferenceSeries& reference,
                              const TestCaseSpec& cfg,
                              const std::string& path) {
    const double deltaP = PressureScale(cfg);
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to write monitor compare csv: " + path);
    out << "sample_id,target_time_s,actual_time_s";
    for (const auto& point : monitorPoints) {
        out << ",p_num_" << point.label << ",p_ref_" << point.label
            << ",abs_err_" << point.label << ",abs_err_over_dP_" << point.label;
    }
    out << ",l1_norm,l2_norm,linf_norm\n";

    double finalL2Norm = 0.0;
    for (const auto& sample : schedule) {
        double sumAbs = 0.0;
        double sumSq = 0.0;
        double maxAbs = 0.0;
        out << sample.sample_id << "," << std::setprecision(12) << sample.target_time_s << "," << sample.actual_time_s;
        for (std::size_t i = 0; i < monitorPoints.size(); ++i) {
            const auto& point = monitorPoints[i];
            if (point.cell_id < 0 || point.cell_id >= static_cast<int>(sample.p_blocks.size())) {
                throw std::runtime_error("[Test_H_CO2_VaryPP] invalid monitor cell id during monitor comparison.");
            }
            const double pNum = sample.p_blocks[static_cast<std::size_t>(point.cell_id)];
            const double pRef = InterpolateMonitorPressure(reference, i, sample.target_time_s);
            const double absErr = std::abs(pNum - pRef);
            sumAbs += absErr;
            sumSq += absErr * absErr;
            maxAbs = std::max(maxAbs, absErr);
            out << "," << pNum << "," << pRef << "," << absErr << "," << (absErr / deltaP);
        }
        const double l1Norm = (sumAbs / static_cast<double>(monitorPoints.size())) / deltaP;
        const double l2Norm = std::sqrt(sumSq / static_cast<double>(monitorPoints.size())) / deltaP;
        const double linfNorm = maxAbs / deltaP;
        out << "," << l1Norm << "," << l2Norm << "," << linfNorm << "\n";
        finalL2Norm = l2Norm;
    }
    return finalL2Norm;
}

std::vector<std::string> RequiredReferenceFiles(const TestCaseSummary& summary, const TestCaseSpec& cfg) {
    std::vector<std::string> files;
    const std::vector<double> fractions = BuildSortedFractions(cfg.report_time_fractions);
    for (double fraction : fractions) {
        files.push_back(summary.case_dir + "/reference/comsol/comsol_profile_" + MakeReportTag(fraction) + ".csv");
    }
    files.push_back(summary.case_dir + "/reference/comsol/comsol_monitor_timeseries.csv");
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
    if (!metrics.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to write metrics csv: " + summary.metrics_csv_path);
    metrics << "case_name,nx,ny,n_cells,h_char,t_end,steps,total_rollbacks,avg_nonlinear_iters,max_nonlinear_iters,"
               "dt_init,dt_min,dt_max,property_mode,solver_route,reference_mode,workflow_mode,validation_status,"
               "reference_ready,validation_performed,validation_passed,grid_convergence_ok,time_sensitivity_ok,"
               "final_profile_l1_norm,final_profile_l2_norm,final_profile_linf_norm,final_monitor_l2_norm,"
               "rho_min,rho_max,mu_min,mu_max\n";
    metrics << cfg.case_name << "," << cfg.nx << "," << cfg.ny << "," << summary.n_cells << ","
            << std::setprecision(12) << summary.h_char << "," << summary.t_end << ","
            << summary.steps << "," << summary.total_rollbacks << "," << summary.avg_iters << "," << summary.max_iters << ","
            << cfg.dt_init << "," << cfg.dt_min << "," << cfg.dt_max << ","
            << "CO2_EOS,RunGenericFIMTransient<1>,"
            << ReferenceModeString(cfg.reference_mode) << "," << WorkflowModeString(cfg.workflow_mode) << ","
            << summary.validation_status << "," << BoolString(summary.reference_ready) << ","
            << BoolString(summary.validation_performed) << "," << BoolString(summary.validation_passed) << ","
            << BoolString(summary.grid_convergence_ok) << "," << BoolString(summary.time_sensitivity_ok) << ","
            << summary.final_profile_l1_norm << "," << summary.final_profile_l2_norm << ","
            << summary.final_profile_linf_norm << "," << summary.final_monitor_l2_norm << ","
            << std::scientific << std::setprecision(8)
            << summary.rho_min << "," << summary.rho_max << "," << summary.mu_min << "," << summary.mu_max << "\n";
}

void WriteValidationSummary(const TestCaseSpec& cfg, const TestCaseSummary& summary) {
    std::ofstream out(summary.validation_summary_path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_CO2_VaryPP] failed to write validation summary: " + summary.validation_summary_path);
    out << "case=" << cfg.case_name << "\n";
    out << "output_dir=" << summary.case_dir << "\n";
    out << "reference_mode=" << ReferenceModeString(cfg.reference_mode) << "\n";
    out << "workflow_mode=" << WorkflowModeString(cfg.workflow_mode) << "\n";
    out << "validation_status=" << summary.validation_status << "\n";
    out << "reference_ready=" << BoolString(summary.reference_ready) << "\n";
    out << "validation_performed=" << BoolString(summary.validation_performed) << "\n";
    out << "validation_passed=" << BoolString(summary.validation_passed) << "\n";
    out << "grid_convergence_ok=" << BoolString(summary.grid_convergence_ok) << "\n";
    out << "time_sensitivity_ok=" << BoolString(summary.time_sensitivity_ok) << "\n";
    out << "gravity_vector=(" << cfg.gravity_vector.m_x << "," << cfg.gravity_vector.m_y << "," << cfg.gravity_vector.m_z << ")\n";
    out << "non_orthogonal_correction=" << BoolString(cfg.enable_non_orthogonal_correction) << "\n";
#ifdef USE_COOLPROP_EOS
    out << "compiled_with_coolprop=true\n";
#else
    out << "compiled_with_coolprop=false\n";
#endif
    out << "final_profile_l1_norm=" << summary.final_profile_l1_norm << "\n";
    out << "final_profile_l2_norm=" << summary.final_profile_l2_norm << "\n";
    out << "final_profile_linf_norm=" << summary.final_profile_linf_norm << "\n";
    out << "final_monitor_l2_norm=" << summary.final_monitor_l2_norm << "\n";
    out << "profile_l2_threshold=" << cfg.profile_l2_threshold << "\n";
    out << "profile_linf_threshold=" << cfg.profile_linf_threshold << "\n";
    out << "coolprop_property_table=" << summary.coolprop_property_table_path << "\n";
    out << "profile_station_definitions=" << summary.profile_station_definitions_path << "\n";
    out << "monitor_point_definitions=" << summary.monitor_point_definitions_path << "\n";
    out << "eng_monitor_timeseries=" << summary.eng_monitor_timeseries_path << "\n";
    out << "comsol_reference_spec=" << summary.comsol_reference_spec_path << "\n";
    if (!summary.grid_convergence_csv_path.empty()) out << "grid_convergence_csv=" << summary.grid_convergence_csv_path << "\n";
    if (!summary.time_sensitivity_csv_path.empty()) out << "time_sensitivity_csv=" << summary.time_sensitivity_csv_path << "\n";
    if (!summary.missing_reference_files.empty()) {
        out << "\n[missing_reference_files]\n";
        for (const auto& path : summary.missing_reference_files) out << path << "\n";
    }
    if (!summary.report_metrics.empty()) {
        out << "\n[profile_report_metrics]\n";
        for (const auto& metric : summary.report_metrics) {
            out << metric.tag
                << " target_time_s=" << metric.target_time_s
                << " actual_time_s=" << metric.actual_time_s
                << " l1_norm=" << metric.l1_norm
                << " l2_norm=" << metric.l2_norm
                << " linf_norm=" << metric.linf_norm
                << " compare_csv=" << metric.compare_csv_path << "\n";
        }
    }
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
    p.max_dT = 1.0;
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
    EnsureDirRecursive(artifacts.summary.case_dir);
    artifacts.summary.convergence_log_path = artifacts.summary.case_dir + "/convergence.log";
    artifacts.summary.metrics_csv_path = artifacts.summary.case_dir + "/metrics.csv";
    artifacts.summary.validation_summary_path = artifacts.summary.case_dir + "/validation_summary.txt";
    artifacts.summary.analytical_feasibility_path = artifacts.summary.case_dir + "/analytical_feasibility.txt";
    artifacts.summary.coolprop_property_table_path = artifacts.summary.case_dir + "/coolprop_property_table.csv";
    artifacts.summary.profile_station_definitions_path = artifacts.summary.case_dir + "/profile_station_definitions.csv";
    artifacts.summary.monitor_point_definitions_path = artifacts.summary.case_dir + "/monitor_point_definitions.csv";
    artifacts.summary.eng_monitor_timeseries_path = artifacts.summary.case_dir + "/eng_monitor_timeseries.csv";
    artifacts.summary.comsol_reference_spec_path = artifacts.summary.case_dir + "/comsol_reference_spec.md";
    artifacts.summary.nx = cfg.nx;
    artifacts.summary.ny = cfg.ny;

    std::ofstream convergenceLog(artifacts.summary.convergence_log_path.c_str(), std::ios::out | std::ios::trunc);
    if (!convergenceLog.good()) {
        throw std::runtime_error("[Test_H_CO2_VaryPP] failed to open convergence log: " + artifacts.summary.convergence_log_path);
    }

    MeshManager mgr(cfg.lx, cfg.ly, 0.0, cfg.nx, cfg.ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OverRelaxed);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(1);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    const std::size_t nCells = mgr.mesh().getCells().size();
    const int totalBlocks = mgr.getTotalDOFCount();
    if (nCells == 0) throw std::runtime_error("[Test_H_CO2_VaryPP] matrix cell count is zero");
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

    const auto pEqCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
    VTKBoundaryVisualizationContext bcVizCtx;
    bcVizCtx.water_family_policy = VTKBCWaterFamilyDerivePolicy::FollowPrimaryFluid;
    bcVizCtx.primary_fluid_model = VTKBCPrimaryFluidModel::CO2;
    bcVizCtx.bindings.push_back(VTKBCVariableBinding{pEqCfg.pressure_field, &bcP, VTKBCTransportKind::Pressure});

    std::vector<double> pBlocksLatest(static_cast<std::size_t>(std::max(totalBlocks, 0)), cfg.p_init);
    SyncPressureFieldsToFM(mgr, fm, pBlocksLatest, cfg.t_init);
    if (cfg.export_vtk) {
        PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(artifacts.summary.case_dir + "/initial.vtk", 0.0);
    }

    bool midExported = false;
    int iterSum = 0;
    int iterCount = 0;
    int maxIters = 0;
    double prevAcceptedTime = 0.0;
    std::vector<double> prevAcceptedPBlocks = pBlocksLatest;

    for (auto& sample : artifacts.monitor_schedule) {
        if (sample.target_time_s <= kPendingTimeTolerance) {
            sample.captured = true;
            sample.actual_time_s = 0.0;
            sample.p_blocks = pBlocksLatest;
        }
    }

    FIM_Engine::TransientOptionalModules<MeshManager, FieldManager_2D> modules;
    modules.pressure_bc = &bcP;
    modules.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::CO2;
    modules.pressure_only_property_mode = FIM_Engine::PressureOnlyPropertyMode::CO2_EOS;
    modules.pressure_only_temperature_k = cfg.t_init;
    modules.disable_default_vtk_output = true;
    modules.property_initializer = [&cfg](MeshManager&, FieldManager_2D& fld) {
        const auto rock = PhysicalProperties_string_op::Rock();
        const auto frac = PhysicalProperties_string_op::Fracture_string();
        const auto water = PhysicalProperties_string_op::Water();
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.k_xx_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.k_yy_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.k_zz_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.phi_tag, cfg.matrix_phi), cfg.matrix_phi);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.c_r_tag, cfg.matrix_ct), cfg.matrix_ct);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.rho_tag, 2600.0), 2600.0);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.cp_tag, 1000.0), 1000.0);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.lambda_tag, 2.0), 2.0);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.k_t_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.k_n_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.phi_tag, cfg.matrix_phi), cfg.matrix_phi);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.c_r_tag, cfg.matrix_ct), cfg.matrix_ct);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(water.k_tag, 0.6), 0.6);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(water.k_tag, 0.6), 0.6);
    };
    modules.on_step_accepted =
        [&](int step, double timeS, double dtUsedS, int newtonIters, double residualInf,
            int totalRollbacks, const std::string& convergeMode,
            const std::vector<double>& pVec, const std::vector<double>&, const std::vector<double>*) {
            if (step <= 0) return;
            if (!pVec.empty()) {
                pBlocksLatest = pVec;
                SyncPressureFieldsToFM(mgr, fm, pBlocksLatest, cfg.t_init);
            }
            artifacts.summary.steps = step;
            artifacts.summary.t_end = timeS;
            artifacts.summary.total_rollbacks = totalRollbacks;
            iterSum += newtonIters;
            iterCount += 1;
            maxIters = std::max(maxIters, newtonIters);
            convergenceLog << "[Step " << step << "] t=" << std::scientific << std::setprecision(8) << timeS
                           << " dt=" << dtUsedS << " iters=" << newtonIters
                           << " residual_inf=" << residualInf
                           << " rollbacks=" << totalRollbacks
                           << " mode=" << convergeMode << "\n";
            for (auto& snapshot : artifacts.snapshots) {
                CaptureSnapshotIfCrossed(prevAcceptedTime, prevAcceptedPBlocks, timeS, pBlocksLatest, snapshot);
            }
            for (auto& sample : artifacts.monitor_schedule) {
                CaptureMonitorIfCrossed(prevAcceptedTime, prevAcceptedPBlocks, timeS, pBlocksLatest, sample);
            }
            prevAcceptedTime = timeS;
            prevAcceptedPBlocks = pBlocksLatest;
            if (cfg.export_vtk && !midExported && timeS >= 0.5 * cfg.target_end_time_s - 1.0e-12) {
                PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(artifacts.summary.case_dir + "/mid.vtk", timeS);
                midExported = true;
            }
        };

    const auto params = BuildSolverParams(cfg);
    FIM_Engine::RunGenericFIMTransient<1>(cfg.case_name, mgr, fm, ic, {}, params, FIM_Engine::SolverRoute::FIM, modules);

    SyncPressureFieldsToFM(mgr, fm, pBlocksLatest, cfg.t_init);
    if (cfg.export_vtk) {
        if (!midExported) PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(artifacts.summary.case_dir + "/mid.vtk", artifacts.summary.t_end);
        PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(artifacts.summary.case_dir + "/final.vtk", artifacts.summary.t_end);
    }

    FinalizeMissingSnapshots(artifacts.summary.t_end, pBlocksLatest, artifacts.snapshots);
    FinalizeMissingMonitorSamples(artifacts.summary.t_end, pBlocksLatest, artifacts.monitor_schedule);

    artifacts.summary.max_iters = maxIters;
    artifacts.summary.avg_iters = (iterCount > 0) ? (static_cast<double>(iterSum) / static_cast<double>(iterCount)) : 0.0;
    artifacts.summary.h_char = ComputeMeshCharLength(mgr);

    {
        const auto water = PhysicalProperties_string_op::Water();
        auto fRho = fm.getMatrixScalar(water.rho_tag);
        auto fMu = fm.getMatrixScalar(water.mu_tag);
        if (fRho && !fRho->data.empty()) {
            double rMin = std::numeric_limits<double>::max();
            double rMax = std::numeric_limits<double>::lowest();
            for (double v : fRho->data) if (v > 0.0) { rMin = std::min(rMin, v); rMax = std::max(rMax, v); }
            if (rMin <= rMax) { artifacts.summary.rho_min = rMin; artifacts.summary.rho_max = rMax; }
        }
        if (fMu && !fMu->data.empty()) {
            double mMin = std::numeric_limits<double>::max();
            double mMax = std::numeric_limits<double>::lowest();
            for (double v : fMu->data) if (v > 0.0) { mMin = std::min(mMin, v); mMax = std::max(mMax, v); }
            if (mMin <= mMax) { artifacts.summary.mu_min = mMin; artifacts.summary.mu_max = mMax; }
        }
    }

    if (cfg.emit_detailed_outputs) {
        const std::string comsolInputDir = artifacts.summary.case_dir + "/reference/comsol_input";
        EnsureDirRecursive(artifacts.summary.case_dir + "/reference/comsol");
        WriteProfileStationDefinitions(artifacts.profile_stations, artifacts.summary.profile_station_definitions_path);
        WriteMonitorPointDefinitions(artifacts.monitor_points, artifacts.summary.monitor_point_definitions_path);
        WriteProfileSchedule(artifacts.snapshots, artifacts.summary.case_dir + "/profile_report_schedule.csv");
        WriteMonitorSchedule(artifacts.monitor_schedule, artifacts.summary.case_dir + "/monitor_sample_schedule.csv");
        for (const auto& snapshot : artifacts.snapshots) {
            WriteEngineeringProfileCsv(artifacts.profile_stations, snapshot,
                                       artifacts.summary.case_dir + "/eng_profile_" + snapshot.tag + ".csv");
        }
        WriteEngineeringMonitorCsv(artifacts.monitor_points, artifacts.monitor_schedule, artifacts.summary.eng_monitor_timeseries_path);
        WriteCoolPropPropertyTables(cfg, artifacts.summary.coolprop_property_table_path, comsolInputDir);
        WriteAnalyticalFeasibility(cfg, artifacts.summary.analytical_feasibility_path);
        WriteComsolReferenceSpec(cfg, artifacts.summary.case_dir, artifacts.summary.comsol_reference_spec_path);
    }

    return artifacts;
}

TestCaseSummary RunCase(const TestCaseSpec& cfg) {
    CaseRunArtifacts mainArtifacts = RunSingleCaseCore(cfg, "");
    TestCaseSummary& summary = mainArtifacts.summary;

    const bool referenceReady = ReferenceFilesReady(summary, cfg);
    const bool wantsValidation = cfg.workflow_mode == WorkflowMode::ValidateAgainstReference
        || (cfg.workflow_mode == WorkflowMode::Full && referenceReady);

    if (cfg.workflow_mode == WorkflowMode::PrepareReference || !wantsValidation) {
        summary.validation_status = referenceReady ? "reference_ready" : "pending_reference";
        summary.validation_performed = false;
        summary.validation_passed = false;
        WriteValidationSummary(cfg, summary);
        WriteMetricsCsv(cfg, summary);
        return summary;
    }

    std::vector<SweepStudyRow> gridRows;
    std::vector<SweepStudyRow> timeRows;
    const std::vector<double> fractions = BuildSortedFractions(cfg.report_time_fractions);
    std::vector<DenseProfileReference> profileReferences;
    profileReferences.reserve(fractions.size());
    for (double fraction : fractions) {
        profileReferences.push_back(LoadProfileReference(
            summary.case_dir + "/reference/comsol/comsol_profile_" + MakeReportTag(fraction) + ".csv",
            mainArtifacts.profile_stations));
    }
    const MonitorReferenceSeries monitorReference =
        LoadMonitorReference(summary.case_dir + "/reference/comsol/comsol_monitor_timeseries.csv",
                             mainArtifacts.monitor_points);

    summary.report_metrics.clear();
    for (std::size_t i = 0; i < mainArtifacts.snapshots.size(); ++i) {
        const SnapshotState& snapshot = mainArtifacts.snapshots[i];
        ProfileCompareMetrics metrics = EvaluateProfileAgainstReference(
            mainArtifacts.profile_stations,
            snapshot,
            profileReferences[i],
            cfg,
            summary.case_dir + "/profile_compare_" + snapshot.tag + ".csv",
            cfg.emit_detailed_outputs);
        summary.report_metrics.push_back(metrics);
    }

    if (!summary.report_metrics.empty()) {
        const ProfileCompareMetrics& finalMetrics = summary.report_metrics.back();
        summary.final_profile_l1_abs = finalMetrics.l1_abs;
        summary.final_profile_l2_abs = finalMetrics.l2_abs;
        summary.final_profile_linf_abs = finalMetrics.linf_abs;
        summary.final_profile_l1_norm = finalMetrics.l1_norm;
        summary.final_profile_l2_norm = finalMetrics.l2_norm;
        summary.final_profile_linf_norm = finalMetrics.linf_norm;
    }

    summary.final_monitor_l2_norm = WriteMonitorCompareCsv(
        mainArtifacts.monitor_points,
        mainArtifacts.monitor_schedule,
        monitorReference,
        cfg,
        summary.case_dir + "/monitor_compare_timeseries.csv");

    if (cfg.enable_grid_convergence_study) {
        const DenseProfileReference& finalReference = profileReferences.back();
        for (const auto& gridCase : cfg.grid_sweep_cases) {
            TestCaseSpec sweepCfg = cfg;
            sweepCfg.nx = gridCase.first;
            sweepCfg.ny = gridCase.second;
            sweepCfg.case_name = cfg.case_name + "_grid_" + std::to_string(sweepCfg.nx) + "x" + std::to_string(sweepCfg.ny);
            sweepCfg.workflow_mode = WorkflowMode::ValidateAgainstReference;
            sweepCfg.enable_grid_convergence_study = false;
            sweepCfg.enable_time_sensitivity_study = false;
            sweepCfg.export_vtk = false;
            sweepCfg.emit_detailed_outputs = false;
            CaseRunArtifacts sweepArtifacts = RunSingleCaseCore(
                sweepCfg,
                summary.case_dir + "/studies/grid_" + std::to_string(sweepCfg.nx) + "x" + std::to_string(sweepCfg.ny));
            const ProfileCompareMetrics sweepMetrics = EvaluateProfileAgainstReference(
                sweepArtifacts.profile_stations,
                sweepArtifacts.snapshots.back(),
                finalReference,
                sweepCfg,
                "",
                false);
            SweepStudyRow row;
            row.label = "grid_" + std::to_string(sweepCfg.nx) + "x" + std::to_string(sweepCfg.ny);
            row.case_dir = sweepArtifacts.summary.case_dir;
            row.nx = sweepCfg.nx;
            row.ny = sweepCfg.ny;
            row.dt_init = sweepCfg.dt_init;
            row.h_char = sweepArtifacts.summary.h_char;
            row.t_end = sweepArtifacts.summary.t_end;
            row.steps = sweepArtifacts.summary.steps;
            row.l1_norm = sweepMetrics.l1_norm;
            row.l2_norm = sweepMetrics.l2_norm;
            row.linf_norm = sweepMetrics.linf_norm;
            gridRows.push_back(row);
        }
        summary.grid_convergence_ok = IsMonotonicNonIncreasingL2(gridRows);
        summary.grid_convergence_csv_path = summary.case_dir + "/grid_convergence.csv";
        WriteStudyCSV(gridRows, summary.grid_convergence_csv_path, "grid");
    }

    if (cfg.enable_time_sensitivity_study) {
        const DenseProfileReference& finalReference = profileReferences.back();
        for (double dtInit : cfg.time_step_sweep) {
            TestCaseSpec sweepCfg = cfg;
            sweepCfg.dt_init = dtInit;
            sweepCfg.case_name = cfg.case_name + "_dt_" + std::to_string(static_cast<int>(std::round(dtInit))) + "s";
            sweepCfg.workflow_mode = WorkflowMode::ValidateAgainstReference;
            sweepCfg.enable_grid_convergence_study = false;
            sweepCfg.enable_time_sensitivity_study = false;
            sweepCfg.export_vtk = false;
            sweepCfg.emit_detailed_outputs = false;
            CaseRunArtifacts sweepArtifacts = RunSingleCaseCore(
                sweepCfg,
                summary.case_dir + "/studies/dt_" + std::to_string(static_cast<int>(std::round(dtInit))) + "s");
            const ProfileCompareMetrics sweepMetrics = EvaluateProfileAgainstReference(
                sweepArtifacts.profile_stations,
                sweepArtifacts.snapshots.back(),
                finalReference,
                sweepCfg,
                "",
                false);
            SweepStudyRow row;
            row.label = "dt_" + std::to_string(static_cast<int>(std::round(dtInit))) + "s";
            row.case_dir = sweepArtifacts.summary.case_dir;
            row.nx = sweepCfg.nx;
            row.ny = sweepCfg.ny;
            row.dt_init = dtInit;
            row.h_char = sweepArtifacts.summary.h_char;
            row.t_end = sweepArtifacts.summary.t_end;
            row.steps = sweepArtifacts.summary.steps;
            row.l1_norm = sweepMetrics.l1_norm;
            row.l2_norm = sweepMetrics.l2_norm;
            row.linf_norm = sweepMetrics.linf_norm;
            timeRows.push_back(row);
        }
        summary.time_sensitivity_ok = IsMonotonicNonIncreasingL2(timeRows);
        summary.time_sensitivity_csv_path = summary.case_dir + "/time_sensitivity.csv";
        WriteStudyCSV(timeRows, summary.time_sensitivity_csv_path, "time");
    }

    summary.validation_performed = true;
    summary.validation_passed =
        summary.final_profile_l2_norm <= cfg.profile_l2_threshold &&
        summary.final_profile_linf_norm <= cfg.profile_linf_threshold &&
        (!cfg.enable_grid_convergence_study || summary.grid_convergence_ok) &&
        (!cfg.enable_time_sensitivity_study || summary.time_sensitivity_ok);
    summary.validation_status = summary.validation_passed ? "passed" : "failed";

    WriteValidationSummary(cfg, summary);
    WriteMetricsCsv(cfg, summary);

    if (!summary.validation_passed) {
        std::ostringstream oss;
        oss << "[Test_H_CO2_VaryPP] COMSOL validation failed: profile L2=" << summary.final_profile_l2_norm
            << " (threshold=" << cfg.profile_l2_threshold << "), profile Linf=" << summary.final_profile_linf_norm
            << " (threshold=" << cfg.profile_linf_threshold << "), grid_ok=" << BoolString(summary.grid_convergence_ok)
            << ", time_ok=" << BoolString(summary.time_sensitivity_ok);
        throw std::runtime_error(oss.str());
    }
    return summary;
}

TestCasePlan BuildDefaultPlan() {
    TestCasePlan plan;
    plan.plan_key = "h_co2_varypp_nofrac_nowell";
    plan.spec.reference_mode = ReferenceMode::Auto;
    plan.spec.workflow_mode = WorkflowMode::Full;
    return plan;
}

TestCasePlan BuildPrepareReferencePlan() {
    TestCasePlan plan = BuildDefaultPlan();
    plan.plan_key = "h_co2_varypp_nofrac_nowell_prepare_reference";
    plan.spec.workflow_mode = WorkflowMode::PrepareReference;
    return plan;
}

TestCasePlan BuildValidateReferencePlan() {
    TestCasePlan plan = BuildDefaultPlan();
    plan.plan_key = "h_co2_varypp_nofrac_nowell_validate_reference";
    plan.spec.reference_mode = ReferenceMode::Comsol;
    plan.spec.workflow_mode = WorkflowMode::ValidateAgainstReference;
    return plan;
}

using BuilderFn = TestCasePlan(*)();

const std::unordered_map<std::string, BuilderFn>& GetRegistry() {
    static const std::unordered_map<std::string, BuilderFn> registry = {
        {"h_co2_varypp_nofrac_nowell", &BuildDefaultPlan},
        {"h_co2_varypp_nofrac_nowell_prepare_reference", &BuildPrepareReferencePlan},
        {"h_co2_varypp_nofrac_nowell_validate_reference", &BuildValidateReferencePlan}
    };
    return registry;
}

void ExecutePlanByKey(const std::string& key) {
    const auto& registry = GetRegistry();
    const auto it = registry.find(key);
    if (it == registry.end()) throw std::runtime_error("[Test_H_CO2_VaryPP] unknown registry key: " + key);
    const TestCasePlan plan = it->second();
    const TestCaseSummary summary = RunCase(plan.spec);
    std::cout << "\n============================================\n";
    std::cout << "[Test_H_CO2_VaryPP] run completed\n";
    std::cout << "  output_dir: " << summary.case_dir << "\n";
    std::cout << "  grid: " << summary.nx << " x " << summary.ny << " (" << summary.n_cells << " cells)\n";
    std::cout << "  steps: " << summary.steps << "  rollbacks: " << summary.total_rollbacks << "\n";
    std::cout << "  Newton iters: avg=" << std::fixed << std::setprecision(2) << summary.avg_iters
              << "  max=" << summary.max_iters << "\n";
    std::cout << "  final_time: " << std::scientific << std::setprecision(4) << summary.t_end << " s\n";
    std::cout << "  rho range: [" << std::setprecision(4) << summary.rho_min << ", " << summary.rho_max << "] kg/m^3\n";
    std::cout << "  mu  range: [" << summary.mu_min << ", " << summary.mu_max << "] Pa*s\n";
    std::cout << "  validation_status: " << summary.validation_status << "\n";
    if (summary.validation_performed) {
        std::cout << "  final profile errors: L1=" << summary.final_profile_l1_norm
                  << "  L2=" << summary.final_profile_l2_norm
                  << "  Linf=" << summary.final_profile_linf_norm << "\n";
    }
    if (!summary.missing_reference_files.empty()) {
        std::cout << "  missing COMSOL reference files:\n";
        for (const auto& path : summary.missing_reference_files) std::cout << "    - " << path << "\n";
    }
    std::cout << "  validation_summary: " << summary.validation_summary_path << "\n";
    std::cout << "============================================\n";
}

} // namespace

void RunTestCase() {
    ExecutePlanByKey("h_co2_varypp_nofrac_nowell");
}

} // namespace Test_H_CO2_VaryPP

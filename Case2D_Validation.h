#pragma once

#include <string>
#include <vector>

class MeshManager;

namespace Case2DValidation {

// Shared analytical-validation surface for 2D pressure-diffusion donors.

struct PressureDiffusionAnalyticalConfig {
    double lx = 0.0;
    double p_init = 0.0;
    double p_left = 0.0;
    double p_right = 0.0;
    double permeability = 0.0;
    double porosity = 0.0;
    double total_compressibility = 0.0;
    double viscosity = 0.0;
    int analytical_terms = 0;
};

struct AnalyticalSnapshot {
    std::string tag;
    double requested_fraction = 0.0;
    double target_time_s = 0.0;
    double actual_time_s = 0.0;
    bool captured = false;
    std::vector<double> p_blocks;
};

struct AnalyticalMetrics {
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
    int max_err_cell = -1;
    std::string compare_csv_path;
    std::string profile_csv_path;
};

struct PressureDiffusionSummaryReport {
    std::string case_name;
    std::string output_dir;
    int nx = 0;
    int ny = 0;
    int n_cells = 0;
    int steps = 0;
    int total_rollbacks = 0;
    double h_char = 0.0;
    double t_end = 0.0;
    double pressure_diffusivity = 0.0;
    double gravity_x = 0.0;
    double gravity_y = 0.0;
    double gravity_z = 0.0;
    bool non_orthogonal_correction = false;
    double final_l1_norm = 0.0;
    double final_l2_norm = 0.0;
    double final_linf_norm = 0.0;
    double l2_threshold = 0.0;
    double linf_threshold = 0.0;
    bool grid_convergence_ok = true;
    bool time_sensitivity_ok = true;
    bool analytical_passed = false;
    std::string grid_convergence_csv_path;
    std::string time_sensitivity_csv_path;
    std::vector<AnalyticalMetrics> report_metrics;
};

double ClampFraction(double fraction);
std::string MakeReportTag(double fraction);
std::vector<double> BuildSortedFractions(const std::vector<double>& rawFractions);

double ComputePressureDiffusivity(const PressureDiffusionAnalyticalConfig& cfg);
double EvaluatePressureDiffusionAnalyticalPressure(const PressureDiffusionAnalyticalConfig& cfg,
                                                   double x,
                                                   double timeS);

std::vector<AnalyticalSnapshot> BuildAnalyticalSnapshots(double targetEndTimeS,
                                                         const std::vector<double>& reportTimeFractions);

AnalyticalMetrics EvaluatePressureDiffusionAnalyticalMetrics(const MeshManager& mgr,
                                                            const std::vector<double>& pBlocks,
                                                            const PressureDiffusionAnalyticalConfig& cfg,
                                                            const AnalyticalSnapshot& snapshot,
                                                            const std::string& csvPath,
                                                            bool writeCsv);

void WritePressureProfileCSV(const MeshManager& mgr,
                             const std::vector<double>& pBlocks,
                             int xStationCount,
                             const PressureDiffusionAnalyticalConfig& cfg,
                             const AnalyticalSnapshot& snapshot,
                             const std::string& csvPath);

void WritePressureDiffusionSummaryReport(const PressureDiffusionSummaryReport& report,
                                         const std::string& path);

} // namespace Case2DValidation

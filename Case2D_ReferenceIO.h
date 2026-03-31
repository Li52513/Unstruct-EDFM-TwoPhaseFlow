#pragma once

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Case2DReferenceIO {

// Shared path-safe CSV/report I/O surface for 2D validation workflows.

struct CsvTable {
    std::vector<std::string> headers;
    std::unordered_map<std::string, std::size_t> index_by_header;
    std::vector<std::vector<std::string> > rows;
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

std::string Trim(std::string value);
std::vector<std::string> SplitCsvLine(const std::string& line);
CsvTable ReadCsvTable(const std::string& path);
double CsvGetDouble(const CsvTable& table, std::size_t row, const std::string& column);
std::string BoolString(bool value);
std::string PathToGenericString(const std::filesystem::path& path);
std::string BuildOutputOpenDiagnostics(const std::filesystem::path& targetPath);
std::filesystem::path BuildShortAsciiStagingPath(const std::string& targetPath);
std::ofstream OpenAsciiStagingStream(const std::string& targetPath,
                                     const char* context,
                                     std::filesystem::path& stagingPath);
void CommitAsciiStagingFile(std::ofstream& out,
                            const std::filesystem::path& stagingPath,
                            const std::string& targetPath,
                            const char* context);

template <typename Writer>
void WriteAsciiFile(const std::string& path, const char* context, Writer&& writer) {
    std::filesystem::path stagingPath;
    std::ofstream out = OpenAsciiStagingStream(path, context, stagingPath);
    writer(out);
    CommitAsciiStagingFile(out, stagingPath, path, context);
}

ProfileReferenceTable LoadProfileReference(const std::string& path);

template <typename StationContainer, typename Snapshot, typename PressureEvaluator>
void WriteAnalyticalProfileReferenceCsv(const StationContainer& stations,
                                        const Snapshot& snapshot,
                                        PressureEvaluator&& evaluatePressure,
                                        const std::string& family,
                                        const std::string& path) {
    WriteAsciiFile(path,
                   "[Case2DReferenceIO] failed to write analytical profile reference",
                   [&](std::ofstream& out) {
                       out << "station_id,label,family,location,target_axis_m,target_x_m,target_y_m,target_time_s,p_ref_pa\n";
                       for (const auto& station : stations) {
                           if (station.family != family) continue;
                           out << station.id << "," << station.label << "," << station.family << "," << station.location << ","
                               << std::setprecision(12) << station.target_axis_m << "," << station.target_x << "," << station.target_y << ","
                               << snapshot.actual_time_s << ","
                               << evaluatePressure(station.actual_x, snapshot.actual_time_s) << "\n";
                       }
                   });
}

template <typename PointContainer, typename ScheduleContainer, typename PressureEvaluator>
void WriteAnalyticalMonitorReferenceCsv(const PointContainer& points,
                                        const ScheduleContainer& schedule,
                                        PressureEvaluator&& evaluatePressure,
                                        const std::string& path) {
    WriteAsciiFile(path,
                   "[Case2DReferenceIO] failed to write analytical monitor reference",
                   [&](std::ofstream& out) {
                       out << "sample_id,target_time_s";
                       for (const auto& point : points) out << ",p_ref_" << point.label;
                       out << "\n";
                       for (const auto& sample : schedule) {
                           out << sample.sample_id << "," << std::setprecision(12) << sample.actual_time_s;
                           for (const auto& point : points) {
                               out << "," << evaluatePressure(point.actual_x, sample.actual_time_s);
                           }
                           out << "\n";
                       }
                   });
}

template <typename StationContainer>
void WriteProfileStationDefinitions(const StationContainer& stations, const std::string& path) {
    WriteAsciiFile(path,
                   "[Case2DReferenceIO] failed to write profile station definitions",
                   [&](std::ofstream& out) {
                       out << "station_id,label,family,location,target_axis_m,target_x_m,target_y_m,block_id,actual_x_m,actual_y_m\n";
                       for (const auto& station : stations) {
                           out << station.id << "," << station.label << "," << station.family << "," << station.location << ","
                               << std::setprecision(12) << station.target_axis_m << "," << station.target_x << "," << station.target_y << ","
                               << station.block_id << "," << station.actual_x << "," << station.actual_y << "\n";
                       }
                   });
}

template <typename PointContainer>
void WriteMonitorPointDefinitions(const PointContainer& points, const std::string& path) {
    WriteAsciiFile(path,
                   "[Case2DReferenceIO] failed to write monitor point definitions",
                   [&](std::ofstream& out) {
                       out << "point_id,label,location,target_axis_m,target_x_m,target_y_m,block_id,actual_x_m,actual_y_m\n";
                       for (const auto& point : points) {
                           out << point.id << "," << point.label << "," << point.location << ","
                               << std::setprecision(12) << point.target_axis_m << "," << point.target_x << "," << point.target_y << ","
                               << point.block_id << "," << point.actual_x << "," << point.actual_y << "\n";
                       }
                   });
}

template <typename SnapshotContainer>
void WriteProfileSchedule(const SnapshotContainer& snapshots, const std::string& path) {
    WriteAsciiFile(path,
                   "[Case2DReferenceIO] failed to write profile schedule",
                   [&](std::ofstream& out) {
                       out << "tag,requested_fraction,target_time_s,actual_time_s\n";
                       for (const auto& snapshot : snapshots) {
                           out << snapshot.tag << ","
                               << std::setprecision(12) << snapshot.requested_fraction << ","
                               << snapshot.target_time_s << "," << snapshot.actual_time_s << "\n";
                       }
                   });
}

template <typename ScheduleContainer>
void WriteMonitorSchedule(const ScheduleContainer& schedule, const std::string& path) {
    WriteAsciiFile(path,
                   "[Case2DReferenceIO] failed to write monitor schedule",
                   [&](std::ofstream& out) {
                       out << "sample_id,target_time_s,actual_time_s\n";
                       for (const auto& sample : schedule) {
                           out << sample.sample_id << ","
                               << std::setprecision(12) << sample.target_time_s << "," << sample.actual_time_s << "\n";
                       }
                   });
}

template <typename StationContainer, typename Snapshot>
void WriteEngineeringProfileCsv(const StationContainer& stations,
                                const std::string& family,
                                const Snapshot& snapshot,
                                const std::string& path) {
    WriteAsciiFile(path,
                   "[Case2DReferenceIO] failed to write engineering profile csv",
                   [&](std::ofstream& out) {
                       out << "station_id,label,family,location,target_axis_m,target_x_m,target_y_m,block_id,actual_x_m,actual_y_m,target_time_s,actual_time_s,p_num_pa,t_num_k\n";
                       for (const auto& station : stations) {
                           if (station.family != family) continue;
                           if (station.block_id < 0 ||
                               station.block_id >= static_cast<int>(snapshot.p_blocks.size()) ||
                               station.block_id >= static_cast<int>(snapshot.t_blocks.size())) {
                               throw std::runtime_error("[Case2DReferenceIO] invalid block id during engineering profile export.");
                           }
                           out << station.id << "," << station.label << "," << station.family << "," << station.location << ","
                               << std::setprecision(12) << station.target_axis_m << "," << station.target_x << "," << station.target_y << ","
                               << station.block_id << "," << station.actual_x << "," << station.actual_y << ","
                               << snapshot.target_time_s << "," << snapshot.actual_time_s << ","
                               << snapshot.p_blocks[static_cast<std::size_t>(station.block_id)] << ","
                               << snapshot.t_blocks[static_cast<std::size_t>(station.block_id)] << "\n";
                       }
                   });
}

template <typename PointContainer, typename ScheduleContainer>
void WriteEngineeringMonitorCsv(const PointContainer& points,
                                const ScheduleContainer& schedule,
                                const std::string& path) {
    WriteAsciiFile(path,
                   "[Case2DReferenceIO] failed to write engineering monitor csv",
                   [&](std::ofstream& out) {
                       out << "sample_id,target_time_s,actual_time_s";
                       for (const auto& point : points) out << ",p_num_" << point.label;
                       for (const auto& point : points) out << ",t_num_" << point.label;
                       out << "\n";
                       for (const auto& sample : schedule) {
                           out << sample.sample_id << "," << std::setprecision(12) << sample.target_time_s << "," << sample.actual_time_s;
                           for (const auto& point : points) {
                               if (point.block_id < 0 || point.block_id >= static_cast<int>(sample.p_blocks.size())) {
                                   throw std::runtime_error("[Case2DReferenceIO] invalid pressure block id during engineering monitor export.");
                               }
                               out << "," << sample.p_blocks[static_cast<std::size_t>(point.block_id)];
                           }
                           for (const auto& point : points) {
                               if (point.block_id < 0 || point.block_id >= static_cast<int>(sample.t_blocks.size())) {
                                   throw std::runtime_error("[Case2DReferenceIO] invalid temperature block id during engineering monitor export.");
                               }
                               out << "," << sample.t_blocks[static_cast<std::size_t>(point.block_id)];
                           }
                           out << "\n";
                       }
                   });
}

template <typename Config>
void WritePropertyTables(const Config& cfg,
                         const std::string& engineeringPath,
                         const std::string& comsolInputPath) {
    const auto writeOne = [&](const std::string& path) {
        WriteAsciiFile(path,
                       "[Case2DReferenceIO] failed to write property table",
                       [&](std::ofstream& out) {
                           out << "key,value,unit\n";
                           out << "lx," << std::setprecision(12) << cfg.lx << ",m\n";
                           out << "ly," << cfg.ly << ",m\n";
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
                           out << "co2_rho_const," << cfg.co2_rho_const << ",kg/m^3\n";
                           out << "co2_mu_const," << cfg.co2_mu_const << ",Pa*s\n";
                           out << "co2_cp_const," << cfg.co2_cp_const << ",J/(kg*K)\n";
                           out << "co2_cv_const," << cfg.co2_cv_const << ",J/(kg*K)\n";
                           out << "co2_k_const," << cfg.co2_k_const << ",W/(m*K)\n";
                           out << "gravity_x," << cfg.gravity_vector.m_x << ",m/s^2\n";
                           out << "gravity_y," << cfg.gravity_vector.m_y << ",m/s^2\n";
                           out << "gravity_z," << cfg.gravity_vector.m_z << ",m/s^2\n";
                       });
    };
    writeOne(engineeringPath);
    writeOne(comsolInputPath);
}

template <typename Config, typename Summary>
void WriteReferenceSpec(const Config& cfg, const Summary& summary) {
    WriteAsciiFile(summary.reference_spec_path,
                   "[Case2DReferenceIO] failed to write reference spec",
                   [&](std::ofstream& out) {
                       out << "# Engineering Reference Spec\n\n";
                       out << "## Case\n";
                       out << "- Case: `" << cfg.case_name << "`\n";
                       out << "- Geometry: `" << cfg.lx << " m x " << cfg.ly << " m`\n";
                       out << "- Mesh: non-orthogonal quadrilateral, `" << cfg.nx << " x " << cfg.ny << "` nominal split\n";
                       out << "- Fluid model: single-phase CO2, constant properties\n";
                       out << "- Fractures / wells: none\n";
                       out << "- Gravity: disabled\n";
                       out << "- Pressure BC: left/right Dirichlet, top/bottom zero-flux\n";
                       out << "- Temperature BC: left/right Dirichlet, top/bottom zero-flux\n\n";
                       out << "## Validation Strategy\n";
                       out << "- Pressure reference: closed-form 1D Fourier analytical solution at cell centers, profiles, and monitor times\n";
                       out << "- Temperature reference: COMSOL 6.3 matrix-only no-fracture reference exported under `" << summary.comsol_output_dir << "`\n";
                       out << "- Coupled acceptance: pressure analytical gate and temperature COMSOL gate must both pass\n\n";
                       out << "## Inputs Shared With COMSOL\n";
                       out << "- `engineering/profile_station_definitions.csv`\n";
                       out << "- `engineering/profile_report_schedule.csv`\n";
                       out << "- `engineering/monitor_point_definitions.csv`\n";
                       out << "- `engineering/monitor_sample_schedule.csv`\n";
                       out << "- `reference/comsol_input/property_table.csv`\n";
                   });
}

template <typename Config>
void WriteAnalyticalNote(const Config& cfg,
                         double pressureDiffusivity,
                         double pressureScale,
                         const std::string& path) {
    WriteAsciiFile(path,
                   "[Case2DReferenceIO] failed to write analytical note",
                   [&](std::ofstream& out) {
                       out << "# Analytical Pressure Reference\n\n";
                       out << "- Pressure admits a 1D Fourier series solution because the target no-fracture case reduces to x-directed diffusion under constant coefficients and zero gravity.\n";
                       out << "- Diffusivity: `" << pressureDiffusivity << " m^2/s`\n";
                       out << "- Pressure normalization: `Delta p = " << pressureScale << " Pa`\n";
                       out << "- Temperature closed form: not enforced in v1; COMSOL is the authority for temperature.\n";
                   });
}

template <typename Config>
void WriteComsolReferenceSpec(const Config& cfg,
                              const std::string& caseDir,
                              const std::string& path,
                              const std::string& primaryRepresentation,
                              const std::string& manualFallbackRepresentation,
                              const std::string& pdeFallbackRepresentation) {
    WriteAsciiFile(path,
                   "[Case2DReferenceIO] failed to write COMSOL reference spec",
                   [&](std::ofstream& out) {
                       out << "# COMSOL Reference Specification\n\n";
                       out << "## Case\n";
                       out << "- Case directory: `" << caseDir << "`\n";
                       out << "- Geometry: `" << cfg.lx << " m x " << cfg.ly << " m`\n";
                       out << "- Pressure window: `" << cfg.p_left / 1.0e6 << " MPa` -> `" << cfg.p_right / 1.0e6 << " MPa`\n";
                       out << "- Temperature window: `" << cfg.t_left << " K` -> `" << cfg.t_right << " K`\n";
                       out << "- Gravity: disabled\n";
                       out << "- Requested route priority: `Nonisothermal Flow in Porous Media -> dl+ht manual coupling -> PDE fallback`\n";
                       out << "- Formal built-in target: `" << primaryRepresentation << "`\n";
                       out << "- Manual built-in fallback: `" << manualFallbackRepresentation << "`\n";
                       out << "- Last-resort fallback: `" << pdeFallbackRepresentation << "`\n\n";
                       out << "## Expected Outputs\n";
                       out << "- `reference/comsol/comsol_profile_matrix_horizontal_*.csv`\n";
                       out << "- `reference/comsol/comsol_profile_matrix_vertical_midline_*.csv`\n";
                       out << "- `reference/comsol/comsol_monitor_timeseries.csv`\n";
                       out << "- `reference/comsol/comsol_model.mph`\n";
                       out << "- `reference/comsol/comsol_progress.log`\n";
                       out << "- `reference/comsol/comsol_run_summary.md`\n";
                       out << "- `reference/comsol/comsol_reference_mesh_check.txt`\n\n";
                       out << "## Commands\n";
                       out << "```powershell\n";
                       out << "powershell -ExecutionPolicy Bypass -File " << cfg.comsol_wrapper_relpath << " -Mode Compile\n";
                       out << "powershell -ExecutionPolicy Bypass -File " << cfg.comsol_wrapper_relpath << " -Mode Run\n";
                       out << "```\n";
                   });
}

template <typename PointContainer>
MonitorReferenceSeries LoadMonitorReference(const std::string& path, const PointContainer& monitorPoints) {
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

template <typename Summary, typename SnapshotContainer, typename FamilyContainer>
bool ReferenceFilesReady(Summary& summary,
                         const SnapshotContainer& snapshots,
                         const FamilyContainer& orderedFamilies) {
    summary.missing_reference_files.clear();
    std::vector<std::string> required;
    required.push_back(summary.comsol_output_dir + "/comsol_monitor_timeseries.csv");
    required.push_back(summary.comsol_output_dir + "/comsol_model.mph");
    required.push_back(summary.comsol_output_dir + "/comsol_progress.log");
    required.push_back(summary.comsol_output_dir + "/comsol_run_summary.md");
    required.push_back(summary.comsol_output_dir + "/comsol_reference_mesh_check.txt");
    for (const auto& snapshot : snapshots) {
        for (const auto& family : orderedFamilies) {
            required.push_back(summary.comsol_output_dir + "/comsol_profile_" + family + "_" + snapshot.tag + ".csv");
        }
    }
    for (const auto& path : required) {
        std::ifstream in(std::filesystem::path(path), std::ios::in | std::ios::binary);
        if (!in.good()) summary.missing_reference_files.push_back(path);
    }
    summary.reference_ready = summary.missing_reference_files.empty();
    return summary.reference_ready;
}

template <typename Config, typename Summary>
void WriteMetricsCsv(const Config& cfg, const Summary& summary) {
    WriteAsciiFile(summary.metrics_csv_path,
                   "[Case2DReferenceIO] failed to write metrics csv",
                   [&](std::ofstream& out) {
                       out << "case_name,nx,ny,n_cells,h_char,t_end,steps,total_rollbacks,avg_nonlinear_iters,max_nonlinear_iters,"
                              "reference_ready,validation_performed,validation_passed,grid_convergence_ok,time_sensitivity_ok,"
                              "final_pressure_cell_l1_norm,final_pressure_cell_l2_norm,final_pressure_cell_linf_norm,"
                              "final_pressure_horizontal_l2_norm,final_pressure_horizontal_linf_norm,"
                              "final_temperature_horizontal_l2_norm,final_temperature_horizontal_linf_norm,"
                              "final_monitor_pressure_l2_norm,final_monitor_pressure_linf_norm,"
                              "final_monitor_temperature_l2_norm,final_monitor_temperature_linf_norm,"
                              "final_pressure_vertical_uniformity_norm_std,final_temperature_vertical_uniformity_norm_std,"
                              "grid_pressure_order_min,grid_temperature_order_min,time_pressure_order_min,time_temperature_order_min,"
                              "comsol_representation,validation_status\n";
                       out << cfg.case_name << "," << cfg.nx << "," << cfg.ny << "," << summary.n_cells << ","
                           << std::setprecision(12) << summary.h_char << "," << summary.t_end << ","
                           << summary.steps << "," << summary.total_rollbacks << "," << summary.avg_iters << "," << summary.max_iters << ","
                           << BoolString(summary.reference_ready) << "," << BoolString(summary.validation_performed) << ","
                           << BoolString(summary.validation_passed) << "," << BoolString(summary.grid_convergence_ok) << ","
                           << BoolString(summary.time_sensitivity_ok) << ","
                           << summary.final_pressure_cell_l1_norm << "," << summary.final_pressure_cell_l2_norm << ","
                           << summary.final_pressure_cell_linf_norm << "," << summary.final_pressure_horizontal_l2_norm << ","
                           << summary.final_pressure_horizontal_linf_norm << "," << summary.final_temperature_horizontal_l2_norm << ","
                           << summary.final_temperature_horizontal_linf_norm << "," << summary.final_monitor_pressure_l2_norm << ","
                           << summary.final_monitor_pressure_linf_norm << "," << summary.final_monitor_temperature_l2_norm << ","
                           << summary.final_monitor_temperature_linf_norm << "," << summary.final_pressure_vertical_uniformity_norm_std << ","
                           << summary.final_temperature_vertical_uniformity_norm_std << "," << summary.grid_pressure_order_min << ","
                           << summary.grid_temperature_order_min << "," << summary.time_pressure_order_min << ","
                           << summary.time_temperature_order_min << "," << summary.comsol_representation << ","
                           << summary.validation_status << "\n";
                   });
}

template <typename Config, typename Summary>
void WriteValidationSummaryCsv(const Config& cfg,
                               const Summary& summary,
                               const std::string& gridOrderSemantics,
                               const std::string& timeOrderSemantics) {
    WriteAsciiFile(summary.validation_summary_csv_path,
                   "[Case2DReferenceIO] failed to write validation summary csv",
                   [&](std::ofstream& out) {
                       out << "case_name,reference_mode,validation_status,reference_ready,validation_performed,validation_passed,"
                              "grid_convergence_ok,time_sensitivity_ok,final_pressure_cell_l2_norm,final_pressure_cell_linf_norm,"
                              "final_temperature_horizontal_l2_norm,final_temperature_horizontal_linf_norm,"
                              "final_monitor_temperature_l2_norm,final_monitor_temperature_linf_norm,"
                              "final_pressure_vertical_uniformity_norm_std,final_temperature_vertical_uniformity_norm_std,"
                              "grid_pressure_order_min,grid_temperature_order_min,time_pressure_order_min,time_temperature_order_min,comsol_representation,"
                              "grid_order_semantics,time_order_semantics\n";
                       out << cfg.case_name << "," << summary.resolved_reference_mode << "," << summary.validation_status << ","
                           << BoolString(summary.reference_ready) << "," << BoolString(summary.validation_performed) << ","
                           << BoolString(summary.validation_passed) << "," << BoolString(summary.grid_convergence_ok) << ","
                           << BoolString(summary.time_sensitivity_ok) << "," << summary.final_pressure_cell_l2_norm << ","
                           << summary.final_pressure_cell_linf_norm << "," << summary.final_temperature_horizontal_l2_norm << ","
                           << summary.final_temperature_horizontal_linf_norm << "," << summary.final_monitor_temperature_l2_norm << ","
                           << summary.final_monitor_temperature_linf_norm << "," << summary.final_pressure_vertical_uniformity_norm_std << ","
                           << summary.final_temperature_vertical_uniformity_norm_std << "," << summary.grid_pressure_order_min << ","
                           << summary.grid_temperature_order_min << "," << summary.time_pressure_order_min << ","
                           << summary.time_temperature_order_min << "," << summary.comsol_representation << ","
                           << gridOrderSemantics << "," << timeOrderSemantics << "\n";
                   });
}

template <typename Config, typename Summary>
void WriteValidationSummary(const Config& cfg,
                            const Summary& summary,
                            const std::string& gridOrderSemantics,
                            const std::string& timeOrderSemantics) {
    WriteAsciiFile(summary.validation_summary_path,
                   "[Case2DReferenceIO] failed to write validation summary",
                   [&](std::ofstream& out) {
                       out << "# Validation Summary\n\n";
                       out << "## Status\n";
                       out << "- Case: `" << cfg.case_name << "`\n";
                       out << "- Output directory: `" << summary.case_dir << "`\n";
                       out << "- Reference mode: `" << summary.resolved_reference_mode << "`\n";
                       out << "- Validation status: `" << summary.validation_status << "`\n";
                       out << "- Reference ready: `" << BoolString(summary.reference_ready) << "`\n";
                       out << "- Validation performed: `" << BoolString(summary.validation_performed) << "`\n";
                       out << "- Validation passed: `" << BoolString(summary.validation_passed) << "`\n";
                       out << "- Grid convergence ok: `" << BoolString(summary.grid_convergence_ok) << "`\n";
                       out << "- Time sensitivity ok: `" << BoolString(summary.time_sensitivity_ok) << "`\n\n";

                       out << "## Physical Alignment\n";
                       out << "- Gravity vector: `(" << cfg.gravity_vector.m_x << ", " << cfg.gravity_vector.m_y << ", " << cfg.gravity_vector.m_z << ")`\n";
                       out << "- Non-orthogonal correction: `" << BoolString(cfg.enable_non_orthogonal_correction) << "`\n";
                       out << "- Pressure reference: analytical 1D Fourier diffusion\n";
                       out << "- Temperature reference: COMSOL 6.3 `" << summary.comsol_representation << "`\n";
                       out << "- COMSOL fine check skipped: `" << BoolString(summary.comsol_fine_check_skipped) << "`\n\n";

                       out << "## Final Acceptance Metrics\n";
                       out << "- Pressure cell `L2_norm`: `" << summary.final_pressure_cell_l2_norm << "`\n";
                       out << "- Pressure cell `Linf_norm`: `" << summary.final_pressure_cell_linf_norm << "`\n";
                       out << "- Pressure horizontal profile `L2_norm`: `" << summary.final_pressure_horizontal_l2_norm << "`\n";
                       out << "- Temperature horizontal profile `L2_norm`: `" << summary.final_temperature_horizontal_l2_norm << "`\n";
                       out << "- Temperature horizontal profile `Linf_norm`: `" << summary.final_temperature_horizontal_linf_norm << "`\n";
                       out << "- Temperature monitor `L2_norm`: `" << summary.final_monitor_temperature_l2_norm << "`\n";
                       out << "- Temperature monitor `Linf_norm`: `" << summary.final_monitor_temperature_linf_norm << "`\n";
                       out << "- Pressure vertical uniformity norm-std: `" << summary.final_pressure_vertical_uniformity_norm_std << "`\n";
                       out << "- Temperature vertical uniformity norm-std: `" << summary.final_temperature_vertical_uniformity_norm_std << "`\n";
                       out << "- Thresholds: pressure `L2 <= " << cfg.pressure_l2_threshold << "`, pressure `Linf <= " << cfg.pressure_linf_threshold
                           << "`, temperature `L2 <= " << cfg.temperature_l2_threshold << "`, temperature `Linf <= " << cfg.temperature_linf_threshold
                           << "`, vertical uniformity `<= " << cfg.vertical_uniformity_threshold << "`\n\n";

                       out << "## Convergence Semantics\n";
                       out << "- Baseline absolute validation: pressure uses analytical reference at cell centers / fixed target points; temperature uses COMSOL fixed-target reference.\n";
                       out << "- Grid order: " << gridOrderSemantics << "\n";
                       out << "- Time order: " << timeOrderSemantics << "\n\n";

                       out << "## Outputs\n";
                       out << "- Engineering dir: `" << summary.engineering_dir << "`\n";
                       out << "- Analytical dir: `" << summary.analytic_dir << "`\n";
                       out << "- COMSOL input dir: `" << summary.comsol_input_dir << "`\n";
                       out << "- COMSOL output dir: `" << summary.comsol_output_dir << "`\n";
                       out << "- Reference spec: `" << summary.reference_spec_path << "`\n";
                       out << "- Analytical note: `" << summary.analytical_note_path << "`\n";
                       out << "- COMSOL spec: `" << summary.comsol_reference_spec_path << "`\n";
                       out << "- MATLAB script: `" << summary.matlab_script_path << "`\n";
                       out << "- Property table: `" << summary.property_table_path << "`\n";
                       out << "- COMSOL property table: `" << summary.comsol_property_table_path << "`\n";
                       if (!summary.grid_convergence_csv_path.empty()) out << "- Grid convergence csv: `" << summary.grid_convergence_csv_path << "`\n";
                       if (!summary.time_sensitivity_csv_path.empty()) out << "- Time sensitivity csv: `" << summary.time_sensitivity_csv_path << "`\n";

                       if (!summary.pressure_cell_metrics.empty()) {
                           out << "\n## Pressure Cell Records\n";
                           for (const auto& metric : summary.pressure_cell_metrics) {
                               out << "- `" << metric.tag << "`: L2=`" << metric.pressure.l2_norm
                                   << "`, Linf=`" << metric.pressure.linf_norm
                                   << "`, csv=`" << metric.compare_csv_path << "`\n";
                           }
                       }

                       if (!summary.report_metrics.empty()) {
                           out << "\n## Profile Compare Records\n";
                           for (const auto& metric : summary.report_metrics) {
                               out << "- `" << metric.family << " / " << metric.tag
                                   << "`: pL2=`" << metric.pressure.l2_norm
                                   << "`, tL2=`" << metric.temperature.l2_norm
                                   << "`, pUniformity=`" << metric.pressure_uniformity_norm_std
                                   << "`, tUniformity=`" << metric.temperature_uniformity_norm_std
                                   << "`, csv=`" << metric.compare_csv_path << "`\n";
                           }
                       }

                       if (!summary.missing_reference_files.empty()) {
                           out << "\n## Missing Reference Files\n";
                           for (const auto& missingPath : summary.missing_reference_files) out << "- `" << missingPath << "`\n";
                       }
                   });
}

} // namespace Case2DReferenceIO

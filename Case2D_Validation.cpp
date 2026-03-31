#include "Case2D_Validation.h"

#include "MeshManager.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace Case2DValidation {
namespace {

constexpr double kPi = 3.1415926535897932384626433832795;

std::string BoolString(bool value) {
    return value ? "true" : "false";
}

} // namespace

double ClampFraction(double fraction) {
    return std::max(0.0, std::min(1.0, fraction));
}

std::string MakeReportTag(double fraction) {
    std::ostringstream oss;
    oss << "t" << std::setw(3) << std::setfill('0')
        << static_cast<int>(std::round(ClampFraction(fraction) * 100.0)) << "pct";
    return oss.str();
}

std::vector<double> BuildSortedFractions(const std::vector<double>& rawFractions) {
    std::vector<double> out;
    for (double fraction : rawFractions) {
        const double clamped = ClampFraction(fraction);
        bool duplicate = false;
        for (double existing : out) {
            if (std::abs(existing - clamped) <= 1.0e-12) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            out.push_back(clamped);
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

double ComputePressureDiffusivity(const PressureDiffusionAnalyticalConfig& cfg) {
    const double denom = cfg.viscosity * cfg.porosity * cfg.total_compressibility;
    if (cfg.permeability <= 0.0 || denom <= 0.0) {
        throw std::runtime_error("[Case2D_Validation] invalid constant-property coefficients for analytical diffusivity.");
    }
    return cfg.permeability / denom;
}

double EvaluatePressureDiffusionAnalyticalPressure(const PressureDiffusionAnalyticalConfig& cfg,
                                                   double x,
                                                   double timeS) {
    if (cfg.lx <= 0.0) {
        throw std::runtime_error("[Case2D_Validation] invalid domain length.");
    }
    if (timeS <= 0.0) {
        return cfg.p_init;
    }

    const double xClamped = std::max(0.0, std::min(cfg.lx, x));
    const double deltaP = cfg.p_left - cfg.p_right;
    const double pSteady = cfg.p_left - deltaP * (xClamped / cfg.lx);
    const double a0 = cfg.p_init - cfg.p_left;
    const double diffusivity = ComputePressureDiffusivity(cfg);

    double transient = 0.0;
    for (int n = 1; n <= cfg.analytical_terms; ++n) {
        const double signN = (n % 2 == 0) ? 1.0 : -1.0;
        const double bn = (2.0 / (static_cast<double>(n) * kPi))
            * (a0 * (1.0 - signN) - deltaP * signN);
        const double lambda = static_cast<double>(n) * kPi / cfg.lx;
        transient += bn * std::sin(lambda * xClamped) * std::exp(-diffusivity * lambda * lambda * timeS);
    }
    return pSteady + transient;
}

std::vector<AnalyticalSnapshot> BuildAnalyticalSnapshots(double targetEndTimeS,
                                                         const std::vector<double>& reportTimeFractions) {
    std::vector<AnalyticalSnapshot> snapshots;
    const std::vector<double> fractions = BuildSortedFractions(reportTimeFractions);
    for (double fraction : fractions) {
        AnalyticalSnapshot snapshot;
        snapshot.tag = MakeReportTag(fraction);
        snapshot.requested_fraction = fraction;
        snapshot.target_time_s = fraction * targetEndTimeS;
        snapshots.push_back(snapshot);
    }
    return snapshots;
}

AnalyticalMetrics EvaluatePressureDiffusionAnalyticalMetrics(const MeshManager& mgr,
                                                            const std::vector<double>& pBlocks,
                                                            const PressureDiffusionAnalyticalConfig& cfg,
                                                            const AnalyticalSnapshot& snapshot,
                                                            const std::string& csvPath,
                                                            bool writeCsv) {
    const auto& cells = mgr.mesh().getCells();
    const std::size_t nUse = std::min(cells.size(), pBlocks.size());
    if (nUse == 0) {
        throw std::runtime_error("[Case2D_Validation] empty state for analytical comparison.");
    }

    const double deltaPAbs = std::max(std::abs(cfg.p_left - cfg.p_right), 1.0);
    double sumAbs = 0.0;
    double sumSq = 0.0;
    double maxAbs = 0.0;
    int maxErrCell = -1;

    std::ofstream csv;
    if (writeCsv) {
        csv.open(csvPath.c_str(), std::ios::out | std::ios::trunc);
        if (!csv.good()) {
            throw std::runtime_error("[Case2D_Validation] failed to open analytical csv: " + csvPath);
        }
        csv << "cell_id,x,y,p_num,p_ana,abs_err,abs_err_over_dP\n";
    }

    for (std::size_t i = 0; i < nUse; ++i) {
        const double x = cells[i].center.m_x;
        const double y = cells[i].center.m_y;
        const double pNum = pBlocks[i];
        const double pAna = EvaluatePressureDiffusionAnalyticalPressure(cfg, x, snapshot.actual_time_s);
        const double absErr = std::abs(pNum - pAna);
        sumAbs += absErr;
        sumSq += absErr * absErr;
        if (absErr > maxAbs) {
            maxAbs = absErr;
            maxErrCell = static_cast<int>(i);
        }

        if (writeCsv) {
            csv << i << "," << std::setprecision(12) << x << "," << y << ","
                << pNum << "," << pAna << "," << absErr << "," << (absErr / deltaPAbs) << "\n";
        }
    }

    AnalyticalMetrics metrics;
    metrics.tag = snapshot.tag;
    metrics.requested_fraction = snapshot.requested_fraction;
    metrics.target_time_s = snapshot.target_time_s;
    metrics.actual_time_s = snapshot.actual_time_s;
    metrics.sample_count = static_cast<int>(nUse);
    metrics.l1_abs = sumAbs / static_cast<double>(nUse);
    metrics.l2_abs = std::sqrt(sumSq / static_cast<double>(nUse));
    metrics.linf_abs = maxAbs;
    metrics.l1_norm = metrics.l1_abs / deltaPAbs;
    metrics.l2_norm = metrics.l2_abs / deltaPAbs;
    metrics.linf_norm = metrics.linf_abs / deltaPAbs;
    metrics.max_err_cell = maxErrCell;
    metrics.compare_csv_path = writeCsv ? csvPath : "";
    return metrics;
}

void WritePressureProfileCSV(const MeshManager& mgr,
                             const std::vector<double>& pBlocks,
                             int xStationCount,
                             const PressureDiffusionAnalyticalConfig& cfg,
                             const AnalyticalSnapshot& snapshot,
                             const std::string& csvPath) {
    if (xStationCount <= 0) {
        return;
    }

    struct Bucket {
        double x_sum = 0.0;
        double p_sum = 0.0;
        int count = 0;
    };

    const auto& cells = mgr.mesh().getCells();
    const std::size_t nUse = std::min(cells.size(), pBlocks.size());
    std::vector<Bucket> buckets(static_cast<std::size_t>(xStationCount));

    for (std::size_t i = 0; i < nUse; ++i) {
        int ix = (cfg.lx > 0.0)
            ? static_cast<int>(std::floor((cells[i].center.m_x / cfg.lx) * static_cast<double>(xStationCount)))
            : 0;
        ix = std::max(0, std::min(xStationCount - 1, ix));
        Bucket& bucket = buckets[static_cast<std::size_t>(ix)];
        bucket.x_sum += cells[i].center.m_x;
        bucket.p_sum += pBlocks[i];
        bucket.count += 1;
    }

    std::ofstream csv(csvPath.c_str(), std::ios::out | std::ios::trunc);
    if (!csv.good()) {
        throw std::runtime_error("[Case2D_Validation] failed to open profile csv: " + csvPath);
    }

    csv << "station_id,x,p_num_avg,p_ana,abs_err,abs_err_over_dP,n_cells\n";
    const double deltaPAbs = std::max(std::abs(cfg.p_left - cfg.p_right), 1.0);
    for (int ix = 0; ix < xStationCount; ++ix) {
        const Bucket& bucket = buckets[static_cast<std::size_t>(ix)];
        if (bucket.count <= 0) {
            continue;
        }
        const double xMean = bucket.x_sum / static_cast<double>(bucket.count);
        const double pNumAvg = bucket.p_sum / static_cast<double>(bucket.count);
        const double pAna = EvaluatePressureDiffusionAnalyticalPressure(cfg, xMean, snapshot.actual_time_s);
        const double absErr = std::abs(pNumAvg - pAna);
        csv << ix << "," << std::setprecision(12) << xMean << "," << pNumAvg << ","
            << pAna << "," << absErr << "," << (absErr / deltaPAbs) << "," << bucket.count << "\n";
    }
}

void WritePressureDiffusionSummaryReport(const PressureDiffusionSummaryReport& report,
                                         const std::string& path) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) {
        throw std::runtime_error("[Case2D_Validation] failed to open analytical summary: " + path);
    }

    out << "case=" << report.case_name << "\n";
    out << "output_dir=" << report.output_dir << "\n";
    out << "mesh=" << report.nx << "x" << report.ny << "\n";
    out << "cells=" << report.n_cells << "\n";
    out << "h_char=" << std::setprecision(12) << report.h_char << "\n";
    out << "t_end=" << report.t_end << "\n";
    out << "steps=" << report.steps << "\n";
    out << "pressure_diffusivity=" << report.pressure_diffusivity << "\n";
    out << "gravity_vector=(" << report.gravity_x << "," << report.gravity_y << "," << report.gravity_z << ")\n";
    out << "non_orthogonal_correction=" << BoolString(report.non_orthogonal_correction) << "\n";
    out << "final_l1_norm=" << report.final_l1_norm << "\n";
    out << "final_l2_norm=" << report.final_l2_norm << "\n";
    out << "final_linf_norm=" << report.final_linf_norm << "\n";
    out << "l2_threshold=" << report.l2_threshold << "\n";
    out << "linf_threshold=" << report.linf_threshold << "\n";
    out << "grid_convergence_ok=" << BoolString(report.grid_convergence_ok) << "\n";
    out << "time_sensitivity_ok=" << BoolString(report.time_sensitivity_ok) << "\n";
    out << "analytical_passed=" << BoolString(report.analytical_passed) << "\n";
    if (!report.grid_convergence_csv_path.empty()) {
        out << "grid_convergence_csv=" << report.grid_convergence_csv_path << "\n";
    }
    if (!report.time_sensitivity_csv_path.empty()) {
        out << "time_sensitivity_csv=" << report.time_sensitivity_csv_path << "\n";
    }
    out << "\n[report_metrics]\n";
    for (const auto& metrics : report.report_metrics) {
        out << metrics.tag
            << " requested_fraction=" << metrics.requested_fraction
            << " target_time_s=" << metrics.target_time_s
            << " actual_time_s=" << metrics.actual_time_s
            << " l1_norm=" << metrics.l1_norm
            << " l2_norm=" << metrics.l2_norm
            << " linf_norm=" << metrics.linf_norm
            << " max_err_cell=" << metrics.max_err_cell << "\n";
        if (!metrics.compare_csv_path.empty()) {
            out << "compare_csv=" << metrics.compare_csv_path << "\n";
        }
        if (!metrics.profile_csv_path.empty()) {
            out << "profile_csv=" << metrics.profile_csv_path << "\n";
        }
    }
}

} // namespace Case2DValidation

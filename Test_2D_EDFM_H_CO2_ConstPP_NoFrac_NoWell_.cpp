/**
 * @file  Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.cpp
 * @brief 独立测试文件�?D 单相 CO2 常物性、无裂缝、无井的压力扩散解析验证
 */

#include "Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.h"

#include "2D_PostProcess.h"
#include "BoundaryConditionManager.h"
#include "FIM_TransientCaseKit.hpp"
#include "MeshDefinitions.h"
#include "MeshManager.h"
#include "SolverContrlStrName_op.h"
#include "Well_WellControlTypes.h"

#include <algorithm>
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

namespace Test_H_CO2_ConstPP {
namespace {

constexpr double kPi = 3.1415926535897932384626433832795;

void EnsureDirRecursive(const std::string& rawPath) {
    if (rawPath.empty()) return;
    std::string path = rawPath;
    for (char& ch : path) if (ch == '\\') ch = '/';
    std::stringstream ss(path);
    std::string token, current;
    while (std::getline(ss, token, '/')) {
        if (token.empty() || token == ".") continue;
        if (!current.empty()) current += "/";
        current += token;
        TEST_MKDIR(current.c_str());
    }
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

std::string BoolString(bool value) { return value ? "true" : "false"; }
double ClampFraction(double f) { return std::max(0.0, std::min(1.0, f)); }

std::string MakeReportTag(double fraction) {
    std::ostringstream oss;
    oss << "t" << std::setw(3) << std::setfill('0')
        << static_cast<int>(std::round(ClampFraction(fraction) * 100.0)) << "pct";
    return oss.str();
}

std::vector<double> BuildSortedFractions(const std::vector<double>& raw) {
    std::vector<double> out;
    for (double f : raw) {
        const double clamped = ClampFraction(f);
        bool duplicate = false;
        for (double existing : out) {
            if (std::abs(existing - clamped) <= 1.0e-12) { duplicate = true; break; }
        }
        if (!duplicate) out.push_back(clamped);
    }
    std::sort(out.begin(), out.end());
    return out;
}

struct TestCaseSpec {
    std::string case_name = "h_co2_constpp_nofrac_nowell";
    std::string output_base_dir = "Test/Transient/FullCaseTest";
    std::string sub_dir = "H_CO2_ConstPP";
    double lx = 400.0, ly = 40.0;
    int nx = 48, ny = 6;
    double matrix_phi = 0.10, matrix_perm = 1.0e-13, matrix_ct = 5.0e-9;
    double mu_const = 6.0e-5, rho_const = 700.0;
    double p_init = 10.0e6, p_left = 12.0e6, p_right = 8.0e6, t_init = 360.0;
    double dt_init = 120.0, dt_min = 1.0, dt_max = 5000.0, target_end_time_s = 1.0e5;
    int max_steps = 12000, max_newton_iter = 12;
    FIM_Engine::LinearSolverType lin_solver = FIM_Engine::LinearSolverType::AMGCL;
    bool amgcl_use_fallback_sparselu = true;
    double amgcl_tol = 1.0e-6;
    int amgcl_maxiter = 500;
    bool enable_non_orthogonal_correction = true;
    bool enable_row_scaling = true;
    double abs_res_tol = 1.0e-10, rel_res_tol = 1.0e-6, rel_update_tol = 1.0e-8;
    double max_dP = 2.0e7;
    bool enable_armijo_line_search = false;
    double rollback_shrink_factor = 0.7, dt_relres_grow_factor = 1.08;
    Vector gravity_vector = Vector(0.0, 0.0, 0.0);
    FIM_Engine::DiagLevel diag_level = FIM_Engine::DiagLevel::Off;
    bool enable_analytical_validation = true;
    int analytical_terms = 200;
    std::vector<double> report_time_fractions = {0.1, 0.5, 1.0};
    double analytical_l2_threshold = 3.0e-2, analytical_linf_threshold = 5.0e-2;
    bool enable_profile_output = true;
    bool enable_grid_convergence_study = true;
    bool enable_time_sensitivity_study = true;
    std::vector<std::pair<int, int> > grid_sweep_cases = {
        std::make_pair(24, 3), std::make_pair(48, 6), std::make_pair(96, 12)
    };
    std::vector<double> time_step_sweep = {480.0, 120.0, 30.0};
    bool export_vtk = true;
    bool emit_detailed_outputs = true;
};

struct AnalyticalSnapshot {
    std::string tag;
    double requested_fraction = 0.0, target_time_s = 0.0, actual_time_s = 0.0;
    bool captured = false;
    std::vector<double> p_blocks;
};

struct AnalyticalMetrics {
    std::string tag;
    double requested_fraction = 0.0, target_time_s = 0.0, actual_time_s = 0.0;
    int sample_count = 0;
    double l1_abs = 0.0, l2_abs = 0.0, linf_abs = 0.0;
    double l1_norm = 0.0, l2_norm = 0.0, linf_norm = 0.0;
    int max_err_cell = -1;
    std::string compare_csv_path, profile_csv_path;
};

struct SweepStudyRow {
    std::string label, case_dir;
    int nx = 0, ny = 0, steps = 0;
    double dt_init = 0.0, h_char = 0.0, t_end = 0.0;
    double l1_norm = 0.0, l2_norm = 0.0, linf_norm = 0.0;
};

struct TestCaseSummary {
    std::string case_dir, convergence_log_path, metrics_csv_path, analytical_summary_path;
    std::string grid_convergence_csv_path, time_sensitivity_csv_path;
    int nx = 0, ny = 0, n_cells = 0, steps = 0, total_rollbacks = 0, max_iters = 0;
    double h_char = 0.0, avg_iters = 0.0, t_end = 0.0;
    double final_l1_abs = 0.0, final_l2_abs = 0.0, final_linf_abs = 0.0;
    double final_l1_norm = 0.0, final_l2_norm = 0.0, final_linf_norm = 0.0;
    bool analytical_passed = false, grid_convergence_ok = true, time_sensitivity_ok = true;
    std::vector<AnalyticalMetrics> report_metrics;
};

struct TestCasePlan { std::string plan_key; TestCaseSpec spec; };

double ComputePressureDiffusivity(const TestCaseSpec& cfg) {
    const double denom = cfg.mu_const * cfg.matrix_phi * cfg.matrix_ct;
    if (cfg.matrix_perm <= 0.0 || denom <= 0.0) {
        throw std::runtime_error("[Test_H_CO2] invalid constant-property coefficients for analytical diffusivity.");
    }
    return cfg.matrix_perm / denom;
}

double EvaluateAnalyticalPressure(const TestCaseSpec& cfg, double x, double timeS) {
    if (cfg.lx <= 0.0) throw std::runtime_error("[Test_H_CO2] invalid domain length.");
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

std::vector<AnalyticalSnapshot> BuildSnapshots(const TestCaseSpec& cfg) {
    std::vector<AnalyticalSnapshot> snapshots;
    if (!cfg.enable_analytical_validation) return snapshots;
    const std::vector<double> fractions = BuildSortedFractions(cfg.report_time_fractions);
    for (double f : fractions) {
        AnalyticalSnapshot s;
        s.tag = MakeReportTag(f);
        s.requested_fraction = f;
        s.target_time_s = f * cfg.target_end_time_s;
        snapshots.push_back(s);
    }
    return snapshots;
}

AnalyticalMetrics EvaluateAnalyticalMetrics(const MeshManager& mgr, const std::vector<double>& pBlocks,
                                            const TestCaseSpec& cfg, const AnalyticalSnapshot& snapshot,
                                            const std::string& csvPath, bool writeCsv) {
    const auto& cells = mgr.mesh().getCells();
    const std::size_t nUse = std::min(cells.size(), pBlocks.size());
    if (nUse == 0) throw std::runtime_error("[Test_H_CO2] empty state for analytical comparison.");
    const double deltaPAbs = std::max(std::abs(cfg.p_left - cfg.p_right), 1.0);
    double sumAbs = 0.0, sumSq = 0.0, maxAbs = 0.0;
    int maxErrCell = -1;
    std::ofstream csv;
    if (writeCsv) {
        csv.open(csvPath, std::ios::out | std::ios::trunc);
        if (!csv.good()) throw std::runtime_error("[Test_H_CO2] failed to open analytical csv: " + csvPath);
        csv << "cell_id,x,y,p_num,p_ana,abs_err,abs_err_over_dP\n";
    }
    for (std::size_t i = 0; i < nUse; ++i) {
        const double x = cells[i].center.m_x;
        const double y = cells[i].center.m_y;
        const double pNum = pBlocks[i];
        const double pAna = EvaluateAnalyticalPressure(cfg, x, snapshot.actual_time_s);
        const double absErr = std::abs(pNum - pAna);
        sumAbs += absErr;
        sumSq += absErr * absErr;
        if (absErr > maxAbs) { maxAbs = absErr; maxErrCell = static_cast<int>(i); }
        if (writeCsv) {
            csv << i << "," << std::setprecision(12) << x << "," << y << ","
                << pNum << "," << pAna << "," << absErr << "," << (absErr / deltaPAbs) << "\n";
        }
    }
    AnalyticalMetrics m;
    m.tag = snapshot.tag;
    m.requested_fraction = snapshot.requested_fraction;
    m.target_time_s = snapshot.target_time_s;
    m.actual_time_s = snapshot.actual_time_s;
    m.sample_count = static_cast<int>(nUse);
    m.l1_abs = sumAbs / static_cast<double>(nUse);
    m.l2_abs = std::sqrt(sumSq / static_cast<double>(nUse));
    m.linf_abs = maxAbs;
    m.l1_norm = m.l1_abs / deltaPAbs;
    m.l2_norm = m.l2_abs / deltaPAbs;
    m.linf_norm = m.linf_abs / deltaPAbs;
    m.max_err_cell = maxErrCell;
    m.compare_csv_path = writeCsv ? csvPath : "";
    return m;
}

void WriteProfileCSV(const MeshManager& mgr, const std::vector<double>& pBlocks,
                     const TestCaseSpec& cfg, const AnalyticalSnapshot& snapshot,
                     const std::string& csvPath) {
    if (cfg.nx <= 0) return;
    struct Bucket { double x_sum = 0.0, p_sum = 0.0; int count = 0; };
    const auto& cells = mgr.mesh().getCells();
    const std::size_t nUse = std::min(cells.size(), pBlocks.size());
    std::vector<Bucket> buckets(static_cast<std::size_t>(cfg.nx));
    for (std::size_t i = 0; i < nUse; ++i) {
        int ix = (cfg.lx > 0.0)
            ? static_cast<int>(std::floor((cells[i].center.m_x / cfg.lx) * static_cast<double>(cfg.nx)))
            : 0;
        ix = std::max(0, std::min(cfg.nx - 1, ix));
        Bucket& b = buckets[static_cast<std::size_t>(ix)];
        b.x_sum += cells[i].center.m_x;
        b.p_sum += pBlocks[i];
        b.count += 1;
    }
    std::ofstream csv(csvPath, std::ios::out | std::ios::trunc);
    if (!csv.good()) throw std::runtime_error("[Test_H_CO2] failed to open profile csv: " + csvPath);
    csv << "station_id,x,p_num_avg,p_ana,abs_err,abs_err_over_dP,n_cells\n";
    const double deltaPAbs = std::max(std::abs(cfg.p_left - cfg.p_right), 1.0);
    for (int ix = 0; ix < cfg.nx; ++ix) {
        const Bucket& b = buckets[static_cast<std::size_t>(ix)];
        if (b.count <= 0) continue;
        const double xMean = b.x_sum / static_cast<double>(b.count);
        const double pNumAvg = b.p_sum / static_cast<double>(b.count);
        const double pAna = EvaluateAnalyticalPressure(cfg, xMean, snapshot.actual_time_s);
        const double absErr = std::abs(pNumAvg - pAna);
        csv << ix << "," << std::setprecision(12) << xMean << "," << pNumAvg << ","
            << pAna << "," << absErr << "," << (absErr / deltaPAbs) << "," << b.count << "\n";
    }
}

void WriteStudyCSV(const std::vector<SweepStudyRow>& rows, const std::string& csvPath, const std::string& study) {
    std::ofstream csv(csvPath, std::ios::out | std::ios::trunc);
    if (!csv.good()) throw std::runtime_error("[Test_H_CO2] failed to open study csv: " + csvPath);
    csv << "study,label,case_dir,nx,ny,dt_init,h_char,t_end,steps,l1_norm,l2_norm,linf_norm\n";
    for (const auto& row : rows) {
        csv << study << "," << row.label << "," << row.case_dir << "," << row.nx << "," << row.ny << ","
            << std::setprecision(12) << row.dt_init << "," << row.h_char << "," << row.t_end << ","
            << row.steps << "," << row.l1_norm << "," << row.l2_norm << "," << row.linf_norm << "\n";
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

void WriteAnalyticalSummary(const TestCaseSpec& cfg, const TestCaseSummary& summary) {
    if (summary.case_dir.empty()) return;
    const std::string path = summary.case_dir + "/analytical_summary.txt";
    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("[Test_H_CO2] failed to open analytical summary: " + path);
    out << "case=" << cfg.case_name << "\n";
    out << "output_dir=" << summary.case_dir << "\n";
    out << "mesh=" << cfg.nx << "x" << cfg.ny << "\n";
    out << "cells=" << summary.n_cells << "\n";
    out << "h_char=" << std::setprecision(12) << summary.h_char << "\n";
    out << "t_end=" << summary.t_end << "\n";
    out << "steps=" << summary.steps << "\n";
    out << "pressure_diffusivity=" << ComputePressureDiffusivity(cfg) << "\n";
    out << "gravity_vector=(" << cfg.gravity_vector.m_x << "," << cfg.gravity_vector.m_y << "," << cfg.gravity_vector.m_z << ")\n";
    out << "non_orthogonal_correction=" << BoolString(cfg.enable_non_orthogonal_correction) << "\n";
    out << "final_l1_norm=" << summary.final_l1_norm << "\n";
    out << "final_l2_norm=" << summary.final_l2_norm << "\n";
    out << "final_linf_norm=" << summary.final_linf_norm << "\n";
    out << "l2_threshold=" << cfg.analytical_l2_threshold << "\n";
    out << "linf_threshold=" << cfg.analytical_linf_threshold << "\n";
    out << "grid_convergence_ok=" << BoolString(summary.grid_convergence_ok) << "\n";
    out << "time_sensitivity_ok=" << BoolString(summary.time_sensitivity_ok) << "\n";
    out << "analytical_passed=" << BoolString(summary.analytical_passed) << "\n";
    if (!summary.grid_convergence_csv_path.empty()) out << "grid_convergence_csv=" << summary.grid_convergence_csv_path << "\n";
    if (!summary.time_sensitivity_csv_path.empty()) out << "time_sensitivity_csv=" << summary.time_sensitivity_csv_path << "\n";
    out << "\n[report_metrics]\n";
    for (const auto& m : summary.report_metrics) {
        out << m.tag << " requested_fraction=" << m.requested_fraction
            << " target_time_s=" << m.target_time_s
            << " actual_time_s=" << m.actual_time_s
            << " l1_norm=" << m.l1_norm
            << " l2_norm=" << m.l2_norm
            << " linf_norm=" << m.linf_norm
            << " max_err_cell=" << m.max_err_cell << "\n";
        if (!m.compare_csv_path.empty()) out << "compare_csv=" << m.compare_csv_path << "\n";
        if (!m.profile_csv_path.empty()) out << "profile_csv=" << m.profile_csv_path << "\n";
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

TestCasePlan BuildDefaultPlan() { TestCasePlan plan; plan.plan_key = "h_co2_constpp_nofrac_nowell"; return plan; }
using BuilderFn = TestCasePlan(*)();

const std::unordered_map<std::string, BuilderFn>& GetRegistry() {
    static const std::unordered_map<std::string, BuilderFn> registry = {
        {"h_co2_constpp_nofrac_nowell", &BuildDefaultPlan}
    };
    return registry;
}

TestCaseSummary RunSingleCaseCore(const TestCaseSpec& cfg, const std::string& outputDirOverride) {
    TestCaseSummary summary;
    summary.nx = cfg.nx;
    summary.ny = cfg.ny;
    summary.case_dir = outputDirOverride.empty()
        ? (cfg.output_base_dir + "/" + cfg.sub_dir + "/" + cfg.case_name)
        : outputDirOverride;
    EnsureDirRecursive(summary.case_dir);
    summary.convergence_log_path = summary.case_dir + "/convergence.log";
    summary.metrics_csv_path = summary.case_dir + "/metrics.csv";
    if (cfg.emit_detailed_outputs && cfg.enable_analytical_validation) {
        summary.analytical_summary_path = summary.case_dir + "/analytical_summary.txt";
    }

    std::ofstream convergenceLog(summary.convergence_log_path, std::ios::out | std::ios::trunc);
    if (!convergenceLog.good()) {
        throw std::runtime_error("[Test_H_CO2] failed to open convergence log: " + summary.convergence_log_path);
    }

    MeshManager mgr(cfg.lx, cfg.ly, 0.0, cfg.nx, cfg.ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OverRelaxed);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(1);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);
    const std::size_t nCells = mgr.mesh().getCells().size();
    if (nCells == 0) throw std::runtime_error("[Test_H_CO2] matrix cell count is zero.");
    summary.n_cells = static_cast<int>(nCells);

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

    std::vector<double> pBlocksLatest(static_cast<std::size_t>(std::max(mgr.getTotalDOFCount(), 0)), cfg.p_init);
    SyncPressureFieldsToFM(mgr, fm, pBlocksLatest, cfg.t_init);
    if (cfg.export_vtk) PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/initial.vtk", 0.0);

    bool midExported = false;
    int iterSum = 0, iterCount = 0, maxIters = 0;
    double prevAcceptedTime = 0.0;
    std::vector<double> prevAcceptedPBlocks;
    std::vector<AnalyticalSnapshot> snapshots = BuildSnapshots(cfg);
    const double captureTol = std::max(1.0e-9 * std::max(cfg.target_end_time_s, 1.0), 1.0e-9);

    FIM_Engine::TransientOptionalModules<MeshManager, FieldManager_2D> modules;
    modules.pressure_bc = &bcP;
    modules.SetFluidModelConfig(FIM_Engine::UnifiedFluidModelConfig::MakePressureOnlyCO2Constant(
        cfg.t_init,
        cfg.rho_const,
        cfg.mu_const));
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
            summary.steps = step;
            summary.t_end = timeS;
            summary.total_rollbacks = totalRollbacks;
            iterSum += newtonIters;
            iterCount += 1;
            maxIters = std::max(maxIters, newtonIters);
            convergenceLog << "[Step " << step << "] t=" << std::scientific << std::setprecision(8) << timeS
                           << " dt=" << dtUsedS << " iters=" << newtonIters
                           << " residual_inf=" << residualInf
                           << " rollbacks=" << totalRollbacks
                           << " mode=" << convergeMode << "\n";
            for (auto& snapshot : snapshots) {
                if (!snapshot.captured && timeS + captureTol >= snapshot.target_time_s) {
                    const bool usePrev =
                        !prevAcceptedPBlocks.empty() &&
                        prevAcceptedTime > 0.0 &&
                        std::abs(prevAcceptedTime - snapshot.target_time_s) < std::abs(timeS - snapshot.target_time_s);
                    snapshot.captured = true;
                    snapshot.actual_time_s = usePrev ? prevAcceptedTime : timeS;
                    snapshot.p_blocks = usePrev ? prevAcceptedPBlocks : pBlocksLatest;
                }
            }
            prevAcceptedTime = timeS;
            prevAcceptedPBlocks = pBlocksLatest;
            if (cfg.export_vtk && !midExported && timeS + captureTol >= 0.5 * cfg.target_end_time_s) {
                PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/mid.vtk", timeS);
                midExported = true;
            }
        };

    const auto params = BuildSolverParams(cfg);
    FIM_Engine::RunGenericFIMTransient<1>(cfg.case_name, mgr, fm, ic, {}, params,
                                          FIM_Engine::SolverRoute::FIM, modules);

    SyncPressureFieldsToFM(mgr, fm, pBlocksLatest, cfg.t_init);
    if (cfg.export_vtk) {
        if (!midExported) PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/mid.vtk", summary.t_end);
        PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/final.vtk", summary.t_end);
    }

    summary.max_iters = maxIters;
    summary.avg_iters = (iterCount > 0) ? (static_cast<double>(iterSum) / static_cast<double>(iterCount)) : 0.0;
    summary.h_char = ComputeMeshCharLength(mgr);

    if (cfg.enable_analytical_validation) {
        for (auto& snapshot : snapshots) {
            if (!snapshot.captured && summary.t_end + captureTol >= snapshot.target_time_s) {
                snapshot.captured = true;
                snapshot.actual_time_s = summary.t_end;
                snapshot.p_blocks = pBlocksLatest;
            }
            if (!snapshot.captured || snapshot.p_blocks.empty()) {
                throw std::runtime_error("[Test_H_CO2] missing analytical snapshot for " + snapshot.tag);
            }
            const std::string compareCsvPath = summary.case_dir + "/analytical_compare_" + snapshot.tag + ".csv";
            const std::string profileCsvPath = summary.case_dir + "/profile_" + snapshot.tag + ".csv";
            AnalyticalMetrics metrics = EvaluateAnalyticalMetrics(
                mgr, snapshot.p_blocks, cfg, snapshot, compareCsvPath, cfg.emit_detailed_outputs);
            if (cfg.enable_profile_output && cfg.emit_detailed_outputs) {
                WriteProfileCSV(mgr, snapshot.p_blocks, cfg, snapshot, profileCsvPath);
                metrics.profile_csv_path = profileCsvPath;
            }
            summary.report_metrics.push_back(metrics);
        }
        if (!summary.report_metrics.empty()) {
            const AnalyticalMetrics& finalMetrics = summary.report_metrics.back();
            summary.final_l1_abs = finalMetrics.l1_abs;
            summary.final_l2_abs = finalMetrics.l2_abs;
            summary.final_linf_abs = finalMetrics.linf_abs;
            summary.final_l1_norm = finalMetrics.l1_norm;
            summary.final_l2_norm = finalMetrics.l2_norm;
            summary.final_linf_norm = finalMetrics.linf_norm;
        }
        summary.analytical_passed =
            summary.final_l2_norm <= cfg.analytical_l2_threshold &&
            summary.final_linf_norm <= cfg.analytical_linf_threshold;
    } else {
        summary.analytical_passed = true;
    }

    std::ofstream metrics(summary.metrics_csv_path, std::ios::out | std::ios::trunc);
    metrics << "case_name,nx,ny,n_cells,h_char,t_end,steps,total_rollbacks,avg_nonlinear_iters,"
               "max_nonlinear_iters,dt_init,dt_min,dt_max,property_mode,solver_route,"
               "final_l1_norm,final_l2_norm,final_linf_norm,analytical_passed\n";
    metrics << cfg.case_name << "," << cfg.nx << "," << cfg.ny << "," << nCells << ","
            << std::setprecision(12) << summary.h_char << "," << summary.t_end << ","
            << summary.steps << "," << summary.total_rollbacks << "," << summary.avg_iters << ","
            << summary.max_iters << "," << cfg.dt_init << "," << cfg.dt_min << "," << cfg.dt_max << ","
            << "ConstantBaseline,RunGenericFIMTransient<1>,"
            << summary.final_l1_norm << "," << summary.final_l2_norm << "," << summary.final_linf_norm << ","
            << BoolString(summary.analytical_passed) << "\n";

    if (cfg.emit_detailed_outputs && cfg.enable_analytical_validation) WriteAnalyticalSummary(cfg, summary);
    return summary;
}

TestCaseSummary RunCase(const TestCaseSpec& cfg) {
    TestCaseSummary summary = RunSingleCaseCore(cfg, "");
    std::vector<SweepStudyRow> gridRows, timeRows;

    if (cfg.enable_grid_convergence_study) {
        for (const auto& gridCase : cfg.grid_sweep_cases) {
            TestCaseSpec sweepCfg = cfg;
            sweepCfg.nx = gridCase.first;
            sweepCfg.ny = gridCase.second;
            sweepCfg.case_name = cfg.case_name + "_grid_" + std::to_string(sweepCfg.nx) + "x" + std::to_string(sweepCfg.ny);
            sweepCfg.enable_grid_convergence_study = false;
            sweepCfg.enable_time_sensitivity_study = false;
            sweepCfg.enable_profile_output = false;
            sweepCfg.export_vtk = false;
            sweepCfg.emit_detailed_outputs = false;
            const TestCaseSummary sweepSummary = RunSingleCaseCore(
                sweepCfg, summary.case_dir + "/studies/grid_" + std::to_string(sweepCfg.nx) + "x" + std::to_string(sweepCfg.ny));
            SweepStudyRow row;
            row.label = "grid_" + std::to_string(sweepCfg.nx) + "x" + std::to_string(sweepCfg.ny);
            row.case_dir = sweepSummary.case_dir;
            row.nx = sweepCfg.nx;
            row.ny = sweepCfg.ny;
            row.dt_init = sweepCfg.dt_init;
            row.h_char = sweepSummary.h_char;
            row.t_end = sweepSummary.t_end;
            row.steps = sweepSummary.steps;
            row.l1_norm = sweepSummary.final_l1_norm;
            row.l2_norm = sweepSummary.final_l2_norm;
            row.linf_norm = sweepSummary.final_linf_norm;
            gridRows.push_back(row);
        }
        summary.grid_convergence_ok = IsMonotonicNonIncreasingL2(gridRows);
        summary.grid_convergence_csv_path = summary.case_dir + "/grid_convergence.csv";
        WriteStudyCSV(gridRows, summary.grid_convergence_csv_path, "grid");
    }

    if (cfg.enable_time_sensitivity_study) {
        for (double dtInit : cfg.time_step_sweep) {
            TestCaseSpec sweepCfg = cfg;
            sweepCfg.dt_init = dtInit;
            sweepCfg.case_name = cfg.case_name + "_dt_" + std::to_string(static_cast<int>(std::round(dtInit))) + "s";
            sweepCfg.enable_grid_convergence_study = false;
            sweepCfg.enable_time_sensitivity_study = false;
            sweepCfg.enable_profile_output = false;
            sweepCfg.export_vtk = false;
            sweepCfg.emit_detailed_outputs = false;
            const TestCaseSummary sweepSummary = RunSingleCaseCore(
                sweepCfg, summary.case_dir + "/studies/dt_" + std::to_string(static_cast<int>(std::round(dtInit))) + "s");
            SweepStudyRow row;
            row.label = "dt_" + std::to_string(static_cast<int>(std::round(dtInit))) + "s";
            row.case_dir = sweepSummary.case_dir;
            row.nx = sweepCfg.nx;
            row.ny = sweepCfg.ny;
            row.dt_init = dtInit;
            row.h_char = sweepSummary.h_char;
            row.t_end = sweepSummary.t_end;
            row.steps = sweepSummary.steps;
            row.l1_norm = sweepSummary.final_l1_norm;
            row.l2_norm = sweepSummary.final_l2_norm;
            row.linf_norm = sweepSummary.final_linf_norm;
            timeRows.push_back(row);
        }
        summary.time_sensitivity_ok = IsMonotonicNonIncreasingL2(timeRows);
        summary.time_sensitivity_csv_path = summary.case_dir + "/time_sensitivity.csv";
        WriteStudyCSV(timeRows, summary.time_sensitivity_csv_path, "time");
    }

    summary.analytical_passed =
        summary.analytical_passed &&
        (!cfg.enable_grid_convergence_study || summary.grid_convergence_ok) &&
        (!cfg.enable_time_sensitivity_study || summary.time_sensitivity_ok);

    if (cfg.emit_detailed_outputs && cfg.enable_analytical_validation) WriteAnalyticalSummary(cfg, summary);

    if (cfg.enable_analytical_validation) {
        if (summary.final_l2_norm > cfg.analytical_l2_threshold || summary.final_linf_norm > cfg.analytical_linf_threshold) {
            std::ostringstream oss;
            oss << "[Test_H_CO2] analytical validation failed: L2=" << summary.final_l2_norm
                << " (threshold=" << cfg.analytical_l2_threshold << "), Linf=" << summary.final_linf_norm
                << " (threshold=" << cfg.analytical_linf_threshold << ")";
            throw std::runtime_error(oss.str());
        }
        if (cfg.enable_grid_convergence_study && !summary.grid_convergence_ok) {
            throw std::runtime_error("[Test_H_CO2] grid convergence check failed: normalized L2 is not monotonically decreasing.");
        }
        if (cfg.enable_time_sensitivity_study && !summary.time_sensitivity_ok) {
            throw std::runtime_error("[Test_H_CO2] time sensitivity check failed: normalized L2 is not monotonically decreasing.");
        }
    }
    return summary;
}

void ExecutePlanByKey(const std::string& key) {
    const auto& registry = GetRegistry();
    const auto it = registry.find(key);
    if (it == registry.end()) throw std::runtime_error("[Test_H_CO2] unknown registry key: " + key);
    const TestCasePlan plan = it->second();
    const TestCaseSummary summary = RunCase(plan.spec);
    std::cout << "\n============================================\n";
    std::cout << "[Test_H_CO2] run completed\n";
    std::cout << "  output_dir: " << summary.case_dir << "\n";
    std::cout << "  grid: " << summary.nx << " x " << summary.ny << " (" << summary.n_cells << " cells)\n";
    std::cout << "  steps: " << summary.steps << "  rollbacks: " << summary.total_rollbacks << "\n";
    std::cout << "  Newton iters: avg=" << std::fixed << std::setprecision(2) << summary.avg_iters
              << "  max=" << summary.max_iters << "\n";
    std::cout << "  final_time: " << std::scientific << std::setprecision(4) << summary.t_end << " s\n";
    std::cout << "  final errors: L1=" << summary.final_l1_norm
              << "  L2=" << summary.final_l2_norm
              << "  Linf=" << summary.final_linf_norm << "\n";
    std::cout << "  analytical_summary: " << summary.analytical_summary_path << "\n";
    if (!summary.grid_convergence_csv_path.empty()) {
        std::cout << "  grid_convergence_csv: " << summary.grid_convergence_csv_path
                  << "  monotonic=" << BoolString(summary.grid_convergence_ok) << "\n";
    }
    if (!summary.time_sensitivity_csv_path.empty()) {
        std::cout << "  time_sensitivity_csv: " << summary.time_sensitivity_csv_path
                  << "  monotonic=" << BoolString(summary.time_sensitivity_ok) << "\n";
    }
    std::cout << "  pass: " << BoolString(summary.analytical_passed) << "\n";
    std::cout << "============================================\n";
}

} // namespace

void RunTestCase() {
    ExecutePlanByKey("h_co2_constpp_nofrac_nowell");
}

} // namespace Test_H_CO2_ConstPP

/**
 * @file Test_Day6_DebugLadder.cpp
 * @brief Day6 debug ladder: pressure-only (N=1), constant-property single-phase CO2, no wells.
 */

#include "Test_Day6_DebugLadder.h"

#include "2D_PostProcess.h"
#include "BoundaryConditionManager.h"
#include "FIM_TransientCaseKit.hpp"
#include "FVM_Grad.h"
#include "MeshDefinitions.h"
#include "MeshManager.h"
#include "SolverContrlStrName_op.h"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#define DAY6L_MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define DAY6L_MKDIR(path) mkdir(path, 0777)
#endif

namespace Test_Day6 {

namespace {

constexpr double kPi = 3.14159265358979323846;

struct PressurePhysicalParams {
    double phi = 0.10;
    double perm = 1.0e-13; // m^2
    double ct = 5.0e-9;    // Pa^-1
    double mu = 6.0e-5;    // Pa*s (constant CO2 viscosity, simplified)
};

struct PressureCaseConfig {
    std::string level_dir;   // L1/L2/L3
    std::string case_name;   // sub-path under level dir

    PressurePhysicalParams physical;

    double lx = 400.0;
    double ly = 40.0;
    int nx = 40;
    int ny = 4;

    double p_init = 10.0e6;
    double p_left = 12.0e6;
    double p_right = 8.0e6;
    double t_init = 360.0;

    double dt_init = 200.0;
    double dt_min = 1.0;
    double dt_max = 2.0e4;
    double target_end_time_s = 2.0e5;
    int max_steps = 20000;

    bool enable_non_orth_correction = true;
    NormalVectorCorrectionMethod normal_corr_method = NormalVectorCorrectionMethod::OverRelaxed;

    int max_nonlinear_iter = 8;         // fixed-point iterations for deferred non-orth term
    double rel_update_tol = 1.0e-8;     // pressure update convergence tolerance
    double linear_res_tol = 1.0e-10;    // linear residual tolerance
    double rollback_shrink_factor = 0.5;

    bool compute_analytical_error = true;
    int analytical_terms = 200;
};

struct StepReport {
    bool converged = false;
    int nonlinear_iters = 0;
    double rel_update = std::numeric_limits<double>::infinity();
    double linear_rel_res = std::numeric_limits<double>::infinity();
};

struct CaseSummary {
    std::string case_dir;
    std::string convergence_log_path;
    std::string metrics_csv_path;
    std::string analytical_csv_path;

    double time_end = 0.0;
    int steps = 0;
    int total_rollbacks = 0;
    int max_iters_observed = 0;
    double avg_iters = 0.0;
    double l2_error_final = std::numeric_limits<double>::quiet_NaN();
    bool initial_vtk_ok = false;
    bool mid_vtk_ok = false;
    bool final_vtk_ok = false;
};

struct GridCaseResult {
    std::string tag;
    int nx = 0;
    int ny = 0;
    int n_cells = 0;
    double h_char = 0.0;
    double l2_error = std::numeric_limits<double>::quiet_NaN();
    int steps = 0;
    int rollbacks = 0;
};

void EnsureDirRecursive(const std::string& raw_path) {
    if (raw_path.empty()) return;
    std::string path = raw_path;
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
        DAY6L_MKDIR(current.c_str());
    }
}

bool FileExistsNonEmpty(const std::string& path) {
    std::ifstream ifs(path, std::ios::binary | std::ios::ate);
    return ifs.good() && ifs.tellg() > 0;
}

void VerifyMandatoryVtkOutputs(const std::string& case_dir) {
    const std::string f0 = case_dir + "/initial.vtk";
    const std::string f1 = case_dir + "/mid.vtk";
    const std::string f2 = case_dir + "/final.vtk";
    if (!FileExistsNonEmpty(f0) || !FileExistsNonEmpty(f1) || !FileExistsNonEmpty(f2)) {
        std::ostringstream oss;
        oss << "[Day6Ladder] Mandatory VTK outputs missing/non-empty check failed in " << case_dir;
        throw std::runtime_error(oss.str());
    }
}

bool CopyBinaryFile(const std::string& src, const std::string& dst) {
    std::ifstream in(src, std::ios::binary);
    if (!in.good()) return false;
    std::ofstream out(dst, std::ios::binary | std::ios::trunc);
    if (!out.good()) return false;
    out << in.rdbuf();
    out.flush();
    return out.good();
}

const char* NormalMethodName(NormalVectorCorrectionMethod method) {
    switch (method) {
    case NormalVectorCorrectionMethod::MinimumCorrection: return "MinimumCorrection";
    case NormalVectorCorrectionMethod::OrthogonalCorrection: return "OrthogonalCorrection";
    case NormalVectorCorrectionMethod::OverRelaxed: return "OverRelaxed";
    default: return "Unknown";
    }
}

double Diffusivity(const PressureCaseConfig& cfg) {
    const double denom = std::max(cfg.physical.mu * cfg.physical.phi * cfg.physical.ct, 1.0e-30);
    return cfg.physical.perm / denom;
}

double AnalyticalPressure1D(const PressureCaseConfig& cfg, double x, double t) {
    const double delta_p = cfg.p_left - cfg.p_right;
    const double a0 = cfg.p_init - cfg.p_left;
    const double d_h = Diffusivity(cfg);
    const double p_ss = cfg.p_left - delta_p * (x / cfg.lx);

    double phi_xt = 0.0;
    for (int n = 1; n <= cfg.analytical_terms; ++n) {
        const double sign_n = (n % 2 == 0) ? 1.0 : -1.0;
        const double bn = (2.0 / (static_cast<double>(n) * kPi))
            * (a0 * (1.0 - sign_n) - delta_p * sign_n);
        const double lambda_n = static_cast<double>(n) * kPi / cfg.lx;
        phi_xt += bn * std::sin(lambda_n * x) * std::exp(-d_h * lambda_n * lambda_n * t);
    }
    return p_ss + phi_xt;
}

double ComputeL2RelativeErrorFinal(const MeshManager& mgr,
                                   const std::vector<double>& p,
                                   const PressureCaseConfig& cfg,
                                   double t_final) {
    const auto& cells = mgr.mesh().getCells();
    if (cells.size() != p.size() || cells.empty()) {
        throw std::runtime_error("[Day6Ladder] pressure/cell size mismatch when computing L2 error.");
    }

    const double dp = std::max(std::abs(cfg.p_left - cfg.p_right), 1.0);
    double sum_e2v = 0.0;
    double sum_v = 0.0;
    for (size_t i = 0; i < cells.size(); ++i) {
        const double p_ana = AnalyticalPressure1D(cfg, cells[i].center.m_x, t_final);
        const double e = p[i] - p_ana;
        const double v = std::max(cells[i].volume, 0.0);
        sum_e2v += e * e * v;
        sum_v += v;
    }

    if (sum_v <= 0.0) return std::numeric_limits<double>::infinity();
    return std::sqrt(sum_e2v / sum_v) / dp;
}

double ComputeMeshCharLength(const MeshManager& mgr) {
    const auto& cells = mgr.mesh().getCells();
    if (cells.empty()) return 0.0;
    double total_v = 0.0;
    for (const auto& c : cells) total_v += std::max(c.volume, 0.0);
    if (total_v <= 0.0) return 0.0;
    return std::sqrt(total_v / static_cast<double>(cells.size()));
}

std::vector<Vector> ComputePressureGradientGG(const MeshManager& mgr,
                                              const std::vector<double>& p,
                                              const BoundarySetting::BoundaryConditionManager& bc_p) {
    volScalarField p_tmp("p_tmp", p.size(), 0.0);
    for (size_t i = 0; i < p.size(); ++i) p_tmp[i] = p[i];

    FVM_Grad grad_solver(mgr.mesh(), nullptr, nullptr, &bc_p);
    auto grad_field = grad_solver.compute(p_tmp, FVM_Grad::Method::GreenGauss);
    std::vector<Vector> grad(p.size(), Vector(0.0));
    if (!grad_field || grad_field->data.size() != p.size()) return grad;
    for (size_t i = 0; i < p.size(); ++i) grad[i] = (*grad_field)[i];
    return grad;
}

void SyncPressureToField(FieldManager_2D& fm,
                         const std::vector<double>& p,
                         double t_const) {
    const auto p_cfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
    const auto t_cfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();

    auto p_field = fm.getOrCreateMatrixScalar(p_cfg.pressure_field, 0.0);
    auto t_field = fm.getOrCreateMatrixScalar(t_cfg.temperatue_field, t_const);
    if (!p_field || p_field->data.size() != p.size()) {
        throw std::runtime_error("[Day6Ladder] failed to sync pressure field to FieldManager.");
    }
    for (size_t i = 0; i < p.size(); ++i) p_field->data[i] = p[i];

    if (t_field) {
        for (size_t i = 0; i < t_field->data.size(); ++i) t_field->data[i] = t_const;
    }
}

void ExportSnapshotVTK(const std::string& case_dir,
                       const std::string& filename,
                       double time_s,
                       MeshManager& mgr,
                       FieldManager_2D& fm,
                       const std::vector<double>& p,
                       double t_const) {
    SyncPressureToField(fm, p, t_const);
    const std::string path = case_dir + "/" + filename;
    PostProcess_2D(mgr, fm).ExportVTK(path, time_s);
}

StepReport SolveSingleImplicitStep(const MeshManager& mgr,
                                   const BoundarySetting::BoundaryConditionManager& bc_p,
                                   const PressureCaseConfig& cfg,
                                   const std::vector<double>& p_old,
                                   double dt,
                                   std::vector<double>& p_new,
                                   std::ostream& log) {
    const auto& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& faces = mesh.getFaces();
    const size_t n = cells.size();

    if (p_old.size() != n || n == 0) {
        throw std::runtime_error("[Day6Ladder] invalid pressure state size for implicit step.");
    }

    StepReport report;
    std::vector<double> p_iter = p_old;
    const double mobility = cfg.physical.perm / std::max(cfg.physical.mu, 1.0e-30);
    constexpr double kEps = 1.0e-14;

    for (int iter = 1; iter <= cfg.max_nonlinear_iter; ++iter) {
        const std::vector<Vector> grad_p = cfg.enable_non_orth_correction
            ? ComputePressureGradientGG(mgr, p_iter, bc_p)
            : std::vector<Vector>(n, Vector(0.0));

        std::vector<double> rhs(n, 0.0);
        std::vector<double> diag(n, 0.0);
        std::vector<Eigen::Triplet<double>> tri;
        tri.reserve(n * 7);

        for (size_t i = 0; i < n; ++i) {
            const double acc = cfg.physical.phi * cfg.physical.ct * std::max(cells[i].volume, 0.0) / std::max(dt, kEps);
            diag[i] += acc;
            rhs[i] += acc * p_old[i];
        }

        for (const auto& face : faces) {
            const int owner = face.ownerCell_index;
            if (owner < 0 || owner >= static_cast<int>(n)) continue;
            const Vector& c_owner = cells[owner].center;
            const double area_e = std::max(face.vectorE.Mag(), kEps);

            if (!face.isBoundary()) {
                const int nei = face.neighborCell_index;
                if (nei < 0 || nei >= static_cast<int>(n)) continue;
                const Vector& c_nei = cells[nei].center;

                const double d = std::max(std::abs((c_nei - c_owner) * face.normal), kEps);
                const double t_orth = mobility * area_e / d;

                diag[owner] += t_orth;
                diag[nei] += t_orth;
                tri.emplace_back(owner, nei, -t_orth);
                tri.emplace_back(nei, owner, -t_orth);

                if (cfg.enable_non_orth_correction) {
                    const double w = std::min(1.0, std::max(0.0, face.f_linearInterpolationCoef));
                    const Vector grad_f = grad_p[owner] * w + grad_p[nei] * (1.0 - w);
                    const double corr = mobility * (grad_f * face.vectorT);
                    rhs[owner] -= corr;
                    rhs[nei] += corr;
                }
                continue;
            }

            if (!bc_p.HasBC(face.physicalGroupId)) continue;
            const auto bc = bc_p.GetBCCoefficients(face.physicalGroupId, face.midpoint);
            const double d_b = std::max(std::abs((face.midpoint - c_owner) * face.normal), kEps);
            const double t_b = mobility * area_e / d_b;

            if (bc.type == BoundarySetting::BoundaryType::Dirichlet) {
                const double p_bc = (std::abs(bc.a) > kEps) ? (bc.c / bc.a) : cfg.p_init;
                diag[owner] += t_b;
                rhs[owner] += t_b * p_bc;
            }
            else if (bc.type == BoundarySetting::BoundaryType::Robin) {
                if (std::abs(bc.b) > kEps && std::abs(bc.a) > kEps) {
                    const double beta = bc.a / bc.b;
                    const double far_field_value = bc.c / bc.a;
                    const double coeff = std::max(0.0, beta) * std::max(face.length, 0.0);
                    diag[owner] += coeff;
                    rhs[owner] += coeff * far_field_value;
                }
            }
            else {
                const double q_out = (std::abs(bc.b) > kEps) ? (bc.c / bc.b) : 0.0;
                rhs[owner] -= q_out * std::max(face.length, 0.0);
            }
        }

        for (size_t i = 0; i < n; ++i) {
            tri.emplace_back(static_cast<int>(i), static_cast<int>(i), diag[i]);
        }

        Eigen::SparseMatrix<double> a(static_cast<int>(n), static_cast<int>(n));
        a.setFromTriplets(tri.begin(), tri.end(), [](double v0, double v1) { return v0 + v1; });
        a.makeCompressed();

        Eigen::VectorXd b(static_cast<int>(n));
        for (size_t i = 0; i < n; ++i) b[static_cast<int>(i)] = rhs[i];

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(a);
        solver.factorize(a);
        if (solver.info() != Eigen::Success) {
            report.converged = false;
            report.nonlinear_iters = iter;
            log << "    [Iter " << iter << "] SparseLU factorization failed.\n";
            return report;
        }

        Eigen::VectorXd x = solver.solve(b);
        if (solver.info() != Eigen::Success) {
            report.converged = false;
            report.nonlinear_iters = iter;
            log << "    [Iter " << iter << "] SparseLU solve failed.\n";
            return report;
        }

        std::vector<double> p_iter_next(n, 0.0);
        for (size_t i = 0; i < n; ++i) p_iter_next[i] = x[static_cast<int>(i)];

        double num = 0.0;
        double den = 0.0;
        for (size_t i = 0; i < n; ++i) {
            const double dval = p_iter_next[i] - p_iter[i];
            num += dval * dval;
            den += p_iter_next[i] * p_iter_next[i];
        }
        report.rel_update = std::sqrt(num) / std::max(std::sqrt(den), 1.0);

        const Eigen::VectorXd r = a * x - b;
        report.linear_rel_res = r.norm() / std::max(b.norm(), 1.0);
        report.nonlinear_iters = iter;

        log << "    [Iter " << std::setw(2) << iter
            << "] rel_update=" << std::scientific << std::setprecision(6) << report.rel_update
            << " lin_rel_res=" << std::scientific << std::setprecision(6) << report.linear_rel_res << "\n";

        p_iter.swap(p_iter_next);
        const bool converged = (report.rel_update <= cfg.rel_update_tol) && (report.linear_rel_res <= cfg.linear_res_tol);
        if (converged) {
            report.converged = true;
            p_new = p_iter;
            return report;
        }
    }

    report.converged = false;
    p_new = p_iter;
    return report;
}

void WriteAnalyticalCSV(const std::string& path,
                        const MeshManager& mgr,
                        const std::vector<double>& p,
                        const PressureCaseConfig& cfg,
                        double t_final) {
    const auto& cells = mgr.mesh().getCells();
    const double dp = std::max(std::abs(cfg.p_left - cfg.p_right), 1.0);

    std::ofstream csv(path, std::ios::out | std::ios::trunc);
    csv << "cell_id,x,y,p_num,p_ana,abs_err,abs_err_over_dP\n";
    for (size_t i = 0; i < cells.size(); ++i) {
        const double x = cells[i].center.m_x;
        const double y = cells[i].center.m_y;
        const double p_ana = AnalyticalPressure1D(cfg, x, t_final);
        const double err = std::abs(p[i] - p_ana);
        csv << i << ","
            << std::setprecision(12) << x << ","
            << std::setprecision(12) << y << ","
            << std::setprecision(12) << p[i] << ","
            << std::setprecision(12) << p_ana << ","
            << std::setprecision(12) << err << ","
            << std::setprecision(12) << (err / dp) << "\n";
    }
}

CaseSummary RunPressureCase(const PressureCaseConfig& cfg) {
    CaseSummary summary;
    summary.case_dir = "Test/Transient/Day6Ladder/" + cfg.level_dir + "/" + cfg.case_name;
    EnsureDirRecursive(summary.case_dir);

    summary.convergence_log_path = summary.case_dir + "/convergence.log";
    summary.metrics_csv_path = summary.case_dir + "/metrics.csv";
    summary.analytical_csv_path = summary.case_dir + "/analytical_compare.csv";

    std::ofstream log(summary.convergence_log_path, std::ios::out | std::ios::trunc);
    if (!log.good()) {
        throw std::runtime_error("[Day6Ladder] failed to open convergence.log: " + summary.convergence_log_path);
    }

    log << "=== Day6 Debug Ladder Pressure-Only Case ===\n";
    log << "case=" << cfg.case_name << "\n";
    log << "level=" << cfg.level_dir << "\n";
    log << "grid=(" << cfg.nx << "," << cfg.ny << ")\n";
    log << "non_orth_enabled=" << (cfg.enable_non_orth_correction ? "true" : "false") << "\n";
    log << "normal_correction_method=" << NormalMethodName(cfg.normal_corr_method) << "\n";
    log << "dt_init=" << cfg.dt_init << " dt_min=" << cfg.dt_min << " dt_max=" << cfg.dt_max << "\n";
    log << "target_end_time_s=" << cfg.target_end_time_s << "\n";
    log << "max_nonlinear_iter=" << cfg.max_nonlinear_iter << " rel_update_tol=" << cfg.rel_update_tol
        << " linear_res_tol=" << cfg.linear_res_tol << "\n";
    log << "Physical(const CO2): phi=" << cfg.physical.phi
        << " k=" << cfg.physical.perm
        << " ct=" << cfg.physical.ct
        << " mu=" << cfg.physical.mu << "\n";

    MeshManager mgr(cfg.lx, cfg.ly, 0.0, cfg.nx, cfg.ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(cfg.normal_corr_method);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(1); // pure pressure route

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    const size_t n_cells = mgr.mesh().getCells().size();
    if (n_cells == 0) {
        throw std::runtime_error("[Day6Ladder] zero matrix cell count.");
    }

    BoundarySetting::BoundaryConditionManager bc_p;
    bc_p.Clear();
    bc_p.SetDirichletBC(MeshTags::LEFT, cfg.p_left);
    bc_p.SetDirichletBC(MeshTags::RIGHT, cfg.p_right);
    bc_p.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bc_p.SetNeumannBC(MeshTags::TOP, 0.0);

    std::vector<double> p(n_cells, cfg.p_init);

    ExportSnapshotVTK(summary.case_dir, "initial.vtk", 0.0, mgr, fm, p, cfg.t_init);

    double t = 0.0;
    double dt = cfg.dt_init;
    int step = 0;
    bool mid_exported = false;
    int total_rollbacks = 0;
    int iters_sum = 0;
    int iters_count = 0;
    int max_iters = 0;

    while (step < cfg.max_steps && t < cfg.target_end_time_s - 1.0e-12) {
        dt = std::max(cfg.dt_min, std::min(dt, cfg.dt_max));
        dt = std::min(dt, cfg.target_end_time_s - t);

        const std::vector<double> p_old = p;
        double dt_try = dt;
        int rollback_this_step = 0;
        StepReport step_report;
        std::vector<double> p_candidate;

        while (true) {
            log << "\n[Step " << (step + 1) << "] t=" << std::scientific << std::setprecision(6) << t
                << " dt_try=" << dt_try << " rollback=" << rollback_this_step << "\n";

            step_report = SolveSingleImplicitStep(mgr, bc_p, cfg, p_old, dt_try, p_candidate, log);
            if (step_report.converged) break;

            rollback_this_step += 1;
            total_rollbacks += 1;
            dt_try *= std::max(0.1, std::min(0.95, cfg.rollback_shrink_factor));
            log << "  [Rollback] nonlinear convergence failed; shrink dt -> " << dt_try << "\n";

            if (dt_try < cfg.dt_min - 1.0e-15) {
                throw std::runtime_error("[Day6Ladder] dt dropped below dt_min during rollback.");
            }
        }

        p = p_candidate;
        t += dt_try;
        ++step;
        summary.steps = step;
        summary.time_end = t;
        summary.total_rollbacks = total_rollbacks;

        iters_sum += step_report.nonlinear_iters;
        iters_count += 1;
        max_iters = std::max(max_iters, step_report.nonlinear_iters);

        log << "  [Step Accepted] t_end=" << std::scientific << std::setprecision(6) << t
            << " nonlinear_iters=" << step_report.nonlinear_iters
            << " rel_update=" << step_report.rel_update
            << " lin_res=" << step_report.linear_rel_res
            << " rollbacks=" << rollback_this_step << "\n";

        if (!mid_exported && t >= 0.5 * cfg.target_end_time_s - 1.0e-12) {
            ExportSnapshotVTK(summary.case_dir, "mid.vtk", t, mgr, fm, p, cfg.t_init);
            mid_exported = true;
        }

        if (step_report.nonlinear_iters <= 2) {
            dt = std::min(cfg.dt_max, dt_try * 1.20);
        }
        else if (step_report.nonlinear_iters <= 4) {
            dt = std::min(cfg.dt_max, dt_try * 1.05);
        }
        else {
            dt = std::max(cfg.dt_min, dt_try * 0.75);
        }
    }

    if (!mid_exported) {
        ExportSnapshotVTK(summary.case_dir, "mid.vtk", t, mgr, fm, p, cfg.t_init);
    }
    ExportSnapshotVTK(summary.case_dir, "final.vtk", t, mgr, fm, p, cfg.t_init);

    VerifyMandatoryVtkOutputs(summary.case_dir);
    summary.initial_vtk_ok = FileExistsNonEmpty(summary.case_dir + "/initial.vtk");
    summary.mid_vtk_ok = FileExistsNonEmpty(summary.case_dir + "/mid.vtk");
    summary.final_vtk_ok = FileExistsNonEmpty(summary.case_dir + "/final.vtk");

    summary.max_iters_observed = max_iters;
    summary.avg_iters = (iters_count > 0) ? (static_cast<double>(iters_sum) / static_cast<double>(iters_count)) : 0.0;

    if (cfg.compute_analytical_error) {
        summary.l2_error_final = ComputeL2RelativeErrorFinal(mgr, p, cfg, t);
        WriteAnalyticalCSV(summary.analytical_csv_path, mgr, p, cfg, t);
    }

    std::ofstream metrics(summary.metrics_csv_path, std::ios::out | std::ios::trunc);
    metrics << "case_name,level,nx,ny,n_cells,h_char,t_end,steps,total_rollbacks,avg_nonlinear_iters,max_nonlinear_iters,dt_init,dt_min,dt_max,enable_non_orth,normal_method,l2_error_final\n";
    metrics << cfg.case_name << ","
            << cfg.level_dir << ","
            << cfg.nx << ","
            << cfg.ny << ","
            << n_cells << ","
            << std::setprecision(12) << ComputeMeshCharLength(mgr) << ","
            << std::setprecision(12) << summary.time_end << ","
            << summary.steps << ","
            << summary.total_rollbacks << ","
            << std::setprecision(8) << summary.avg_iters << ","
            << summary.max_iters_observed << ","
            << std::setprecision(12) << cfg.dt_init << ","
            << std::setprecision(12) << cfg.dt_min << ","
            << std::setprecision(12) << cfg.dt_max << ","
            << (cfg.enable_non_orth_correction ? 1 : 0) << ","
            << NormalMethodName(cfg.normal_corr_method) << ","
            << std::setprecision(12) << summary.l2_error_final << "\n";

    log << "\n[Case PASS] " << cfg.case_name
        << " steps=" << summary.steps
        << " rollbacks=" << summary.total_rollbacks
        << " avg_iters=" << summary.avg_iters
        << " max_iters=" << summary.max_iters_observed
        << " l2_error_final=" << summary.l2_error_final << "\n";

    return summary;
}

PressureCaseConfig MakeBaseConfig(const std::string& level_dir, const std::string& case_name, int nx, int ny) {
    PressureCaseConfig cfg;
    cfg.level_dir = level_dir;
    cfg.case_name = case_name;
    cfg.nx = nx;
    cfg.ny = ny;

    cfg.enable_non_orth_correction = true;
    cfg.normal_corr_method = NormalVectorCorrectionMethod::OverRelaxed;
    cfg.compute_analytical_error = true;
    cfg.analytical_terms = 200;
    return cfg;
}

void CopyCaseFramesToRoot(const std::string& src_case_dir, const std::string& dst_case_dir) {
    EnsureDirRecursive(dst_case_dir);
    const std::vector<std::string> names = { "initial.vtk", "mid.vtk", "final.vtk" };
    for (const auto& n : names) {
        const std::string src = src_case_dir + "/" + n;
        const std::string dst = dst_case_dir + "/" + n;
        if (!CopyBinaryFile(src, dst)) {
            throw std::runtime_error("[Day6Ladder] failed to copy VTK frame from " + src + " to " + dst);
        }
    }
    VerifyMandatoryVtkOutputs(dst_case_dir);
}

} // namespace

void Run_Day6L1_2D_SP_CO2Const_NoWell_Analytical() {
    std::cout << "\n=== Run_Day6L1_2D_SP_CO2Const_NoWell_Analytical ===\n";

    PressureCaseConfig cfg = MakeBaseConfig(
        "L1",
        "day6l1_2d_sp_co2_const_nowell_analytical",
        48, 6);

    cfg.dt_init = 120.0;
    cfg.dt_min = 1.0;
    cfg.dt_max = 5000.0;
    cfg.target_end_time_s = 1.0e5;
    cfg.max_steps = 12000;
    cfg.rel_update_tol = 5.0e-9;
    cfg.linear_res_tol = 1.0e-10;
    cfg.max_nonlinear_iter = 8;

    const CaseSummary summary = RunPressureCase(cfg);
    constexpr double kL1ErrorThreshold = 5.0e-2;

    std::cout << "[L1] case_dir = " << summary.case_dir << "\n";
    std::cout << "[L1] l2_error_final = " << summary.l2_error_final
              << " (threshold=" << kL1ErrorThreshold << ")\n";

    if (!(summary.l2_error_final <= kL1ErrorThreshold)) {
        std::ostringstream oss;
        oss << "[L1] analytical validation failed: l2_error_final=" << summary.l2_error_final
            << " > " << kL1ErrorThreshold;
        throw std::runtime_error(oss.str());
    }
}

void Run_Day6L2_2D_SP_CO2Const_NoWell_GridConvergence() {
    std::cout << "\n=== Run_Day6L2_2D_SP_CO2Const_NoWell_GridConvergence ===\n";

    const std::string root_case_name = "day6l2_2d_sp_co2_const_nowell_grid";
    const std::string root_case_dir = "Test/Transient/Day6Ladder/L2/" + root_case_name;
    EnsureDirRecursive(root_case_dir);

    struct GridSpec { int nx; int ny; double dt_init; };
    const std::vector<GridSpec> grids = {
        { 24, 3, 200.0 },
        { 48, 6, 120.0 },
        { 96, 12, 60.0 }
    };

    std::vector<GridCaseResult> results;
    results.reserve(grids.size());
    std::string finest_case_dir;

    std::ofstream root_log(root_case_dir + "/convergence.log", std::ios::out | std::ios::trunc);
    root_log << "=== L2 Grid Convergence ===\n";

    for (size_t i = 0; i < grids.size(); ++i) {
        const auto& g = grids[i];
        std::ostringstream tag;
        tag << "nx" << g.nx << "_ny" << g.ny;

        PressureCaseConfig cfg = MakeBaseConfig(
            "L2",
            root_case_name + "/" + tag.str(),
            g.nx, g.ny);

        cfg.dt_init = g.dt_init;
        cfg.dt_min = 1.0;
        cfg.dt_max = 4000.0;
        cfg.target_end_time_s = 8.0e4;
        cfg.max_steps = 20000;
        cfg.rel_update_tol = 5.0e-9;
        cfg.linear_res_tol = 1.0e-10;
        cfg.max_nonlinear_iter = 10;

        const CaseSummary summary = RunPressureCase(cfg);
        finest_case_dir = summary.case_dir;

        MeshManager mesh_probe(cfg.lx, cfg.ly, 0.0, cfg.nx, cfg.ny, 0, true, false);
        mesh_probe.BuildSolidMatrixGrid_2D(cfg.normal_corr_method);
        const int n_cells = static_cast<int>(mesh_probe.mesh().getCells().size());
        const double h_char = ComputeMeshCharLength(mesh_probe);

        GridCaseResult r;
        r.tag = tag.str();
        r.nx = cfg.nx;
        r.ny = cfg.ny;
        r.n_cells = n_cells;
        r.h_char = h_char;
        r.l2_error = summary.l2_error_final;
        r.steps = summary.steps;
        r.rollbacks = summary.total_rollbacks;
        results.push_back(r);

        root_log << "[" << r.tag << "] h=" << r.h_char
                 << " l2_error=" << r.l2_error
                 << " steps=" << r.steps
                 << " rollbacks=" << r.rollbacks << "\n";
    }

    std::ofstream metrics(root_case_dir + "/metrics.csv", std::ios::out | std::ios::trunc);
    metrics << "grid_tag,nx,ny,n_cells,h_char,l2_error_final,order_vs_prev,steps,total_rollbacks\n";
    for (size_t i = 0; i < results.size(); ++i) {
        double order = std::numeric_limits<double>::quiet_NaN();
        if (i > 0 && results[i].l2_error > 0.0 && results[i - 1].l2_error > 0.0
            && results[i].h_char > 0.0 && results[i - 1].h_char > 0.0) {
            order = std::log(results[i - 1].l2_error / results[i].l2_error)
                / std::log(results[i - 1].h_char / results[i].h_char);
        }

        metrics << results[i].tag << ","
                << results[i].nx << ","
                << results[i].ny << ","
                << results[i].n_cells << ","
                << std::setprecision(12) << results[i].h_char << ","
                << std::setprecision(12) << results[i].l2_error << ","
                << std::setprecision(12) << order << ","
                << results[i].steps << ","
                << results[i].rollbacks << "\n";

        root_log << "[order] " << results[i].tag << " order_vs_prev=" << order << "\n";
    }

    if (results.size() >= 2) {
        for (size_t i = 1; i < results.size(); ++i) {
            if (!(results[i].l2_error < results[i - 1].l2_error)) {
                std::ostringstream oss;
                oss << "[L2] error is not decreasing with refinement at " << results[i].tag
                    << ": prev=" << results[i - 1].l2_error
                    << " curr=" << results[i].l2_error;
                throw std::runtime_error(oss.str());
            }
        }
    }

    if (!finest_case_dir.empty()) {
        CopyCaseFramesToRoot(finest_case_dir, root_case_dir);
    }

    std::cout << "[L2] root_case_dir = " << root_case_dir << "\n";
    std::cout << "[L2] grid-convergence metrics written: " << (root_case_dir + "/metrics.csv") << "\n";
}

void Run_Day6L3_2D_SP_CO2Const_NoWell_SolverRobustness() {
    std::cout << "\n=== Run_Day6L3_2D_SP_CO2Const_NoWell_SolverRobustness ===\n";

    PressureCaseConfig cfg = MakeBaseConfig(
        "L3",
        "day6l3_2d_sp_co2_const_nowell_solver",
        48, 6);

    cfg.dt_init = 2.0e4;     // stress large dt
    cfg.dt_min = 50.0;
    cfg.dt_max = 3.0e4;
    cfg.target_end_time_s = 1.2e5;
    cfg.max_steps = 5000;
    cfg.max_nonlinear_iter = 3;  // intentionally tight to trigger rollback behaviour
    cfg.rel_update_tol = 1.0e-10;
    cfg.linear_res_tol = 1.0e-11;
    cfg.rollback_shrink_factor = 0.5;
    cfg.compute_analytical_error = false;

    const CaseSummary summary = RunPressureCase(cfg);

    std::cout << "[L3] case_dir = " << summary.case_dir << "\n";
    std::cout << "[L3] steps=" << summary.steps
              << " total_rollbacks=" << summary.total_rollbacks
              << " avg_iters=" << summary.avg_iters
              << " max_iters=" << summary.max_iters_observed << "\n";
}

void Run_Day6Ladder_2D_SP_CO2Const_NoWell_All() {
    std::cout << "\n=== Run_Day6Ladder_2D_SP_CO2Const_NoWell_All ===\n";
    Run_Day6L1_2D_SP_CO2Const_NoWell_Analytical();
    Run_Day6L2_2D_SP_CO2Const_NoWell_GridConvergence();
    Run_Day6L3_2D_SP_CO2Const_NoWell_SolverRobustness();
    std::cout << "[Day6Ladder] ALL PASS\n";
}

} // namespace Test_Day6

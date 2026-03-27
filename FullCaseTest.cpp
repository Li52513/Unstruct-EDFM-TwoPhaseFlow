#include "FullCaseTest.h"

#include "2D_PostProcess.h"
#include "BoundaryConditionManager.h"
#include "FIM_TransientCaseKit.hpp"
#include "MeshDefinitions.h"
#include "MeshManager.h"
#include "SolverContrlStrName_op.h"
#include "Test_Day6_TransientSolver.h"
#include "Well_WellControlTypes.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#define FULLCASE_MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define FULLCASE_MKDIR(path) mkdir(path, 0777)
#endif

namespace FullCaseTest {
namespace {

constexpr char kAuditSolverRoute[] = "RunGenericFIMTransient<1>";
constexpr char kAuditGradRoute[] = "FVM_Grad+BC(via RunGeneric_impl non-orth path)";
constexpr char kAuditBcGeomRoute[] = "BoundaryAssembler(dist=projection,area=vectorE.Mag)";
constexpr bool kAuditLocalSolverImpl = false;
constexpr bool kAuditLocalFluxImpl = false;
constexpr bool kAuditLocalAssemblyImpl = false;

enum class N1PlanPolicy {
    SingleRun = 0,
    GridConvergence = 1
};

struct N1CaseSpec {
    std::string case_name;
    std::string level_dir = "TemplateN1";
    TopologyVariant topology = TopologyVariant::NoFrac;

    double lx = 400.0;
    double ly = 40.0;
    int nx = 48;
    int ny = 6;

    // Fracture geometry (single-fracture default is vertical barrier).
    double frac_x_ratio = 0.50;
    double frac_y0_ratio = 0.10;
    double frac_y1_ratio = 0.90;

    // Matrix properties
    double matrix_phi = 0.10;
    double matrix_perm = 1.0e-13;
    double matrix_ct = 5.0e-9;

    // Fracture properties
    double fracture_phi = 0.10;
    double fracture_kt = 1.0e-13;
    double fracture_kn = 1.0e-13;
    double fracture_ct = 5.0e-9;

    // Fluid / initial-boundary
    double p_init = 10.0e6;
    double p_left = 12.0e6;
    double p_right = 8.0e6;
    double t_init = 360.0;
    double mu_const = 6.0e-5;
    double rho_const = 700.0;

    // Time controls
    double dt_init = 120.0;
    double dt_min = 1.0;
    double dt_max = 5000.0;
    double target_end_time_s = 1.0e5;
    int max_steps = 12000;
    int max_newton_iter = 12;

    // Property mode
    bool use_co2_eos = false;

    // Analytical validation
    bool analytical_check = false;
    int analytical_terms = 200;
    double analytical_l2_threshold = 5.0e-2;

    // Barrier probe control
    bool enable_barrier_probe = false;
    double probe_y_min_ratio = 0.15;
    double probe_y_max_ratio = 0.85;
};

struct N1CaseSummary {
    std::string case_dir;
    std::string convergence_log_path;
    std::string metrics_csv_path;

    int nx = 0;
    int ny = 0;
    int n_cells = 0;
    double h_char = 0.0;

    int steps = 0;
    int total_rollbacks = 0;
    double avg_iters = 0.0;
    int max_iters = 0;
    double t_end = 0.0;

    double l2_error_final = std::numeric_limits<double>::quiet_NaN();
    double delta_p_cross_fracture = std::numeric_limits<double>::quiet_NaN();
    double delta_p_cross_norm = std::numeric_limits<double>::quiet_NaN();
    int probe_left_count = 0;
    int probe_right_count = 0;
};

struct N1CasePlan {
    std::string plan_key;
    std::vector<N1CaseSpec> cases;
    N1PlanPolicy policy = N1PlanPolicy::SingleRun;
    std::string aggregate_dir;
    std::string aggregate_csv_name;
};

using BuilderFn = N1CasePlan(*)();

constexpr std::array<TopologyVariant, 3> kTopologyOrder = {
    TopologyVariant::NoFrac,
    TopologyVariant::SingleFrac,
    TopologyVariant::CrossFrac
};

void ValidateScenarioIdOrThrow(int scenarioId) {
    if (scenarioId < 1 || scenarioId > 8) {
        throw std::runtime_error("[FullCaseTest] scenarioId must be in [1, 8], got " + std::to_string(scenarioId));
    }
}

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
        FULLCASE_MKDIR(current.c_str());
    }
}

bool FileExistsNonEmpty(const std::string& path) {
    std::ifstream ifs(path, std::ios::binary | std::ios::ate);
    return ifs.good() && ifs.tellg() > 0;
}

void VerifyMandatoryVtkOutputs(const std::string& caseDir) {
    const std::string initialVtk = caseDir + "/initial.vtk";
    const std::string midVtk = caseDir + "/mid.vtk";
    const std::string finalVtk = caseDir + "/final.vtk";
    if (!FileExistsNonEmpty(initialVtk) || !FileExistsNonEmpty(midVtk) || !FileExistsNonEmpty(finalVtk)) {
        throw std::runtime_error("[FullCaseTest-N1] mandatory vtk outputs missing/non-empty check failed: " + caseDir);
    }
}

void ApplyUniformScalarField(const std::shared_ptr<volScalarField>& field, double value) {
    if (!field) return;
    for (double& v : field->data) v = value;
}

void WriteAuditLog(std::ostream& os) {
    os << "[AUDIT] solver_route=" << kAuditSolverRoute << "\n";
    os << "[AUDIT] grad_route=" << kAuditGradRoute << "\n";
    os << "[AUDIT] bc_geom_route=" << kAuditBcGeomRoute << "\n";
    os << "[AUDIT] local_solver_impl=" << (kAuditLocalSolverImpl ? 1 : 0) << "\n";
    os << "[AUDIT] local_flux_impl=" << (kAuditLocalFluxImpl ? 1 : 0) << "\n";
    os << "[AUDIT] local_assembly_impl=" << (kAuditLocalAssemblyImpl ? 1 : 0) << "\n";
}

void AssertNoLocalSolverPathStatic() {
    static bool checked = false;
    if (checked) return;
    checked = true;

    std::ifstream ifs(__FILE__, std::ios::in);
    if (!ifs.good()) {
        throw std::runtime_error("[FullCaseTest-N1] static purity scan failed: cannot open source file.");
    }

    const std::string src((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    const std::array<std::string, 5> forbidden = {
        std::string("SolveSingleImplicitStep") + "(",
        std::string("AssemblePressureMatrixLocal") + "(",
        std::string("Eigen::SparseLU") + "<",
        std::string("Eigen::Triplet") + "<double>",
        std::string("ManualNewtonPressureStep") + "("
    };
    for (const std::string& token : forbidden) {
        if (src.find(token) != std::string::npos) {
            throw std::runtime_error(std::string("[FullCaseTest-N1] static purity scan failed: forbidden token: ") + token);
        }
    }
}

void SyncN1SnapshotPressureFields(MeshManager& mgr,
                                  FieldManager_2D& fm,
                                  const std::vector<double>& pBlocks,
                                  double tConst) {
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
    const int nTotal = mgr.getTotalDOFCount();
    const int nUse = std::min(static_cast<int>(pBlocks.size()), nTotal);

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

void AddTemplateFractures(MeshManager& mgr, const N1CaseSpec& cfg) {
    if (cfg.topology == TopologyVariant::NoFrac) return;

    const Vector primaryStart(cfg.frac_x_ratio * cfg.lx, cfg.frac_y0_ratio * cfg.ly, 0.0);
    const Vector primaryEnd(cfg.frac_x_ratio * cfg.lx, cfg.frac_y1_ratio * cfg.ly, 0.0);
    mgr.addFracture(primaryStart, primaryEnd);

    if (cfg.topology == TopologyVariant::CrossFrac) {
        const double yCross = 0.5 * (cfg.frac_y0_ratio + cfg.frac_y1_ratio) * cfg.ly;
        mgr.addFracture(Vector(0.20 * cfg.lx, yCross, 0.0), Vector(0.80 * cfg.lx, yCross, 0.0));
    }

    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
}

double Diffusivity(const N1CaseSpec& cfg) {
    const double denom = std::max(cfg.mu_const * cfg.matrix_phi * cfg.matrix_ct, 1.0e-30);
    return cfg.matrix_perm / denom;
}

double AnalyticalPressure1D(const N1CaseSpec& cfg, double x, double t) {
    constexpr double kPi = 3.14159265358979323846;
    const double deltaP = cfg.p_left - cfg.p_right;
    const double a0 = cfg.p_init - cfg.p_left;
    const double dH = Diffusivity(cfg);
    const double pSS = cfg.p_left - deltaP * (x / cfg.lx);

    double phiXT = 0.0;
    for (int n = 1; n <= cfg.analytical_terms; ++n) {
        const double signN = (n % 2 == 0) ? 1.0 : -1.0;
        const double bn = (2.0 / (static_cast<double>(n) * kPi))
            * (a0 * (1.0 - signN) - deltaP * signN);
        const double lambdaN = static_cast<double>(n) * kPi / cfg.lx;
        phiXT += bn * std::sin(lambdaN * x) * std::exp(-dH * lambdaN * lambdaN * t);
    }

    return pSS + phiXT;
}

double ComputeL2RelativeErrorFinal(const MeshManager& mgr,
                                   const std::vector<double>& p,
                                   const N1CaseSpec& cfg,
                                   double tFinal) {
    const auto& cells = mgr.mesh().getCells();
    if (cells.empty() || p.size() != cells.size()) {
        throw std::runtime_error("[FullCaseTest-N1] invalid pressure/cell size for analytical L2 error.");
    }

    const double dp = std::max(std::abs(cfg.p_left - cfg.p_right), 1.0);
    double sumE2V = 0.0;
    double sumV = 0.0;
    for (size_t i = 0; i < cells.size(); ++i) {
        const double pana = AnalyticalPressure1D(cfg, cells[i].center.m_x, tFinal);
        const double e = p[i] - pana;
        const double v = std::max(cells[i].volume, 0.0);
        sumE2V += e * e * v;
        sumV += v;
    }

    if (sumV <= 0.0) return std::numeric_limits<double>::infinity();
    return std::sqrt(sumE2V / sumV) / dp;
}

double ComputeMeshCharLength(const MeshManager& mgr) {
    const auto& cells = mgr.mesh().getCells();
    if (cells.empty()) return 0.0;
    double totalV = 0.0;
    for (const auto& c : cells) totalV += std::max(c.volume, 0.0);
    if (totalV <= 0.0) return 0.0;
    return std::sqrt(totalV / static_cast<double>(cells.size()));
}

struct BarrierProbeResult {
    bool valid = false;
    double left_avg = std::numeric_limits<double>::quiet_NaN();
    double right_avg = std::numeric_limits<double>::quiet_NaN();
    double delta = std::numeric_limits<double>::quiet_NaN();
    double delta_norm = std::numeric_limits<double>::quiet_NaN();
    int left_count = 0;
    int right_count = 0;
};

BarrierProbeResult ComputeBarrierProbe(const MeshManager& mgr, const std::vector<double>& pMatrix, const N1CaseSpec& cfg) {
    BarrierProbeResult out;
    if (!cfg.enable_barrier_probe) return out;

    const auto& cells = mgr.mesh().getCells();
    if (cells.size() != pMatrix.size() || cells.empty()) return out;

    const double h = std::max(ComputeMeshCharLength(mgr), 1.0e-12);
    const double fx = cfg.frac_x_ratio * cfg.lx;
    const double yMin = cfg.probe_y_min_ratio * cfg.ly;
    const double yMax = cfg.probe_y_max_ratio * cfg.ly;
    const double probeOffset = std::max(2.0 * h, 0.01 * cfg.lx);
    const double probeBand = std::max(1.2 * h, 0.005 * cfg.lx);
    const double leftTargetX = fx - probeOffset;
    const double rightTargetX = fx + probeOffset;

    double leftSum = 0.0;
    double rightSum = 0.0;
    int leftCnt = 0;
    int rightCnt = 0;

    for (size_t i = 0; i < cells.size(); ++i) {
        const Vector& c = cells[i].center;
        if (c.m_y < yMin || c.m_y > yMax) continue;
        if (std::abs(c.m_x - leftTargetX) <= probeBand) {
            leftSum += pMatrix[i];
            ++leftCnt;
        }
        if (std::abs(c.m_x - rightTargetX) <= probeBand) {
            rightSum += pMatrix[i];
            ++rightCnt;
        }
    }

    if (leftCnt <= 0 || rightCnt <= 0) return out;

    out.valid = true;
    out.left_avg = leftSum / static_cast<double>(leftCnt);
    out.right_avg = rightSum / static_cast<double>(rightCnt);
    out.delta = out.left_avg - out.right_avg;
    out.delta_norm = out.delta / std::max(std::abs(cfg.p_left - cfg.p_right), 1.0);
    out.left_count = leftCnt;
    out.right_count = rightCnt;
    return out;
}

FIM_Engine::TransientSolverParams BuildN1Params(const N1CaseSpec& cfg) {
    FIM_Engine::TransientSolverParams p;
    p.max_steps = cfg.max_steps;
    p.dt_init = cfg.dt_init;
    p.dt_min = cfg.dt_min;
    p.dt_max = cfg.dt_max;
    p.target_end_time_s = cfg.target_end_time_s;

    p.max_newton_iter = cfg.max_newton_iter;
    p.abs_res_tol = 1.0e-10;
    p.rel_res_tol = 1.0e-6;
    p.rel_update_tol = 1.0e-8;

    p.enable_non_orthogonal_correction = true;
    p.lin_solver = FIM_Engine::LinearSolverType::AMGCL;
    p.amgcl_use_fallback_sparselu = true;
    p.enable_row_scaling = true;
    p.enable_armijo_line_search = false;
    p.diag_level = FIM_Engine::DiagLevel::Off;

    p.max_dP = 2.0e7;
    p.max_dT = 1.0;
    p.max_dSw = 0.1;
    p.min_alpha = 1.0e-8;
    return p;
}

N1CaseSummary RunN1Common(const N1CaseSpec& cfg) {
    AssertNoLocalSolverPathStatic();

    N1CaseSummary summary;
    summary.case_dir = "Test/Transient/FullCaseTest/N1/" + cfg.level_dir + "/" + cfg.case_name;
    EnsureDirRecursive(summary.case_dir);
    summary.convergence_log_path = summary.case_dir + "/convergence.log";
    summary.metrics_csv_path = summary.case_dir + "/metrics.csv";
    summary.nx = cfg.nx;
    summary.ny = cfg.ny;

    std::ofstream log(summary.convergence_log_path, std::ios::out | std::ios::trunc);
    if (!log.good()) {
        throw std::runtime_error("[FullCaseTest-N1] failed to open convergence.log: " + summary.convergence_log_path);
    }
    WriteAuditLog(log);

    MeshManager mgr(cfg.lx, cfg.ly, 0.0, cfg.nx, cfg.ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OverRelaxed);
    AddTemplateFractures(mgr, cfg);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(1);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    const size_t nCells = mgr.mesh().getCells().size();
    const int totalBlocks = mgr.getTotalDOFCount();
    if (nCells == 0) {
        throw std::runtime_error("[FullCaseTest-N1] zero matrix cell count.");
    }
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

    const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
    std::vector<double> pBlocksLatest(static_cast<std::size_t>(std::max(totalBlocks, 0)), cfg.p_init);
    SyncN1SnapshotPressureFields(mgr, fm, pBlocksLatest, cfg.t_init);
    PostProcess_2D(mgr, fm).ExportVTK(summary.case_dir + "/initial.vtk", 0.0);

    bool midExported = false;
    int iterSum = 0;
    int iterCount = 0;
    int maxIters = 0;

    FIM_Engine::TransientOptionalModules<MeshManager, FieldManager_2D> modules;
    modules.pressure_bc = &bcP;
    if (cfg.use_co2_eos) {
        modules.SetFluidModelConfig(FIM_Engine::UnifiedFluidModelConfig::MakePressureOnlyCO2EOS(cfg.t_init));
    }
    else {
        modules.SetFluidModelConfig(FIM_Engine::UnifiedFluidModelConfig::MakePressureOnlyCO2Constant(
            cfg.t_init,
            cfg.rho_const,
            cfg.mu_const));
    }
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

        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.k_t_tag, cfg.fracture_kt), cfg.fracture_kt);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.k_n_tag, cfg.fracture_kn), cfg.fracture_kn);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.phi_tag, cfg.fracture_phi), cfg.fracture_phi);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.c_r_tag, cfg.fracture_ct), cfg.fracture_ct);

        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(water.k_tag, 0.6), 0.6);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(water.k_tag, 0.6), 0.6);
    };

    modules.on_step_accepted =
        [&](int step, double timeS, double dtUsedS, int newtonIters, double residualInf, int totalRollbacks,
            const std::string& convergeMode, const std::vector<double>& pVec, const std::vector<double>&, const std::vector<double>*) {
            if (step <= 0) return;

            if (!pVec.empty()) {
                pBlocksLatest = pVec;
                SyncN1SnapshotPressureFields(mgr, fm, pBlocksLatest, cfg.t_init);
            }

            summary.steps = step;
            summary.t_end = timeS;
            summary.total_rollbacks = totalRollbacks;

            iterSum += newtonIters;
            iterCount += 1;
            maxIters = std::max(maxIters, newtonIters);

            log << "[Step " << step << "] t=" << std::scientific << std::setprecision(8) << timeS
                << " dt=" << dtUsedS
                << " iters=" << newtonIters
                << " residual_inf=" << residualInf
                << " rollbacks=" << totalRollbacks
                << " mode=" << convergeMode << "\n";

            if (!midExported && timeS >= 0.5 * cfg.target_end_time_s - 1.0e-12) {
                PostProcess_2D(mgr, fm).ExportVTK(summary.case_dir + "/mid.vtk", timeS);
                midExported = true;
            }
        };

    const auto params = BuildN1Params(cfg);
    FIM_Engine::RunGenericFIMTransient<1>(
        cfg.case_name,
        mgr,
        fm,
        ic,
        {},
        params,
        FIM_Engine::SolverRoute::FIM,
        modules);

    SyncN1SnapshotPressureFields(mgr, fm, pBlocksLatest, cfg.t_init);
    if (!midExported) {
        PostProcess_2D(mgr, fm).ExportVTK(summary.case_dir + "/mid.vtk", summary.t_end);
    }
    PostProcess_2D(mgr, fm).ExportVTK(summary.case_dir + "/final.vtk", summary.t_end);
    VerifyMandatoryVtkOutputs(summary.case_dir);

    summary.max_iters = maxIters;
    summary.avg_iters = (iterCount > 0) ? (static_cast<double>(iterSum) / static_cast<double>(iterCount)) : 0.0;
    summary.h_char = ComputeMeshCharLength(mgr);

    std::vector<double> pFinal;
    auto pFinalField = fm.getMatrixScalar(pCfg.pressure_field);
    if (pFinalField) pFinal = pFinalField->data;
    if (cfg.analytical_check && pFinal.size() == nCells) {
        summary.l2_error_final = ComputeL2RelativeErrorFinal(mgr, pFinal, cfg, summary.t_end);
    }

    const BarrierProbeResult probe = ComputeBarrierProbe(mgr, pFinal, cfg);
    summary.delta_p_cross_fracture = probe.delta;
    summary.delta_p_cross_norm = probe.delta_norm;
    summary.probe_left_count = probe.left_count;
    summary.probe_right_count = probe.right_count;

    std::ofstream metrics(summary.metrics_csv_path, std::ios::out | std::ios::trunc);
    metrics << "case_name,level,topology,nx,ny,n_cells,h_char,t_end,steps,total_rollbacks,avg_nonlinear_iters,max_nonlinear_iters,dt_init,dt_min,dt_max,enable_non_orth,normal_method,property_mode,solver_route,grad_route,bc_geom_route,local_solver_impl,local_flux_impl,local_assembly_impl,l2_error_final,delta_p_cross_fracture,delta_p_cross_norm,probe_left_count,probe_right_count\n";
    metrics << cfg.case_name << ","
            << cfg.level_dir << ","
            << static_cast<int>(cfg.topology) << ","
            << cfg.nx << ","
            << cfg.ny << ","
            << nCells << ","
            << std::setprecision(12) << summary.h_char << ","
            << summary.t_end << ","
            << summary.steps << ","
            << summary.total_rollbacks << ","
            << summary.avg_iters << ","
            << summary.max_iters << ","
            << cfg.dt_init << ","
            << cfg.dt_min << ","
            << cfg.dt_max << ","
            << 1 << ","
            << "OverRelaxed" << ","
            << (cfg.use_co2_eos ? "CO2_EOS" : "ConstantBaseline") << ","
            << kAuditSolverRoute << ","
            << kAuditGradRoute << ","
            << kAuditBcGeomRoute << ","
            << (kAuditLocalSolverImpl ? 1 : 0) << ","
            << (kAuditLocalFluxImpl ? 1 : 0) << ","
            << (kAuditLocalAssemblyImpl ? 1 : 0) << ","
            << summary.l2_error_final << ","
            << summary.delta_p_cross_fracture << ","
            << summary.delta_p_cross_norm << ","
            << summary.probe_left_count << ","
            << summary.probe_right_count
            << "\n";

    return summary;
}

N1CaseSpec BuildBaseSpec(const std::string& caseName, const std::string& levelDir, TopologyVariant topology) {
    N1CaseSpec cfg;
    cfg.case_name = caseName;
    cfg.level_dir = levelDir;
    cfg.topology = topology;
    return cfg;
}

N1CasePlan BuildPlanTemplateNoFrac() {
    N1CasePlan plan;
    plan.plan_key = "full_n1_template_const_nowell_nofrac";
    N1CaseSpec cfg = BuildBaseSpec("n1_template_2d_sp_co2_const_nowell_nofrac", "TemplateN1", TopologyVariant::NoFrac);
    cfg.analytical_check = true;
    plan.cases.push_back(cfg);
    return plan;
}

N1CasePlan BuildPlanTemplateSingleFracBarrier() {
    N1CasePlan plan;
    plan.plan_key = "full_n1_template_const_nowell_singlefrac";
    N1CaseSpec cfg = BuildBaseSpec("n1_template_2d_sp_co2_const_nowell_singlefrac", "TemplateN1", TopologyVariant::SingleFrac);
    cfg.analytical_check = false;
    cfg.enable_barrier_probe = true;

    // Barrier-enhanced setup: very low normal conductivity across fracture.
    cfg.matrix_perm = 5.0e-13;
    cfg.matrix_phi = 0.10;
    cfg.matrix_ct = 5.0e-9;

    cfg.fracture_kt = 5.0e-12;
    cfg.fracture_kn = 1.0e-19;
    cfg.fracture_phi = 0.30;
    cfg.fracture_ct = 5.0e-9;

    cfg.p_left = 16.0e6;
    cfg.p_right = 8.0e6;
    cfg.target_end_time_s = 3.0e4;
    cfg.dt_init = 60.0;
    cfg.dt_max = 2000.0;
    plan.cases.push_back(cfg);
    return plan;
}

N1CasePlan BuildPlanTemplateCrossFrac() {
    N1CasePlan plan;
    plan.plan_key = "full_n1_template_const_nowell_crossfrac";
    N1CaseSpec cfg = BuildBaseSpec("n1_template_2d_sp_co2_const_nowell_crossfrac", "TemplateN1", TopologyVariant::CrossFrac);
    cfg.analytical_check = false;
    cfg.enable_barrier_probe = true;
    cfg.fracture_kn = 5.0e-19;
    cfg.fracture_kt = 5.0e-12;
    plan.cases.push_back(cfg);
    return plan;
}

N1CasePlan BuildPlanDay6L1() {
    N1CasePlan plan;
    plan.plan_key = "day6l1_2d_sp_co2_const_nowell_analytical";
    N1CaseSpec cfg = BuildBaseSpec("day6l1_2d_sp_co2_const_nowell_analytical", "L1", TopologyVariant::NoFrac);
    cfg.analytical_check = true;
    cfg.analytical_l2_threshold = 5.0e-2;
    plan.cases.push_back(cfg);
    return plan;
}

N1CasePlan BuildPlanDay6L2() {
    N1CasePlan plan;
    plan.plan_key = "day6l2_2d_sp_co2_const_nowell_grid";
    plan.policy = N1PlanPolicy::GridConvergence;
    plan.aggregate_dir = "Test/Transient/FullCaseTest/N1/L2/day6l2_2d_sp_co2_const_nowell_grid";
    plan.aggregate_csv_name = "grid_convergence.csv";

    struct GridLevel { int nx; int ny; };
    const std::array<GridLevel, 3> grids = { GridLevel{24, 3}, GridLevel{36, 5}, GridLevel{48, 6} };
    for (const auto& g : grids) {
        std::ostringstream name;
        name << "day6l2_2d_sp_co2_const_nowell_grid_" << g.nx << "x" << g.ny;
        N1CaseSpec cfg = BuildBaseSpec(name.str(), "L2", TopologyVariant::NoFrac);
        cfg.nx = g.nx;
        cfg.ny = g.ny;
        cfg.analytical_check = true;
        cfg.analytical_l2_threshold = 5.0e-2;
        plan.cases.push_back(cfg);
    }
    return plan;
}

N1CasePlan BuildPlanDay6L3() {
    N1CasePlan plan;
    plan.plan_key = "day6l3_2d_sp_co2_const_nowell_solver";
    N1CaseSpec cfg = BuildBaseSpec("day6l3_2d_sp_co2_const_nowell_solver", "L3", TopologyVariant::NoFrac);
    cfg.target_end_time_s = 1.2e5;
    cfg.dt_init = 8000.0;
    cfg.dt_min = 1.0;
    cfg.dt_max = 2.0e4;
    cfg.max_newton_iter = 14;
    plan.cases.push_back(cfg);
    return plan;
}

N1CasePlan BuildPlanDay6L4() {
    N1CasePlan plan;
    plan.plan_key = "day6l4_2d_sp_co2_varprop_nowell";
    N1CaseSpec cfg = BuildBaseSpec("day6l4_2d_sp_co2_varprop_nowell", "L4", TopologyVariant::NoFrac);
    cfg.use_co2_eos = true;
    cfg.target_end_time_s = 1.2e5;
    cfg.dt_init = 5000.0;
    cfg.dt_min = 1.0;
    cfg.dt_max = 2.0e4;
    cfg.max_newton_iter = 16;
    plan.cases.push_back(cfg);
    return plan;
}

N1CasePlan BuildPlanDay6L4SingleFrac() {
    N1CasePlan plan;
    plan.plan_key = "day6l4_2d_sp_co2_varprop_nowell_singlefrac";
    N1CaseSpec cfg = BuildBaseSpec("day6l4_2d_sp_co2_varprop_nowell_singlefrac", "L4", TopologyVariant::SingleFrac);
    cfg.use_co2_eos = true;
    cfg.target_end_time_s = 1.2e5;
    cfg.dt_init = 5000.0;
    cfg.dt_min = 1.0;
    cfg.dt_max = 2.0e4;
    cfg.max_newton_iter = 16;
    cfg.enable_barrier_probe = true;
    plan.cases.push_back(cfg);
    return plan;
}

N1CasePlan BuildPlanDay6L1Legacy() {
    N1CasePlan plan = BuildPlanDay6L1();
    plan.plan_key = "day6l1_2d_sp_co2_const_nowell_analytical_legacy";
    plan.cases[0].case_name = "day6l1_2d_sp_co2_const_nowell_analytical_legacy";
    plan.cases[0].level_dir = "Legacy/L1";
    return plan;
}

N1CasePlan BuildPlanDay6L2Legacy() {
    N1CasePlan plan = BuildPlanDay6L2();
    plan.plan_key = "day6l2_2d_sp_co2_const_nowell_grid_legacy";
    plan.aggregate_dir = "Test/Transient/FullCaseTest/N1/Legacy/L2/day6l2_2d_sp_co2_const_nowell_grid_legacy";
    for (auto& c : plan.cases) {
        c.case_name += "_legacy";
        c.level_dir = "Legacy/L2";
    }
    return plan;
}

N1CasePlan BuildPlanDay6L3Legacy() {
    N1CasePlan plan = BuildPlanDay6L3();
    plan.plan_key = "day6l3_2d_sp_co2_const_nowell_solver_legacy";
    plan.cases[0].case_name = "day6l3_2d_sp_co2_const_nowell_solver_legacy";
    plan.cases[0].level_dir = "Legacy/L3";
    return plan;
}

N1CasePlan BuildPlanDay6L4Legacy() {
    N1CasePlan plan = BuildPlanDay6L4();
    plan.plan_key = "day6l4_2d_sp_co2_varprop_nowell_legacy";
    plan.cases[0].case_name = "day6l4_2d_sp_co2_varprop_nowell_legacy";
    plan.cases[0].level_dir = "Legacy/L4";
    return plan;
}

const std::unordered_map<std::string, BuilderFn>& GetN1Registry() {
    static const std::unordered_map<std::string, BuilderFn> registry = {
        {"full_n1_template_const_nowell_nofrac", &BuildPlanTemplateNoFrac},
        {"full_n1_template_const_nowell_singlefrac", &BuildPlanTemplateSingleFracBarrier},
        {"full_n1_template_const_nowell_crossfrac", &BuildPlanTemplateCrossFrac},
        {"day6l1_2d_sp_co2_const_nowell_analytical", &BuildPlanDay6L1},
        {"day6l2_2d_sp_co2_const_nowell_grid", &BuildPlanDay6L2},
        {"day6l3_2d_sp_co2_const_nowell_solver", &BuildPlanDay6L3},
        {"day6l4_2d_sp_co2_varprop_nowell", &BuildPlanDay6L4},
        {"day6l4_2d_sp_co2_varprop_nowell_singlefrac", &BuildPlanDay6L4SingleFrac},
        {"day6l1_2d_sp_co2_const_nowell_analytical_legacy", &BuildPlanDay6L1Legacy},
        {"day6l2_2d_sp_co2_const_nowell_grid_legacy", &BuildPlanDay6L2Legacy},
        {"day6l3_2d_sp_co2_const_nowell_solver_legacy", &BuildPlanDay6L3Legacy},
        {"day6l4_2d_sp_co2_varprop_nowell_legacy", &BuildPlanDay6L4Legacy}
    };
    return registry;
}

void WriteGridConvergenceCsv(const N1CasePlan& plan, const std::vector<N1CaseSummary>& summaries) {
    if (plan.aggregate_dir.empty() || plan.aggregate_csv_name.empty()) return;
    EnsureDirRecursive(plan.aggregate_dir);
    const std::string csvPath = plan.aggregate_dir + "/" + plan.aggregate_csv_name;

    std::ofstream ofs(csvPath, std::ios::out | std::ios::trunc);
    ofs << "case_name,nx,ny,h_char,l2_error,order_vs_prev\n";
    for (size_t i = 0; i < summaries.size(); ++i) {
        double order = std::numeric_limits<double>::quiet_NaN();
        if (i > 0) {
            const double ePrev = summaries[i - 1].l2_error_final;
            const double eNow = summaries[i].l2_error_final;
            const double hPrev = summaries[i - 1].h_char;
            const double hNow = summaries[i].h_char;
            if (ePrev > 0.0 && eNow > 0.0 && hPrev > 0.0 && hNow > 0.0 && hPrev != hNow) {
                order = std::log(ePrev / eNow) / std::log(hPrev / hNow);
            }
        }
        ofs << plan.cases[i].case_name << ","
            << summaries[i].nx << ","
            << summaries[i].ny << ","
            << summaries[i].h_char << ","
            << summaries[i].l2_error_final << ","
            << order << "\n";
    }
}

void ValidatePlanAcceptance(const N1CasePlan& plan, const std::vector<N1CaseSummary>& summaries) {
    for (size_t i = 0; i < plan.cases.size(); ++i) {
        const N1CaseSpec& spec = plan.cases[i];
        const N1CaseSummary& sum = summaries[i];
        if (spec.analytical_check && !(sum.l2_error_final <= spec.analytical_l2_threshold)) {
            std::ostringstream oss;
            oss << "[FullCase-N1] analytical validation failed for " << spec.case_name
                << ": l2_error_final=" << sum.l2_error_final
                << " > " << spec.analytical_l2_threshold;
            throw std::runtime_error(oss.str());
        }
    }

    if (plan.policy == N1PlanPolicy::GridConvergence) {
        if (summaries.size() < 2) {
            throw std::runtime_error("[FullCase-N1] L2 grid convergence requires at least 2 grid levels.");
        }
        for (size_t i = 1; i < summaries.size(); ++i) {
            const double ePrev = summaries[i - 1].l2_error_final;
            const double eNow = summaries[i].l2_error_final;
            if (!(eNow <= ePrev + 1.0e-14)) {
                std::ostringstream oss;
                oss << "[FullCase-N1] L2 non-monotonic error: level " << i
                    << " error=" << eNow << " > prev=" << ePrev;
                throw std::runtime_error(oss.str());
            }
        }
    }
}

std::vector<N1CaseSummary> ExecutePlanByKey(const std::string& registryKey) {
    const auto& registry = GetN1Registry();
    const auto it = registry.find(registryKey);
    if (it == registry.end()) {
        throw std::runtime_error("[FullCase-N1] unknown registry key: " + registryKey);
    }

    const N1CasePlan plan = it->second();
    if (plan.cases.empty()) {
        throw std::runtime_error("[FullCase-N1] empty case plan: " + registryKey);
    }

    std::vector<N1CaseSummary> summaries;
    summaries.reserve(plan.cases.size());

    for (const auto& spec : plan.cases) {
        const N1CaseSummary summary = RunN1Common(spec);
        summaries.push_back(summary);
        std::cout << "[FullCase-N1] case_dir=" << summary.case_dir << "\n";
        std::cout << "[FullCase-N1] steps=" << summary.steps
                  << " rollbacks=" << summary.total_rollbacks
                  << " avg_iters=" << summary.avg_iters
                  << " max_iters=" << summary.max_iters
                  << " l2_error_final=" << summary.l2_error_final
                  << " delta_p_cross=" << summary.delta_p_cross_fracture
                  << " delta_p_cross_norm=" << summary.delta_p_cross_norm << "\n";
    }

    WriteGridConvergenceCsv(plan, summaries);
    ValidatePlanAcceptance(plan, summaries);
    return summaries;
}

// ===== Migrated from Test_ValidationSuite.cpp (N=2/N=3) =====
namespace ValMigrated {

namespace {

// ============================================================
// Constants
// ============================================================
constexpr double kPi = 3.14159265358979323846;

// Fluid (water)
constexpr double kMu_w   = 1.0e-3;    // Pa*s
constexpr double kRho_w  = 1000.0;    // kg/m^3
constexpr double kCp_w   = 4200.0;    // J/(kg*K)
constexpr double kCt_w   = 4.5e-10;   // Pa^-1 (fluid compressibility)
constexpr double kLam_w  = 0.6;       // W/(m*K)

// Rock
constexpr double kRho_r  = 2650.0;    // kg/m^3
constexpr double kCp_r   = 920.0;     // J/(kg*K)
constexpr double kLam_r  = 3.0;       // W/(m*K)

// Two-phase (CO2-water Corey params for analytical BL)
constexpr double kKrw_max = 0.5;
constexpr double kKro_max = 0.8;
constexpr double kNw      = 2.0;
constexpr double kNo      = 2.0;
constexpr double kSwi     = 0.2;
constexpr double kSor     = 0.2;
constexpr double kMu_o    = 5.0e-5;   // CO2 viscosity [Pa*s]

// ============================================================
// Directory helper
// ============================================================
static void EnsureDir(const std::string& raw) {
    if (raw.empty()) return;
    std::string p = raw;
    for (char& c : p) if (c == '\\') c = '/';
    std::stringstream ss(p); std::string tok, cur;
    while (std::getline(ss, tok, '/')) {
        if (tok.empty() || tok == ".") continue;
        if (!cur.empty()) cur += "/";
        cur += tok;
        FULLCASE_MKDIR(cur.c_str());
    }
}

static bool FileNonEmpty(const std::string& p) {
    std::ifstream f(p, std::ios::binary | std::ios::ate);
    return f.good() && f.tellg() > 0;
}

// ============================================================
// Analytical solutions
// ============================================================

/**
 * @brief 1D diffusion: Fourier series solution P(x,t) for [0,L] with
 *   Dirichlet BC at x=0 (P_left), x=L (P_right), IC P(x,0)=P0.
 *   D = hydraulic or thermal diffusivity.
 */
static double Analytic1D(double x, double t, double L,
                          double P_left, double P_right, double P0,
                          double D, int nterms = 100)
{
    const double dP     = P_right - P_left;
    const double P_ss   = P_left + dP * (x / L);          // steady-state
    const double phi0   = P0 - P_ss;                       // IC deviation from steady-state

    // phi0 = P0 - P_left - (P_right-P_left)*x/L
    // bn = (2/L) * integral_0^L phi0 * sin(n*pi*x/L) dx
    //    = 2*(P0-P_left)/(n*pi)*(1-(-1)^n) - 2*(P_right-P_left)*(-1)^n/(n*pi)
    // simplified: phi0 = a - b*x/L  with a=P0-P_left, b=P_right-P_left
    const double a = P0 - P_left;
    const double b = dP;
    double sum = 0.0;
    for (int n = 1; n <= nterms; ++n) {
        const double np = static_cast<double>(n) * kPi;
        const double sign_n = (n % 2 == 0) ? 1.0 : -1.0;
        // bn = (2/np)*(a*(1-sign_n) - b*sign_n)
        const double bn  = (2.0 / np) * (a * (1.0 - sign_n) - b * sign_n);
        const double lam = np / L;
        sum += bn * std::sin(lam * x) * std::exp(-D * lam * lam * t);
    }
    return P_ss + sum;
}

/**
 * @brief Theis well function W(u) = -Ei(-u), via series expansion.
 * W(u) = -gamma - ln(u) + u - u^2/8 + ...  for u < 1
 * For u >= 10: use exponential integral upper bound approximation
 */
static double TheisW(double u) {
    if (u <= 0.0) return 1.0e30;
    if (u >= 60.0) return 0.0;

    // Use series expansion for all u (converges well for u<10)
    // W(u) = -Ei(-u) = integral_u^inf exp(-t)/t dt
    // For u<1: power series is accurate
    // For u>=1: use asymptotic series (careful)
    if (u < 10.0) {
        // Series: W(u) = -EulerGamma - ln(u) + u - u^2/(2*2!) + u^3/(3*3!) - ...
        const double gamma = 0.5772156649015328606;
        double sum = 0.0;
        double term = 1.0;
        double sign = 1.0;
        for (int k = 1; k <= 50; ++k) {
            term *= -u / static_cast<double>(k);
            sum += sign * term / static_cast<double>(k);
            sign = -sign;
        }
        // Equivalent form: W(u) = -EulerGamma - ln(u) - sum_{k=1}^inf (-u)^k / (k*k!)
        // Re-implement cleanly:
        double s2 = 0.0;
        double uk = u;
        for (int k = 1; k <= 80; ++k) {
            double fac = 1.0;
            for (int j = 1; j <= k; ++j) fac *= j;
            double contrib = (k % 2 == 0 ? 1.0 : -1.0) * uk / (k * fac);
            s2 += contrib;
            uk *= u;
            if (std::abs(contrib) < 1e-15 * std::abs(s2)) break;
        }
        return -gamma - std::log(u) - s2;
    }
    else {
        // Asymptotic series: W(u) ~ exp(-u)/u * (1 - 1/u + 2/u^2 - 6/u^3 + ...)
        // Only first few terms are useful for moderate u
        double emu = std::exp(-u);
        if (emu < 1e-30) return 0.0;
        double s = 1.0;
        double term2 = 1.0;
        for (int k = 1; k <= 8; ++k) {
            term2 *= -static_cast<double>(k) / u;
            if (std::abs(term2) > std::abs(s)) break;
            s += term2;
        }
        return emu / u * s;
    }
}

// Corey relative permeability (for analytical BL)
static double Corey_krw(double Sw) {
    double Se = std::max(0.0, std::min(1.0, (Sw - kSwi) / (1.0 - kSwi - kSor)));
    return kKrw_max * std::pow(Se, kNw);
}
static double Corey_kro(double Sw) {
    double Se = std::max(0.0, std::min(1.0, (1.0 - Sw - kSor) / (1.0 - kSwi - kSor)));
    return kKro_max * std::pow(Se, kNo);
}
static double BL_fw(double Sw) {
    double lam_w = Corey_krw(Sw) / kMu_w;
    double lam_o = Corey_kro(Sw) / kMu_o;
    double denom = lam_w + lam_o;
    return (denom > 1e-30) ? lam_w / denom : 0.0;
}
static double BL_dfw_dSw(double Sw) {
    const double eps = 1e-5;
    return (BL_fw(Sw + eps) - BL_fw(Sw - eps)) / (2.0 * eps);
}

/**
 * @brief Welge tangent method: find frontal saturation Sw_f.
 * Returns Sw_f such that tangent from (Swi, fw(Swi)) touches fw curve.
 */
static double BL_FrontalSw() {
    const double fw_swi = BL_fw(kSwi);
    const int N = 2000;
    double Sw_f = kSwi + 0.01;
    double min_err = 1e30;
    for (int i = 1; i <= N; ++i) {
        double Sw = kSwi + (1.0 - kSwi - kSor - 1e-5) * i / N;
        double fw_sw = BL_fw(Sw);
        double chord = (Sw > kSwi + 1e-10) ? (fw_sw - fw_swi) / (Sw - kSwi) : 0.0;
        double slope = BL_dfw_dSw(Sw);
        double err = std::abs(slope - chord);
        if (err < min_err) { min_err = err; Sw_f = Sw; }
    }
    return Sw_f;
}

// ============================================================
// Error computation helpers
// ============================================================
struct ErrorMetrics {
    double L_inf = 0.0;
    double L2    = 0.0;
};

static ErrorMetrics ComputeError(const std::vector<double>& num,
                                  const std::vector<double>& ana,
                                  double scale) {
    if (num.size() != ana.size() || num.empty()) return {};
    ErrorMetrics em;
    double sum2 = 0.0;
    for (size_t i = 0; i < num.size(); ++i) {
        double e = std::abs(num[i] - ana[i]) / std::max(std::abs(scale), 1e-30);
        em.L_inf = std::max(em.L_inf, e);
        sum2 += e * e;
    }
    em.L2 = std::sqrt(sum2 / static_cast<double>(num.size()));
    return em;
}

static void PrintRelError(const std::string& tag,
                           const ErrorMetrics& em,
                           double tol_warn)
{
    const bool over = (em.L_inf > tol_warn);
    std::cout << "[VAL] " << tag
              << "  L_inf=" << std::scientific << std::setprecision(3) << em.L_inf
              << "  L2="    << std::scientific << std::setprecision(3) << em.L2
              << "  tol="   << tol_warn;
    if (over) std::cout << "  *** WARNING: L_inf > tol ***";
    std::cout << "\n";
}

static void PrintStability(const std::string& tag, bool ok) {
    std::cout << "[VAL] " << tag << (ok ? "  PASS (converged)" : "  *** WARNING: run failed ***") << "\n";
}

// ============================================================
// Field reading helpers
// ============================================================
static std::vector<double> ReadField(FieldManager_2D& fm,
                                      const std::string& name,
                                      int nCells) {
    auto f = fm.getOrCreateMatrixScalar(name, 0.0);
    if (!f || static_cast<int>(f->data.size()) < nCells) return {};
    return std::vector<double>(f->data.begin(), f->data.begin() + nCells);
}

// ============================================================
// Mesh / topology builders
// ============================================================
static void BuildMatrixOnly(MeshManager& mgr) {
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
}

static void AddSingleFrac(MeshManager& mgr, double Lx, double Ly) {
    // Diagonal fracture from (0.2Lx, 0.1Ly) to (0.8Lx, 0.9Ly)
    mgr.addFracture(Vector(0.2 * Lx, 0.1 * Ly, 0.0),
                    Vector(0.8 * Lx, 0.9 * Ly, 0.0));
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
}

static void AddCrossFrac(MeshManager& mgr, double Lx, double Ly) {
    // Horizontal frac: (0.1Lx, 0.5Ly) -> (0.9Lx, 0.5Ly)
    mgr.addFracture(Vector(0.1 * Lx, 0.5 * Ly, 0.0),
                    Vector(0.9 * Lx, 0.5 * Ly, 0.0));
    // Vertical frac: (0.5Lx, 0.1Ly) -> (0.5Lx, 0.9Ly)
    mgr.addFracture(Vector(0.5 * Lx, 0.1 * Ly, 0.0),
                    Vector(0.5 * Lx, 0.9 * Ly, 0.0));
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
}

// ============================================================
// Preset builders
// ============================================================
enum class FracConfig { None, Single, Cross };
enum class FluidConfig { PressureOnly, ConstantWater, VarWater };

static FIM_CaseKit::PropertyPreset2D BuildPreset_SP(
    double phi, double k, double c_r,
    double rho_r, double cp_r, double k_r,
    FluidConfig fc)
{
    auto p = FIM_CaseKit::MakeDefaultPropertyPreset2D();
    p.enable_rock_region = false;
    p.rock_bg.phi_r  = phi;
    p.rock_bg.kxx    = k;
    p.rock_bg.kyy    = k;
    p.rock_bg.kzz    = k;
    p.rock_bg.compressibility = c_r;
    p.rock_bg.rho_r  = rho_r;
    p.rock_bg.cp_r   = cp_r;
    p.rock_bg.k_r    = k_r;
    p.frac.permeability   = 1.0e-11;
    p.frac.aperture       = 1.0e-3;
    p.frac.phi_f          = 0.6;
    p.frac.compressibility = c_r;
    p.frac.rho_f = rho_r;
    p.frac.cp_f  = cp_r;
    p.frac.k_f   = k_r;
    p.use_unified_fluid_model = true;
    FluidConstantProperties water_props;
    if (fc == FluidConfig::PressureOnly) {
        p.fluid_model = FIM_Engine::UnifiedFluidModelConfig::MakeSinglePhaseWaterConstant(water_props, true);
        p.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::ConstantWaterNoConvection;
    }
    else if (fc == FluidConfig::ConstantWater) {
        p.fluid_model = FIM_Engine::UnifiedFluidModelConfig::MakeSinglePhaseWaterConstant(water_props, false);
        p.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::ConstantWater;
    }
    else {
        p.fluid_model = FIM_Engine::UnifiedFluidModelConfig::MakeSinglePhaseWaterEOS();
        p.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::Water;
    }
    return p;
}

static FIM_CaseKit::PropertyPreset2D BuildPreset_TP(
    double phi, double k, double c_r,
    double rho_r, double cp_r, double k_r,
    FIM_Engine::FluidPropertyMode tp_mode = FIM_Engine::FluidPropertyMode::EOS)
{
    // Two-phase public switch: Water+CO2 EOS or Water+CO2 constant properties.
    auto p = FIM_CaseKit::MakeDefaultPropertyPreset2D();
    p.enable_rock_region = false;
    p.rock_bg.phi_r  = phi;
    p.rock_bg.kxx    = k;
    p.rock_bg.kyy    = k;
    p.rock_bg.kzz    = k;
    p.rock_bg.compressibility = c_r;
    p.rock_bg.rho_r  = rho_r;
    p.rock_bg.cp_r   = cp_r;
    p.rock_bg.k_r    = k_r;
    p.frac.permeability   = 1.0e-11;
    p.frac.aperture       = 1.0e-3;
    p.frac.phi_f          = 0.6;
    p.frac.compressibility = c_r;
    p.frac.rho_f = rho_r;
    p.frac.cp_f  = cp_r;
    p.frac.k_f   = k_r;
    p.use_unified_fluid_model = true;
    if (tp_mode == FIM_Engine::FluidPropertyMode::Constant) {
        p.fluid_model = FIM_Engine::UnifiedFluidModelConfig::MakeTwoPhaseWaterCO2Constant(
            FluidConstantProperties{},
            FluidConstantProperties{ 700.0, 5.0e-5, 1200.0, 900.0, 0.08 });
    }
    else {
        p.fluid_model = FIM_Engine::UnifiedFluidModelConfig::MakeTwoPhaseWaterCO2EOS();
    }

    // vG params: Swr=0.2, Sgr=0.2 matching Swi/Sor
    p.vg.Swr    = kSwi;
    p.vg.Sgr    = kSor;
    p.vg.n      = 2.0;
    p.vg.alpha  = 5.0e-6;
    p.vg.Pc_max = 1.0e6;
    p.vg.Se_eps = 1.0e-3;
    p.rp.L      = 0.5;
    return p;
}

// ============================================================
// Solver parameters builder for validation
// ============================================================
static FIM_Engine::TransientSolverParams BuildValParams(
    bool twoPhase, double dt_init, double dt_max,
    double t_final, int max_steps)
{
    auto params = FIM_CaseKit::BuildSolverParams(twoPhase, max_steps, dt_init);
    params.gravity_vector = Vector(0.0, 0.0, 0.0);
    params.enable_non_orthogonal_correction = true;
    params.target_end_time_s = t_final;
    params.dt_max    = dt_max;
    params.dt_min    = dt_init * 0.01;
    params.diag_level = FIM_Engine::DiagLevel::Off;
    params.enable_ls_trace = false;
    params.lin_solver = twoPhase
        ? FIM_Engine::LinearSolverType::AMGCL_CPR
        : FIM_Engine::LinearSolverType::AMGCL;
    if (twoPhase) {
        params.max_dSw = 0.08;
        params.max_newton_iter = 20;
    }
    else {
        params.max_newton_iter = 14;
    }
    return params;
}

// ============================================================
// T1: No-well single-phase pressure diffusion
//     N=2, 40脳4 grid, 400m脳40m
//     ConstantWater model, T frozen at T_init
// ============================================================
struct T1Config {
    std::string tag;
    FracConfig  frac;
    FluidConfig fluid;
    bool        compare_analytic;
};

static void RunT1(const T1Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 400.0, Ly = 40.0;
    const int nx = 40, ny = 4;
    const double phi  = 0.10;
    const double k    = 1.0e-13;   // 100 mD
    const double c_r  = kCt_w;    // rock compressibility = fluid compressibility (ConstantWater)
    const double P0   = 10.0e6;
    const double P_L  = 12.0e6;
    const double P_R  = 8.0e6;
    const double T0   = 380.0;
    const double t_final = 2.0e5;
    const double D_h  = k / (kMu_w * phi * kCt_w);

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_SP(phi, k, c_r, kRho_r, kCp_r, kLam_r, cfg.fluid);

    BoundarySetting::BoundaryConditionManager bcP, bcT;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P_L);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P_R);
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcP.SetNeumannBC(MeshTags::TOP,    0.0);

    // Freeze T: all sides Dirichlet at T0
    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T0);
    bcT.SetDirichletBC(MeshTags::RIGHT,  T0);
    bcT.SetDirichletBC(MeshTags::BOTTOM, T0);
    bcT.SetDirichletBC(MeshTags::TOP,    T0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, nullptr);

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = 1.0;

    auto params = BuildValParams(false, 100.0, 2.0e4, t_final, 50000);
    // ConstantWaterNoConvection: P equation is linear (const rho, const mu).
    // (a) PTC: row_floor=1.0 >> real Jacobian (~4e-7) 鈫?Newton step O(1 Pa) instead of O(1 MPa).
    // (b) Armijo: probe lacks non-orthogonal correction 鈫?probe_res >> conv_res 鈫?LS-BASE-CHECK
    //     sets use_ref=probe_res; at alpha鈫? trial_res鈮坧robe_res鈮se_ref 鈫?Armijo never passes.
    // (c) dt-aware alpha damping: max_dP_eff = 2e5*(dt/dt_ref); at dt=1s 鈫?2000 Pa cap 鈫?alpha=0.001.
    // Fix: disable PTC+Armijo; raise max_dP so damp_scale never limits below full step.
    params.enable_ptc = false;
    params.enable_armijo_line_search = false;
    params.max_dP = 2.0e8;  // 2e8*(dt=1/dt_ref=100)=2e6 Pa 鈮?any expected 螖P 鈫?alpha=1 always
    params.lin_solver = FIM_Engine::LinearSolverType::SparseLU;

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<2>(
            cfg.tag, mgr, fm, ic, {}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const auto& cells = mgr.mesh().getCells();
        const int nMat = mgr.getMatrixDOFCount();
        auto p_num = ReadField(fm, "P", nMat);
        if (p_num.empty()) {
            std::cout << "[VAL] " << cfg.tag << " WARNING: could not read P field.\n";
            return;
        }

        const double DP = P_L - P_R;
        std::vector<double> p_ana(nMat);
        for (int i = 0; i < nMat; ++i) {
            p_ana[i] = Analytic1D(cells[i].center.m_x, t_final,
                                   Lx, P_L, P_R, P0, D_h, 100);
        }
        auto em = ComputeError(p_num, p_ana, DP);
        PrintRelError(cfg.tag + " P_error/螖P", em, 0.03);

        // Write CSV
        std::ofstream csv(caseDir + "/analytical_compare.csv");
        csv << "cell_id,x,y,p_num,p_ana,rel_err\n";
        for (int i = 0; i < nMat; ++i) {
            double e = std::abs(p_num[i] - p_ana[i]) / std::max(std::abs(DP), 1.0);
            csv << i << "," << cells[i].center.m_x << "," << cells[i].center.m_y
                << "," << p_num[i] << "," << p_ana[i] << "," << e << "\n";
        }
    }
}


// ============================================================
// T2: No-well single-phase thermal diffusion
//     N=2, k=1e-20 (no-flow), 40脳4 grid, 400m脳40m
//     ConstantWater model
// ============================================================
struct T2Config {
    std::string tag;
    FracConfig  frac;
    FluidConfig fluid;
    bool        compare_analytic;
};

static void RunT2(const T2Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 400.0, Ly = 40.0;
    const int nx = 40, ny = 4;
    const double phi  = 0.10;
    const double k    = 1.0e-13;   // still need some k to avoid solver issues
    const double c_r  = kCt_w;
    const double P0   = 10.0e6;
    const double T0   = 380.0;
    const double T_L  = 400.0;
    const double T_R  = 360.0;
    const double t_final = 1.0e7;  // ~115 days

    // Effective thermal diffusivity
    const double lam_eff = phi * kLam_w + (1.0 - phi) * kLam_r;
    const double C_eff   = phi * kRho_w * kCp_w + (1.0 - phi) * kRho_r * kCp_r;
    const double D_T     = lam_eff / C_eff;

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_SP(phi, k, c_r, kRho_r, kCp_r, kLam_r, cfg.fluid);

    BoundarySetting::BoundaryConditionManager bcP, bcT;
    // Same pressure left/right = no Darcy flow
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P0);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P0);
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcP.SetNeumannBC(MeshTags::TOP,    0.0);

    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T_L);
    bcT.SetDirichletBC(MeshTags::RIGHT,  T_R);
    bcT.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcT.SetNeumannBC(MeshTags::TOP,    0.0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, nullptr);

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = 1.0;

    auto params = BuildValParams(false, 1000.0, 5.0e5, t_final, 100000);

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<2>(
            cfg.tag, mgr, fm, ic, {}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const auto& cells = mgr.mesh().getCells();
        const int nMat = mgr.getMatrixDOFCount();
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        auto t_num = ReadField(fm, tCfg.temperatue_field, nMat);
        if (t_num.empty()) {
            std::cout << "[VAL] " << cfg.tag << " WARNING: could not read T field.\n";
            return;
        }

        const double DT = T_L - T_R;
        std::vector<double> t_ana(nMat);
        for (int i = 0; i < nMat; ++i) {
            t_ana[i] = Analytic1D(cells[i].center.m_x, t_final,
                                   Lx, T_L, T_R, T0, D_T, 100);
        }
        auto em = ComputeError(t_num, t_ana, DT);
        PrintRelError(cfg.tag + " T_error/螖T", em, 0.03);

        std::ofstream csv(caseDir + "/analytical_compare.csv");
        csv << "cell_id,x,y,t_num,t_ana,rel_err\n";
        for (int i = 0; i < nMat; ++i) {
            double e = std::abs(t_num[i] - t_ana[i]) / std::max(std::abs(DT), 1.0);
            csv << i << "," << cells[i].center.m_x << "," << cells[i].center.m_y
                << "," << t_num[i] << "," << t_ana[i] << "," << e << "\n";
        }
    }
}


// ============================================================
// T3: No-well two-phase Buckley-Leverett
//     N=3, 40脳4 grid, 400m脳40m
// ============================================================
struct T3Config {
    std::string tag;
    FracConfig  frac;
    bool        compare_analytic;
};

static void RunT3(const T3Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 400.0, Ly = 40.0;
    const int nx = 40, ny = 4;
    const double phi = 0.20;
    const double k   = 1.0e-13;
    const double c_r = kCt_w;
    const double P_L = 13.0e6, P_R = 11.0e6, P0 = 12.0e6;
    const double T0  = 380.0;
    const double Sw0 = kSwi;
    const double t_final = 5.0e6;

    // Compute BL frontal saturation and position
    const double Sw_f = BL_FrontalSw();
    const double fw_f = BL_fw(Sw_f);
    const double dfw_f = BL_dfw_dSw(Sw_f);
    // Total flux approx: q_darcy = k * A * 螖P / (mu_w * L)
    const double A_cross = Ly * 1.0;  // unit depth
    const double q_vol   = k * A_cross * (P_L - P_R) / (kMu_w * Lx);  // rough average
    const double x_f_ana = dfw_f * q_vol * t_final / (phi * A_cross);

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(3);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_TP(phi, k, c_r, kRho_r, kCp_r, kLam_r);

    BoundarySetting::BoundaryConditionManager bcP, bcT, bcS;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P_L);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P_R);
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcP.SetNeumannBC(MeshTags::TOP,    0.0);

    // Freeze T
    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T0);
    bcT.SetDirichletBC(MeshTags::RIGHT,  T0);
    bcT.SetDirichletBC(MeshTags::BOTTOM, T0);
    bcT.SetDirichletBC(MeshTags::TOP,    T0);

    bcS.Clear();
    bcS.SetDirichletBC(MeshTags::LEFT,  1.0);   // inject water Sw=1
    bcS.SetDirichletBC(MeshTags::RIGHT, Sw0);
    bcS.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcS.SetNeumannBC(MeshTags::TOP,    0.0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, &bcS);

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = Sw0;

    auto params = BuildValParams(true, 1000.0, 1.0e5, t_final, 200000);

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<3>(
            cfg.tag, mgr, fm, ic, {}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const auto& cells = mgr.mesh().getCells();
        const int nMat = mgr.getMatrixDOFCount();
        auto sw_num = ReadField(fm, "S_w", nMat);
        if (sw_num.empty()) {
            std::cout << "[VAL] " << cfg.tag << " WARNING: could not read S_w field.\n";
            return;
        }

        // Find numerical frontal position: where Sw transitions from high to low
        // Use x where Sw = (1.0 + Sw0) / 2
        double Sw_mid = 0.5 * (1.0 + Sw0);
        double x_f_num = -1.0;
        for (int i = 0; i < nMat; ++i) {
            if (sw_num[i] < Sw_mid) {
                x_f_num = cells[i].center.m_x;
                break;
            }
        }

        if (x_f_num > 0.0) {
            double err = std::abs(x_f_num - x_f_ana) / Lx;
            std::cout << "[VAL] " << cfg.tag
                      << "  Sw_f=" << Sw_f << "  x_f_ana=" << x_f_ana
                      << "  x_f_num=" << x_f_num
                      << "  rel_err/L=" << err;
            if (err > 0.05) std::cout << "  *** WARNING: front position error > 5% ***";
            std::cout << "\n";
        }
        else {
            std::cout << "[VAL] " << cfg.tag << "  WARNING: front not detected in domain.\n";
        }
    }
}


// ============================================================
// T4: No-well two-phase + heat (BL + thermal front coupling)
//     N=3, 40脳4, 400m脳40m
// ============================================================
struct T4Config {
    std::string tag;
    FracConfig  frac;
    bool        compare_analytic;
};

static void RunT4(const T4Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 400.0, Ly = 40.0;
    const int nx = 40, ny = 4;
    const double phi = 0.20;
    const double k   = 1.0e-13;
    const double c_r = kCt_w;
    const double P_L = 13.0e6, P_R = 11.0e6, P0 = 12.0e6;
    const double T0  = 380.0, T_inj = 350.0;
    const double Sw0 = kSwi;
    const double t_final = 5.0e6;

    // Adiabatic thermal front velocity (piston approximation)
    // v_T = q / (phi * A) * (rho_w * cp_w) / (rho_w * cp_w + (1-phi)/phi * rho_r * cp_r)
    // Total flux q = k * A * 螖P / (mu * L)  [approximate with water mobility]
    const double A_cross = Ly * 1.0;
    const double q_vol   = k * A_cross * (P_L - P_R) / (kMu_w * Lx);
    const double v_q     = q_vol / (phi * A_cross);
    const double Rch     = (kRho_w * kCp_w) /
                           (kRho_w * kCp_w + (1.0 - phi) / phi * kRho_r * kCp_r);
    const double x_T_ana = v_q * Rch * t_final;

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(3);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_TP(phi, k, c_r, kRho_r, kCp_r, kLam_r);

    BoundarySetting::BoundaryConditionManager bcP, bcT, bcS;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P_L);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P_R);
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcP.SetNeumannBC(MeshTags::TOP,    0.0);

    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T_inj);   // cold injection
    bcT.SetDirichletBC(MeshTags::RIGHT,  T0);
    bcT.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcT.SetNeumannBC(MeshTags::TOP,    0.0);

    bcS.Clear();
    bcS.SetDirichletBC(MeshTags::LEFT,  1.0);
    bcS.SetDirichletBC(MeshTags::RIGHT, Sw0);
    bcS.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcS.SetNeumannBC(MeshTags::TOP,    0.0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, &bcS);

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = Sw0;

    auto params = BuildValParams(true, 1000.0, 1.0e5, t_final, 200000);

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<3>(
            cfg.tag, mgr, fm, ic, {}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const auto& cells = mgr.mesh().getCells();
        const int nMat = mgr.getMatrixDOFCount();
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        auto t_num = ReadField(fm, tCfg.temperatue_field, nMat);
        if (t_num.empty()) { std::cout << "[VAL] " << cfg.tag << " WARNING: T field missing.\n"; return; }

        // Find numerical thermal front: where T = (T_inj + T0)/2
        double T_mid = 0.5 * (T_inj + T0);
        double x_T_num = -1.0;
        for (int i = 0; i < nMat; ++i) {
            if (t_num[i] < T_mid) { x_T_num = cells[i].center.m_x; break; }
        }

        if (x_T_num > 0.0) {
            double err = std::abs(x_T_num - x_T_ana) / Lx;
            std::cout << "[VAL] " << cfg.tag
                      << "  x_T_ana=" << x_T_ana << "  x_T_num=" << x_T_num
                      << "  rel_err/L=" << err;
            if (err > 0.10) std::cout << "  *** WARNING: thermal front error > 10% ***";
            std::cout << "\n";
        }
        else {
            std::cout << "[VAL] " << cfg.tag << "  INFO: thermal front not reached domain (ok if x_T_ana > L).\n";
        }
    }
}


// ============================================================
// T5: With-well single-phase Theis (pressure drawdown)
//     N=2, 50脳50 grid, 2000m脳2000m
//     Rate-controlled production well at center
// ============================================================
struct T5Config {
    std::string tag;
    FracConfig  frac;
    FluidConfig fluid;
    bool        compare_analytic;
};

static void RunT5(const T5Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 2000.0, Ly = 2000.0;
    const int nx = 50, ny = 50;
    const double phi = 0.10;
    const double k   = 1.0e-14;   // 10 mD
    const double c_r = kCt_w;
    const double h   = 1.0;       // unit thickness (2D, Lz=1)
    const double P0  = 30.0e6;
    const double T0  = 380.0;
    const double t_final = 3.0e5;
    const double Q_mass  = 0.1;   // kg/s production rate

    // Theis parameters
    // T_darcy = k * h / mu  [m^2/(Pa*s)]
    const double T_darcy = k * h / kMu_w;
    // S = phi * ct * h  (storativity, dimensionless)
    const double S = phi * kCt_w * h;
    const double Q_vol  = Q_mass / kRho_w;  // m^3/s

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_SP(phi, k, c_r, kRho_r, kCp_r, kLam_r, cfg.fluid);

    // All boundaries: Dirichlet P=P0 (constant pressure outer boundary)
    BoundarySetting::BoundaryConditionManager bcP, bcT;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P0);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P0);
    bcP.SetDirichletBC(MeshTags::BOTTOM, P0);
    bcP.SetDirichletBC(MeshTags::TOP,    P0);

    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T0);
    bcT.SetDirichletBC(MeshTags::RIGHT,  T0);
    bcT.SetDirichletBC(MeshTags::BOTTOM, T0);
    bcT.SetDirichletBC(MeshTags::TOP,    T0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, nullptr);

    // Find center cell for production well
    const auto& cells = mgr.mesh().getCells();
    int well_cell = 0;
    double best_d2 = 1e30;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        double dx = cells[i].center.m_x - 0.5 * Lx;
        double dy = cells[i].center.m_y - 0.5 * Ly;
        double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) { best_d2 = d2; well_cell = i; }
    }

    WellScheduleStep prod;
    prod.well_name     = "PROD_T5";
    prod.domain        = WellTargetDomain::Matrix;
    prod.control_mode  = WellControlMode::Rate;
    prod.target_value  = -Q_mass;   // negative = production
    prod.completion_id = well_cell;
    prod.component_mode = WellComponentMode::Water;
    prod.frac_w = 1.0;
    prod.frac_g = 0.0;

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = 1.0;

    auto params = BuildValParams(false, 100.0, 3.0e4, t_final, 50000);
    // Same issues as T1: PTC destroys Newton direction, Armijo probe lacks non-orth correction,
    // and dt-aware damping caps max_dP_eff at 2 kPa for dt=1s.
    params.enable_ptc = false;
    params.enable_armijo_line_search = false;
    params.max_dP = 2.0e8;
    params.lin_solver = FIM_Engine::LinearSolverType::SparseLU;

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<2>(
            cfg.tag, mgr, fm, ic, {prod}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const int nMat = mgr.getMatrixDOFCount();
        auto p_num = ReadField(fm, "P", nMat);
        if (p_num.empty()) { std::cout << "[VAL] " << cfg.tag << " WARNING: P field missing.\n"; return; }

        const double wx = cells[well_cell].center.m_x;
        const double wy = cells[well_cell].center.m_y;

        // Theis comparison at several radii
        std::array<double,4> r_vals = { 56.0, 112.0, 200.0, 400.0 };
        double max_rel_err = 0.0;
        std::cout << "[VAL] " << cfg.tag << " Theis comparison at t=" << t_final << "s:\n";
        for (double r : r_vals) {
            // Find nearest cell at distance r from well
            int best_i = 0; double best_dr = 1e30;
            for (int i = 0; i < nMat; ++i) {
                double dist = std::hypot(cells[i].center.m_x - wx, cells[i].center.m_y - wy);
                if (std::abs(dist - r) < best_dr) { best_dr = std::abs(dist - r); best_i = i; }
            }
            double r_actual = std::hypot(cells[best_i].center.m_x - wx, cells[best_i].center.m_y - wy);
            double u  = r_actual * r_actual * S / (4.0 * T_darcy * t_final);
            // Drawdown: s(r,t) = Q_vol/(4*pi*T_darcy) * W(u)
            // Note: Q_vol in [m^3/s], T_darcy in [m^2/(Pa*s)], so:
            // s = Q_vol/(4*pi*T_darcy) * W(u) has units [m^3/s * Pa*s/m^2] = Pa
            double s_ana = Q_vol / (4.0 * kPi * T_darcy) * TheisW(u);
            double P_ana = P0 - s_ana;
            double dP_max = Q_vol / (4.0 * kPi * T_darcy) * TheisW(r_vals[0] * r_vals[0] * S / (4.0 * T_darcy * t_final));
            double rel_err = (dP_max > 1.0) ? std::abs(p_num[best_i] - P_ana) / dP_max : 0.0;
            max_rel_err = std::max(max_rel_err, rel_err);
            std::cout << "    r=" << r_actual << "  u=" << u << "  W(u)=" << TheisW(u)
                      << "  P_ana=" << P_ana << "  P_num=" << p_num[best_i]
                      << "  rel_err=" << rel_err << "\n";
        }
        if (max_rel_err > 0.05) std::cout << "  *** WARNING: max Theis relative error > 5% ***\n";
    }
}


// ============================================================
// T6: With-well single-phase thermal front
//     N=2, 50脳50, 2000m脳2000m, injection well at center
// ============================================================
struct T6Config {
    std::string tag;
    FracConfig  frac;
    FluidConfig fluid;
    bool        compare_analytic;
};

static void RunT6(const T6Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 2000.0, Ly = 2000.0;
    const int nx = 50, ny = 50;
    const double phi = 0.10;
    const double k   = 1.0e-13;
    const double c_r = kCt_w;
    const double P0  = 30.0e6;
    const double T0  = 380.0, T_inj = 340.0;
    const double t_final = 5.0e6;
    const double Q_mass  = 0.1;   // kg/s injection

    const double Q_vol  = Q_mass / kRho_w;
    const double Rch    = kRho_w * kCp_w /
                          (phi * kRho_w * kCp_w + (1.0 - phi) * kRho_r * kCp_r);
    // Thermal front radius: r_T = sqrt(Q_vol * Rch * t / pi)
    const double r_T_ana = std::sqrt(Q_vol * Rch * t_final / kPi);

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_SP(phi, k, c_r, kRho_r, kCp_r, kLam_r, cfg.fluid);

    BoundarySetting::BoundaryConditionManager bcP, bcT;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P0);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P0);
    bcP.SetDirichletBC(MeshTags::BOTTOM, P0);
    bcP.SetDirichletBC(MeshTags::TOP,    P0);

    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T0);
    bcT.SetDirichletBC(MeshTags::RIGHT,  T0);
    bcT.SetDirichletBC(MeshTags::BOTTOM, T0);
    bcT.SetDirichletBC(MeshTags::TOP,    T0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, nullptr);

    const auto& cells = mgr.mesh().getCells();
    int well_cell = 0; double best_d2 = 1e30;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        double dx = cells[i].center.m_x - 0.5 * Lx;
        double dy = cells[i].center.m_y - 0.5 * Ly;
        double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) { best_d2 = d2; well_cell = i; }
    }

    WellScheduleStep inj;
    inj.well_name          = "INJ_T6";
    inj.domain             = WellTargetDomain::Matrix;
    inj.control_mode       = WellControlMode::Rate;
    inj.target_value       = Q_mass;
    inj.completion_id      = well_cell;
    inj.component_mode     = WellComponentMode::Water;
    inj.frac_w             = 1.0;
    inj.frac_g             = 0.0;
    inj.injection_temperature = T_inj;

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = 1.0;

    auto params = BuildValParams(false, 100.0, 1.0e5, t_final, 100000);

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<2>(
            cfg.tag, mgr, fm, ic, {inj}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const int nMat = mgr.getMatrixDOFCount();
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        auto t_num = ReadField(fm, tCfg.temperatue_field, nMat);
        if (t_num.empty()) { std::cout << "[VAL] " << cfg.tag << " WARNING: T field missing.\n"; return; }

        const double wx = cells[well_cell].center.m_x;
        const double wy = cells[well_cell].center.m_y;
        double T_mid = 0.5 * (T_inj + T0);

        // Find radius where T = T_mid
        double r_T_num = -1.0;
        double max_r_seen = 0.0;
        for (int i = 0; i < nMat; ++i) {
            double r = std::hypot(cells[i].center.m_x - wx, cells[i].center.m_y - wy);
            if (r > max_r_seen) max_r_seen = r;
            if (t_num[i] < T_mid && r > r_T_num) r_T_num = r;
        }

        std::cout << "[VAL] " << cfg.tag
                  << "  r_T_ana=" << r_T_ana
                  << "  r_T_num=" << r_T_num;
        if (r_T_num > 0.0) {
            double err = std::abs(r_T_num - r_T_ana) / r_T_ana;
            std::cout << "  rel_err=" << err;
            if (err > 0.10) std::cout << "  *** WARNING: thermal radius error > 10% ***";
        }
        std::cout << "\n";
    }
}


// ============================================================
// T7: With-well two-phase radial BL
//     N=3, 50脳50, 2000m脳2000m, water injection at center
// ============================================================
struct T7Config {
    std::string tag;
    FracConfig  frac;
    bool        compare_analytic;
};

static void RunT7(const T7Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 2000.0, Ly = 2000.0;
    const int nx = 50, ny = 50;
    const double phi = 0.20;
    const double k   = 1.0e-13;
    const double c_r = kCt_w;
    const double P0  = 25.0e6;
    const double T0  = 380.0;
    const double Sw0 = kSwi;
    const double t_final = 3.0e6;
    const double Q_mass  = 0.1;

    const double Q_vol  = Q_mass / kRho_w;
    // Radial BL front: r_f虏 鈮?Q_vol * fw(Sw_f) * t / (pi * phi)
    const double Sw_f   = BL_FrontalSw();
    const double fw_f   = BL_fw(Sw_f);
    const double r_f_ana = std::sqrt(Q_vol * fw_f * t_final / (kPi * phi));

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(3);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_TP(phi, k, c_r, kRho_r, kCp_r, kLam_r);

    BoundarySetting::BoundaryConditionManager bcP, bcT, bcS;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P0);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P0);
    bcP.SetDirichletBC(MeshTags::BOTTOM, P0);
    bcP.SetDirichletBC(MeshTags::TOP,    P0);

    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T0);
    bcT.SetDirichletBC(MeshTags::RIGHT,  T0);
    bcT.SetDirichletBC(MeshTags::BOTTOM, T0);
    bcT.SetDirichletBC(MeshTags::TOP,    T0);

    bcS.Clear();
    bcS.SetNeumannBC(MeshTags::LEFT,   0.0);
    bcS.SetNeumannBC(MeshTags::RIGHT,  0.0);
    bcS.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcS.SetNeumannBC(MeshTags::TOP,    0.0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, &bcS);

    const auto& cells = mgr.mesh().getCells();
    int well_cell = 0; double best_d2 = 1e30;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        double dx = cells[i].center.m_x - 0.5 * Lx;
        double dy = cells[i].center.m_y - 0.5 * Ly;
        double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) { best_d2 = d2; well_cell = i; }
    }

    WellScheduleStep inj;
    inj.well_name          = "INJ_T7";
    inj.domain             = WellTargetDomain::Matrix;
    inj.control_mode       = WellControlMode::Rate;
    inj.target_value       = Q_mass;
    inj.completion_id      = well_cell;
    inj.component_mode     = WellComponentMode::Water;
    inj.frac_w             = 1.0;
    inj.frac_g             = 0.0;
    inj.injection_is_co2   = false;

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = Sw0;

    auto params = BuildValParams(true, 1000.0, 1.0e5, t_final, 200000);

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<3>(
            cfg.tag, mgr, fm, ic, {inj}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const int nMat = mgr.getMatrixDOFCount();
        auto sw_num = ReadField(fm, "S_w", nMat);
        if (sw_num.empty()) { std::cout << "[VAL] " << cfg.tag << " WARNING: S_w field missing.\n"; return; }

        const double wx = cells[well_cell].center.m_x;
        const double wy = cells[well_cell].center.m_y;
        double Sw_mid = 0.5 * (1.0 + Sw0);

        // Find numerical frontal radius
        double r_f_num = -1.0;
        for (int i = 0; i < nMat; ++i) {
            double r = std::hypot(cells[i].center.m_x - wx, cells[i].center.m_y - wy);
            if (sw_num[i] > Sw_mid && r > r_f_num) r_f_num = r;
        }

        std::cout << "[VAL] " << cfg.tag
                  << "  Sw_f=" << Sw_f << "  fw_f=" << fw_f
                  << "  r_f_ana=" << r_f_ana << "  r_f_num=" << r_f_num;
        if (r_f_num > 0.0 && r_f_ana > 0.0) {
            double err = std::abs(r_f_num - r_f_ana) / r_f_ana;
            std::cout << "  rel_err=" << err;
            if (err > 0.10) std::cout << "  *** WARNING: radial BL front error > 10% ***";
        }
        std::cout << "\n";
    }
}


// ============================================================
// T8: With-well two-phase radial BL + heat
//     N=3, 50脳50, 2000m脳2000m, cold water injection
// ============================================================
struct T8Config {
    std::string tag;
    FracConfig  frac;
    bool        compare_analytic;
};

static void RunT8(const T8Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 2000.0, Ly = 2000.0;
    const int nx = 50, ny = 50;
    const double phi = 0.20;
    const double k   = 1.0e-13;
    const double c_r = kCt_w;
    const double P0  = 25.0e6;
    const double T0  = 380.0, T_inj = 340.0;
    const double Sw0 = kSwi;
    const double t_final = 3.0e6;
    const double Q_mass  = 0.1;

    const double Q_vol  = Q_mass / kRho_w;
    const double Sw_f   = BL_FrontalSw();
    const double fw_f   = BL_fw(Sw_f);
    const double r_f_ana = std::sqrt(Q_vol * fw_f * t_final / (kPi * phi));
    // Thermal front: two-phase heat capacity
    const double C_tp   = phi * (Sw0 * kRho_w * kCp_w + (1.0 - Sw0) * kRho_r * kCp_r * 0.0)
                        + (1.0 - phi) * kRho_r * kCp_r;  // approx Sw=Swi region
    // Simpler: use single-fluid capacity for approximate thermal front
    const double Rch_tp = kRho_w * kCp_w /
                          (phi * kRho_w * kCp_w + (1.0 - phi) * kRho_r * kCp_r);
    const double r_T_ana = std::sqrt(Q_vol * Rch_tp * t_final / kPi);

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(3);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_TP(phi, k, c_r, kRho_r, kCp_r, kLam_r);

    BoundarySetting::BoundaryConditionManager bcP, bcT, bcS;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P0);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P0);
    bcP.SetDirichletBC(MeshTags::BOTTOM, P0);
    bcP.SetDirichletBC(MeshTags::TOP,    P0);

    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T0);
    bcT.SetDirichletBC(MeshTags::RIGHT,  T0);
    bcT.SetDirichletBC(MeshTags::BOTTOM, T0);
    bcT.SetDirichletBC(MeshTags::TOP,    T0);

    bcS.Clear();
    bcS.SetNeumannBC(MeshTags::LEFT,   0.0);
    bcS.SetNeumannBC(MeshTags::RIGHT,  0.0);
    bcS.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcS.SetNeumannBC(MeshTags::TOP,    0.0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, &bcS);

    const auto& cells = mgr.mesh().getCells();
    int well_cell = 0; double best_d2 = 1e30;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        double dx = cells[i].center.m_x - 0.5 * Lx;
        double dy = cells[i].center.m_y - 0.5 * Ly;
        double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) { best_d2 = d2; well_cell = i; }
    }

    WellScheduleStep inj;
    inj.well_name             = "INJ_T8";
    inj.domain                = WellTargetDomain::Matrix;
    inj.control_mode          = WellControlMode::Rate;
    inj.target_value          = Q_mass;
    inj.completion_id         = well_cell;
    inj.component_mode        = WellComponentMode::Water;
    inj.frac_w                = 1.0;
    inj.frac_g                = 0.0;
    inj.injection_is_co2      = false;
    inj.injection_temperature = T_inj;

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = Sw0;

    auto params = BuildValParams(true, 1000.0, 1.0e5, t_final, 200000);

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<3>(
            cfg.tag, mgr, fm, ic, {inj}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const int nMat = mgr.getMatrixDOFCount();
        auto sw_num = ReadField(fm, "S_w", nMat);
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        auto t_num = ReadField(fm, tCfg.temperatue_field, nMat);

        const double wx = cells[well_cell].center.m_x;
        const double wy = cells[well_cell].center.m_y;

        // BL front
        if (!sw_num.empty()) {
            double Sw_mid = 0.5 * (1.0 + Sw0);
            double r_f_num = -1.0;
            for (int i = 0; i < nMat; ++i) {
                double r = std::hypot(cells[i].center.m_x - wx, cells[i].center.m_y - wy);
                if (sw_num[i] > Sw_mid && r > r_f_num) r_f_num = r;
            }
            std::cout << "[VAL] " << cfg.tag
                      << "  BL: r_f_ana=" << r_f_ana << "  r_f_num=" << r_f_num;
            if (r_f_num > 0.0 && r_f_ana > 0.0) {
                double err = std::abs(r_f_num - r_f_ana) / r_f_ana;
                std::cout << "  rel_err=" << err;
                if (err > 0.10) std::cout << "  *** WARNING ***";
            }
            std::cout << "\n";
        }

        // Thermal front
        if (!t_num.empty()) {
            double T_mid = 0.5 * (T_inj + T0);
            double r_T_num = -1.0;
            for (int i = 0; i < nMat; ++i) {
                double r = std::hypot(cells[i].center.m_x - wx, cells[i].center.m_y - wy);
                if (t_num[i] < T_mid && r > r_T_num) r_T_num = r;
            }
            std::cout << "[VAL] " << cfg.tag
                      << "  TH: r_T_ana=" << r_T_ana << "  r_T_num=" << r_T_num;
            if (r_T_num > 0.0 && r_T_ana > 0.0) {
                double err = std::abs(r_T_num - r_T_ana) / r_T_ana;
                std::cout << "  rel_err=" << err;
                if (err > 0.15) std::cout << "  *** WARNING ***";
            }
            std::cout << "\n";
        }
    }
}

} // namespace (anonymous)

// ============================================================
// Public wrapper definitions (outside anonymous namespace,
// matching declarations in Test_ValidationSuite.h)
// ============================================================

// T1
void Val_T1_Const_NoFrac()    { RunT1({"Val_T1_Const_NoFrac",    FracConfig::None,   FluidConfig::PressureOnly, true});  }
void Val_T1_Const_SingleFrac(){ RunT1({"Val_T1_Const_SingleFrac",FracConfig::Single, FluidConfig::PressureOnly, false}); }
void Val_T1_Const_CrossFrac() { RunT1({"Val_T1_Const_CrossFrac", FracConfig::Cross,  FluidConfig::PressureOnly, false}); }
void Val_T1_Var_NoFrac()      { RunT1({"Val_T1_Var_NoFrac",      FracConfig::None,   FluidConfig::VarWater,      false}); }
void Val_T1_Var_SingleFrac()  { RunT1({"Val_T1_Var_SingleFrac",  FracConfig::Single, FluidConfig::VarWater,      false}); }
void Val_T1_Var_CrossFrac()   { RunT1({"Val_T1_Var_CrossFrac",   FracConfig::Cross,  FluidConfig::VarWater,      false}); }

// T2
void Val_T2_Const_NoFrac()    { RunT2({"Val_T2_Const_NoFrac",    FracConfig::None,   FluidConfig::ConstantWater, true});  }
void Val_T2_Const_SingleFrac(){ RunT2({"Val_T2_Const_SingleFrac",FracConfig::Single, FluidConfig::ConstantWater, false}); }
void Val_T2_Const_CrossFrac() { RunT2({"Val_T2_Const_CrossFrac", FracConfig::Cross,  FluidConfig::ConstantWater, false}); }
void Val_T2_Var_NoFrac()      { RunT2({"Val_T2_Var_NoFrac",      FracConfig::None,   FluidConfig::VarWater,      false}); }
void Val_T2_Var_SingleFrac()  { RunT2({"Val_T2_Var_SingleFrac",  FracConfig::Single, FluidConfig::VarWater,      false}); }
void Val_T2_Var_CrossFrac()   { RunT2({"Val_T2_Var_CrossFrac",   FracConfig::Cross,  FluidConfig::VarWater,      false}); }

// T3
void Val_T3_Const_NoFrac()    { RunT3({"Val_T3_Const_NoFrac",    FracConfig::None,   true});  }
void Val_T3_Const_SingleFrac(){ RunT3({"Val_T3_Const_SingleFrac",FracConfig::Single, false}); }
void Val_T3_Const_CrossFrac() { RunT3({"Val_T3_Const_CrossFrac", FracConfig::Cross,  false}); }
void Val_T3_Var_NoFrac()      { RunT3({"Val_T3_Var_NoFrac",      FracConfig::None,   false}); }
void Val_T3_Var_SingleFrac()  { RunT3({"Val_T3_Var_SingleFrac",  FracConfig::Single, false}); }
void Val_T3_Var_CrossFrac()   { RunT3({"Val_T3_Var_CrossFrac",   FracConfig::Cross,  false}); }

// T4
void Val_T4_Const_NoFrac()    { RunT4({"Val_T4_Const_NoFrac",    FracConfig::None,   true});  }
void Val_T4_Const_SingleFrac(){ RunT4({"Val_T4_Const_SingleFrac",FracConfig::Single, false}); }
void Val_T4_Const_CrossFrac() { RunT4({"Val_T4_Const_CrossFrac", FracConfig::Cross,  false}); }
void Val_T4_Var_NoFrac()      { RunT4({"Val_T4_Var_NoFrac",      FracConfig::None,   false}); }
void Val_T4_Var_SingleFrac()  { RunT4({"Val_T4_Var_SingleFrac",  FracConfig::Single, false}); }
void Val_T4_Var_CrossFrac()   { RunT4({"Val_T4_Var_CrossFrac",   FracConfig::Cross,  false}); }

// T5
void Val_T5_Const_NoFrac()    { RunT5({"Val_T5_Const_NoFrac",    FracConfig::None,   FluidConfig::PressureOnly, true});  }
void Val_T5_Const_SingleFrac(){ RunT5({"Val_T5_Const_SingleFrac",FracConfig::Single, FluidConfig::PressureOnly, false}); }
void Val_T5_Const_CrossFrac() { RunT5({"Val_T5_Const_CrossFrac", FracConfig::Cross,  FluidConfig::PressureOnly, false}); }
void Val_T5_Var_NoFrac()      { RunT5({"Val_T5_Var_NoFrac",      FracConfig::None,   FluidConfig::VarWater,      false}); }
void Val_T5_Var_SingleFrac()  { RunT5({"Val_T5_Var_SingleFrac",  FracConfig::Single, FluidConfig::VarWater,      false}); }
void Val_T5_Var_CrossFrac()   { RunT5({"Val_T5_Var_CrossFrac",   FracConfig::Cross,  FluidConfig::VarWater,      false}); }

// T6
void Val_T6_Const_NoFrac()    { RunT6({"Val_T6_Const_NoFrac",    FracConfig::None,   FluidConfig::ConstantWater, true});  }
void Val_T6_Const_SingleFrac(){ RunT6({"Val_T6_Const_SingleFrac",FracConfig::Single, FluidConfig::ConstantWater, false}); }
void Val_T6_Const_CrossFrac() { RunT6({"Val_T6_Const_CrossFrac", FracConfig::Cross,  FluidConfig::ConstantWater, false}); }
void Val_T6_Var_NoFrac()      { RunT6({"Val_T6_Var_NoFrac",      FracConfig::None,   FluidConfig::VarWater,      false}); }
void Val_T6_Var_SingleFrac()  { RunT6({"Val_T6_Var_SingleFrac",  FracConfig::Single, FluidConfig::VarWater,      false}); }
void Val_T6_Var_CrossFrac()   { RunT6({"Val_T6_Var_CrossFrac",   FracConfig::Cross,  FluidConfig::VarWater,      false}); }

// T7
void Val_T7_Const_NoFrac()    { RunT7({"Val_T7_Const_NoFrac",    FracConfig::None,   true});  }
void Val_T7_Const_SingleFrac(){ RunT7({"Val_T7_Const_SingleFrac",FracConfig::Single, false}); }
void Val_T7_Const_CrossFrac() { RunT7({"Val_T7_Const_CrossFrac", FracConfig::Cross,  false}); }
void Val_T7_Var_NoFrac()      { RunT7({"Val_T7_Var_NoFrac",      FracConfig::None,   false}); }
void Val_T7_Var_SingleFrac()  { RunT7({"Val_T7_Var_SingleFrac",  FracConfig::Single, false}); }
void Val_T7_Var_CrossFrac()   { RunT7({"Val_T7_Var_CrossFrac",   FracConfig::Cross,  false}); }

// T8
void Val_T8_Const_NoFrac()    { RunT8({"Val_T8_Const_NoFrac",    FracConfig::None,   true});  }
void Val_T8_Const_SingleFrac(){ RunT8({"Val_T8_Const_SingleFrac",FracConfig::Single, false}); }
void Val_T8_Const_CrossFrac() { RunT8({"Val_T8_Const_CrossFrac", FracConfig::Cross,  false}); }
void Val_T8_Var_NoFrac()      { RunT8({"Val_T8_Var_NoFrac",      FracConfig::None,   false}); }
void Val_T8_Var_SingleFrac()  { RunT8({"Val_T8_Var_SingleFrac",  FracConfig::Single, false}); }
void Val_T8_Var_CrossFrac()   { RunT8({"Val_T8_Var_CrossFrac",   FracConfig::Cross,  false}); }

// ============================================================
// Grouped runners
// ============================================================
void Run_All_T1() {
    std::cout << "\n========== ValidationSuite T1: Pressure Diffusion ==========\n";
    Val_T1_Const_NoFrac();    Val_T1_Const_SingleFrac(); Val_T1_Const_CrossFrac();
    Val_T1_Var_NoFrac();      Val_T1_Var_SingleFrac();   Val_T1_Var_CrossFrac();
}
void Run_All_T2() {
    std::cout << "\n========== ValidationSuite T2: Thermal Diffusion ==========\n";
    Val_T2_Const_NoFrac();    Val_T2_Const_SingleFrac(); Val_T2_Const_CrossFrac();
    Val_T2_Var_NoFrac();      Val_T2_Var_SingleFrac();   Val_T2_Var_CrossFrac();
}
void Run_All_T3() {
    std::cout << "\n========== ValidationSuite T3: Buckley-Leverett ==========\n";
    Val_T3_Const_NoFrac();    Val_T3_Const_SingleFrac(); Val_T3_Const_CrossFrac();
    Val_T3_Var_NoFrac();      Val_T3_Var_SingleFrac();   Val_T3_Var_CrossFrac();
}
void Run_All_T4() {
    std::cout << "\n========== ValidationSuite T4: BL + Heat ==========\n";
    Val_T4_Const_NoFrac();    Val_T4_Const_SingleFrac(); Val_T4_Const_CrossFrac();
    Val_T4_Var_NoFrac();      Val_T4_Var_SingleFrac();   Val_T4_Var_CrossFrac();
}
void Run_All_T5() {
    std::cout << "\n========== ValidationSuite T5: Theis Well ==========\n";
    Val_T5_Const_NoFrac();    Val_T5_Const_SingleFrac(); Val_T5_Const_CrossFrac();
    Val_T5_Var_NoFrac();      Val_T5_Var_SingleFrac();   Val_T5_Var_CrossFrac();
}
void Run_All_T6() {
    std::cout << "\n========== ValidationSuite T6: Thermal Front Well ==========\n";
    Val_T6_Const_NoFrac();    Val_T6_Const_SingleFrac(); Val_T6_Const_CrossFrac();
    Val_T6_Var_NoFrac();      Val_T6_Var_SingleFrac();   Val_T6_Var_CrossFrac();
}
void Run_All_T7() {
    std::cout << "\n========== ValidationSuite T7: Radial BL ==========\n";
    Val_T7_Const_NoFrac();    Val_T7_Const_SingleFrac(); Val_T7_Const_CrossFrac();
    Val_T7_Var_NoFrac();      Val_T7_Var_SingleFrac();   Val_T7_Var_CrossFrac();
}
void Run_All_T8() {
    std::cout << "\n========== ValidationSuite T8: Radial BL + Heat ==========\n";
    Val_T8_Const_NoFrac();    Val_T8_Const_SingleFrac(); Val_T8_Const_CrossFrac();
    Val_T8_Var_NoFrac();      Val_T8_Var_SingleFrac();   Val_T8_Var_CrossFrac();
}

void Run_ValidationSuite_All() {
    std::cout << "\n#################### ValidationSuite ALL 48 CASES ####################\n";
    Run_All_T1(); Run_All_T2(); Run_All_T3(); Run_All_T4();
    Run_All_T5(); Run_All_T6(); Run_All_T7(); Run_All_T8();
    std::cout << "\n#################### ValidationSuite COMPLETE ####################\n";
}

} // namespace ValMigrated
// ===== End migrated validation block =====

using LegacyValFn = void(*)();

struct LegacyValCase {
    const char* key;
    int t_id; // 1..8
    LegacyValFn fn;
};

const std::vector<LegacyValCase>& GetLegacyValCaseRegistry() {
    static const std::vector<LegacyValCase> kCases = {
        {"val_t1_const_nofrac", 1, &ValMigrated::Val_T1_Const_NoFrac},
        {"val_t1_const_singlefrac", 1, &ValMigrated::Val_T1_Const_SingleFrac},
        {"val_t1_const_crossfrac", 1, &ValMigrated::Val_T1_Const_CrossFrac},
        {"val_t1_var_nofrac", 1, &ValMigrated::Val_T1_Var_NoFrac},
        {"val_t1_var_singlefrac", 1, &ValMigrated::Val_T1_Var_SingleFrac},
        {"val_t1_var_crossfrac", 1, &ValMigrated::Val_T1_Var_CrossFrac},

        {"val_t2_const_nofrac", 2, &ValMigrated::Val_T2_Const_NoFrac},
        {"val_t2_const_singlefrac", 2, &ValMigrated::Val_T2_Const_SingleFrac},
        {"val_t2_const_crossfrac", 2, &ValMigrated::Val_T2_Const_CrossFrac},
        {"val_t2_var_nofrac", 2, &ValMigrated::Val_T2_Var_NoFrac},
        {"val_t2_var_singlefrac", 2, &ValMigrated::Val_T2_Var_SingleFrac},
        {"val_t2_var_crossfrac", 2, &ValMigrated::Val_T2_Var_CrossFrac},

        {"val_t3_const_nofrac", 3, &ValMigrated::Val_T3_Const_NoFrac},
        {"val_t3_const_singlefrac", 3, &ValMigrated::Val_T3_Const_SingleFrac},
        {"val_t3_const_crossfrac", 3, &ValMigrated::Val_T3_Const_CrossFrac},
        {"val_t3_var_nofrac", 3, &ValMigrated::Val_T3_Var_NoFrac},
        {"val_t3_var_singlefrac", 3, &ValMigrated::Val_T3_Var_SingleFrac},
        {"val_t3_var_crossfrac", 3, &ValMigrated::Val_T3_Var_CrossFrac},

        {"val_t4_const_nofrac", 4, &ValMigrated::Val_T4_Const_NoFrac},
        {"val_t4_const_singlefrac", 4, &ValMigrated::Val_T4_Const_SingleFrac},
        {"val_t4_const_crossfrac", 4, &ValMigrated::Val_T4_Const_CrossFrac},
        {"val_t4_var_nofrac", 4, &ValMigrated::Val_T4_Var_NoFrac},
        {"val_t4_var_singlefrac", 4, &ValMigrated::Val_T4_Var_SingleFrac},
        {"val_t4_var_crossfrac", 4, &ValMigrated::Val_T4_Var_CrossFrac},

        {"val_t5_const_nofrac", 5, &ValMigrated::Val_T5_Const_NoFrac},
        {"val_t5_const_singlefrac", 5, &ValMigrated::Val_T5_Const_SingleFrac},
        {"val_t5_const_crossfrac", 5, &ValMigrated::Val_T5_Const_CrossFrac},
        {"val_t5_var_nofrac", 5, &ValMigrated::Val_T5_Var_NoFrac},
        {"val_t5_var_singlefrac", 5, &ValMigrated::Val_T5_Var_SingleFrac},
        {"val_t5_var_crossfrac", 5, &ValMigrated::Val_T5_Var_CrossFrac},

        {"val_t6_const_nofrac", 6, &ValMigrated::Val_T6_Const_NoFrac},
        {"val_t6_const_singlefrac", 6, &ValMigrated::Val_T6_Const_SingleFrac},
        {"val_t6_const_crossfrac", 6, &ValMigrated::Val_T6_Const_CrossFrac},
        {"val_t6_var_nofrac", 6, &ValMigrated::Val_T6_Var_NoFrac},
        {"val_t6_var_singlefrac", 6, &ValMigrated::Val_T6_Var_SingleFrac},
        {"val_t6_var_crossfrac", 6, &ValMigrated::Val_T6_Var_CrossFrac},

        {"val_t7_const_nofrac", 7, &ValMigrated::Val_T7_Const_NoFrac},
        {"val_t7_const_singlefrac", 7, &ValMigrated::Val_T7_Const_SingleFrac},
        {"val_t7_const_crossfrac", 7, &ValMigrated::Val_T7_Const_CrossFrac},
        {"val_t7_var_nofrac", 7, &ValMigrated::Val_T7_Var_NoFrac},
        {"val_t7_var_singlefrac", 7, &ValMigrated::Val_T7_Var_SingleFrac},
        {"val_t7_var_crossfrac", 7, &ValMigrated::Val_T7_Var_CrossFrac},

        {"val_t8_const_nofrac", 8, &ValMigrated::Val_T8_Const_NoFrac},
        {"val_t8_const_singlefrac", 8, &ValMigrated::Val_T8_Const_SingleFrac},
        {"val_t8_const_crossfrac", 8, &ValMigrated::Val_T8_Const_CrossFrac},
        {"val_t8_var_nofrac", 8, &ValMigrated::Val_T8_Var_NoFrac},
        {"val_t8_var_singlefrac", 8, &ValMigrated::Val_T8_Var_SingleFrac},
        {"val_t8_var_crossfrac", 8, &ValMigrated::Val_T8_Var_CrossFrac}
    };
    return kCases;
}

const LegacyValCase& FindLegacyValCaseByKeyOrThrow(const std::string& key) {
    const auto& cases = GetLegacyValCaseRegistry();
    const auto it = std::find_if(
        cases.begin(), cases.end(),
        [&key](const LegacyValCase& c) { return key == c.key; });
    if (it == cases.end()) {
        throw std::runtime_error("[FullCaseTest] unknown migrated validation key: " + key);
    }
    return *it;
}

void RunLegacyValCaseByKeyInternal(const std::string& key) {
    const LegacyValCase& c = FindLegacyValCaseByKeyOrThrow(key);
    std::cout << "[FullCase-ValidationMigrated] run key=" << key << " t=" << c.t_id << "\n";
    c.fn();
}

std::vector<const LegacyValCase*> CollectLegacyValCasesForT(int t_id) {
    std::vector<const LegacyValCase*> out;
    const auto& cases = GetLegacyValCaseRegistry();
    for (const auto& c : cases) {
        if (c.t_id == t_id) out.push_back(&c);
    }
    return out;
}

std::string BuildConstCaseKey(int scenarioId, TopologyVariant topology) {
    const std::string t = "val_t" + std::to_string(scenarioId) + "_const_";
    switch (topology) {
    case TopologyVariant::NoFrac: return t + "nofrac";
    case TopologyVariant::SingleFrac: return t + "singlefrac";
    case TopologyVariant::CrossFrac: return t + "crossfrac";
    default:
        throw std::runtime_error("[FullCaseTest] unknown topology variant.");
    }
}

Test_Day6::CampaignTopologyAxis ToDay6Topology(TopologyVariant topology) {
    switch (topology) {
    case TopologyVariant::NoFrac:
        return Test_Day6::CampaignTopologyAxis::F0;
    case TopologyVariant::SingleFrac:
        return Test_Day6::CampaignTopologyAxis::F1;
    case TopologyVariant::CrossFrac:
        return Test_Day6::CampaignTopologyAxis::F2;
    default:
        throw std::runtime_error("[FullCaseTest] unknown topology variant.");
    }
}

} // namespace

void Run2DCase(int scenarioId, TopologyVariant topology) {
    ValidateScenarioIdOrThrow(scenarioId);
    RunLegacyValCaseByKeyInternal(BuildConstCaseKey(scenarioId, topology));
}

void Run3DCase(int scenarioId, TopologyVariant topology) {
    ValidateScenarioIdOrThrow(scenarioId);
    Test_Day6::Run_Day6_Campaign_3D_Single(scenarioId, ToDay6Topology(topology));
}

void Run2DAll() {
    for (int scenarioId = 1; scenarioId <= 8; ++scenarioId) {
        for (const auto topology : kTopologyOrder) {
            Run2DCase(scenarioId, topology);
        }
    }
}

void Run3DAll() {
    for (int scenarioId = 1; scenarioId <= 8; ++scenarioId) {
        for (const auto topology : kTopologyOrder) {
            Run3DCase(scenarioId, topology);
        }
    }
}

void RunAll() {
    Run2DAll();
    Run3DAll();
}

void RunN1L1ConstNoWellAnalytical() {
    ExecutePlanByKey("day6l1_2d_sp_co2_const_nowell_analytical");
}

void RunN1L2ConstNoWellGrid() {
    ExecutePlanByKey("day6l2_2d_sp_co2_const_nowell_grid");
}

void RunN1L3ConstNoWellSolver() {
    ExecutePlanByKey("day6l3_2d_sp_co2_const_nowell_solver");
}

void RunN1L4VarPropNoWell() {
    ExecutePlanByKey("day6l4_2d_sp_co2_varprop_nowell");
}

void RunN1L4VarPropNoWellSingleFrac() {
    ExecutePlanByKey("day6l4_2d_sp_co2_varprop_nowell_singlefrac");
}

void RunN1LadderAll() {
    RunN1L1ConstNoWellAnalytical();
    RunN1L2ConstNoWellGrid();
    RunN1L3ConstNoWellSolver();
    RunN1L4VarPropNoWell();
    std::cout << "[Day6Ladder] ALL PASS\n";
}

void RunN1L1ConstNoWellAnalyticalLegacy() {
    ExecutePlanByKey("day6l1_2d_sp_co2_const_nowell_analytical_legacy");
}

void RunN1L2ConstNoWellGridLegacy() {
    ExecutePlanByKey("day6l2_2d_sp_co2_const_nowell_grid_legacy");
}

void RunN1L3ConstNoWellSolverLegacy() {
    ExecutePlanByKey("day6l3_2d_sp_co2_const_nowell_solver_legacy");
}

void RunN1L4VarPropNoWellLegacy() {
    ExecutePlanByKey("day6l4_2d_sp_co2_varprop_nowell_legacy");
}

void RunN1LadderAllLegacy() {
    RunN1L1ConstNoWellAnalyticalLegacy();
    RunN1L2ConstNoWellGridLegacy();
    RunN1L3ConstNoWellSolverLegacy();
    RunN1L4VarPropNoWellLegacy();
    std::cout << "[Day6Ladder-LegacyCompat] ALL PASS\n";
}

void RunN1TemplateConstNoWellNoFrac() {
    ExecutePlanByKey("full_n1_template_const_nowell_nofrac");
}

void RunN1TemplateConstNoWellSingleFrac() {
    ExecutePlanByKey("full_n1_template_const_nowell_singlefrac");
}

void RunN1TemplateConstNoWellCrossFrac() {
    ExecutePlanByKey("full_n1_template_const_nowell_crossfrac");
}

void RunValidationCaseByKey(const std::string& case_key) {
    RunLegacyValCaseByKeyInternal(case_key);
}

void RunValidationGroup(int t_id) {
    ValidateScenarioIdOrThrow(t_id);
    const auto group = CollectLegacyValCasesForT(t_id);
    if (group.empty()) {
        throw std::runtime_error("[FullCaseTest] no migrated validation cases for T" + std::to_string(t_id));
    }
    for (const auto* c : group) c->fn();
}

void RunValidationAllMigrated() {
    for (int t = 1; t <= 8; ++t) {
        RunValidationGroup(t);
    }
}

} // namespace FullCaseTest




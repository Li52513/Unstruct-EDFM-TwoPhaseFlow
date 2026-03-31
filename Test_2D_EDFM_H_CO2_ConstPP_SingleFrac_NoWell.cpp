/**
 * @file Test_2D_EDFM_H_CO2_ConstPP_SingleFrac_NoWell.cpp
 * @brief Standalone test: 2D single-phase CO2 constant-property, single-fracture, no-well.
 */

#include "Test_2D_EDFM_H_CO2_ConstPP_SingleFrac_NoWell.h"

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
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#define TEST_MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define TEST_MKDIR(path) mkdir(path, 0777)
#endif

namespace Test_H_CO2_ConstPP_SingleFrac {
namespace {

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

void ApplyUniformScalarField(const std::shared_ptr<volScalarField>& field, double value) {
    if (!field) return;
    for (double& v : field->data) v = value;
}

void SyncPressureFieldsToFM(MeshManager& mgr,
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
        }
        else {
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
    if (totalV <= 0.0) return 0.0;
    return std::sqrt(totalV / static_cast<double>(cells.size()));
}

struct TestCaseSpec {
    std::string case_name = "h_co2_constpp_singlefrac_nowell";
    std::string output_base_dir = "Test/Transient/FullCaseTest";
    std::string sub_dir = "H_CO2_ConstPP";
    bool use_complex_fractures = false;

    double lx = 400.0;
    double ly = 40.0;
    int nx = 48;
    int ny = 6;

    double matrix_phi = 0.10;
    double matrix_perm = 1.0e-13;
    double matrix_ct = 5.0e-9;

    double fracture_kt = 1.0e-17;
    double fracture_kn = 1.0e-20;
    double fracture_phi = 0.15;
    double fracture_ct = 5.0e-9;

    double mu_const = 6.0e-5;
    double rho_const = 700.0;

    double p_init = 10.0e6;
    double p_left = 12.0e6;
    double p_right = 8.0e6;
    double t_init = 360.0;

    double dt_init = 120.0;
    double dt_min = 1.0;
    double dt_max = 5000.0;
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

    bool enable_armijo_line_search = false;

    double rollback_shrink_factor = 0.7;
    double dt_relres_grow_factor = 1.08;

    bool enable_barrier_probe = true;
    double probe_y_min_ratio = 0.15;
    double probe_y_max_ratio = 0.85;

    Vector gravity_vector = Vector(0.0, -9.81, 0.0);
    FIM_Engine::DiagLevel diag_level = FIM_Engine::DiagLevel::Off;
};

struct TestCaseSummary {
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

    double delta_p_cross_fracture = std::numeric_limits<double>::quiet_NaN();
    double delta_p_cross_norm = std::numeric_limits<double>::quiet_NaN();
    int probe_left_count = 0;
    int probe_right_count = 0;
};

struct BarrierProbeResult {
    bool valid = false;
    double left_avg = std::numeric_limits<double>::quiet_NaN();
    double right_avg = std::numeric_limits<double>::quiet_NaN();
    double delta = std::numeric_limits<double>::quiet_NaN();
    double delta_norm = std::numeric_limits<double>::quiet_NaN();
    int left_count = 0;
    int right_count = 0;
};

struct TestCasePlan {
    std::string plan_key;
    TestCaseSpec spec;
};

void AddConfiguredFractures(MeshManager& mgr, const TestCaseSpec& cfg) {
    mgr.addFracture(Vector(0.2 * cfg.lx, 0.1 * cfg.ly, 0.0),
                    Vector(0.8 * cfg.lx, 0.9 * cfg.ly, 0.0));
    if (!cfg.use_complex_fractures) {
        return;
    }

    // Keep the baseline diagonal and add a deterministic multi-fracture network
    // so A5 stays reproducible and solver-stable inside the current worktree.
    mgr.addFracture(Vector(0.16 * cfg.lx, 0.82 * cfg.ly, 0.0),
                    Vector(0.76 * cfg.lx, 0.18 * cfg.ly, 0.0));
    mgr.addFracture(Vector(0.12 * cfg.lx, 0.54 * cfg.ly, 0.0),
                    Vector(0.88 * cfg.lx, 0.58 * cfg.ly, 0.0));
    mgr.addFracture(Vector(0.44 * cfg.lx, 0.12 * cfg.ly, 0.0),
                    Vector(0.64 * cfg.lx, 0.90 * cfg.ly, 0.0));
}

BarrierProbeResult ComputeBarrierProbe(const MeshManager& mgr,
                                       const std::vector<double>& pMatrix,
                                       const TestCaseSpec& cfg) {
    BarrierProbeResult out;
    if (!cfg.enable_barrier_probe) return out;

    const auto& cells = mgr.mesh().getCells();
    if (cells.size() != pMatrix.size() || cells.empty()) return out;

    const double h = std::max(ComputeMeshCharLength(mgr), 1.0e-12);
    const double fx = 0.5 * cfg.lx;
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

TestCasePlan BuildDefaultPlan() {
    TestCasePlan plan;
    plan.plan_key = "h_co2_constpp_singlefrac_nowell";
    return plan;
}

TestCasePlan BuildComplexPlan() {
    TestCasePlan plan = BuildDefaultPlan();
    plan.plan_key = "h_co2_constpp_complexfrac_nowell";
    plan.spec.case_name = "h_co2_constpp_complexfrac_nowell";
    plan.spec.use_complex_fractures = true;
    return plan;
}

using BuilderFn = TestCasePlan(*)();

const std::unordered_map<std::string, BuilderFn>& GetRegistry() {
    static const std::unordered_map<std::string, BuilderFn> registry = {
        {"h_co2_constpp_singlefrac_nowell", &BuildDefaultPlan},
        {"h_co2_constpp_complexfrac_nowell", &BuildComplexPlan}
    };
    return registry;
}

TestCaseSummary RunCase(const TestCaseSpec& cfg) {
    TestCaseSummary summary;

    summary.case_dir = cfg.output_base_dir + "/" + cfg.sub_dir + "/" + cfg.case_name;
    EnsureDirRecursive(summary.case_dir);
    summary.convergence_log_path = summary.case_dir + "/convergence.log";
    summary.metrics_csv_path = summary.case_dir + "/metrics.csv";
    summary.nx = cfg.nx;
    summary.ny = cfg.ny;

    std::ofstream convergenceLog(summary.convergence_log_path, std::ios::out | std::ios::trunc);
    if (!convergenceLog.good()) {
        throw std::runtime_error("[Test_H_CO2_ConstPP_SingleFrac] failed to open convergence log: " + summary.convergence_log_path);
    }

    MeshManager mgr(cfg.lx, cfg.ly, 0.0, cfg.nx, cfg.ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OverRelaxed);
    AddConfiguredFractures(mgr, cfg);
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);

    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(1);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    const size_t nCells = mgr.mesh().getCells().size();
    const int totalBlocks = mgr.getTotalDOFCount();
    if (nCells == 0) {
        throw std::runtime_error("[Test_H_CO2_ConstPP_SingleFrac] matrix cell count is zero");
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

    const auto pEqCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
    VTKBoundaryVisualizationContext bcVizCtx;
    bcVizCtx.water_family_policy = VTKBCWaterFamilyDerivePolicy::FollowPrimaryFluid;
    bcVizCtx.primary_fluid_model = VTKBCPrimaryFluidModel::CO2;
    bcVizCtx.bindings.push_back(
        VTKBCVariableBinding{ pEqCfg.pressure_field, &bcP, VTKBCTransportKind::Pressure });

    std::vector<double> pBlocksLatest(static_cast<std::size_t>(std::max(totalBlocks, 0)), cfg.p_init);
    SyncPressureFieldsToFM(mgr, fm, pBlocksLatest, cfg.t_init);

    PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/initial.vtk", 0.0);

    bool midExported = false;
    int iterSum = 0;
    int iterCount = 0;
    int maxIters = 0;

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

        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.k_t_tag, cfg.fracture_kt), cfg.fracture_kt);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.k_n_tag, cfg.fracture_kn), cfg.fracture_kn);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.phi_tag, cfg.fracture_phi), cfg.fracture_phi);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.c_r_tag, cfg.fracture_ct), cfg.fracture_ct);

        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(water.k_tag, 0.6), 0.6);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(water.k_tag, 0.6), 0.6);
    };

    modules.on_step_accepted =
        [&](int step, double timeS, double dtUsedS, int newtonIters, double residualInf,
            int totalRollbacks, const std::string& convergeMode,
            const std::vector<double>& pVec, const std::vector<double>&,
            const std::vector<double>*) {
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
                           << " dt=" << dtUsedS
                           << " iters=" << newtonIters
                           << " residual_inf=" << residualInf
                           << " rollbacks=" << totalRollbacks
                           << " mode=" << convergeMode << "\n";

            if (!midExported && timeS >= 0.5 * cfg.target_end_time_s - 1.0e-12) {
                PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/mid.vtk", timeS);
                midExported = true;
            }
        };

    const auto params = BuildSolverParams(cfg);
    FIM_Engine::RunGenericFIMTransient<1>(
        cfg.case_name,
        mgr,
        fm,
        ic,
        {},
        params,
        FIM_Engine::SolverRoute::FIM,
        modules);

    SyncPressureFieldsToFM(mgr, fm, pBlocksLatest, cfg.t_init);
    if (!midExported) {
        PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/mid.vtk", summary.t_end);
    }
    PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/final.vtk", summary.t_end);

    summary.max_iters = maxIters;
    summary.avg_iters = (iterCount > 0) ? (static_cast<double>(iterSum) / static_cast<double>(iterCount)) : 0.0;
    summary.h_char = ComputeMeshCharLength(mgr);

    std::vector<double> pFinal;
    auto pFinalField = fm.getMatrixScalar(pEqCfg.pressure_field);
    if (pFinalField) pFinal = pFinalField->data;

    const BarrierProbeResult probe = ComputeBarrierProbe(mgr, pFinal, cfg);
    summary.delta_p_cross_fracture = probe.delta;
    summary.delta_p_cross_norm = probe.delta_norm;
    summary.probe_left_count = probe.left_count;
    summary.probe_right_count = probe.right_count;

    std::ofstream metrics(summary.metrics_csv_path, std::ios::out | std::ios::trunc);
    metrics << "case_name,nx,ny,n_cells,h_char,t_end,steps,total_rollbacks,"
               "avg_nonlinear_iters,max_nonlinear_iters,"
               "dt_init,dt_min,dt_max,property_mode,solver_route,"
               "delta_p_cross_fracture,delta_p_cross_norm,probe_left_count,probe_right_count\n";

    metrics << cfg.case_name << ","
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
            << "ConstantBaseline" << ","
            << "RunGenericFIMTransient<1>" << ","
            << std::scientific << std::setprecision(8)
            << summary.delta_p_cross_fracture << ","
            << summary.delta_p_cross_norm << ","
            << summary.probe_left_count << ","
            << summary.probe_right_count << "\n";

    return summary;
}

void ExecutePlanByKeyImpl(const std::string& key) {
    const auto& registry = GetRegistry();
    auto it = registry.find(key);
    if (it == registry.end()) {
        throw std::runtime_error("[Test_H_CO2_ConstPP_SingleFrac] unknown registry key: " + key);
    }

    const TestCasePlan plan = it->second();
    const TestCaseSummary summary = RunCase(plan.spec);

    std::cout << "\n============================================\n";
    std::cout << "[Test_H_CO2_ConstPP_SingleFrac] run completed\n";
    std::cout << "  output_dir: " << summary.case_dir << "\n";
    std::cout << "  grid: " << summary.nx << " x " << summary.ny
              << " (" << summary.n_cells << " cells)\n";
    std::cout << "  steps: " << summary.steps
              << "  rollbacks: " << summary.total_rollbacks << "\n";
    std::cout << "  Newton iters: avg=" << std::fixed << std::setprecision(2) << summary.avg_iters
              << "  max=" << summary.max_iters << "\n";
    std::cout << "  final_time: " << std::scientific << std::setprecision(4) << summary.t_end << " s\n";
    std::cout << "  delta_p_cross: " << summary.delta_p_cross_fracture
              << "  norm=" << summary.delta_p_cross_norm << "\n";
    std::cout << "  probes: left=" << summary.probe_left_count
              << "  right=" << summary.probe_right_count << "\n";
    std::cout << "============================================\n";
}

} // namespace

void RunTestCase() {
    ExecutePlanByKeyImpl("h_co2_constpp_singlefrac_nowell");
}

void ExecutePlanByKey(const std::string& key) {
    ExecutePlanByKeyImpl(key);
}

void RunSolveOnly() {
    ExecutePlanByKeyImpl("h_co2_constpp_singlefrac_nowell");
}

void RunPrepareReference() {
    ExecutePlanByKeyImpl("h_co2_constpp_singlefrac_nowell");
}

void RunValidateOnly() {
    ExecutePlanByKeyImpl("h_co2_constpp_singlefrac_nowell");
}

void RunFullWorkflow() {
    ExecutePlanByKeyImpl("h_co2_constpp_singlefrac_nowell");
}

} // namespace Test_H_CO2_ConstPP_SingleFrac

/**
 * @file Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.cpp
 * @brief Standalone test: 2D single-phase CO2 constant-property coupled pressure-temperature, single-fracture, no-well.
 */

#include "Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.h"

#include "2D_PostProcess.h"
#include "BoundaryConditionManager.h"
#include "FIM_TransientCaseKit.hpp"
#include "FractureElement.h"
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

double ComputeMeshCharLength(const MeshManager& mgr) {
    const auto& cells = mgr.mesh().getCells();
    if (cells.empty()) return 0.0;
    double totalV = 0.0;
    for (const auto& c : cells) totalV += std::max(c.volume, 0.0);
    if (totalV <= 0.0) return 0.0;
    return std::sqrt(totalV / static_cast<double>(cells.size()));
}

double Clamp01(double v) {
    return std::max(0.0, std::min(1.0, v));
}

double PointToSegmentDistance(const Vector& point, const Vector& segStart, const Vector& segEnd) {
    const Vector seg = segEnd - segStart;
    const double lenSq = std::max(seg.Mag2(), 1.0e-20);
    const double t = Clamp01(((point - segStart) * seg) / lenSq);
    const Vector projection = segStart + t * seg;
    return (point - projection).Mag();
}

struct LeftRightAverageResult {
    bool valid = false;
    double left_avg = std::numeric_limits<double>::quiet_NaN();
    double right_avg = std::numeric_limits<double>::quiet_NaN();
    int left_count = 0;
    int right_count = 0;
};

LeftRightAverageResult ComputeLeftRightAverage(const MeshManager& mgr,
                                               const std::vector<double>& values,
                                               double xSplit) {
    LeftRightAverageResult out;
    const auto& cells = mgr.mesh().getCells();
    if (cells.size() != values.size() || cells.empty()) return out;

    double leftSum = 0.0;
    double rightSum = 0.0;
    int leftCnt = 0;
    int rightCnt = 0;

    for (size_t i = 0; i < cells.size(); ++i) {
        if (cells[i].center.m_x < xSplit) {
            leftSum += values[i];
            ++leftCnt;
        } else {
            rightSum += values[i];
            ++rightCnt;
        }
    }

    if (leftCnt <= 0 || rightCnt <= 0) return out;
    out.valid = true;
    out.left_avg = leftSum / static_cast<double>(leftCnt);
    out.right_avg = rightSum / static_cast<double>(rightCnt);
    out.left_count = leftCnt;
    out.right_count = rightCnt;
    return out;
}

double ComputeMaxAbsDelta(const std::vector<double>& values, double ref) {
    double maxAbs = 0.0;
    for (double v : values) {
        maxAbs = std::max(maxAbs, std::abs(v - ref));
    }
    return maxAbs;
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
    const int nTotal = mgr.getTotalDOFCount();
    const int nUse = std::min(static_cast<int>(pBlocks.size()), nTotal);

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
    double matrix_cp_r = 1000.0;
    double matrix_lambda_r = 2.0;

    double fracture_phi = 0.15;
    double fracture_kt = 5.0e-12;
    double fracture_kn = 5.0e-13;
    double fracture_ct = 5.0e-9;
    double fracture_rho_r = 2600.0;
    double fracture_cp_r = 1000.0;
    double fracture_lambda_r = 4.0;

    double co2_rho_const = 700.0;
    double co2_mu_const = 6.0e-5;
    double co2_cp_const = 1100.0;
    double co2_cv_const = 850.0;
    double co2_k_const = 0.03;

    double p_init = 10.0e6;
    double p_left = 12.0e6;
    double p_right = 8.0e6;
    double t_init = 380.0;
    double t_left = 400.0;
    double t_right = 360.0;

    double dt_init = 120.0;
    double dt_min = 1.0;
    double dt_max = 5000.0;
    double target_end_time_s = 1.0e5;
    int max_steps = 12000;
    int max_newton_iter = 14;

    FIM_Engine::LinearSolverType lin_solver = FIM_Engine::LinearSolverType::AMGCL;
    bool amgcl_use_fallback_sparselu = true;
    double amgcl_tol = 1.0e-6;
    int amgcl_maxiter = 500;

    bool enable_non_orthogonal_correction = true;
    bool enable_row_scaling = true;

    double abs_res_tol = 1.0e-8;
    double rel_res_tol = 1.0e-6;
    double rel_update_tol = 1.0e-8;

    double max_dP = 2.0e7;
    double max_dT = 10.0;
    bool enable_armijo_line_search = false;

    double rollback_shrink_factor = 0.7;
    double dt_relres_grow_factor = 1.08;

    Vector gravity_vector = Vector(0.0, -9.81, 0.0);
    FIM_Engine::DiagLevel diag_level = FIM_Engine::DiagLevel::Off;

    double min_max_abs_dp_for_evolution = 1.0e3;
    double min_max_abs_dt_for_evolution = 1.0e-1;
};

struct TestCaseSummary {
    std::string case_dir;
    std::string convergence_log_path;
    std::string metrics_csv_path;

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

    double p_left_avg = std::numeric_limits<double>::quiet_NaN();
    double p_right_avg = std::numeric_limits<double>::quiet_NaN();
    int p_left_count = 0;
    int p_right_count = 0;
    int p_lr_pass = 0;

    double t_left_avg = std::numeric_limits<double>::quiet_NaN();
    double t_right_avg = std::numeric_limits<double>::quiet_NaN();
    int t_left_count = 0;
    int t_right_count = 0;
    int t_lr_pass = 0;

    double max_abs_dp = 0.0;
    double max_abs_dt = 0.0;
    int p_evolution_pass = 0;
    int t_evolution_pass = 0;

    int near_fracture_band_count = 0;
    double frac_mean_abs_dp = std::numeric_limits<double>::quiet_NaN();
    double near_band_mean_abs_dp = std::numeric_limits<double>::quiet_NaN();
    int frac_dp_advantage_pass = 0;

    double frac_mean_abs_dt = std::numeric_limits<double>::quiet_NaN();
    double near_band_mean_abs_dt = std::numeric_limits<double>::quiet_NaN();
    int frac_dt_advantage_pass = 0;

    double frac_p_upstream_avg = std::numeric_limits<double>::quiet_NaN();
    double frac_p_downstream_avg = std::numeric_limits<double>::quiet_NaN();
    int frac_p_upstream_count = 0;
    int frac_p_downstream_count = 0;
    int frac_p_along_pass = 0;

    double frac_t_upstream_avg = std::numeric_limits<double>::quiet_NaN();
    double frac_t_downstream_avg = std::numeric_limits<double>::quiet_NaN();
    int frac_t_upstream_count = 0;
    int frac_t_downstream_count = 0;
    int frac_t_along_pass = 0;

    int physics_checks_pass = 0;
};

struct FractureChannelResult {
    bool valid = false;
    int fracture_count = 0;
    int near_matrix_count = 0;
    int upstream_count = 0;
    int downstream_count = 0;
    double fracture_mean_abs_delta = std::numeric_limits<double>::quiet_NaN();
    double near_matrix_mean_abs_delta = std::numeric_limits<double>::quiet_NaN();
    double upstream_avg = std::numeric_limits<double>::quiet_NaN();
    double downstream_avg = std::numeric_limits<double>::quiet_NaN();
};

struct TestCasePlan {
    std::string plan_key;
    TestCaseSpec spec;
};

FractureChannelResult ComputeFractureChannelResult(const MeshManager& mgr,
                                                   const std::vector<double>& blockValues,
                                                   double initialValue,
                                                   const TestCaseSpec& cfg) {
    FractureChannelResult out;
    const int nMat = mgr.getMatrixDOFCount();
    const int nTotal = mgr.getTotalDOFCount();
    if (nMat <= 0 || nTotal <= nMat) return out;
    if (static_cast<int>(blockValues.size()) < nTotal) return out;

    const auto& cells = mgr.mesh().getCells();
    if (static_cast<int>(cells.size()) < nMat) return out;

    const Vector fracStart(cfg.frac_x0_ratio * cfg.lx, cfg.frac_y0_ratio * cfg.ly, 0.0);
    const Vector fracEnd(cfg.frac_x1_ratio * cfg.lx, cfg.frac_y1_ratio * cfg.ly, 0.0);
    const double bandWidth = std::max(ComputeMeshCharLength(mgr), 1.0e-12);

    double fracSumAbs = 0.0;
    double nearSumAbs = 0.0;
    double upSum = 0.0;
    double downSum = 0.0;
    int fracCnt = 0;
    int nearCnt = 0;
    int upCnt = 0;
    int downCnt = 0;

    for (int i = 0; i < nMat; ++i) {
        if (PointToSegmentDistance(cells[static_cast<std::size_t>(i)].center, fracStart, fracEnd) <= bandWidth) {
            nearSumAbs += std::abs(blockValues[static_cast<std::size_t>(i)] - initialValue);
            ++nearCnt;
        }
    }

    for (int solverIdx = nMat; solverIdx < nTotal; ++solverIdx) {
        const FractureElement* elem = mgr.getFractureElementBySolverIndex(solverIdx);
        if (!elem) continue;

        const double value = blockValues[static_cast<std::size_t>(solverIdx)];
        const double tau = Clamp01(0.5 * (elem->param0 + elem->param1));

        fracSumAbs += std::abs(value - initialValue);
        ++fracCnt;

        if (tau < 0.5) {
            upSum += value;
            ++upCnt;
        } else {
            downSum += value;
            ++downCnt;
        }
    }

    if (fracCnt <= 0 || nearCnt <= 0 || upCnt <= 0 || downCnt <= 0) return out;

    out.valid = true;
    out.fracture_count = fracCnt;
    out.near_matrix_count = nearCnt;
    out.upstream_count = upCnt;
    out.downstream_count = downCnt;
    out.fracture_mean_abs_delta = fracSumAbs / static_cast<double>(fracCnt);
    out.near_matrix_mean_abs_delta = nearSumAbs / static_cast<double>(nearCnt);
    out.upstream_avg = upSum / static_cast<double>(upCnt);
    out.downstream_avg = downSum / static_cast<double>(downCnt);
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

TestCasePlan BuildDefaultPlan() {
    TestCasePlan plan;
    plan.plan_key = "h_t_co2_constpp_singlefrac_nowell";
    return plan;
}

using BuilderFn = TestCasePlan(*)();

const std::unordered_map<std::string, BuilderFn>& GetRegistry() {
    static const std::unordered_map<std::string, BuilderFn> registry = {
        {"h_t_co2_constpp_singlefrac_nowell", &BuildDefaultPlan}
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
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] failed to open convergence log: " + summary.convergence_log_path);
    }

    MeshManager mgr(cfg.lx, cfg.ly, 0.0, cfg.nx, cfg.ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OverRelaxed);
    mgr.addFracture(Vector(cfg.frac_x0_ratio * cfg.lx, cfg.frac_y0_ratio * cfg.ly, 0.0),
                    Vector(cfg.frac_x1_ratio * cfg.lx, cfg.frac_y1_ratio * cfg.ly, 0.0));
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    const size_t nCells = mgr.mesh().getCells().size();
    const int nMat = mgr.getMatrixDOFCount();
    const int totalBlocks = mgr.getTotalDOFCount();
    if (nCells == 0) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] matrix cell count is zero");
    }
    if (nMat <= 0 || totalBlocks <= nMat) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] fracture DOF count is zero");
    }
    summary.n_cells = static_cast<int>(nCells);
    summary.n_fracture_dofs = totalBlocks - nMat;

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
    bcVizCtx.bindings.push_back(
        VTKBCVariableBinding{ pEqCfg.pressure_field, &bcP, VTKBCTransportKind::Pressure });
    bcVizCtx.bindings.push_back(
        VTKBCVariableBinding{ tEqCfg.temperatue_field, &bcT, VTKBCTransportKind::Temperature });

    std::vector<double> pBlocksLatest(static_cast<std::size_t>(std::max(totalBlocks, 0)), cfg.p_init);
    std::vector<double> tBlocksLatest(static_cast<std::size_t>(std::max(totalBlocks, 0)), cfg.t_init);
    SyncPTFieldsToFM(mgr, fm, pBlocksLatest, tBlocksLatest, cfg.p_init, cfg.t_init);
    PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/initial.vtk", 0.0);

    bool midExported = false;
    int iterSum = 0;
    int iterCount = 0;
    int maxIters = 0;

    FIM_Engine::TransientOptionalModules<MeshManager, FieldManager_2D> modules;
    modules.pressure_bc = &bcP;
    modules.temperature_bc = &bcT;
    modules.disable_default_vtk_output = true;

    FluidConstantProperties co2Props;
    co2Props.rho = cfg.co2_rho_const;
    co2Props.mu = cfg.co2_mu_const;
    co2Props.cp = cfg.co2_cp_const;
    co2Props.cv = cfg.co2_cv_const;
    co2Props.k = cfg.co2_k_const;
    modules.SetFluidModelConfig(FIM_Engine::UnifiedFluidModelConfig::MakeSinglePhaseCO2Constant(co2Props));

    modules.property_initializer = [&cfg](MeshManager&, FieldManager_2D& fld) {
        const auto rock = PhysicalProperties_string_op::Rock();
        const auto frac = PhysicalProperties_string_op::Fracture_string();
        const auto co2 = PhysicalProperties_string_op::CO2();

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

        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(co2.rho_tag, cfg.co2_rho_const), cfg.co2_rho_const);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(co2.mu_tag, cfg.co2_mu_const), cfg.co2_mu_const);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(co2.cp_tag, cfg.co2_cp_const), cfg.co2_cp_const);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(co2.cv_tag, cfg.co2_cv_const), cfg.co2_cv_const);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(co2.k_tag, cfg.co2_k_const), cfg.co2_k_const);

        ApplyUniformScalarField(fld.getOrCreateFractureScalar(co2.rho_tag, cfg.co2_rho_const), cfg.co2_rho_const);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(co2.mu_tag, cfg.co2_mu_const), cfg.co2_mu_const);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(co2.cp_tag, cfg.co2_cp_const), cfg.co2_cp_const);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(co2.cv_tag, cfg.co2_cv_const), cfg.co2_cv_const);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(co2.k_tag, cfg.co2_k_const), cfg.co2_k_const);
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
    FIM_Engine::RunGenericFIMTransient<2>(
        cfg.case_name,
        mgr,
        fm,
        ic,
        {},
        params,
        FIM_Engine::SolverRoute::FIM,
        modules);

    SyncPTFieldsToFM(mgr, fm, pBlocksLatest, tBlocksLatest, cfg.p_init, cfg.t_init);
    if (!midExported) {
        PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/mid.vtk", summary.t_end);
    }
    PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/final.vtk", summary.t_end);

    summary.max_iters = maxIters;
    summary.avg_iters = (iterCount > 0) ? (static_cast<double>(iterSum) / static_cast<double>(iterCount)) : 0.0;
    summary.h_char = ComputeMeshCharLength(mgr);

    std::vector<double> pFinal(pBlocksLatest.begin(), pBlocksLatest.begin() + static_cast<std::size_t>(nMat));
    std::vector<double> tFinal(tBlocksLatest.begin(), tBlocksLatest.begin() + static_cast<std::size_t>(nMat));

    const double xSplit = 0.5 * cfg.lx;
    const auto pLR = ComputeLeftRightAverage(mgr, pFinal, xSplit);
    const auto tLR = ComputeLeftRightAverage(mgr, tFinal, xSplit);
    if (pLR.valid) {
        summary.p_left_avg = pLR.left_avg;
        summary.p_right_avg = pLR.right_avg;
        summary.p_left_count = pLR.left_count;
        summary.p_right_count = pLR.right_count;
    }
    if (tLR.valid) {
        summary.t_left_avg = tLR.left_avg;
        summary.t_right_avg = tLR.right_avg;
        summary.t_left_count = tLR.left_count;
        summary.t_right_count = tLR.right_count;
    }
    summary.p_lr_pass = (pLR.valid && pLR.left_avg > pLR.right_avg) ? 1 : 0;
    summary.t_lr_pass = (tLR.valid && tLR.left_avg > tLR.right_avg) ? 1 : 0;

    summary.max_abs_dp = ComputeMaxAbsDelta(pFinal, cfg.p_init);
    summary.max_abs_dt = ComputeMaxAbsDelta(tFinal, cfg.t_init);
    summary.p_evolution_pass = (summary.max_abs_dp >= cfg.min_max_abs_dp_for_evolution) ? 1 : 0;
    summary.t_evolution_pass = (summary.max_abs_dt >= cfg.min_max_abs_dt_for_evolution) ? 1 : 0;

    const FractureChannelResult pFrac = ComputeFractureChannelResult(mgr, pBlocksLatest, cfg.p_init, cfg);
    if (!pFrac.valid) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] invalid fracture pressure metrics: empty fracture DOFs or near-fracture matrix band");
    }
    const FractureChannelResult tFrac = ComputeFractureChannelResult(mgr, tBlocksLatest, cfg.t_init, cfg);
    if (!tFrac.valid) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] invalid fracture temperature metrics: empty fracture DOFs or near-fracture matrix band");
    }

    summary.near_fracture_band_count = pFrac.near_matrix_count;
    summary.frac_mean_abs_dp = pFrac.fracture_mean_abs_delta;
    summary.near_band_mean_abs_dp = pFrac.near_matrix_mean_abs_delta;
    summary.frac_dp_advantage_pass = (pFrac.fracture_mean_abs_delta > pFrac.near_matrix_mean_abs_delta) ? 1 : 0;

    summary.frac_mean_abs_dt = tFrac.fracture_mean_abs_delta;
    summary.near_band_mean_abs_dt = tFrac.near_matrix_mean_abs_delta;
    summary.frac_dt_advantage_pass = (tFrac.fracture_mean_abs_delta > tFrac.near_matrix_mean_abs_delta) ? 1 : 0;

    summary.frac_p_upstream_avg = pFrac.upstream_avg;
    summary.frac_p_downstream_avg = pFrac.downstream_avg;
    summary.frac_p_upstream_count = pFrac.upstream_count;
    summary.frac_p_downstream_count = pFrac.downstream_count;
    summary.frac_p_along_pass = (pFrac.upstream_avg > pFrac.downstream_avg) ? 1 : 0;

    summary.frac_t_upstream_avg = tFrac.upstream_avg;
    summary.frac_t_downstream_avg = tFrac.downstream_avg;
    summary.frac_t_upstream_count = tFrac.upstream_count;
    summary.frac_t_downstream_count = tFrac.downstream_count;
    summary.frac_t_along_pass = (tFrac.upstream_avg > tFrac.downstream_avg) ? 1 : 0;

    summary.physics_checks_pass = (summary.p_lr_pass && summary.t_lr_pass &&
                                   summary.p_evolution_pass && summary.t_evolution_pass &&
                                   summary.frac_dp_advantage_pass && summary.frac_dt_advantage_pass &&
                                   summary.frac_p_along_pass && summary.frac_t_along_pass) ? 1 : 0;

    std::ofstream metrics(summary.metrics_csv_path, std::ios::out | std::ios::trunc);
    metrics << "case_name,nx,ny,n_cells,n_fracture_dofs,h_char,t_end,steps,total_rollbacks,"
               "avg_nonlinear_iters,max_nonlinear_iters,"
               "dt_init,dt_min,dt_max,property_mode,solver_route,"
               "p_left_avg,p_right_avg,p_left_count,p_right_count,p_lr_pass,"
               "t_left_avg,t_right_avg,t_left_count,t_right_count,t_lr_pass,"
               "max_abs_dp,max_abs_dt,dp_evolution_threshold,dt_evolution_threshold,"
               "p_evolution_pass,t_evolution_pass,"
               "near_fracture_band_count,"
               "frac_mean_abs_dp,near_band_mean_abs_dp,frac_dp_advantage_pass,"
               "frac_mean_abs_dt,near_band_mean_abs_dt,frac_dt_advantage_pass,"
               "frac_p_upstream_avg,frac_p_downstream_avg,frac_p_upstream_count,frac_p_downstream_count,frac_p_along_pass,"
               "frac_t_upstream_avg,frac_t_downstream_avg,frac_t_upstream_count,frac_t_downstream_count,frac_t_along_pass,"
               "physics_checks_pass\n";

    metrics << cfg.case_name << ","
            << cfg.nx << ","
            << cfg.ny << ","
            << nCells << ","
            << summary.n_fracture_dofs << ","
            << std::setprecision(12) << summary.h_char << ","
            << summary.t_end << ","
            << summary.steps << ","
            << summary.total_rollbacks << ","
            << summary.avg_iters << ","
            << summary.max_iters << ","
            << cfg.dt_init << ","
            << cfg.dt_min << ","
            << cfg.dt_max << ","
            << "ConstantSinglePhaseCO2" << ","
            << "RunGenericFIMTransient<2>" << ","
            << std::scientific << std::setprecision(8)
            << summary.p_left_avg << ","
            << summary.p_right_avg << ","
            << summary.p_left_count << ","
            << summary.p_right_count << ","
            << summary.p_lr_pass << ","
            << summary.t_left_avg << ","
            << summary.t_right_avg << ","
            << summary.t_left_count << ","
            << summary.t_right_count << ","
            << summary.t_lr_pass << ","
            << summary.max_abs_dp << ","
            << summary.max_abs_dt << ","
            << cfg.min_max_abs_dp_for_evolution << ","
            << cfg.min_max_abs_dt_for_evolution << ","
            << summary.p_evolution_pass << ","
            << summary.t_evolution_pass << ","
            << summary.near_fracture_band_count << ","
            << summary.frac_mean_abs_dp << ","
            << summary.near_band_mean_abs_dp << ","
            << summary.frac_dp_advantage_pass << ","
            << summary.frac_mean_abs_dt << ","
            << summary.near_band_mean_abs_dt << ","
            << summary.frac_dt_advantage_pass << ","
            << summary.frac_p_upstream_avg << ","
            << summary.frac_p_downstream_avg << ","
            << summary.frac_p_upstream_count << ","
            << summary.frac_p_downstream_count << ","
            << summary.frac_p_along_pass << ","
            << summary.frac_t_upstream_avg << ","
            << summary.frac_t_downstream_avg << ","
            << summary.frac_t_upstream_count << ","
            << summary.frac_t_downstream_count << ","
            << summary.frac_t_along_pass << ","
            << summary.physics_checks_pass << "\n";

    return summary;
}

void ExecutePlanByKey(const std::string& key) {
    const auto& registry = GetRegistry();
    auto it = registry.find(key);
    if (it == registry.end()) {
        throw std::runtime_error("[Test_H_T_CO2_ConstPP_SingleFrac] unknown registry key: " + key);
    }

    const TestCasePlan plan = it->second();
    const TestCaseSummary summary = RunCase(plan.spec);

    std::cout << "\n============================================\n";
    std::cout << "[Test_H_T_CO2_ConstPP_SingleFrac] run completed\n";
    std::cout << "  output_dir: " << summary.case_dir << "\n";
    std::cout << "  grid: " << summary.nx << " x " << summary.ny
              << " (" << summary.n_cells << " cells, fracture_dofs=" << summary.n_fracture_dofs << ")\n";
    std::cout << "  steps: " << summary.steps
              << "  rollbacks: " << summary.total_rollbacks << "\n";
    std::cout << "  Newton iters: avg=" << std::fixed << std::setprecision(2) << summary.avg_iters
              << "  max=" << summary.max_iters << "\n";
    std::cout << "  final_time: " << std::scientific << std::setprecision(4) << summary.t_end << " s\n";
    std::cout << "  checks: p_lr=" << summary.p_lr_pass
              << " t_lr=" << summary.t_lr_pass
              << " p_evo=" << summary.p_evolution_pass
              << " t_evo=" << summary.t_evolution_pass
              << " frac_dp=" << summary.frac_dp_advantage_pass
              << " frac_dt=" << summary.frac_dt_advantage_pass
              << " frac_p_along=" << summary.frac_p_along_pass
              << " frac_t_along=" << summary.frac_t_along_pass
              << " overall=" << summary.physics_checks_pass << "\n";
    std::cout << "============================================\n";
}

} // namespace

void RunTestCase() {
    ExecutePlanByKey("h_t_co2_constpp_singlefrac_nowell");
}

} // namespace Test_H_T_CO2_ConstPP_SingleFrac

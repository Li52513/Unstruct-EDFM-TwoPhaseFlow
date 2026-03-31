/**
 * @file  Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.cpp
 * @brief 
 */

#include "Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.h"

#include "2D_PostProcess.h"
#include "BoundaryConditionManager.h"
#include "Case2D_ReferenceIO.h"
#include "Case2D_Studies.h"
#include "Case2D_Validation.h"
#include "CaseCommon_Artifacts.h"
#include "CaseCommon_Catalog.h"
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

namespace Test_H_CO2_ConstPP {
namespace {

void EnsureDirRecursive(const std::string& rawPath) {
    CaseCommon::EnsureDirRecursive(rawPath);
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

using AnalyticalSnapshot = Case2DValidation::AnalyticalSnapshot;
using AnalyticalMetrics = Case2DValidation::AnalyticalMetrics;
using SweepStudyRow = Case2DStudies::SweepStudyRow;

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

Case2DValidation::PressureDiffusionAnalyticalConfig BuildAnalyticalConfig(const TestCaseSpec& cfg) {
    Case2DValidation::PressureDiffusionAnalyticalConfig out;
    out.lx = cfg.lx;
    out.p_init = cfg.p_init;
    out.p_left = cfg.p_left;
    out.p_right = cfg.p_right;
    out.permeability = cfg.matrix_perm;
    out.porosity = cfg.matrix_phi;
    out.total_compressibility = cfg.matrix_ct;
    out.viscosity = cfg.mu_const;
    out.analytical_terms = cfg.analytical_terms;
    return out;
}

double ComputePressureDiffusivity(const TestCaseSpec& cfg) {
    return Case2DValidation::ComputePressureDiffusivity(BuildAnalyticalConfig(cfg));
}

double EvaluateAnalyticalPressure(const TestCaseSpec& cfg, double x, double timeS) {
    return Case2DValidation::EvaluatePressureDiffusionAnalyticalPressure(BuildAnalyticalConfig(cfg), x, timeS);
}

std::vector<AnalyticalSnapshot> BuildSnapshots(const TestCaseSpec& cfg) {
    if (!cfg.enable_analytical_validation) return {};
    return Case2DValidation::BuildAnalyticalSnapshots(cfg.target_end_time_s, cfg.report_time_fractions);
}

AnalyticalMetrics EvaluateAnalyticalMetrics(const MeshManager& mgr, const std::vector<double>& pBlocks,
                                            const TestCaseSpec& cfg, const AnalyticalSnapshot& snapshot,
                                            const std::string& csvPath, bool writeCsv) {
    return Case2DValidation::EvaluatePressureDiffusionAnalyticalMetrics(
        mgr, pBlocks, BuildAnalyticalConfig(cfg), snapshot, csvPath, writeCsv);
}

void WriteProfileCSV(const MeshManager& mgr, const std::vector<double>& pBlocks,
                     const TestCaseSpec& cfg, const AnalyticalSnapshot& snapshot,
                     const std::string& csvPath) {
    Case2DValidation::WritePressureProfileCSV(
        mgr, pBlocks, cfg.nx, BuildAnalyticalConfig(cfg), snapshot, csvPath);
}

void WriteAnalyticalSummary(const TestCaseSpec& cfg, const TestCaseSummary& summary) {
    if (summary.case_dir.empty()) return;
    Case2DValidation::PressureDiffusionSummaryReport report;
    report.case_name = cfg.case_name;
    report.output_dir = summary.case_dir;
    report.nx = cfg.nx;
    report.ny = cfg.ny;
    report.n_cells = summary.n_cells;
    report.steps = summary.steps;
    report.total_rollbacks = summary.total_rollbacks;
    report.h_char = summary.h_char;
    report.t_end = summary.t_end;
    report.pressure_diffusivity = Case2DValidation::ComputePressureDiffusivity(BuildAnalyticalConfig(cfg));
    report.gravity_x = cfg.gravity_vector.m_x;
    report.gravity_y = cfg.gravity_vector.m_y;
    report.gravity_z = cfg.gravity_vector.m_z;
    report.non_orthogonal_correction = cfg.enable_non_orthogonal_correction;
    report.final_l1_norm = summary.final_l1_norm;
    report.final_l2_norm = summary.final_l2_norm;
    report.final_linf_norm = summary.final_linf_norm;
    report.l2_threshold = cfg.analytical_l2_threshold;
    report.linf_threshold = cfg.analytical_linf_threshold;
    report.grid_convergence_ok = summary.grid_convergence_ok;
    report.time_sensitivity_ok = summary.time_sensitivity_ok;
    report.analytical_passed = summary.analytical_passed;
    report.grid_convergence_csv_path = summary.grid_convergence_csv_path;
    report.time_sensitivity_csv_path = summary.time_sensitivity_csv_path;
    report.report_metrics = summary.report_metrics;
    Case2DValidation::WritePressureDiffusionSummaryReport(report, summary.case_dir + "/analytical_summary.txt");
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

const CaseCommon::CaseCatalogEntry& GetA1CatalogEntryOrThrow() {
    const CaseCommon::CaseCatalogEntry* entry = CaseCommon::FindCaseCatalogEntry("A1");
    if (!entry) throw std::runtime_error("[Test_H_CO2] failed to resolve A1 catalog entry.");
    return *entry;
}

CaseCommon::CaseArtifactPaths BuildA1ArtifactPaths() {
    const CaseCommon::CaseCatalogEntry& entry = GetA1CatalogEntryOrThrow();
    return CaseCommon::BuildArtifactPaths(
        entry.metadata.output_root,
        entry.metadata.case_code,
        entry.metadata.case_slug);
}

void EnsureA1ArtifactContractDirs(const CaseCommon::CaseArtifactPaths& artifacts) {
    EnsureDirRecursive(artifacts.case_dir);
    EnsureDirRecursive(artifacts.studies_dir);
    EnsureDirRecursive(artifacts.figures_dir);
    EnsureDirRecursive(artifacts.engineering_dir);
    EnsureDirRecursive(artifacts.reference_dir);
    EnsureDirRecursive(artifacts.report_dir);
    EnsureDirRecursive(artifacts.report_scripts_dir);
}

std::string JoinFractions(const std::vector<double>& values) {
    std::ostringstream oss;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) oss << ",";
        oss << std::setprecision(12) << values[i];
    }
    return oss.str();
}

void WriteA1StageManifest(const CaseCommon::CaseArtifactPaths& artifacts,
                          CaseCommon::CaseStage stage,
                          const std::string& status,
                          const std::string& primaryOutputDir) {
    Case2DReferenceIO::WriteAsciiFile(
        artifacts.engineering_stage_manifest_path,
        "[Test_H_CO2] failed to write A1 stage manifest",
        [&](std::ofstream& out) {
            out << "case_code=A1\n";
            out << "case_slug=" << GetA1CatalogEntryOrThrow().metadata.case_slug << "\n";
            out << "stage=" << CaseCommon::ToString(stage) << "\n";
            out << "status=" << status << "\n";
            out << "reference_mode=analytical\n";
            out << "primary_output_dir=" << primaryOutputDir << "\n";
        });
}

void WriteA1ReferenceContract(const TestCaseSpec& cfg,
                              const CaseCommon::CaseArtifactPaths& artifacts,
                              CaseCommon::CaseStage stage) {
    Case2DReferenceIO::WriteAsciiFile(
        artifacts.reference_contract_path,
        "[Test_H_CO2] failed to write A1 reference contract",
        [&](std::ofstream& out) {
            out << "reference_mode=analytical\n";
            out << "case_code=A1\n";
            out << "case_slug=" << GetA1CatalogEntryOrThrow().metadata.case_slug << "\n";
            out << "stage=" << CaseCommon::ToString(stage) << "\n";
            out << "family=pressure_diffusion\n";
            out << "variables=pressure\n";
            out << "report_time_fractions=" << JoinFractions(cfg.report_time_fractions) << "\n";
            out << "target_end_time_s=" << std::setprecision(12) << cfg.target_end_time_s << "\n";
            out << "analytical_terms=" << cfg.analytical_terms << "\n";
            out << "pressure_diffusivity=" << ComputePressureDiffusivity(cfg) << "\n";
            out << "prepare_reference_behavior=artifact_only\n";
            out << "validate_only_behavior=rerun_solver_until_engineering_snapshot_cache_exists\n";
        });
}

void WriteA1StageStatus(const TestCaseSpec& cfg,
                        const CaseCommon::CaseArtifactPaths& artifacts,
                        CaseCommon::CaseStage stage,
                        const std::string& status,
                        const std::string& primaryOutputDir) {
    Case2DReferenceIO::WriteAsciiFile(
        artifacts.report_status_markdown_path,
        "[Test_H_CO2] failed to write A1 status markdown",
        [&](std::ofstream& out) {
            out << "# A1 Thin Template Status\n\n";
            out << "- Stage: `" << CaseCommon::ToString(stage) << "`\n";
            out << "- Status: `" << status << "`\n";
            out << "- Reference mode: `analytical`\n";
            out << "- Primary output dir: `" << primaryOutputDir << "`\n";
            out << "- Case dir: `" << artifacts.case_dir << "`\n\n";
            out << "## Stage semantics\n";
            out << "- `solve_only`: run core solve without grid/time studies.\n";
            out << "- `prepare_reference`: materialize analytical-reference contract only.\n";
            out << "- `validate_only`: rerun validation workflow until engineering snapshot persistence exists.\n";
            out << "- `full_workflow`: run the full legacy analytical workflow under the template-system case root.\n\n";
            out << "## Numerical settings\n";
            out << "- Grid: `" << cfg.nx << "x" << cfg.ny << "`\n";
            out << "- End time: `" << std::setprecision(12) << cfg.target_end_time_s << " s`\n";
            out << "- Report fractions: `" << JoinFractions(cfg.report_time_fractions) << "`\n";
        });
    Case2DReferenceIO::WriteAsciiFile(
        artifacts.report_scripts_dir + "/README.m",
        "[Test_H_CO2] failed to write A1 Matlab stub",
        [&](std::ofstream& out) {
            out << "% A1 Matlab script placeholder\n";
            out << "% Stage: " << CaseCommon::ToString(stage) << "\n";
            out << "% Output root: " << artifacts.case_dir << "\n";
        });
}

TestCaseSpec BuildStageSpec(const TestCaseSpec& base, CaseCommon::CaseStage stage) {
    TestCaseSpec cfg = base;
    switch (stage) {
    case CaseCommon::CaseStage::SolveOnly:
        cfg.enable_grid_convergence_study = false;
        cfg.enable_time_sensitivity_study = false;
        cfg.enable_profile_output = false;
        cfg.emit_detailed_outputs = false;
        break;
    case CaseCommon::CaseStage::ValidateOnly:
        cfg.export_vtk = false;
        cfg.enable_grid_convergence_study = true;
        cfg.enable_time_sensitivity_study = true;
        cfg.enable_profile_output = true;
        cfg.emit_detailed_outputs = true;
        break;
    case CaseCommon::CaseStage::PrepareReference:
    case CaseCommon::CaseStage::FullWorkflow:
        break;
    default:
        throw std::runtime_error("[Test_H_CO2] unsupported stage in BuildStageSpec.");
    }
    return cfg;
}

void PrintRunSummary(const std::string& banner, const TestCaseSummary& summary) {
    std::cout << "\n============================================\n";
    std::cout << "[Test_H_CO2] " << banner << "\n";
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

void PrintPrepareReferenceSummary(const CaseCommon::CaseArtifactPaths& artifacts) {
    std::cout << "\n============================================\n";
    std::cout << "[Test_H_CO2] prepared analytical reference artifacts\n";
    std::cout << "  case_dir: " << artifacts.case_dir << "\n";
    std::cout << "  stage_manifest: " << artifacts.engineering_stage_manifest_path << "\n";
    std::cout << "  reference_contract: " << artifacts.reference_contract_path << "\n";
    std::cout << "  status_report: " << artifacts.report_status_markdown_path << "\n";
    std::cout << "============================================\n";
}

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

TestCaseSummary RunCase(const TestCaseSpec& cfg, const std::string& outputDirOverride) {
    TestCaseSummary summary = RunSingleCaseCore(cfg, outputDirOverride);
    std::vector<SweepStudyRow> gridRows, timeRows;
    const std::string studiesDir = summary.case_dir + "/studies";
    EnsureDirRecursive(studiesDir);

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
        summary.grid_convergence_ok = Case2DStudies::IsMonotonicNonIncreasingL2(gridRows);
        summary.grid_convergence_csv_path = studiesDir + "/grid_convergence.csv";
        Case2DStudies::WriteStudyCSV(gridRows, summary.grid_convergence_csv_path, "grid");
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
        summary.time_sensitivity_ok = Case2DStudies::IsMonotonicNonIncreasingL2(timeRows);
        summary.time_sensitivity_csv_path = studiesDir + "/time_sensitivity.csv";
        Case2DStudies::WriteStudyCSV(timeRows, summary.time_sensitivity_csv_path, "time");
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

void RunStageByKeyImpl(const std::string& key, CaseCommon::CaseStage stage) {
    const auto& registry = GetRegistry();
    const auto it = registry.find(key);
    if (it == registry.end()) throw std::runtime_error("[Test_H_CO2] unknown registry key: " + key);
    const TestCasePlan plan = it->second();
    const CaseCommon::CaseArtifactPaths artifacts = BuildA1ArtifactPaths();
    EnsureA1ArtifactContractDirs(artifacts);
    WriteA1StageManifest(artifacts, stage, "started", artifacts.case_dir);
    WriteA1ReferenceContract(plan.spec, artifacts, stage);
    WriteA1StageStatus(plan.spec, artifacts, stage, "started", artifacts.case_dir);

    if (stage == CaseCommon::CaseStage::PrepareReference) {
        WriteA1StageManifest(artifacts, stage, "prepared_reference_inputs", artifacts.reference_dir);
        WriteA1StageStatus(plan.spec, artifacts, stage, "prepared_reference_inputs", artifacts.reference_dir);
        PrintPrepareReferenceSummary(artifacts);
        return;
    }

    const TestCaseSpec stageSpec = BuildStageSpec(plan.spec, stage);
    TestCaseSummary summary;
    switch (stage) {
    case CaseCommon::CaseStage::SolveOnly:
        summary = RunSingleCaseCore(stageSpec, artifacts.engineering_dir);
        WriteA1StageManifest(artifacts, stage, "completed", artifacts.engineering_dir);
        WriteA1StageStatus(stageSpec, artifacts, stage, "completed", artifacts.engineering_dir);
        PrintRunSummary("solve_only completed", summary);
        return;
    case CaseCommon::CaseStage::ValidateOnly:
        summary = RunCase(stageSpec, artifacts.case_dir);
        WriteA1StageManifest(artifacts, stage, "completed", artifacts.case_dir);
        WriteA1StageStatus(stageSpec, artifacts, stage, "completed", artifacts.case_dir);
        PrintRunSummary("validate_only completed", summary);
        return;
    case CaseCommon::CaseStage::FullWorkflow:
        summary = RunCase(stageSpec, artifacts.case_dir);
        WriteA1StageManifest(artifacts, stage, "completed", artifacts.case_dir);
        WriteA1StageStatus(stageSpec, artifacts, stage, "completed", artifacts.case_dir);
        PrintRunSummary("full_workflow completed", summary);
        return;
    default:
        throw std::runtime_error("[Test_H_CO2] unsupported stage in RunStageByKeyImpl.");
    }
}

void ExecutePlanByKeyImpl(const std::string& key) {
    RunStageByKeyImpl(key, CaseCommon::CaseStage::FullWorkflow);
}

} // namespace

void RunTestCase() {
    ExecutePlanByKeyImpl("h_co2_constpp_nofrac_nowell");
}

void ExecutePlanByKey(const std::string& key) {
    ExecutePlanByKeyImpl(key);
}

void RunSolveOnly() {
    RunStageByKeyImpl("h_co2_constpp_nofrac_nowell", CaseCommon::CaseStage::SolveOnly);
}

void RunPrepareReference() {
    RunStageByKeyImpl("h_co2_constpp_nofrac_nowell", CaseCommon::CaseStage::PrepareReference);
}

void RunValidateOnly() {
    RunStageByKeyImpl("h_co2_constpp_nofrac_nowell", CaseCommon::CaseStage::ValidateOnly);
}

void RunFullWorkflow() {
    RunStageByKeyImpl("h_co2_constpp_nofrac_nowell", CaseCommon::CaseStage::FullWorkflow);
}

} // namespace Test_H_CO2_ConstPP

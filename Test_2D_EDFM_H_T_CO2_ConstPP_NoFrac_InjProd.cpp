#include "Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_InjProd.h"

#include "2D_PostProcess.h"
#include "BoundaryAssembler.h"
#include "BoundaryConditionManager.h"
#include "Case2D_ReferenceIO.h"
#include "CaseCommon_Artifacts.h"
#include "CaseCommon_Catalog.h"
#include "FIM_TransientCaseKit.hpp"
#include "FIM_TransientEngine/Types.hpp"
#include "MeshDefinitions.h"
#include "MeshManager.h"
#include "SolverContrlStrName_op.h"
#include "WellDOFManager.h"
#include "Well_WellControlTypes.h"

#include <algorithm>
#include <array>
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

namespace Test_H_T_CO2_ConstPP_NoFrac_InjProd {
namespace {

constexpr const char* kPlanKey = "h_t_co2_constpp_nofrac_injprod";
constexpr const char* kFamilyMatrixHorizontal = "matrix_horizontal";
constexpr const char* kFamilyMatrixVerticalMidline = "matrix_vertical_midline";
constexpr double kCellSearchTieTolerance = 1.0e-12;

struct TestCaseSpec {
    std::string case_name = kPlanKey;
    double lx = 400.0;
    double ly = 40.0;
    int nx = 48;
    int ny = 6;
    double matrix_phi = 0.10;
    double matrix_perm = 1.0e-13;
    double matrix_ct = 5.0e-9;
    double matrix_rho_r = 2600.0;
    double matrix_cp_r = 800.0;
    double matrix_lambda_r = 3.0;
    double co2_rho_const = 700.0;
    double co2_mu_const = 6.0e-5;
    double co2_cp_const = 1100.0;
    double co2_cv_const = 850.0;
    double co2_k_const = 0.03;
    double p_init = 10.0e6;
    double t_init = 380.0;
    double dt_init = 120.0;
    double dt_min = 1.0;
    double dt_max = 5000.0;
    double target_end_time_s = 1.0e5;
    int max_steps = 6000;
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
    double max_dT = 25.0;
    bool enable_armijo_line_search = false;
    double rollback_shrink_factor = 0.7;
    double dt_relres_grow_factor = 1.08;
    Vector gravity_vector = Vector(0.0, 0.0, 0.0);
    FIM_Engine::DiagLevel diag_level = FIM_Engine::DiagLevel::Off;
    bool export_vtk = true;
    CaseCommon::WellTemplateSpec well_template;
    double well_rw = 0.1;
    double well_skin = 0.0;
    std::vector<double> report_time_fractions = {0.1, 0.5, 1.0};
};

struct TestCasePlan {
    std::string plan_key;
    TestCaseSpec spec;
};

struct ResolvedWellSpec {
    std::string label;
    std::string role;
    int completion_id = -1;
    double target_x = 0.0;
    double target_y = 0.0;
    double actual_x = 0.0;
    double actual_y = 0.0;
    WellScheduleStep schedule;
};

struct WellTimeSeriesRow {
    int sample_id = -1;
    double time_s = 0.0;
    std::string well_name;
    std::string role;
    std::string control_mode;
    std::string component_mode;
    double target_value = 0.0;
    double actual_bhp_pa = std::numeric_limits<double>::quiet_NaN();
    double actual_rate_kg_per_s = std::numeric_limits<double>::quiet_NaN();
    double completion_pressure_pa = std::numeric_limits<double>::quiet_NaN();
    double completion_temperature_k = std::numeric_limits<double>::quiet_NaN();
    double reported_temperature_k = std::numeric_limits<double>::quiet_NaN();
    double target_x = 0.0;
    double target_y = 0.0;
    double actual_x = 0.0;
    double actual_y = 0.0;
};

struct B7OutputPaths {
    std::string property_table_path;
    std::string well_schedule_path;
    std::string well_timeseries_path;
    std::string reference_spec_path;
    std::string validation_summary_path;
    std::string matlab_script_path;
    std::string run_summary_path;
    std::string final_vtk_path;
};

struct CaseRunSummary {
    std::string output_dir;
    int n_cells = 0;
    int steps = 0;
    int total_rollbacks = 0;
    int max_iters = 0;
    double avg_iters = 0.0;
    double t_end = 0.0;
    double h_char = 0.0;
    std::string validation_status = "not_run";
    bool reference_ready = false;
    std::vector<std::string> missing_reference_files;
};

struct CaseRunArtifacts {
    CaseRunSummary summary;
    std::vector<WellTimeSeriesRow> well_rows;
};

std::string BoolString(bool value) {
    return value ? "true" : "false";
}

std::string JoinFractions(const std::vector<double>& values) {
    std::ostringstream oss;
    oss << std::setprecision(6);
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) oss << ", ";
        oss << values[i];
    }
    return oss.str();
}

std::string BuildSnapshotTag(double fraction) {
    const int pct = static_cast<int>(std::round(100.0 * fraction));
    std::ostringstream oss;
    oss << "t" << std::setfill('0') << std::setw(3) << pct << "pct";
    return oss.str();
}

void ApplyUniformScalarField(const std::shared_ptr<volScalarField>& field, double value) {
    if (!field) return;
    for (double& v : field->data) v = value;
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
    const int nUse = std::min(static_cast<int>(pBlocks.size()), mgr.getTotalDOFCount());
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

double ComputeMeshCharLength(const MeshManager& mgr) {
    const auto& cells = mgr.mesh().getCells();
    if (cells.empty()) return 0.0;
    double totalV = 0.0;
    for (const auto& c : cells) totalV += std::max(c.volume, 0.0);
    return (totalV > 0.0) ? std::sqrt(totalV / static_cast<double>(cells.size())) : 0.0;
}

const char* ToString(WellTargetDomain domain) {
    switch (domain) {
    case WellTargetDomain::Matrix: return "matrix";
    case WellTargetDomain::Fracture: return "fracture";
    default: return "unknown";
    }
}

const char* ToString(WellControlMode mode) {
    switch (mode) {
    case WellControlMode::BHP: return "bhp";
    case WellControlMode::Rate: return "rate";
    default: return "unknown";
    }
}

const char* ToString(WellComponentMode mode) {
    switch (mode) {
    case WellComponentMode::Water: return "water";
    case WellComponentMode::Gas: return "gas";
    case WellComponentMode::Total: return "total";
    default: return "unknown";
    }
}

const CaseCommon::CaseCatalogEntry& GetB7CatalogEntryOrThrow() {
    const CaseCommon::CaseCatalogEntry* entry = CaseCommon::FindCaseCatalogEntry("B7");
    if (!entry) throw std::runtime_error("[Test_H_T_CO2_InjProd] missing B7 catalog entry.");
    return *entry;
}

CaseCommon::CaseArtifactPaths BuildB7ArtifactPaths() {
    const CaseCommon::CaseCatalogEntry& entry = GetB7CatalogEntryOrThrow();
    return CaseCommon::BuildArtifactPaths(
        entry.metadata.output_root,
        entry.metadata.case_code,
        entry.metadata.case_slug);
}

void EnsureB7ArtifactContractDirs(const CaseCommon::CaseArtifactPaths& artifacts) {
    CaseCommon::EnsureDirRecursive(artifacts.root_dir);
    CaseCommon::EnsureDirRecursive(artifacts.case_dir);
    CaseCommon::EnsureDirRecursive(artifacts.studies_dir);
    CaseCommon::EnsureDirRecursive(artifacts.figures_dir);
    CaseCommon::EnsureDirRecursive(artifacts.engineering_dir);
    CaseCommon::EnsureDirRecursive(artifacts.reference_dir);
    CaseCommon::EnsureDirRecursive(artifacts.report_dir);
    CaseCommon::EnsureDirRecursive(artifacts.report_scripts_dir);
    CaseCommon::EnsureDirRecursive(artifacts.reference_dir + "/comsol");
}

B7OutputPaths BuildB7OutputPaths(const CaseCommon::CaseArtifactPaths& artifacts,
                                 const std::string& outputDir) {
    B7OutputPaths paths;
    paths.property_table_path = artifacts.engineering_dir + "/property_table.csv";
    paths.well_schedule_path = artifacts.engineering_dir + "/well_schedule.csv";
    paths.well_timeseries_path = artifacts.engineering_dir + "/well_timeseries.csv";
    paths.reference_spec_path = artifacts.engineering_dir + "/reference_spec.md";
    paths.validation_summary_path = artifacts.report_dir + "/validation_summary.md";
    paths.matlab_script_path = artifacts.report_scripts_dir + "/plot_validation_results.m";
    paths.run_summary_path = outputDir + "/run_summary.txt";
    paths.final_vtk_path = outputDir + "/final.vtk";
    return paths;
}

TestCaseSpec BuildBaseSpec() {
    TestCaseSpec cfg;
    cfg.well_template.injector_target_value = -1.0;
    cfg.well_template.producer_target_value = 9.5e6;
    cfg.well_template.injector_component_mode = WellComponentMode::Gas;
    cfg.well_template.producer_component_mode = WellComponentMode::Total;
    cfg.well_template.injection_fluid = "co2";
    cfg.well_template.thermal_policy = "cold_injection";
    cfg.well_template.injector_temperature = 320.0;
    return cfg;
}

TestCasePlan BuildPlanByKey(const std::string& key) {
    if (key != kPlanKey) {
        throw std::runtime_error("[Test_H_T_CO2_InjProd] unknown plan key: " + key);
    }
    TestCasePlan plan;
    plan.plan_key = key;
    plan.spec = BuildBaseSpec();
    return plan;
}

TestCaseSpec BuildStageSpec(const TestCaseSpec& base, CaseCommon::CaseStage stage) {
    TestCaseSpec cfg = base;
    if (stage == CaseCommon::CaseStage::PrepareReference || stage == CaseCommon::CaseStage::ValidateOnly) {
        cfg.export_vtk = false;
    }
    return cfg;
}

int FindNearestMatrixCell(const MeshManager& mgr, double targetX, double targetY) {
    const auto& cells = mgr.mesh().getCells();
    if (cells.empty()) {
        throw std::runtime_error("[Test_H_T_CO2_InjProd] mesh has no matrix cells.");
    }
    int bestIdx = -1;
    double bestDist2 = std::numeric_limits<double>::max();
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        const double dx = cells[static_cast<std::size_t>(i)].center.m_x - targetX;
        const double dy = cells[static_cast<std::size_t>(i)].center.m_y - targetY;
        const double dist2 = dx * dx + dy * dy;
        if (dist2 + kCellSearchTieTolerance < bestDist2) {
            bestDist2 = dist2;
            bestIdx = i;
        }
    }
    if (bestIdx < 0) {
        throw std::runtime_error("[Test_H_T_CO2_InjProd] failed to resolve nearest matrix completion.");
    }
    return bestIdx;
}

ResolvedWellSpec BuildResolvedWell(const TestCaseSpec& cfg,
                                   const MeshManager& mgr,
                                   const CaseCommon::WellTemplateLocation& location,
                                   const std::string& role,
                                   const std::string& wellName,
                                   WellControlMode controlMode,
                                   double targetValue,
                                   WellComponentMode componentMode,
                                   bool injectionIsCO2) {
    ResolvedWellSpec out;
    out.role = role;
    out.label = wellName;
    out.target_x = location.x_fraction * cfg.lx;
    out.target_y = location.y_fraction * cfg.ly;
    out.completion_id = FindNearestMatrixCell(mgr, out.target_x, out.target_y);
    const auto& cell = mgr.mesh().getCells()[static_cast<std::size_t>(out.completion_id)];
    out.actual_x = cell.center.m_x;
    out.actual_y = cell.center.m_y;

    out.schedule.t_start = 0.0;
    out.schedule.t_end = cfg.target_end_time_s;
    out.schedule.well_name = wellName;
    out.schedule.domain = location.domain;
    out.schedule.control_mode = controlMode;
    out.schedule.target_value = targetValue;
    out.schedule.component_mode = componentMode;
    out.schedule.rw = cfg.well_rw;
    out.schedule.skin = cfg.well_skin;
    out.schedule.well_axis = location.axis;
    out.schedule.completion_id = out.completion_id;
    out.schedule.wi_override = -1.0;
    out.schedule.L_override = -1.0;
    out.schedule.frac_w = injectionIsCO2 ? 0.0 : 1.0;
    out.schedule.frac_g = injectionIsCO2 ? 1.0 : 0.0;
    out.schedule.injection_is_co2 = injectionIsCO2;
    out.schedule.injection_temperature = injectionIsCO2 ? cfg.well_template.injector_temperature : -1.0;
    return out;
}

std::vector<ResolvedWellSpec> BuildWellSchedule(const TestCaseSpec& cfg, const MeshManager& mgr) {
    std::vector<ResolvedWellSpec> wells;
    wells.push_back(BuildResolvedWell(
        cfg,
        mgr,
        cfg.well_template.injector_location,
        "injector",
        cfg.well_template.injector_name,
        cfg.well_template.injector_control_mode,
        cfg.well_template.injector_target_value,
        cfg.well_template.injector_component_mode,
        true));
    wells.push_back(BuildResolvedWell(
        cfg,
        mgr,
        cfg.well_template.producer_location,
        "producer",
        cfg.well_template.producer_name,
        cfg.well_template.producer_control_mode,
        cfg.well_template.producer_target_value,
        cfg.well_template.producer_component_mode,
        false));
    return wells;
}

std::vector<std::string> OrderedProfileFamilies() {
    return {kFamilyMatrixHorizontal, kFamilyMatrixVerticalMidline};
}

std::vector<std::string> BuildExpectedReferenceFiles(const TestCaseSpec& cfg) {
    std::vector<std::string> files;
    for (double fraction : cfg.report_time_fractions) {
        const std::string tag = BuildSnapshotTag(fraction);
        for (const std::string& family : OrderedProfileFamilies()) {
            files.push_back("profile_" + family + "_" + tag + ".csv");
        }
    }
    files.push_back("monitor_timeseries.csv");
    files.push_back("well_timeseries.csv");
    return files;
}
void WriteB7StageManifest(const CaseCommon::CaseArtifactPaths& artifacts,
                          CaseCommon::CaseStage stage,
                          const std::string& status,
                          const std::string& outputDir) {
    Case2DReferenceIO::WriteAsciiFile(
        artifacts.engineering_stage_manifest_path,
        "[Test_H_T_CO2_InjProd] failed to write B7 stage manifest",
        [&](std::ofstream& out) {
            out << "case_code=B7\n";
            out << "stage=" << CaseCommon::ToString(stage) << "\n";
            out << "status=" << status << "\n";
            out << "output_dir=" << outputDir << "\n";
            out << "case_dir=" << artifacts.case_dir << "\n";
            out << "engineering_dir=" << artifacts.engineering_dir << "\n";
            out << "reference_dir=" << artifacts.reference_dir << "\n";
            out << "report_dir=" << artifacts.report_dir << "\n";
        });
}

void WriteB7ReferenceContract(const TestCaseSpec& cfg,
                              const CaseCommon::CaseArtifactPaths& artifacts,
                              CaseCommon::CaseStage stage) {
    Case2DReferenceIO::WriteAsciiFile(
        artifacts.reference_contract_path,
        "[Test_H_T_CO2_InjProd] failed to write B7 reference contract",
        [&](std::ofstream& out) {
            out << "# B7 Reference Contract\n\n";
            out << "- Stage: `" << CaseCommon::ToString(stage) << "`\n";
            out << "- Reference mode: `comsol`\n";
            out << "- Validation variables: `pressure`, `temperature`, `well_bhp`, `well_rate`, `production_temperature`\n";
            out << "- Default well policy: `injector rate + producer BHP`\n";
            out << "- Injection fluid: `co2`\n";
            out << "- Thermal policy: `cold_injection`\n";
            out << "- Injector temperature: `" << std::setprecision(12) << cfg.well_template.injector_temperature << " K`\n";
            out << "- Output end time: `" << std::setprecision(12) << cfg.target_end_time_s << " s`\n";
            out << "- Expected reference files:\n";
            for (const auto& file : BuildExpectedReferenceFiles(cfg)) {
                out << "  - `reference/comsol/" << file << "`\n";
            }
        });
}

void WriteB7StageStatus(const TestCaseSpec& cfg,
                        const CaseCommon::CaseArtifactPaths& artifacts,
                        CaseCommon::CaseStage stage,
                        const std::string& status,
                        const std::string& outputDir) {
    Case2DReferenceIO::WriteAsciiFile(
        artifacts.report_status_markdown_path,
        "[Test_H_T_CO2_InjProd] failed to write B7 stage status",
        [&](std::ofstream& out) {
            out << "# B7 Template Status\n\n";
            out << "- Stage: `" << CaseCommon::ToString(stage) << "`\n";
            out << "- Status: `" << status << "`\n";
            out << "- Output dir: `" << outputDir << "`\n";
            out << "- Grid: `" << cfg.nx << "x" << cfg.ny << "`\n";
            out << "- End time: `" << std::setprecision(12) << cfg.target_end_time_s << " s`\n";
            out << "- Well policy: `injector rate + producer BHP`\n";
            out << "- Injection fluid: `co2`\n";
            out << "- Thermal policy: `cold_injection`\n";
            out << "- Injector temperature: `" << cfg.well_template.injector_temperature << " K`\n\n";
            out << "## Stage semantics\n";
            out << "- `solve_only`: run the coupled N=2 injector-producer solve and write engineering outputs plus well time series.\n";
            out << "- `prepare_reference`: materialize weak-coupling reference inputs without running the solve.\n";
            out << "- `validate_only`: rerun the case under the case root and stop at `missing_reference` until COMSOL data exists.\n";
            out << "- `full_workflow`: run solve + reference check under the case root.\n";
        });
}

void WriteB7PropertyTable(const TestCaseSpec& cfg, const B7OutputPaths& paths) {
    Case2DReferenceIO::WriteAsciiFile(
        paths.property_table_path,
        "[Test_H_T_CO2_InjProd] failed to write B7 property table",
        [&](std::ofstream& out) {
            out << "key,value,unit\n";
            out << "lx," << std::setprecision(12) << cfg.lx << ",m\n";
            out << "ly," << cfg.ly << ",m\n";
            out << "nx," << cfg.nx << ",count\n";
            out << "ny," << cfg.ny << ",count\n";
            out << "matrix_perm," << cfg.matrix_perm << ",m2\n";
            out << "matrix_phi," << cfg.matrix_phi << ",1\n";
            out << "matrix_ct," << cfg.matrix_ct << ",1/Pa\n";
            out << "matrix_lambda_r," << cfg.matrix_lambda_r << ",W/m/K\n";
            out << "co2_rho_const," << cfg.co2_rho_const << ",kg/m3\n";
            out << "co2_mu_const," << cfg.co2_mu_const << ",Pa*s\n";
            out << "co2_cp_const," << cfg.co2_cp_const << ",J/kg/K\n";
            out << "co2_k_const," << cfg.co2_k_const << ",W/m/K\n";
            out << "p_init," << cfg.p_init << ",Pa\n";
            out << "t_init," << cfg.t_init << ",K\n";
            out << "injector_temperature," << cfg.well_template.injector_temperature << ",K\n";
            out << "dt_init," << cfg.dt_init << ",s\n";
            out << "target_end_time," << cfg.target_end_time_s << ",s\n";
        });
}

void WriteB7WellSchedule(const std::vector<ResolvedWellSpec>& wells, const B7OutputPaths& paths) {
    Case2DReferenceIO::WriteAsciiFile(
        paths.well_schedule_path,
        "[Test_H_T_CO2_InjProd] failed to write B7 well schedule",
        [&](std::ofstream& out) {
            out << "role,well_name,domain,control_mode,target_value,component_mode,completion_id,target_x_m,target_y_m,actual_x_m,actual_y_m,injection_is_co2,injection_temperature_k\n";
            for (const auto& well : wells) {
                out << well.role << "," << well.schedule.well_name << "," << ToString(well.schedule.domain) << ","
                    << ToString(well.schedule.control_mode) << "," << std::setprecision(12)
                    << well.schedule.target_value << "," << ToString(well.schedule.component_mode) << ","
                    << well.completion_id << "," << well.target_x << "," << well.target_y << ","
                    << well.actual_x << "," << well.actual_y << "," << BoolString(well.schedule.injection_is_co2) << ","
                    << well.schedule.injection_temperature << "\n";
            }
        });
}

void WriteB7ReferenceSpec(const TestCaseSpec& cfg,
                          const CaseCommon::CaseArtifactPaths& artifacts,
                          const B7OutputPaths& paths,
                          const std::vector<ResolvedWellSpec>& wells) {
    Case2DReferenceIO::WriteAsciiFile(
        paths.reference_spec_path,
        "[Test_H_T_CO2_InjProd] failed to write B7 reference spec",
        [&](std::ofstream& out) {
            out << "# B7 Reference Input Spec\n\n";
            out << "- Case dir: `" << artifacts.case_dir << "`\n";
            out << "- Grid: `" << cfg.nx << "x" << cfg.ny << "`\n";
            out << "- Domain: `" << cfg.lx << " m x " << cfg.ly << " m`\n";
            out << "- Variables: `pressure`, `temperature`, `well_bhp`, `well_rate`, `production_temperature`\n";
            out << "- Report fractions: `" << JoinFractions(cfg.report_time_fractions) << "`\n";
            out << "- Wells:\n";
            for (const auto& well : wells) {
                out << "  - `" << well.label << "` (" << well.role << "), cell `" << well.completion_id
                    << "`, actual location `" << std::setprecision(12) << well.actual_x << ", " << well.actual_y << "`\n";
            }
            out << "- Expected reference payloads:\n";
            for (const auto& file : BuildExpectedReferenceFiles(cfg)) {
                out << "  - `reference/comsol/" << file << "`\n";
            }
        });
}

void WriteB7MatlabStub(const B7OutputPaths& paths,
                       CaseCommon::CaseStage stage,
                       const CaseCommon::CaseArtifactPaths& artifacts) {
    Case2DReferenceIO::WriteAsciiFile(
        paths.matlab_script_path,
        "[Test_H_T_CO2_InjProd] failed to write B7 matlab stub",
        [&](std::ofstream& out) {
            out << "% B7 Matlab plotting placeholder\n";
            out << "% Stage: " << CaseCommon::ToString(stage) << "\n";
            out << "% Case root: " << artifacts.case_dir << "\n";
            out << "% Pending shared B/E well-aware plotting extraction.\n";
        });
}

void WriteRunSummary(const CaseRunSummary& summary, const B7OutputPaths& paths) {
    Case2DReferenceIO::WriteAsciiFile(
        paths.run_summary_path,
        "[Test_H_T_CO2_InjProd] failed to write B7 run summary",
        [&](std::ofstream& out) {
            out << "output_dir=" << summary.output_dir << "\n";
            out << "n_cells=" << summary.n_cells << "\n";
            out << "steps=" << summary.steps << "\n";
            out << "total_rollbacks=" << summary.total_rollbacks << "\n";
            out << "max_iters=" << summary.max_iters << "\n";
            out << std::setprecision(12);
            out << "avg_iters=" << summary.avg_iters << "\n";
            out << "t_end=" << summary.t_end << "\n";
            out << "h_char=" << summary.h_char << "\n";
            out << "validation_status=" << summary.validation_status << "\n";
            out << "reference_ready=" << BoolString(summary.reference_ready) << "\n";
        });
}

void WriteB7WellTimeSeries(const std::vector<WellTimeSeriesRow>& rows, const B7OutputPaths& paths) {
    Case2DReferenceIO::WriteAsciiFile(
        paths.well_timeseries_path,
        "[Test_H_T_CO2_InjProd] failed to write B7 well time series",
        [&](std::ofstream& out) {
            out << "sample_id,time_s,well_name,role,control_mode,component_mode,target_value,actual_bhp_pa,actual_rate_kg_per_s,completion_pressure_pa,completion_temperature_k,reported_temperature_k,target_x_m,target_y_m,actual_x_m,actual_y_m\n";
            out << std::setprecision(12);
            for (const auto& row : rows) {
                out << row.sample_id << "," << row.time_s << "," << row.well_name << "," << row.role << ","
                    << row.control_mode << "," << row.component_mode << "," << row.target_value << ","
                    << row.actual_bhp_pa << "," << row.actual_rate_kg_per_s << ","
                    << row.completion_pressure_pa << "," << row.completion_temperature_k << ","
                    << row.reported_temperature_k << "," << row.target_x << "," << row.target_y << ","
                    << row.actual_x << "," << row.actual_y << "\n";
            }
        });
}

void ValidateReferencePayloadOrThrow(const TestCaseSpec& cfg,
                                     const CaseCommon::CaseArtifactPaths& artifacts,
                                     const B7OutputPaths& paths,
                                     CaseRunSummary& summary) {
    std::vector<std::string> missingFiles;
    const std::string comsolDir = artifacts.reference_dir + "/comsol";
    for (const auto& file : BuildExpectedReferenceFiles(cfg)) {
        std::ifstream in((comsolDir + "/" + file).c_str(), std::ios::in);
        if (!in.good()) missingFiles.push_back("reference/comsol/" + file);
    }

    if (!missingFiles.empty()) {
        summary.reference_ready = false;
        summary.validation_status = "missing_reference";
        summary.missing_reference_files = missingFiles;
        Case2DReferenceIO::WriteAsciiFile(
            paths.validation_summary_path,
            "[Test_H_T_CO2_InjProd] failed to write B7 validation summary",
            [&](std::ofstream& out) {
                out << "# B7 Validation Summary\n\n";
                out << "- Status: `missing_reference`\n";
                out << "- Missing files:\n";
                for (const auto& file : missingFiles) out << "  - `" << file << "`\n";
            });
        throw std::runtime_error("[Test_H_T_CO2_InjProd] B7 reference files are missing.");
    }

    summary.reference_ready = true;
    summary.validation_status = "reference_payload_ready_metrics_pending";
    Case2DReferenceIO::WriteAsciiFile(
        paths.validation_summary_path,
        "[Test_H_T_CO2_InjProd] failed to write B7 validation summary",
        [&](std::ofstream& out) {
            out << "# B7 Validation Summary\n\n";
            out << "- Status: `reference_payload_ready_metrics_pending`\n";
            out << "- Note: quantitative B7 profile/monitor/well comparison remains the next B/E extraction target.\n";
        });
}

void MaterializeB7ReferenceInputs(const TestCaseSpec& cfg,
                                  const CaseCommon::CaseArtifactPaths& artifacts,
                                  CaseCommon::CaseStage stage) {
    MeshManager mgr(cfg.lx, cfg.ly, 0.0, cfg.nx, cfg.ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OverRelaxed);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(2);

    const std::vector<ResolvedWellSpec> wells = BuildWellSchedule(cfg, mgr);
    const B7OutputPaths paths = BuildB7OutputPaths(artifacts, artifacts.case_dir);
    WriteB7PropertyTable(cfg, paths);
    WriteB7WellSchedule(wells, paths);
    WriteB7ReferenceSpec(cfg, artifacts, paths, wells);
    WriteB7MatlabStub(paths, stage, artifacts);
}

FIM_Engine::TransientSolverParams BuildSolverParams(const TestCaseSpec& cfg) {
    FIM_Engine::TransientSolverParams params;
    params.max_steps = cfg.max_steps;
    params.dt_init = cfg.dt_init;
    params.dt_min = cfg.dt_min;
    params.dt_max = cfg.dt_max;
    params.target_end_time_s = cfg.target_end_time_s;
    params.max_newton_iter = cfg.max_newton_iter;
    params.abs_res_tol = cfg.abs_res_tol;
    params.rel_res_tol = cfg.rel_res_tol;
    params.rel_update_tol = cfg.rel_update_tol;
    params.lin_solver = cfg.lin_solver;
    params.amgcl_use_fallback_sparselu = cfg.amgcl_use_fallback_sparselu;
    params.amgcl_tol = cfg.amgcl_tol;
    params.amgcl_maxiter = cfg.amgcl_maxiter;
    params.enable_non_orthogonal_correction = cfg.enable_non_orthogonal_correction;
    params.enable_row_scaling = cfg.enable_row_scaling;
    params.max_dP = cfg.max_dP;
    params.max_dT = cfg.max_dT;
    params.max_dSw = 0.1;
    params.min_alpha = 1.0e-8;
    params.enable_armijo_line_search = cfg.enable_armijo_line_search;
    params.rollback_shrink_factor = cfg.rollback_shrink_factor;
    params.dt_relres_grow_factor = cfg.dt_relres_grow_factor;
    params.gravity_vector = cfg.gravity_vector;
    params.diag_level = cfg.diag_level;
    return params;
}

void AppendWellSnapshotRows(const TestCaseSpec& cfg,
                            double timeS,
                            int sampleId,
                            const std::vector<ResolvedWellSpec>& resolvedWells,
                            const std::vector<WellScheduleStep>& rawWells,
                            const FIM_Engine::WellDOFManager<2>& wellMgr,
                            MeshManager& mgr,
                            FieldManager_2D& fm,
                            const FluidPropertyEvalConfig& fluidCfg,
                            const CapRelPerm::VGParams& vg,
                            const CapRelPerm::RelPermParams& rp,
                            const std::vector<double>& pBlocks,
                            const std::vector<double>& tBlocks,
                            std::vector<WellTimeSeriesRow>& rows) {
    const std::vector<WellScheduleStep> activeWells =
        FIM_Engine::detail::SelectActiveAndNormalizeWells(rawWells, mgr, timeS, true);

    FIM_StateMap<2> state;
    state.InitSizes(static_cast<std::size_t>(wellMgr.TotalBlocksWithWells()));
    for (std::size_t i = 0; i < state.P.size(); ++i) {
        state.P[i] = (i < pBlocks.size()) ? pBlocks[i] : cfg.p_init;
        state.T[i] = (i < tBlocks.size()) ? tBlocks[i] : cfg.t_init;
    }

    const auto wellBhpByName = wellMgr.BuildWellBhpMap(state);
    const int totalEq = wellMgr.TotalBlocksWithWells() * 2;
    std::vector<double> residual(static_cast<std::size_t>(std::max(totalEq, 0)), 0.0);
    std::vector<std::array<double, 3>> jacobian(
        static_cast<std::size_t>(std::max(totalEq, 0)),
        std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<WellCompletionLinearization> completionLins;

    BoundaryAssembler::Assemble_Wells_2D_FullJac(
        mgr,
        fm,
        activeWells,
        0,
        -1,
        0,
        1,
        residual,
        jacobian,
        fluidCfg,
        vg,
        rp,
        &completionLins,
        &wellBhpByName);

    std::unordered_map<std::string, double> rateByName;
    rateByName.reserve(completionLins.size());
    for (const auto& lin : completionLins) {
        rateByName[lin.well_name] += lin.q_total;
    }

    for (const auto& well : resolvedWells) {
        WellTimeSeriesRow row;
        row.sample_id = sampleId;
        row.time_s = timeS;
        row.well_name = well.schedule.well_name;
        row.role = well.role;
        row.control_mode = ToString(well.schedule.control_mode);
        row.component_mode = ToString(well.schedule.component_mode);
        row.target_value = well.schedule.target_value;
        row.target_x = well.target_x;
        row.target_y = well.target_y;
        row.actual_x = well.actual_x;
        row.actual_y = well.actual_y;

        const int completionIdx =
            (well.schedule.completion_solver_index >= 0) ? well.schedule.completion_solver_index : well.completion_id;
        if (completionIdx >= 0 && completionIdx < static_cast<int>(pBlocks.size())) {
            row.completion_pressure_pa = pBlocks[static_cast<std::size_t>(completionIdx)];
        }
        if (completionIdx >= 0 && completionIdx < static_cast<int>(tBlocks.size())) {
            row.completion_temperature_k = tBlocks[static_cast<std::size_t>(completionIdx)];
        }
        row.reported_temperature_k =
            (well.schedule.injection_temperature > 0.0 && well.role == "injector")
                ? well.schedule.injection_temperature
                : row.completion_temperature_k;

        const auto bhpIt = wellBhpByName.find(well.schedule.well_name);
        if (bhpIt != wellBhpByName.end()) row.actual_bhp_pa = bhpIt->second;

        const auto rateIt = rateByName.find(well.schedule.well_name);
        if (rateIt != rateByName.end()) row.actual_rate_kg_per_s = rateIt->second;

        rows.push_back(row);
    }
}

CaseRunArtifacts RunSingleCaseCore(const TestCaseSpec& cfg, const std::string& outputDir) {
    CaseRunArtifacts artifacts;
    artifacts.summary.output_dir = outputDir;
    CaseCommon::EnsureDirRecursive(outputDir);

    MeshManager mgr(cfg.lx, cfg.ly, 0.0, cfg.nx, cfg.ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OverRelaxed);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);
    artifacts.summary.n_cells = static_cast<int>(mgr.mesh().getCells().size());
    if (artifacts.summary.n_cells <= 0) {
        throw std::runtime_error("[Test_H_T_CO2_InjProd] matrix cell count is zero.");
    }

    FIM_Engine::InitialConditions ic;
    ic.P_init = cfg.p_init;
    ic.T_init = cfg.t_init;
    ic.Sw_init = 1.0;

    BoundarySetting::BoundaryConditionManager bcP;
    bcP.Clear();
    bcP.SetNeumannBC(MeshTags::LEFT, 0.0);
    bcP.SetNeumannBC(MeshTags::RIGHT, 0.0);
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcP.SetNeumannBC(MeshTags::TOP, 0.0);

    BoundarySetting::BoundaryConditionManager bcT;
    bcT.Clear();
    bcT.SetNeumannBC(MeshTags::LEFT, 0.0);
    bcT.SetNeumannBC(MeshTags::RIGHT, 0.0);
    bcT.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcT.SetNeumannBC(MeshTags::TOP, 0.0);

    const auto pEqCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
    const auto tEqCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
    VTKBoundaryVisualizationContext bcVizCtx;
    bcVizCtx.water_family_policy = VTKBCWaterFamilyDerivePolicy::FollowPrimaryFluid;
    bcVizCtx.primary_fluid_model = VTKBCPrimaryFluidModel::CO2;
    bcVizCtx.bindings.push_back(VTKBCVariableBinding{pEqCfg.pressure_field, &bcP, VTKBCTransportKind::Pressure});
    bcVizCtx.bindings.push_back(VTKBCVariableBinding{tEqCfg.temperatue_field, &bcT, VTKBCTransportKind::Temperature});

    const std::vector<ResolvedWellSpec> resolvedWells = BuildWellSchedule(cfg, mgr);
    std::vector<WellScheduleStep> rawWells;
    rawWells.reserve(resolvedWells.size());
    for (const auto& well : resolvedWells) rawWells.push_back(well.schedule);

    const std::vector<WellScheduleStep> normalizedWellsSuperset =
        FIM_Engine::detail::SelectActiveAndNormalizeWells(rawWells, mgr, 0.0, false);
    FIM_Engine::WellDOFManager<2> wellMgr;
    wellMgr.Setup(normalizedWellsSuperset, mgr.getTotalDOFCount());

    FIM_StateMap<2> initialState;
    initialState.InitSizes(static_cast<std::size_t>(wellMgr.TotalBlocksWithWells()));
    for (std::size_t i = 0; i < initialState.P.size(); ++i) {
        initialState.P[i] = cfg.p_init;
        initialState.T[i] = cfg.t_init;
    }
    const std::vector<WellScheduleStep> initialActiveWells =
        FIM_Engine::detail::SelectActiveAndNormalizeWells(rawWells, mgr, 0.0, true);
    wellMgr.InitWellState(initialState, cfg.p_init, cfg.t_init, 1.0, initialActiveWells);

    std::vector<double> pBlocksLatest = initialState.P;
    std::vector<double> tBlocksLatest = initialState.T;
    SyncPTFieldsToFM(mgr, fm, pBlocksLatest, tBlocksLatest, cfg.p_init, cfg.t_init);
    if (cfg.export_vtk) {
        PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(outputDir + "/initial.vtk", 0.0);
    }

    int sampleId = 0;
    int iterSum = 0;
    int iterCount = 0;
    int maxIters = 0;

    FIM_Engine::TransientOptionalModules<MeshManager, FieldManager_2D> modules;
    modules.pressure_bc = &bcP;
    modules.temperature_bc = &bcT;
    modules.disable_default_vtk_output = true;
    modules.SetFluidModelConfig(FIM_Engine::UnifiedFluidModelConfig::MakeSinglePhaseCO2Constant(
        FluidConstantProperties{
            cfg.co2_rho_const,
            cfg.co2_mu_const,
            cfg.co2_cp_const,
            cfg.co2_cv_const,
            cfg.co2_k_const
        }));
    modules.property_initializer = [&cfg](MeshManager&, FieldManager_2D& fld) {
        const auto rock = PhysicalProperties_string_op::Rock();
        const auto frac = PhysicalProperties_string_op::Fracture_string();
        const auto gas = PhysicalProperties_string_op::CO2();
        const auto water = PhysicalProperties_string_op::Water();

        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.k_xx_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.k_yy_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.k_zz_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.phi_tag, cfg.matrix_phi), cfg.matrix_phi);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.c_r_tag, cfg.matrix_ct), cfg.matrix_ct);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.rho_tag, cfg.matrix_rho_r), cfg.matrix_rho_r);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.cp_tag, cfg.matrix_cp_r), cfg.matrix_cp_r);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.lambda_tag, cfg.matrix_lambda_r), cfg.matrix_lambda_r);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.k_t_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.k_n_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.phi_tag, cfg.matrix_phi), cfg.matrix_phi);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.c_r_tag, cfg.matrix_ct), cfg.matrix_ct);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(gas.k_tag, cfg.co2_k_const), cfg.co2_k_const);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(gas.k_tag, cfg.co2_k_const), cfg.co2_k_const);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(water.k_tag, cfg.co2_k_const), cfg.co2_k_const);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(water.k_tag, cfg.co2_k_const), cfg.co2_k_const);
    };

    const auto resolvedModules = FIM_Engine::ResolveTransientFluidModelConfig(modules);
    AppendWellSnapshotRows(
        cfg,
        0.0,
        sampleId++,
        resolvedWells,
        rawWells,
        wellMgr,
        mgr,
        fm,
        resolvedModules.fluid_property_eval,
        resolvedModules.vg_params,
        resolvedModules.rp_params,
        pBlocksLatest,
        tBlocksLatest,
        artifacts.well_rows);

    modules.on_step_accepted =
        [&](int step, double timeS, double, int newtonIters, double,
            int totalRollbacks, const std::string&,
            const std::vector<double>& pVec, const std::vector<double>& tVec,
            const std::vector<double>*) {
            if (step <= 0) return;
            if (!pVec.empty()) pBlocksLatest = pVec;
            if (!tVec.empty()) tBlocksLatest = tVec;
            SyncPTFieldsToFM(mgr, fm, pBlocksLatest, tBlocksLatest, cfg.p_init, cfg.t_init);

            artifacts.summary.steps = step;
            artifacts.summary.t_end = timeS;
            artifacts.summary.total_rollbacks = totalRollbacks;
            iterSum += newtonIters;
            iterCount += 1;
            maxIters = std::max(maxIters, newtonIters);

            AppendWellSnapshotRows(
                cfg,
                timeS,
                sampleId++,
                resolvedWells,
                rawWells,
                wellMgr,
                mgr,
                fm,
                resolvedModules.fluid_property_eval,
                resolvedModules.vg_params,
                resolvedModules.rp_params,
                pBlocksLatest,
                tBlocksLatest,
                artifacts.well_rows);
        };

    const auto params = BuildSolverParams(cfg);
    FIM_Engine::RunGenericFIMTransient<2>(
        cfg.case_name,
        mgr,
        fm,
        ic,
        rawWells,
        params,
        FIM_Engine::SolverRoute::FIM,
        modules);

    SyncPTFieldsToFM(mgr, fm, pBlocksLatest, tBlocksLatest, cfg.p_init, cfg.t_init);
    if (cfg.export_vtk) {
        PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(outputDir + "/final.vtk", artifacts.summary.t_end);
    }

    artifacts.summary.max_iters = maxIters;
    artifacts.summary.avg_iters = (iterCount > 0) ? (static_cast<double>(iterSum) / static_cast<double>(iterCount)) : 0.0;
    artifacts.summary.h_char = ComputeMeshCharLength(mgr);
    artifacts.summary.validation_status = "not_run";
    return artifacts;
}

void PrintPrepareReferenceSummary(const CaseCommon::CaseArtifactPaths& artifacts) {
    std::cout << "[B7] prepare_reference completed\n";
    std::cout << "  case_dir      : " << artifacts.case_dir << "\n";
    std::cout << "  engineering   : " << artifacts.engineering_dir << "\n";
    std::cout << "  reference_dir : " << artifacts.reference_dir << "\n";
}

void PrintRunSummary(const std::string& banner, const CaseRunSummary& summary) {
    std::cout << "[B7] " << banner << "\n";
    std::cout << "  output_dir : " << summary.output_dir << "\n";
    std::cout << "  n_cells    : " << summary.n_cells << "\n";
    std::cout << "  steps      : " << summary.steps << "\n";
    std::cout << "  t_end      : " << std::scientific << std::setprecision(8) << summary.t_end << " s\n";
    std::cout << "  h_char     : " << std::scientific << std::setprecision(8) << summary.h_char << " m\n";
    std::cout << "  status     : " << summary.validation_status << "\n";
}

void RunStageByKeyImpl(const std::string& key, CaseCommon::CaseStage stage) {
    const TestCasePlan plan = BuildPlanByKey(key);
    const CaseCommon::CaseArtifactPaths artifacts = BuildB7ArtifactPaths();
    EnsureB7ArtifactContractDirs(artifacts);
    WriteB7StageManifest(artifacts, stage, "started", artifacts.case_dir);
    WriteB7ReferenceContract(plan.spec, artifacts, stage);
    WriteB7StageStatus(plan.spec, artifacts, stage, "started", artifacts.case_dir);
    MaterializeB7ReferenceInputs(plan.spec, artifacts, stage);

    if (stage == CaseCommon::CaseStage::PrepareReference) {
        WriteB7StageManifest(artifacts, stage, "prepared_reference_inputs", artifacts.reference_dir);
        WriteB7StageStatus(plan.spec, artifacts, stage, "prepared_reference_inputs", artifacts.reference_dir);
        PrintPrepareReferenceSummary(artifacts);
        return;
    }

    const TestCaseSpec stageSpec = BuildStageSpec(plan.spec, stage);
    const std::string defaultOutputDir =
        (stage == CaseCommon::CaseStage::SolveOnly) ? artifacts.engineering_dir : artifacts.case_dir;

    try {
        CaseRunArtifacts runArtifacts;
        const B7OutputPaths paths = BuildB7OutputPaths(artifacts, defaultOutputDir);
        switch (stage) {
        case CaseCommon::CaseStage::SolveOnly:
            runArtifacts = RunSingleCaseCore(stageSpec, artifacts.engineering_dir);
            WriteB7WellTimeSeries(runArtifacts.well_rows, paths);
            WriteRunSummary(runArtifacts.summary, BuildB7OutputPaths(artifacts, artifacts.engineering_dir));
            WriteB7StageManifest(artifacts, stage, "completed", artifacts.engineering_dir);
            WriteB7StageStatus(stageSpec, artifacts, stage, "completed", artifacts.engineering_dir);
            PrintRunSummary("solve_only completed", runArtifacts.summary);
            return;
        case CaseCommon::CaseStage::ValidateOnly:
            runArtifacts = RunSingleCaseCore(stageSpec, artifacts.case_dir);
            WriteB7WellTimeSeries(runArtifacts.well_rows, paths);
            ValidateReferencePayloadOrThrow(stageSpec, artifacts, paths, runArtifacts.summary);
            WriteRunSummary(runArtifacts.summary, paths);
            WriteB7StageManifest(artifacts, stage, "completed", artifacts.case_dir);
            WriteB7StageStatus(stageSpec, artifacts, stage, "completed", artifacts.case_dir);
            PrintRunSummary("validate_only completed", runArtifacts.summary);
            return;
        case CaseCommon::CaseStage::FullWorkflow:
            runArtifacts = RunSingleCaseCore(stageSpec, artifacts.case_dir);
            WriteB7WellTimeSeries(runArtifacts.well_rows, paths);
            ValidateReferencePayloadOrThrow(stageSpec, artifacts, paths, runArtifacts.summary);
            WriteRunSummary(runArtifacts.summary, paths);
            WriteB7StageManifest(artifacts, stage, "completed", artifacts.case_dir);
            WriteB7StageStatus(stageSpec, artifacts, stage, "completed", artifacts.case_dir);
            PrintRunSummary("full_workflow completed", runArtifacts.summary);
            return;
        default:
            throw std::runtime_error("[Test_H_T_CO2_InjProd] unsupported stage in RunStageByKeyImpl.");
        }
    } catch (const std::exception& ex) {
        const std::string message = ex.what();
        const bool missingReference =
            message.find("missing") != std::string::npos &&
            message.find("reference") != std::string::npos;
        const std::string failureStatus = missingReference ? "missing_reference" : "failed";
        WriteB7StageManifest(artifacts, stage, failureStatus, defaultOutputDir);
        WriteB7StageStatus(stageSpec, artifacts, stage, failureStatus, defaultOutputDir);
        throw;
    }
}

void ExecutePlanByKeyImpl(const std::string& key) {
    RunStageByKeyImpl(key, CaseCommon::CaseStage::FullWorkflow);
}

} // namespace

void RunTestCase() {
    ExecutePlanByKeyImpl(kPlanKey);
}

void RunSolveOnly() {
    RunStageByKeyImpl(kPlanKey, CaseCommon::CaseStage::SolveOnly);
}

void RunPrepareReference() {
    RunStageByKeyImpl(kPlanKey, CaseCommon::CaseStage::PrepareReference);
}

void RunValidateOnly() {
    RunStageByKeyImpl(kPlanKey, CaseCommon::CaseStage::ValidateOnly);
}

void RunFullWorkflow() {
    RunStageByKeyImpl(kPlanKey, CaseCommon::CaseStage::FullWorkflow);
}

} // namespace Test_H_T_CO2_ConstPP_NoFrac_InjProd

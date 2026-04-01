#include "Test_2D_EDFM_H_TP_CO2H2O_ConstPP_NoFrac_NoWell.h"

#include "2D_PostProcess.h"
#include "Case2D_ReferenceIO.h"
#include "CaseCommon_Artifacts.h"
#include "CaseCommon_Catalog.h"
#include "FIM_TransientCaseKit.hpp"
#include "FVM_Grad.h"
#include "MeshManager.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Test_H_TP_CO2H2O_ConstPP_NoFrac {
namespace {

constexpr double kPendingTimeTolerance = 1.0e-9;
constexpr const char* kCaseCodeC1 = "C1";
constexpr const char* kCaseCodeC2 = "C2";
constexpr const char* kCaseCodeC3 = "C3";
constexpr const char* kCaseCodeC4 = "C4";
constexpr const char* kCaseCodeC5 = "C5";
constexpr const char* kCaseCodeC6 = "C6";
constexpr const char* kPlanKeyC1 = "h_tp_co2h2o_constpp_nofrac_nowell";
constexpr const char* kPlanKeyC2 = "h_tp_co2h2o_varypp_nofrac_nowell";
constexpr const char* kPlanKeyC2Smoke = "h_tp_co2h2o_varypp_nofrac_nowell_smoke";
constexpr const char* kPlanKeyC3 = "h_tp_co2h2o_constpp_singlefrac_nowell";
constexpr const char* kPlanKeyC3Smoke = "h_tp_co2h2o_constpp_singlefrac_nowell_smoke";
constexpr const char* kPlanKeyC4 = "h_tp_co2h2o_varypp_singlefrac_nowell";
constexpr const char* kPlanKeyC4Smoke = "h_tp_co2h2o_varypp_singlefrac_nowell_smoke";
constexpr const char* kPlanKeyC5 = "h_tp_co2h2o_constpp_complexfrac_nowell";
constexpr const char* kPlanKeyC5Smoke = "h_tp_co2h2o_constpp_complexfrac_nowell_smoke";
constexpr const char* kPlanKeyC6 = "h_tp_co2h2o_varypp_complexfrac_nowell";
constexpr const char* kPlanKeyC6Smoke = "h_tp_co2h2o_varypp_complexfrac_nowell_smoke";
constexpr const char* kFamilyMatrixHorizontal = "matrix_horizontal";
constexpr const char* kFamilyMatrixVerticalMidline = "matrix_vertical_midline";
constexpr const char* kLocationMatrix = "matrix";

enum class FractureVariant { None, Single, Complex };

struct SpatialSamplePoint {
    int id = -1;
    std::string label;
    std::string family;
    std::string location;
    double target_axis_m = 0.0;
    double target_x = 0.0;
    double target_y = 0.0;
    int block_id = -1;
    double actual_x = 0.0;
    double actual_y = 0.0;
};

struct SnapshotState {
    std::string tag;
    double requested_fraction = 0.0;
    double target_time_s = 0.0;
    double actual_time_s = 0.0;
    bool captured = false;
    std::vector<double> p_blocks;
    std::vector<double> t_blocks;
    std::vector<double> sw_blocks;
};

struct MonitorScheduleState {
    int sample_id = -1;
    double target_time_s = 0.0;
    double actual_time_s = 0.0;
    bool captured = false;
    std::vector<double> p_blocks;
    std::vector<double> t_blocks;
    std::vector<double> sw_blocks;
};

struct TestCaseSpec {
    std::string case_code = kCaseCodeC1;
    std::string case_name = kPlanKeyC1;
    bool use_variable_properties = false;
    FractureVariant fracture_variant = FractureVariant::None;
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
    double matrix_cp_r = 800.0;
    double matrix_lambda_r = 3.0;
    double fracture_phi = 0.15;
    double fracture_kt = 5.0e-12;
    double fracture_kn = 5.0e-13;
    double fracture_ct = 5.0e-9;
    double fracture_rho_r = 2600.0;
    double fracture_cp_r = 800.0;
    double fracture_lambda_r = 6.0;
    FluidConstantProperties water_props = FluidConstantProperties{1000.0, 1.0e-3, 4200.0, 4180.0, 0.6};
    FluidConstantProperties gas_props = FluidConstantProperties{700.0, 6.0e-5, 1100.0, 850.0, 0.03};
    double p_init = 10.0e6;
    double p_left = 12.0e6;
    double p_right = 8.0e6;
    double t_init = 380.0;
    double t_left = 430.0;
    double t_right = 330.0;
    double sw_init = 0.90;
    double sw_left = 0.25;
    double sw_right = 0.95;
    double dt_init = 1.0e3;
    double dt_min = 10.0;
    double dt_max = 2.0e4;
    double target_end_time_s = 2.0e5;
    int max_steps = 600;
    int max_newton_iter = 12;
    FIM_Engine::LinearSolverType lin_solver = FIM_Engine::LinearSolverType::SparseLU;
    bool amgcl_use_fallback_sparselu = true;
    double amgcl_tol = 1.0e-6;
    int amgcl_maxiter = 300;
    bool enable_non_orthogonal_correction = true;
    bool enable_row_scaling = true;
    double abs_res_tol = 1.0e-7;
    double rel_res_tol = 1.0e-5;
    double rel_update_tol = 1.0e-7;
    double max_dP = 2.0e6;
    double max_dT = 10.0;
    double max_dSw = 0.05;
    bool export_vtk = true;
    bool write_engineering_profiles = true;
    std::vector<double> report_time_fractions = {0.1, 0.5, 1.0};
    std::vector<double> monitor_time_fractions = {0.0, 0.25, 0.5, 0.75, 1.0};
    int matrix_horizontal_station_count = 81;
    int matrix_vertical_station_count = 25;
};

struct TestCasePlan {
    std::string plan_key;
    TestCaseSpec spec;
};

struct CaseOutputPaths {
    std::string output_dir;
    std::string property_table_path;
    std::string reference_spec_path;
    std::string validation_summary_path;
    std::string matlab_script_path;
    std::string profile_station_definitions_path;
    std::string monitor_point_definitions_path;
    std::string profile_schedule_path;
    std::string monitor_schedule_path;
    std::string engineering_monitor_path;
    std::string engineering_profile_prefix;
    std::string run_summary_path;
    std::string final_vtk_path;
};

struct CaseRunSummary {
    std::string case_code;
    std::string output_dir;
    int n_cells = 0;
    int steps = 0;
    double t_end = 0.0;
    double h_char = 0.0;
    std::string validation_status = "not_run";
    bool reference_ready = false;
    std::vector<std::string> missing_reference_files;
};

struct CaseRunArtifacts {
    CaseRunSummary summary;
    std::vector<SpatialSamplePoint> profile_stations;
    std::vector<SpatialSamplePoint> monitor_points;
    std::vector<SnapshotState> snapshots;
    std::vector<MonitorScheduleState> monitor_schedule;
};

std::string BoolString(bool value) {
    return value ? "true" : "false";
}

std::string FractureVariantString(FractureVariant value) {
    switch (value) {
    case FractureVariant::None: return "nofrac";
    case FractureVariant::Single: return "singlefrac";
    case FractureVariant::Complex: return "complexfrac";
    default: return "unknown";
    }
}

std::string PropertyModeString(const TestCaseSpec& cfg) {
    return cfg.use_variable_properties ? "varypp" : "constpp";
}

std::string CaseBanner(const std::string& caseCode) {
    return "[" + caseCode + "]";
}

std::string FailurePrefix(const std::string& caseCode) {
    return "[Test_H_TP_CO2H2O_" + caseCode + "]";
}

const CaseCommon::CaseCatalogEntry& GetCatalogEntryOrThrow(const std::string& caseCode) {
    const CaseCommon::CaseCatalogEntry* entry = CaseCommon::FindCaseCatalogEntry(caseCode);
    if (!entry) throw std::runtime_error(FailurePrefix(caseCode) + " missing catalog entry.");
    return *entry;
}

CaseCommon::CaseArtifactPaths BuildArtifactPathsForCase(const std::string& caseCode) {
    const CaseCommon::CaseCatalogEntry& entry = GetCatalogEntryOrThrow(caseCode);
    return CaseCommon::BuildArtifactPaths(entry.metadata.output_root, entry.metadata.case_code, entry.metadata.case_slug);
}

void EnsureArtifactContractDirs(const CaseCommon::CaseArtifactPaths& artifacts) {
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

CaseOutputPaths BuildOutputPaths(const CaseCommon::CaseArtifactPaths& artifacts,
                                 const std::string& outputDir) {
    CaseOutputPaths paths;
    paths.output_dir = outputDir;
    paths.property_table_path = artifacts.engineering_dir + "/property_table.csv";
    paths.reference_spec_path = artifacts.engineering_dir + "/reference_spec.md";
    paths.validation_summary_path = artifacts.report_dir + "/validation_summary.md";
    paths.matlab_script_path = artifacts.report_scripts_dir + "/plot_validation_results.m";
    paths.profile_station_definitions_path = artifacts.engineering_dir + "/profile_station_definitions.csv";
    paths.monitor_point_definitions_path = artifacts.engineering_dir + "/monitor_point_definitions.csv";
    paths.profile_schedule_path = artifacts.engineering_dir + "/profile_report_schedule.csv";
    paths.monitor_schedule_path = artifacts.engineering_dir + "/monitor_sample_schedule.csv";
    paths.engineering_monitor_path = artifacts.engineering_dir + "/eng_monitor_timeseries.csv";
    paths.engineering_profile_prefix = artifacts.engineering_dir + "/eng_profile_";
    paths.run_summary_path = outputDir + "/run_summary.txt";
    paths.final_vtk_path = outputDir + "/final.vtk";
    return paths;
}

TestCaseSpec BuildBaseSpec() {
    return TestCaseSpec();
}

TestCasePlan BuildVariantPlan(const std::string& planKey,
                              const std::string& caseCode,
                              const std::string& caseName,
                              bool useVariableProperties,
                              FractureVariant fractureVariant) {
    TestCasePlan plan;
    plan.plan_key = planKey;
    plan.spec = BuildBaseSpec();
    plan.spec.case_code = caseCode;
    plan.spec.case_name = caseName;
    plan.spec.use_variable_properties = useVariableProperties;
    plan.spec.fracture_variant = fractureVariant;
    return plan;
}

TestCasePlan BuildC1Plan() {
    return BuildVariantPlan(kPlanKeyC1, kCaseCodeC1, kPlanKeyC1, false, FractureVariant::None);
}

TestCasePlan BuildC2Plan() {
    return BuildVariantPlan(kPlanKeyC2, kCaseCodeC2, "c2_varypp_nofrac", true, FractureVariant::None);
}

TestCasePlan BuildC2SmokePlan() {
    TestCasePlan plan = BuildC2Plan();
    plan.plan_key = kPlanKeyC2Smoke;
    plan.spec.case_name = "c2_smoke";
    plan.spec.nx = 24;
    plan.spec.ny = 3;
    plan.spec.dt_init = 5.0e2;
    plan.spec.dt_max = 5.0e3;
    plan.spec.target_end_time_s = 2.0e4;
    plan.spec.max_steps = 300;
    return plan;
}

TestCasePlan BuildC3Plan() {
    return BuildVariantPlan(kPlanKeyC3, kCaseCodeC3, "c3_constpp_singlefrac", false, FractureVariant::Single);
}

TestCasePlan BuildC3SmokePlan() {
    TestCasePlan plan = BuildC3Plan();
    plan.plan_key = kPlanKeyC3Smoke;
    plan.spec.case_name = "c3_smoke";
    plan.spec.nx = 24;
    plan.spec.ny = 3;
    plan.spec.dt_init = 5.0e2;
    plan.spec.dt_max = 5.0e3;
    plan.spec.target_end_time_s = 2.0e4;
    plan.spec.max_steps = 300;
    return plan;
}

TestCasePlan BuildC4Plan() {
    return BuildVariantPlan(kPlanKeyC4, kCaseCodeC4, "c4_varypp_singlefrac", true, FractureVariant::Single);
}

TestCasePlan BuildC4SmokePlan() {
    TestCasePlan plan = BuildC4Plan();
    plan.plan_key = kPlanKeyC4Smoke;
    plan.spec.case_name = "c4_smoke";
    plan.spec.nx = 16;
    plan.spec.ny = 3;
    plan.spec.dt_init = 1.0e2;
    plan.spec.dt_max = 1.0e3;
    plan.spec.target_end_time_s = 5.0e3;
    plan.spec.max_steps = 160;
    return plan;
}

TestCasePlan BuildC5Plan() {
    return BuildVariantPlan(kPlanKeyC5, kCaseCodeC5, "c5_constpp_complexfrac", false, FractureVariant::Complex);
}

TestCasePlan BuildC5SmokePlan() {
    TestCasePlan plan = BuildC5Plan();
    plan.plan_key = kPlanKeyC5Smoke;
    plan.spec.case_name = "c5_smoke";
    plan.spec.nx = 16;
    plan.spec.ny = 3;
    plan.spec.dt_init = 1.0e2;
    plan.spec.dt_max = 1.0e3;
    plan.spec.target_end_time_s = 5.0e3;
    plan.spec.max_steps = 160;
    return plan;
}

TestCasePlan BuildC6Plan() {
    return BuildVariantPlan(kPlanKeyC6, kCaseCodeC6, "c6_varypp_complexfrac", true, FractureVariant::Complex);
}

TestCasePlan BuildC6SmokePlan() {
    TestCasePlan plan = BuildC6Plan();
    plan.plan_key = kPlanKeyC6Smoke;
    plan.spec.case_name = "c6_smoke";
    plan.spec.nx = 8;
    plan.spec.ny = 2;
    plan.spec.dt_init = 5.0e1;
    plan.spec.dt_max = 2.0e2;
    plan.spec.target_end_time_s = 1.0e3;
    plan.spec.max_steps = 120;
    return plan;
}

TestCasePlan BuildPlanByKey(const std::string& key) {
    using BuilderFn = TestCasePlan(*)();
    static const std::unordered_map<std::string, BuilderFn> registry = {
        {kPlanKeyC1, &BuildC1Plan},
        {kPlanKeyC2, &BuildC2Plan},
        {kPlanKeyC2Smoke, &BuildC2SmokePlan},
        {kPlanKeyC3, &BuildC3Plan},
        {kPlanKeyC3Smoke, &BuildC3SmokePlan},
        {kPlanKeyC4, &BuildC4Plan},
        {kPlanKeyC4Smoke, &BuildC4SmokePlan},
        {kPlanKeyC5, &BuildC5Plan},
        {kPlanKeyC5Smoke, &BuildC5SmokePlan},
        {kPlanKeyC6, &BuildC6Plan},
        {kPlanKeyC6Smoke, &BuildC6SmokePlan}
    };
    const auto it = registry.find(key);
    if (it == registry.end()) {
        throw std::runtime_error("[Test_H_TP_CO2H2O] unknown plan key: " + key);
    }
    return it->second();
}

TestCaseSpec BuildStageSpec(const TestCaseSpec& base, CaseCommon::CaseStage stage) {
    TestCaseSpec cfg = base;
    if (stage == CaseCommon::CaseStage::PrepareReference) {
        cfg.export_vtk = false;
        cfg.write_engineering_profiles = false;
    } else if (stage == CaseCommon::CaseStage::ValidateOnly) {
        cfg.export_vtk = false;
    }
    return cfg;
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

template <typename SampleType>
void CaptureIfReached(double timeS,
                      const std::vector<double>& p,
                      const std::vector<double>& t,
                      const std::vector<double>& sw,
                      SampleType& sample) {
    if (sample.captured) return;
    if (timeS + kPendingTimeTolerance < sample.target_time_s) return;
    sample.captured = true;
    sample.actual_time_s = timeS;
    sample.p_blocks = p;
    sample.t_blocks = t;
    sample.sw_blocks = sw;
}

template <typename SampleContainer>
void FinalizeMissingSamples(double timeS,
                            const std::vector<double>& p,
                            const std::vector<double>& t,
                            const std::vector<double>& sw,
                            SampleContainer& samples) {
    for (auto& sample : samples) {
        if (sample.captured) continue;
        sample.captured = true;
        sample.actual_time_s = timeS;
        sample.p_blocks = p;
        sample.t_blocks = t;
        sample.sw_blocks = sw;
    }
}

std::vector<SnapshotState> BuildSnapshots(const TestCaseSpec& cfg) {
    std::vector<SnapshotState> snapshots;
    snapshots.reserve(cfg.report_time_fractions.size());
    for (double fraction : cfg.report_time_fractions) {
        SnapshotState snapshot;
        snapshot.tag = BuildSnapshotTag(fraction);
        snapshot.requested_fraction = fraction;
        snapshot.target_time_s = fraction * cfg.target_end_time_s;
        snapshots.push_back(snapshot);
    }
    return snapshots;
}

std::vector<MonitorScheduleState> BuildMonitorSchedule(const TestCaseSpec& cfg) {
    std::vector<MonitorScheduleState> schedule;
    schedule.reserve(cfg.monitor_time_fractions.size());
    for (std::size_t i = 0; i < cfg.monitor_time_fractions.size(); ++i) {
        MonitorScheduleState sample;
        sample.sample_id = static_cast<int>(i);
        sample.target_time_s = cfg.monitor_time_fractions[i] * cfg.target_end_time_s;
        schedule.push_back(sample);
    }
    return schedule;
}

Vector BuildPrimaryFractureStart(const TestCaseSpec& cfg) {
    return Vector(cfg.frac_x0_ratio * cfg.lx, cfg.frac_y0_ratio * cfg.ly, 0.0);
}

Vector BuildPrimaryFractureEnd(const TestCaseSpec& cfg) {
    return Vector(cfg.frac_x1_ratio * cfg.lx, cfg.frac_y1_ratio * cfg.ly, 0.0);
}

void AddConfiguredFractures(MeshManager& mgr, const TestCaseSpec& cfg) {
    if (cfg.fracture_variant == FractureVariant::None) return;

    mgr.addFracture(BuildPrimaryFractureStart(cfg), BuildPrimaryFractureEnd(cfg));
    if (cfg.fracture_variant != FractureVariant::Complex) return;

    mgr.addFracture(
        Vector(0.18 * cfg.lx, 0.80 * cfg.ly, 0.0),
        Vector(0.52 * cfg.lx, 0.38 * cfg.ly, 0.0));
    mgr.addFracture(
        Vector(0.58 * cfg.lx, 0.82 * cfg.ly, 0.0),
        Vector(0.88 * cfg.lx, 0.58 * cfg.ly, 0.0));
}

void ConfigureMeshForCase(MeshManager& mgr, const TestCaseSpec& cfg) {
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OverRelaxed);
    if (cfg.fracture_variant != FractureVariant::None) {
        AddConfiguredFractures(mgr, cfg);
        mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    }
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(3);
}

double DistanceSq(const Vector& a, double x, double y) {
    const double dx = a.m_x - x;
    const double dy = a.m_y - y;
    return dx * dx + dy * dy;
}

SpatialSamplePoint BuildNearestSamplePoint(const MeshManager& mgr,
                                          int pointId,
                                          const std::string& label,
                                          const std::string& family,
                                          const std::string& location,
                                          double targetAxis,
                                          double targetX,
                                          double targetY) {
    SpatialSamplePoint point;
    point.id = pointId;
    point.label = label;
    point.family = family;
    point.location = location;
    point.target_axis_m = targetAxis;
    point.target_x = targetX;
    point.target_y = targetY;

    const auto& cells = mgr.mesh().getCells();
    if (cells.empty()) {
        throw std::runtime_error("[Test_H_TP_CO2H2O_ConstPP_NoFrac] mesh has no cells.");
    }

    double bestDist = DistanceSq(cells.front().center, targetX, targetY);
    int bestIdx = 0;
    for (int i = 1; i < static_cast<int>(cells.size()); ++i) {
        const double dist = DistanceSq(cells[static_cast<std::size_t>(i)].center, targetX, targetY);
        if (dist < bestDist) {
            bestDist = dist;
            bestIdx = i;
        }
    }

    point.block_id = bestIdx;
    point.actual_x = cells[static_cast<std::size_t>(bestIdx)].center.m_x;
    point.actual_y = cells[static_cast<std::size_t>(bestIdx)].center.m_y;
    return point;
}

std::vector<SpatialSamplePoint> BuildProfileStations(const MeshManager& mgr, const TestCaseSpec& cfg) {
    std::vector<SpatialSamplePoint> stations;
    int pointId = 0;

    const double yMid = 0.5 * cfg.ly;
    const double xMid = 0.5 * cfg.lx;
    const int horizontalCount = std::max(cfg.matrix_horizontal_station_count, 2);
    const int verticalCount = std::max(cfg.matrix_vertical_station_count, 2);

    for (int i = 0; i < horizontalCount; ++i) {
        const double fraction = static_cast<double>(i) / static_cast<double>(horizontalCount - 1);
        const double x = fraction * cfg.lx;
        std::ostringstream label;
        label << "mh_" << i;
        stations.push_back(
            BuildNearestSamplePoint(
                mgr,
                pointId++,
                label.str(),
                kFamilyMatrixHorizontal,
                kLocationMatrix,
                x,
                x,
                yMid));
    }

    for (int i = 0; i < verticalCount; ++i) {
        const double fraction = static_cast<double>(i) / static_cast<double>(verticalCount - 1);
        const double y = fraction * cfg.ly;
        std::ostringstream label;
        label << "mv_" << i;
        stations.push_back(
            BuildNearestSamplePoint(
                mgr,
                pointId++,
                label.str(),
                kFamilyMatrixVerticalMidline,
                kLocationMatrix,
                y,
                xMid,
                y));
    }

    return stations;
}

std::vector<SpatialSamplePoint> BuildMonitorPoints(const MeshManager& mgr, const TestCaseSpec& cfg) {
    std::vector<SpatialSamplePoint> points;
    points.push_back(BuildNearestSamplePoint(mgr, 0, "p25", kLocationMatrix, kLocationMatrix, 0.25 * cfg.lx, 0.25 * cfg.lx, 0.5 * cfg.ly));
    points.push_back(BuildNearestSamplePoint(mgr, 1, "p50", kLocationMatrix, kLocationMatrix, 0.50 * cfg.lx, 0.50 * cfg.lx, 0.5 * cfg.ly));
    points.push_back(BuildNearestSamplePoint(mgr, 2, "p75", kLocationMatrix, kLocationMatrix, 0.75 * cfg.lx, 0.75 * cfg.lx, 0.5 * cfg.ly));
    return points;
}

void WriteStageManifest(const std::string& caseCode,
                        const CaseCommon::CaseArtifactPaths& artifacts,
                        CaseCommon::CaseStage stage,
                        const std::string& status,
                        const std::string& outputDir) {
    const std::string context = FailurePrefix(caseCode) + " failed to write stage manifest";
    Case2DReferenceIO::WriteAsciiFile(
        artifacts.engineering_stage_manifest_path,
        context.c_str(),
        [&](std::ofstream& out) {
            out << "case_code=" << caseCode << "\n";
            out << "stage=" << CaseCommon::ToString(stage) << "\n";
            out << "status=" << status << "\n";
            out << "output_dir=" << outputDir << "\n";
        });
}

void WriteReferenceContract(const TestCaseSpec& cfg,
                            const CaseCommon::CaseArtifactPaths& artifacts,
                            CaseCommon::CaseStage stage) {
    const std::string context = FailurePrefix(cfg.case_code) + " failed to write reference contract";
    Case2DReferenceIO::WriteAsciiFile(
        artifacts.reference_contract_path,
        context.c_str(),
        [&](std::ofstream& out) {
            out << "# " << cfg.case_code << " Reference Contract\n\n";
            out << "- Stage: `" << CaseCommon::ToString(stage) << "`\n";
            out << "- Case dir: `" << artifacts.case_dir << "`\n";
            out << "- Reference mode: `comsol`\n";
            out << "- Property mode: `" << PropertyModeString(cfg) << "`\n";
            out << "- Fracture mode: `" << FractureVariantString(cfg.fracture_variant) << "`\n";
            out << "- Validation variables: `pressure`, `temperature`, `co2_saturation`\n";
            out << "- End time: `" << std::setprecision(12) << cfg.target_end_time_s << " s`\n";
        });
}

void WriteStageStatus(const TestCaseSpec& cfg,
                      const CaseCommon::CaseArtifactPaths& artifacts,
                      CaseCommon::CaseStage stage,
                      const std::string& status,
                      const std::string& outputDir) {
    const std::string context = FailurePrefix(cfg.case_code) + " failed to write stage status";
    Case2DReferenceIO::WriteAsciiFile(
        artifacts.report_status_markdown_path,
        context.c_str(),
        [&](std::ofstream& out) {
            out << "# " << cfg.case_code << " Template Status\n\n";
            out << "- Stage: `" << CaseCommon::ToString(stage) << "`\n";
            out << "- Status: `" << status << "`\n";
            out << "- Output dir: `" << outputDir << "`\n";
            out << "- Grid: `" << cfg.nx << "x" << cfg.ny << "`\n";
            out << "- Property mode: `" << PropertyModeString(cfg) << "`\n";
            out << "- Fracture mode: `" << FractureVariantString(cfg.fracture_variant) << "`\n";
            out << "- End time: `" << std::setprecision(12) << cfg.target_end_time_s << " s`\n";
        });
}

void WritePropertyTable(const TestCaseSpec& cfg, const CaseOutputPaths& paths) {
    const std::string context = FailurePrefix(cfg.case_code) + " failed to write property table";
    Case2DReferenceIO::WriteAsciiFile(
        paths.property_table_path,
        context.c_str(),
        [&](std::ofstream& out) {
            out << "phase,rho,mu,cp,cv,k\n";
            out << "water," << cfg.water_props.rho << "," << cfg.water_props.mu << "," << cfg.water_props.cp << ","
                << cfg.water_props.cv << "," << cfg.water_props.k << "\n";
            out << "co2," << cfg.gas_props.rho << "," << cfg.gas_props.mu << "," << cfg.gas_props.cp << ","
                << cfg.gas_props.cv << "," << cfg.gas_props.k << "\n";
        });
}

void WriteReferenceSpec(const TestCaseSpec& cfg,
                        const CaseOutputPaths& paths,
                        const std::vector<SpatialSamplePoint>& profileStations,
                        const std::vector<SpatialSamplePoint>& monitorPoints,
                        const std::vector<SnapshotState>& snapshots,
                        const std::vector<MonitorScheduleState>& monitorSchedule) {
    const std::string context = FailurePrefix(cfg.case_code) + " failed to write reference spec";
    Case2DReferenceIO::WriteAsciiFile(
        paths.reference_spec_path,
        context.c_str(),
        [&](std::ofstream& out) {
            out << "# " << cfg.case_code << " Reference Input Spec\n\n";
            out << "- Grid: `" << cfg.nx << "x" << cfg.ny << "`\n";
            out << "- Domain: `" << cfg.lx << " m x " << cfg.ly << " m`\n";
            out << "- Property mode: `" << PropertyModeString(cfg) << "`\n";
            out << "- Fracture mode: `" << FractureVariantString(cfg.fracture_variant) << "`\n";
            out << "- Variables: `pressure`, `temperature`, `co2_saturation`\n";
            out << "- Stored transported saturation: `water_saturation`\n";
            out << "- Engineering profile stations: `" << profileStations.size() << "`\n";
            out << "- Engineering monitor points: `" << monitorPoints.size() << "`\n";
            out << "- Profile snapshots: `" << snapshots.size() << "`\n";
            out << "- Monitor samples: `" << monitorSchedule.size() << "`\n";
        });
}

void WriteMatlabStub(const TestCaseSpec& cfg,
                     const CaseOutputPaths& paths,
                     const CaseCommon::CaseStage stage,
                     const CaseCommon::CaseArtifactPaths& artifacts) {
    const std::string context = FailurePrefix(cfg.case_code) + " failed to write Matlab stub";
    Case2DReferenceIO::WriteAsciiFile(
        paths.matlab_script_path,
        context.c_str(),
        [&](std::ofstream& out) {
            out << "% " << cfg.case_code << " Matlab plotting placeholder\n";
            out << "% Stage: " << CaseCommon::ToString(stage) << "\n";
            out << "% Case root: " << artifacts.case_dir << "\n";
            out << "% Shared C/F two-phase plotting module is pending.\n";
        });
}

void WriteEngineeringProfileCsvN3(const std::vector<SpatialSamplePoint>& stations,
                                  const std::string& family,
                                  const SnapshotState& snapshot,
                                  const std::string& path,
                                  const std::string& caseCode) {
    const std::string context = FailurePrefix(caseCode) + " failed to write engineering profile csv";
    Case2DReferenceIO::WriteAsciiFile(
        path,
        context.c_str(),
        [&](std::ofstream& out) {
            out << "station_id,label,family,location,target_axis_m,target_x_m,target_y_m,block_id,actual_x_m,actual_y_m,target_time_s,actual_time_s,p_num_pa,t_num_k,sw_num\n";
            for (const auto& station : stations) {
                if (station.family != family) continue;
                const std::size_t block = static_cast<std::size_t>(station.block_id);
                out << station.id << "," << station.label << "," << station.family << "," << station.location << ","
                    << std::setprecision(12) << station.target_axis_m << "," << station.target_x << "," << station.target_y << ","
                    << station.block_id << "," << station.actual_x << "," << station.actual_y << ","
                    << snapshot.target_time_s << "," << snapshot.actual_time_s << ","
                    << snapshot.p_blocks[block] << "," << snapshot.t_blocks[block] << "," << snapshot.sw_blocks[block] << "\n";
            }
        });
}

void WriteEngineeringMonitorCsvN3(const std::vector<SpatialSamplePoint>& points,
                                  const std::vector<MonitorScheduleState>& schedule,
                                  const std::string& path,
                                  const std::string& caseCode) {
    const std::string context = FailurePrefix(caseCode) + " failed to write engineering monitor csv";
    Case2DReferenceIO::WriteAsciiFile(
        path,
        context.c_str(),
        [&](std::ofstream& out) {
            out << "sample_id,target_time_s,actual_time_s";
            for (const auto& point : points) out << ",p_num_" << point.label;
            for (const auto& point : points) out << ",t_num_" << point.label;
            for (const auto& point : points) out << ",sw_num_" << point.label;
            out << "\n";
            for (const auto& sample : schedule) {
                out << sample.sample_id << "," << std::setprecision(12) << sample.target_time_s << "," << sample.actual_time_s;
                for (const auto& point : points) out << "," << sample.p_blocks[static_cast<std::size_t>(point.block_id)];
                for (const auto& point : points) out << "," << sample.t_blocks[static_cast<std::size_t>(point.block_id)];
                for (const auto& point : points) out << "," << sample.sw_blocks[static_cast<std::size_t>(point.block_id)];
                out << "\n";
            }
        });
}

std::vector<std::string> BuildExpectedReferenceFiles(const std::vector<SnapshotState>& snapshots) {
    std::vector<std::string> files;
    for (const auto& snapshot : snapshots) {
        files.push_back(std::string("profile_") + kFamilyMatrixHorizontal + "_" + snapshot.tag + ".csv");
        files.push_back(std::string("profile_") + kFamilyMatrixVerticalMidline + "_" + snapshot.tag + ".csv");
    }
    files.push_back("monitor_timeseries.csv");
    return files;
}

void ValidateReferencePayloadOrThrow(const CaseCommon::CaseArtifactPaths& artifacts,
                                     const CaseOutputPaths& paths,
                                     CaseRunSummary& summary,
                                     const std::vector<SnapshotState>& snapshots) {
    std::vector<std::string> missingFiles;
    const std::string comsolDir = artifacts.reference_dir + "/comsol";
    for (const auto& name : BuildExpectedReferenceFiles(snapshots)) {
        std::ifstream in((comsolDir + "/" + name).c_str(), std::ios::in);
        if (!in.good()) missingFiles.push_back("reference/comsol/" + name);
    }

    if (!missingFiles.empty()) {
        summary.reference_ready = false;
        summary.validation_status = "missing_reference";
        summary.missing_reference_files = missingFiles;
        const std::string context = FailurePrefix(summary.case_code) + " failed to write validation summary";
        Case2DReferenceIO::WriteAsciiFile(
            paths.validation_summary_path,
            context.c_str(),
            [&](std::ofstream& out) {
                out << "# " << summary.case_code << " Validation Summary\n\n";
                out << "- Status: `missing_reference`\n";
                out << "- Missing files:\n";
                for (const auto& file : missingFiles) out << "  - `" << file << "`\n";
            });
        throw std::runtime_error(FailurePrefix(summary.case_code) + " reference files are missing.");
    }

    summary.reference_ready = true;
    summary.validation_status = "reference_payload_ready_metrics_pending";
    const std::string context = FailurePrefix(summary.case_code) + " failed to write validation summary";
    Case2DReferenceIO::WriteAsciiFile(
        paths.validation_summary_path,
        context.c_str(),
        [&](std::ofstream& out) {
            out << "# " << summary.case_code << " Validation Summary\n\n";
            out << "- Status: `reference_payload_ready_metrics_pending`\n";
            out << "- Note: quantitative `C/F` reference comparison is the next shared-module extraction target.\n";
        });
}

void WriteRunSummary(const CaseRunSummary& summary, const CaseOutputPaths& paths) {
    const std::string context = FailurePrefix(summary.case_code) + " failed to write run summary";
    Case2DReferenceIO::WriteAsciiFile(
        paths.run_summary_path,
        context.c_str(),
        [&](std::ofstream& out) {
            out << "case_code=" << summary.case_code << "\n";
            out << "output_dir=" << summary.output_dir << "\n";
            out << "n_cells=" << summary.n_cells << "\n";
            out << "steps=" << summary.steps << "\n";
            out << std::setprecision(12);
            out << "t_end=" << summary.t_end << "\n";
            out << "h_char=" << summary.h_char << "\n";
            out << "validation_status=" << summary.validation_status << "\n";
            out << "reference_ready=" << BoolString(summary.reference_ready) << "\n";
        });
}

void MaterializeReferenceInputs(const TestCaseSpec& cfg,
                                const CaseCommon::CaseArtifactPaths& artifacts,
                                CaseCommon::CaseStage stage) {
    const CaseOutputPaths paths = BuildOutputPaths(artifacts, artifacts.case_dir);

    MeshManager mgr(cfg.lx, cfg.ly, 0.0, cfg.nx, cfg.ny, 0, true, false);
    ConfigureMeshForCase(mgr, cfg);

    const std::vector<SpatialSamplePoint> profileStations = BuildProfileStations(mgr, cfg);
    const std::vector<SpatialSamplePoint> monitorPoints = BuildMonitorPoints(mgr, cfg);
    const std::vector<SnapshotState> snapshots = BuildSnapshots(cfg);
    const std::vector<MonitorScheduleState> monitorSchedule = BuildMonitorSchedule(cfg);

    WritePropertyTable(cfg, paths);
    WriteReferenceSpec(cfg, paths, profileStations, monitorPoints, snapshots, monitorSchedule);
    WriteMatlabStub(cfg, paths, stage, artifacts);
    Case2DReferenceIO::WriteProfileStationDefinitions(profileStations, paths.profile_station_definitions_path);
    Case2DReferenceIO::WriteMonitorPointDefinitions(monitorPoints, paths.monitor_point_definitions_path);
    Case2DReferenceIO::WriteProfileSchedule(snapshots, paths.profile_schedule_path);
    Case2DReferenceIO::WriteMonitorSchedule(monitorSchedule, paths.monitor_schedule_path);
}

CaseRunArtifacts RunCase(const TestCaseSpec& cfg,
                         const CaseCommon::CaseArtifactPaths& artifacts,
                         const std::string& outputDir) {
    CaseRunArtifacts result;
    result.summary.case_code = cfg.case_code;
    result.summary.output_dir = outputDir;

    const CaseOutputPaths paths = BuildOutputPaths(artifacts, outputDir);
    CaseCommon::EnsureDirRecursive(outputDir);

    MeshManager mgr(cfg.lx, cfg.ly, 0.0, cfg.nx, cfg.ny, 0, true, false);
    ConfigureMeshForCase(mgr, cfg);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    const int totalBlocks = mgr.getTotalDOFCount();
    result.summary.n_cells = static_cast<int>(mgr.mesh().getCells().size());
    result.summary.h_char = std::sqrt(std::pow(cfg.lx / static_cast<double>(cfg.nx), 2.0) +
                                      std::pow(cfg.ly / static_cast<double>(cfg.ny), 2.0));

    result.profile_stations = BuildProfileStations(mgr, cfg);
    result.monitor_points = BuildMonitorPoints(mgr, cfg);
    result.snapshots = BuildSnapshots(cfg);
    result.monitor_schedule = BuildMonitorSchedule(cfg);

    std::vector<double> pLatest(static_cast<std::size_t>(std::max(totalBlocks, 0)), cfg.p_init);
    std::vector<double> tLatest(static_cast<std::size_t>(std::max(totalBlocks, 0)), cfg.t_init);
    std::vector<double> swLatest(static_cast<std::size_t>(std::max(totalBlocks, 0)), cfg.sw_init);

    for (auto& sample : result.monitor_schedule) {
        if (sample.target_time_s <= kPendingTimeTolerance) {
            sample.captured = true;
            sample.actual_time_s = 0.0;
            sample.p_blocks = pLatest;
            sample.t_blocks = tLatest;
            sample.sw_blocks = swLatest;
        }
    }

    FIM_Engine::InitialConditions ic;
    ic.P_init = cfg.p_init;
    ic.T_init = cfg.t_init;
    ic.Sw_init = cfg.sw_init;

    BoundarySetting::BoundaryConditionManager bcP;
    BoundarySetting::BoundaryConditionManager bcT;
    BoundarySetting::BoundaryConditionManager bcS;
    FIM_CaseKit::ConfigurePressureBC2D(bcP, cfg.p_left, cfg.p_right, cfg.p_right);
    FIM_CaseKit::ConfigureTemperatureBC2D(bcT, cfg.t_right);
    bcT.SetDirichletBC(MeshTags::LEFT, cfg.t_left);
    bcT.SetDirichletBC(MeshTags::RIGHT, cfg.t_right);
    FIM_CaseKit::ConfigureSaturationBC2D(bcS, cfg.sw_left, cfg.sw_right);

    FIM_CaseKit::PropertyPreset2D preset = FIM_CaseKit::MakeDefaultPropertyPreset2D();
    preset.enable_rock_region = false;
    preset.rock_bg.phi_r = cfg.matrix_phi;
    preset.rock_bg.kxx = cfg.matrix_perm;
    preset.rock_bg.kyy = cfg.matrix_perm;
    preset.rock_bg.kzz = cfg.matrix_perm;
    preset.rock_bg.compressibility = cfg.matrix_ct;
    preset.rock_bg.rho_r = cfg.matrix_rho_r;
    preset.rock_bg.cp_r = cfg.matrix_cp_r;
    preset.rock_bg.k_r = cfg.matrix_lambda_r;
    preset.frac.phi_f = cfg.fracture_phi;
    preset.frac.permeability = cfg.fracture_kt;
    preset.frac.compressibility = cfg.fracture_ct;
    preset.frac.k_f = cfg.fracture_lambda_r;
    preset.use_unified_fluid_model = true;
    if (cfg.use_variable_properties) {
        preset.fluid_model = FIM_Engine::UnifiedFluidModelConfig::MakeTwoPhaseWaterCO2EOS();
    } else {
        preset.fluid_model = FIM_Engine::UnifiedFluidModelConfig::MakeTwoPhaseWaterCO2Constant(
            cfg.water_props,
            cfg.gas_props);
    }

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, &bcS);
    modules.disable_default_vtk_output = true;
    modules.on_step_accepted =
        [&](int step,
            double timeS,
            double,
            int,
            double,
            int,
            const std::string&,
            const std::vector<double>& p,
            const std::vector<double>& t,
            const std::vector<double>* sw) {
            if (!p.empty()) pLatest = p;
            if (!t.empty()) tLatest = t;
            if (sw && !sw->empty()) swLatest = *sw;

            result.summary.steps = step;
            result.summary.t_end = timeS;
            for (auto& snapshot : result.snapshots) CaptureIfReached(timeS, pLatest, tLatest, swLatest, snapshot);
            for (auto& sample : result.monitor_schedule) CaptureIfReached(timeS, pLatest, tLatest, swLatest, sample);
        };

    FIM_Engine::TransientSolverParams params;
    params.output_root_dir = outputDir;
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
    params.max_dSw = cfg.max_dSw;
    params.diag_level = FIM_Engine::DiagLevel::Off;

    FIM_Engine::RunGenericFIMTransient<3>(
        cfg.case_name,
        mgr,
        fm,
        ic,
        {},
        params,
        FIM_Engine::SolverRoute::FIM,
        modules);

    FinalizeMissingSamples(result.summary.t_end, pLatest, tLatest, swLatest, result.snapshots);
    FinalizeMissingSamples(result.summary.t_end, pLatest, tLatest, swLatest, result.monitor_schedule);

    if (cfg.export_vtk) {
        PostProcess_2D(mgr, fm).ExportVTK(paths.final_vtk_path, result.summary.t_end);
    }

    if (cfg.write_engineering_profiles) {
        for (const auto& snapshot : result.snapshots) {
            WriteEngineeringProfileCsvN3(
                result.profile_stations,
                kFamilyMatrixHorizontal,
                snapshot,
                paths.engineering_profile_prefix + kFamilyMatrixHorizontal + "_" + snapshot.tag + ".csv",
                cfg.case_code);
            WriteEngineeringProfileCsvN3(
                result.profile_stations,
                kFamilyMatrixVerticalMidline,
                snapshot,
                paths.engineering_profile_prefix + kFamilyMatrixVerticalMidline + "_" + snapshot.tag + ".csv",
                cfg.case_code);
        }
        WriteEngineeringMonitorCsvN3(
            result.monitor_points,
            result.monitor_schedule,
            paths.engineering_monitor_path,
            cfg.case_code);
    }

    WriteRunSummary(result.summary, paths);
    return result;
}

void PrintPrepareReferenceSummary(const std::string& caseCode,
                                  const CaseCommon::CaseArtifactPaths& artifacts) {
    std::cout << CaseBanner(caseCode) << " prepare_reference completed\n";
    std::cout << "  case_dir      : " << artifacts.case_dir << "\n";
    std::cout << "  engineering   : " << artifacts.engineering_dir << "\n";
    std::cout << "  reference_dir : " << artifacts.reference_dir << "\n";
}

void PrintRunSummary(const std::string& banner, const CaseRunSummary& summary) {
    std::cout << CaseBanner(summary.case_code) << " " << banner << "\n";
    std::cout << "  output_dir : " << summary.output_dir << "\n";
    std::cout << "  n_cells    : " << summary.n_cells << "\n";
    std::cout << "  steps      : " << summary.steps << "\n";
    std::cout << "  t_end      : " << std::scientific << std::setprecision(8) << summary.t_end << " s\n";
    std::cout << "  validation : " << summary.validation_status << "\n";
}

void RunStageByKeyImpl(const std::string& key, CaseCommon::CaseStage stage) {
    const TestCasePlan plan = BuildPlanByKey(key);
    const CaseCommon::CaseArtifactPaths artifacts = BuildArtifactPathsForCase(plan.spec.case_code);
    EnsureArtifactContractDirs(artifacts);
    WriteStageManifest(plan.spec.case_code, artifacts, stage, "started", artifacts.case_dir);
    WriteReferenceContract(plan.spec, artifacts, stage);
    WriteStageStatus(plan.spec, artifacts, stage, "started", artifacts.case_dir);
    MaterializeReferenceInputs(plan.spec, artifacts, stage);

    if (stage == CaseCommon::CaseStage::PrepareReference) {
        WriteStageManifest(plan.spec.case_code, artifacts, stage, "prepared_reference_inputs", artifacts.reference_dir);
        WriteStageStatus(plan.spec, artifacts, stage, "prepared_reference_inputs", artifacts.reference_dir);
        PrintPrepareReferenceSummary(plan.spec.case_code, artifacts);
        return;
    }

    const TestCaseSpec stageSpec = BuildStageSpec(plan.spec, stage);
    const std::string defaultOutputDir =
        (stage == CaseCommon::CaseStage::SolveOnly) ? artifacts.engineering_dir : artifacts.case_dir;
    try {
        CaseRunArtifacts result;
        switch (stage) {
        case CaseCommon::CaseStage::SolveOnly:
            result = RunCase(stageSpec, artifacts, artifacts.engineering_dir);
            result.summary.validation_status = "not_run";
            WriteRunSummary(result.summary, BuildOutputPaths(artifacts, artifacts.engineering_dir));
            WriteStageManifest(plan.spec.case_code, artifacts, stage, "completed", artifacts.engineering_dir);
            WriteStageStatus(stageSpec, artifacts, stage, "completed", artifacts.engineering_dir);
            PrintRunSummary("solve_only completed", result.summary);
            return;
        case CaseCommon::CaseStage::ValidateOnly:
            result = RunCase(stageSpec, artifacts, artifacts.case_dir);
            ValidateReferencePayloadOrThrow(
                artifacts,
                BuildOutputPaths(artifacts, artifacts.case_dir),
                result.summary,
                result.snapshots);
            WriteRunSummary(result.summary, BuildOutputPaths(artifacts, artifacts.case_dir));
            WriteStageManifest(plan.spec.case_code, artifacts, stage, "completed", artifacts.case_dir);
            WriteStageStatus(stageSpec, artifacts, stage, "completed", artifacts.case_dir);
            PrintRunSummary("validate_only completed", result.summary);
            return;
        case CaseCommon::CaseStage::FullWorkflow:
            result = RunCase(stageSpec, artifacts, artifacts.case_dir);
            ValidateReferencePayloadOrThrow(
                artifacts,
                BuildOutputPaths(artifacts, artifacts.case_dir),
                result.summary,
                result.snapshots);
            WriteRunSummary(result.summary, BuildOutputPaths(artifacts, artifacts.case_dir));
            WriteStageManifest(plan.spec.case_code, artifacts, stage, "completed", artifacts.case_dir);
            WriteStageStatus(stageSpec, artifacts, stage, "completed", artifacts.case_dir);
            PrintRunSummary("full_workflow completed", result.summary);
            return;
        default:
            throw std::runtime_error(FailurePrefix(plan.spec.case_code) + " unsupported stage in RunStageByKeyImpl.");
        }
    } catch (const std::exception& ex) {
        const std::string message = ex.what();
        const bool missingReference =
            message.find("missing") != std::string::npos &&
            message.find("reference") != std::string::npos;
        const std::string failureStatus = missingReference ? "missing_reference" : "failed";
        WriteStageManifest(plan.spec.case_code, artifacts, stage, failureStatus, defaultOutputDir);
        WriteStageStatus(stageSpec, artifacts, stage, failureStatus, defaultOutputDir);
        throw;
    }
}

void ExecutePlanByKeyImpl(const std::string& key) {
    RunStageByKeyImpl(key, CaseCommon::CaseStage::FullWorkflow);
}

} // namespace

void RunTestCase() {
    ExecutePlanByKeyImpl(kPlanKeyC1);
}

void ExecutePlanByKey(const std::string& key) {
    ExecutePlanByKeyImpl(key);
}

void RunStageByKey(const std::string& key, CaseCommon::CaseStage stage) {
    RunStageByKeyImpl(key, stage);
}

void RunSolveOnly() {
    RunStageByKeyImpl(kPlanKeyC1, CaseCommon::CaseStage::SolveOnly);
}

void RunPrepareReference() {
    RunStageByKeyImpl(kPlanKeyC1, CaseCommon::CaseStage::PrepareReference);
}

void RunValidateOnly() {
    RunStageByKeyImpl(kPlanKeyC1, CaseCommon::CaseStage::ValidateOnly);
}

void RunFullWorkflow() {
    RunStageByKeyImpl(kPlanKeyC1, CaseCommon::CaseStage::FullWorkflow);
}

} // namespace Test_H_TP_CO2H2O_ConstPP_NoFrac

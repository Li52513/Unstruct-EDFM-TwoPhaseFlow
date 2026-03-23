/**
 * @file Test_Day6_TransientSolver.cpp
 * @brief Day6 transient scenarios (thin wrappers over FIM_TransientCaseKit)
 *
 * Integrated features exercised by Run_Day6_Transient_2D_SP_InjProd():
 *   - CoolProp EOS (Step 1): Span-Wagner CO2 + IAPWS-95 Water via USE_COOLPROP_EOS
 *   - Non-orthogonal correction (Step 2): deferred T-vector gradient correction
 *   - Independent well DOF (Step 3): WellDOFManager adds P_wbh block per well
 *   - VTU/PVD post-processing (Step 4): ParaView XML Unstructured Grid export
 */

#include "Test_Day6_TransientSolver.h"
#include "FIM_TransientCaseKit.hpp"
#include "2D_PostProcess.h"
#include "Test_Issue11_FrozenMatrix.h"
#include "Test_Issue12_LinearSolverMemory.h"

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#define DAY6_MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define DAY6_MKDIR(path) mkdir(path, 0777)
#endif

namespace Test_Day6 {

namespace {

enum class TopologyAxis {
    F0 = 0,
    F1 = 1,
    F2 = 2
};

struct PhysicalScenarioDef {
    int id = 0;
    bool variable_properties = false;
    bool with_wells = false;
    bool two_phase = false;
    bool heat_enabled = false;
    const char* description = "";
};

constexpr std::array<PhysicalScenarioDef, 16> kPhysicalScenarios = {{
    {  1, false, false, false, false, "Const / no-well / single-phase flow" },
    {  2, false, false, false, true,  "Const / no-well / single-phase flow+heat" },
    {  3, false, false, true,  false, "Const / no-well / two-phase flow" },
    {  4, false, false, true,  true,  "Const / no-well / two-phase flow+heat" },

    {  5, false, true,  false, false, "Const / with-well / single-phase flow" },
    {  6, false, true,  false, true,  "Const / with-well / single-phase flow+heat" },
    {  7, false, true,  true,  false, "Const / with-well / two-phase flow" },
    {  8, false, true,  true,  true,  "Const / with-well / two-phase flow+heat" },

    {  9, true,  false, false, false, "Variable / no-well / single-phase flow" },
    { 10, true,  false, false, true,  "Variable / no-well / single-phase flow+heat" },
    { 11, true,  false, true,  false, "Variable / no-well / two-phase flow" },
    { 12, true,  false, true,  true,  "Variable / no-well / two-phase flow+heat" },

    { 13, true,  true,  false, false, "Variable / with-well / single-phase flow" },
    { 14, true,  true,  false, true,  "Variable / with-well / single-phase flow+heat" },
    { 15, true,  true,  true,  false, "Variable / with-well / two-phase flow" },
    { 16, true,  true,  true,  true,  "Variable / with-well / two-phase flow+heat" }
}};

constexpr std::array<TopologyAxis, 3> kTopologyOrder = {
    TopologyAxis::F0,
    TopologyAxis::F1,
    TopologyAxis::F2
};

std::string ToTxx(int id) {
    std::ostringstream oss;
    oss << "T" << std::setw(2) << std::setfill('0') << id;
    return oss.str();
}

std::string ToFy(TopologyAxis topology) {
    switch (topology) {
    case TopologyAxis::F0: return "F0";
    case TopologyAxis::F1: return "F1";
    case TopologyAxis::F2: return "F2";
    default:               return "F0";
    }
}

std::string DimLabel(int dim) {
    return (dim == 3) ? "3D" : "2D";
}

std::string BuildScenarioKey(const PhysicalScenarioDef& scenario, TopologyAxis topology) {
    return ToTxx(scenario.id) + "_" + ToFy(topology);
}

std::string BuildCaseFolderName(int dim, const std::string& scenarioKey) {
    return DimLabel(dim) + "_" + scenarioKey;
}

bool IsMilestoneScenario(int scenarioId) {
    return scenarioId == 4 || scenarioId == 8 || scenarioId == 12 || scenarioId == 16;
}

void EnsureDirRecursive(const std::string& rawPath) {
    if (rawPath.empty()) {
        return;
    }
    std::string path = rawPath;
    for (char& ch : path) {
        if (ch == '\\') {
            ch = '/';
        }
    }

    std::stringstream ss(path);
    std::string token;
    std::string current;
    while (std::getline(ss, token, '/')) {
        if (token.empty() || token == ".") {
            continue;
        }
        if (!current.empty()) {
            current += "/";
        }
        current += token;
        DAY6_MKDIR(current.c_str());
    }
}

bool FileExistsNonEmpty(const std::string& path) {
    std::ifstream ifs(path, std::ios::binary | std::ios::ate);
    if (!ifs.good()) {
        return false;
    }
    return ifs.tellg() > 0;
}

std::string CsvEscape(const std::string& input) {
    const bool needQuote = (input.find(',') != std::string::npos) ||
                           (input.find('"') != std::string::npos) ||
                           (input.find('\n') != std::string::npos) ||
                           (input.find('\r') != std::string::npos);
    if (!needQuote) {
        return input;
    }
    std::string escaped = "\"";
    for (char ch : input) {
        if (ch == '"') {
            escaped += "\"\"";
        }
        else {
            escaped += ch;
        }
    }
    escaped += "\"";
    return escaped;
}

std::string NowUtcIso8601() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t nowT = std::chrono::system_clock::to_time_t(now);
    std::tm gmt;
#ifdef _WIN32
    gmtime_s(&gmt, &nowT);
#else
    gmtime_r(&nowT, &gmt);
#endif
    std::ostringstream oss;
    oss << std::put_time(&gmt, "%Y-%m-%dT%H:%M:%SZ");
    return oss.str();
}

class DualStreamBuf : public std::streambuf {
public:
    DualStreamBuf() = default;

    void Reset(std::streambuf* left, std::streambuf* right) {
        left_ = left;
        right_ = right;
    }

protected:
    int overflow(int c) override {
        if (c == EOF) {
            return !EOF;
        }
        const int leftResult = left_ ? left_->sputc(static_cast<char>(c)) : c;
        const int rightResult = right_ ? right_->sputc(static_cast<char>(c)) : c;
        if (leftResult == EOF || rightResult == EOF) {
            return EOF;
        }
        return c;
    }

    int sync() override {
        const int leftSync = left_ ? left_->pubsync() : 0;
        const int rightSync = right_ ? right_->pubsync() : 0;
        return (leftSync == 0 && rightSync == 0) ? 0 : -1;
    }

private:
    std::streambuf* left_ = nullptr;
    std::streambuf* right_ = nullptr;
};

class ScopedCoutLog {
public:
    explicit ScopedCoutLog(const std::string& logPath) {
        log_.open(logPath, std::ios::out | std::ios::trunc);
        previous_ = std::cout.rdbuf();
        if (log_.is_open()) {
            tee_.Reset(previous_, log_.rdbuf());
            std::cout.rdbuf(&tee_);
            active_ = true;
        }
    }

    ~ScopedCoutLog() {
        if (active_) {
            std::cout.rdbuf(previous_);
        }
    }

    ScopedCoutLog(const ScopedCoutLog&) = delete;
    ScopedCoutLog& operator=(const ScopedCoutLog&) = delete;

private:
    std::ofstream log_;
    std::streambuf* previous_ = nullptr;
    DualStreamBuf tee_;
    bool active_ = false;
};

int FindNearestCell2D(const MeshManager& mgr, double x, double y) {
    const auto& cells = mgr.mesh().getCells();
    int bestIdx = -1;
    double bestD2 = std::numeric_limits<double>::max();
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        const double dx = cells[i].center.m_x - x;
        const double dy = cells[i].center.m_y - y;
        const double d2 = dx * dx + dy * dy;
        if (d2 < bestD2) {
            bestD2 = d2;
            bestIdx = i;
        }
    }
    return (bestIdx >= 0) ? bestIdx : 0;
}

int FindNearestCell3D(const MeshManager_3D& mgr, double x, double y, double z) {
    const auto& cells = mgr.mesh().getCells();
    int bestIdx = -1;
    double bestD2 = std::numeric_limits<double>::max();
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        const double dx = cells[i].center.m_x - x;
        const double dy = cells[i].center.m_y - y;
        const double dz = cells[i].center.m_z - z;
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < bestD2) {
            bestD2 = d2;
            bestIdx = i;
        }
    }
    return (bestIdx >= 0) ? bestIdx : 0;
}

void BuildTopology2D(MeshManager& mgr, TopologyAxis topology) {
    if (topology == TopologyAxis::F1 || topology == TopologyAxis::F2) {
        mgr.addFracture(Vector(20.0, 4.0, 0.0), Vector(380.0, 36.0, 0.0));
    }
    if (topology == TopologyAxis::F2) {
        mgr.addFracture(Vector(20.0, 36.0, 0.0), Vector(380.0, 4.0, 0.0));
    }

    if (topology != TopologyAxis::F0) {
        mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    }
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
}

void BuildTopology3D(MeshManager_3D& mgr, TopologyAxis topology) {
    if (topology == TopologyAxis::F1 || topology == TopologyAxis::F2) {
        std::vector<Vector> pts1 = {
            Vector(10.0, 4.0, 4.0),
            Vector(110.0, 4.0, 4.0),
            Vector(110.0, 36.0, 36.0),
            Vector(10.0, 36.0, 36.0)
        };
        mgr.addFracturetoFractureNetwork(Fracture_2D(0, pts1));
    }

    if (topology == TopologyAxis::F2) {
        std::vector<Vector> pts2 = {
            Vector(10.0, 36.0, 4.0),
            Vector(110.0, 36.0, 4.0),
            Vector(110.0, 4.0, 36.0),
            Vector(10.0, 4.0, 36.0)
        };
        mgr.addFracturetoFractureNetwork(Fracture_2D(1, pts2));
    }

    if (topology != TopologyAxis::F0) {
        mgr.meshAllFracturesinNetwork(3, 3);
    }

    mgr.setupGlobalIndices();
    if (topology != TopologyAxis::F0) {
        mgr.DetectFractureFractureIntersectionsInNetwork();
        mgr.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);
        mgr.fracture_network().rebuildEdgeProperties();
        mgr.removeDuplicateInteractions();
        mgr.resolveCoplanarInteractions();
    }
    mgr.buildTopologyMaps();
}

void ConfigurePressureBC2D(BoundarySetting::BoundaryConditionManager& bc, bool withWells, double pInit) {
    bc.Clear();
    if (withWells) {
        bc.SetNeumannBC(MeshTags::LEFT, 0.0);
        bc.SetNeumannBC(MeshTags::RIGHT, 0.0);
        bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
        bc.SetNeumannBC(MeshTags::TOP, 0.0);
        return;
    }

    bc.SetDirichletBC(MeshTags::LEFT, pInit + 2.0e6);
    bc.SetDirichletBC(MeshTags::RIGHT, pInit - 2.0e6);
    bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bc.SetNeumannBC(MeshTags::TOP, 0.0);
}

void ConfigureTemperatureBC2D(BoundarySetting::BoundaryConditionManager& bc, bool heatEnabled, double tInit) {
    bc.Clear();
    if (heatEnabled) {
        bc.SetDirichletBC(MeshTags::LEFT, tInit + 12.0);
        bc.SetDirichletBC(MeshTags::RIGHT, tInit - 12.0);
        bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
        bc.SetNeumannBC(MeshTags::TOP, 0.0);
        return;
    }

    bc.SetDirichletBC(MeshTags::LEFT, tInit);
    bc.SetDirichletBC(MeshTags::RIGHT, tInit);
    bc.SetDirichletBC(MeshTags::BOTTOM, tInit);
    bc.SetDirichletBC(MeshTags::TOP, tInit);
}

void ConfigureSaturationBC2D(BoundarySetting::BoundaryConditionManager& bc, bool withWells) {
    bc.Clear();
    if (withWells) {
        bc.SetNeumannBC(MeshTags::LEFT, 0.0);
        bc.SetNeumannBC(MeshTags::RIGHT, 0.0);
        bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
        bc.SetNeumannBC(MeshTags::TOP, 0.0);
        return;
    }

    bc.SetDirichletBC(MeshTags::LEFT, 0.85);
    bc.SetDirichletBC(MeshTags::RIGHT, 0.25);
    bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bc.SetNeumannBC(MeshTags::TOP, 0.0);
}

void ConfigurePressureBC3D(BoundarySetting::BoundaryConditionManager& bc, bool withWells, double pInit) {
    bc.Clear();
    if (withWells) {
        bc.SetNeumannBC(MeshTags::LEFT, 0.0);
        bc.SetNeumannBC(MeshTags::RIGHT, 0.0);
        bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
        bc.SetNeumannBC(MeshTags::TOP, 0.0);
        bc.SetNeumannBC(MeshTags::TAG_FRONT, 0.0);
        bc.SetNeumannBC(MeshTags::TAG_BACK, 0.0);
        return;
    }

    bc.SetDirichletBC(MeshTags::LEFT, pInit + 2.0e6);
    bc.SetDirichletBC(MeshTags::RIGHT, pInit - 2.0e6);
    bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bc.SetNeumannBC(MeshTags::TOP, 0.0);
    bc.SetNeumannBC(MeshTags::TAG_FRONT, 0.0);
    bc.SetNeumannBC(MeshTags::TAG_BACK, 0.0);
}

void ConfigureTemperatureBC3D(BoundarySetting::BoundaryConditionManager& bc, bool heatEnabled, double tInit) {
    bc.Clear();
    if (heatEnabled) {
        bc.SetDirichletBC(MeshTags::LEFT, tInit + 12.0);
        bc.SetDirichletBC(MeshTags::RIGHT, tInit - 12.0);
        bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
        bc.SetNeumannBC(MeshTags::TOP, 0.0);
        bc.SetNeumannBC(MeshTags::TAG_FRONT, 0.0);
        bc.SetNeumannBC(MeshTags::TAG_BACK, 0.0);
        return;
    }

    bc.SetDirichletBC(MeshTags::LEFT, tInit);
    bc.SetDirichletBC(MeshTags::RIGHT, tInit);
    bc.SetDirichletBC(MeshTags::BOTTOM, tInit);
    bc.SetDirichletBC(MeshTags::TOP, tInit);
    bc.SetDirichletBC(MeshTags::TAG_FRONT, tInit);
    bc.SetDirichletBC(MeshTags::TAG_BACK, tInit);
}

void ConfigureSaturationBC3D(BoundarySetting::BoundaryConditionManager& bc, bool withWells) {
    bc.Clear();
    if (withWells) {
        bc.SetNeumannBC(MeshTags::LEFT, 0.0);
        bc.SetNeumannBC(MeshTags::RIGHT, 0.0);
        bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
        bc.SetNeumannBC(MeshTags::TOP, 0.0);
        bc.SetNeumannBC(MeshTags::TAG_FRONT, 0.0);
        bc.SetNeumannBC(MeshTags::TAG_BACK, 0.0);
        return;
    }

    bc.SetDirichletBC(MeshTags::LEFT, 0.85);
    bc.SetDirichletBC(MeshTags::RIGHT, 0.25);
    bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bc.SetNeumannBC(MeshTags::TOP, 0.0);
    bc.SetNeumannBC(MeshTags::TAG_FRONT, 0.0);
    bc.SetNeumannBC(MeshTags::TAG_BACK, 0.0);
}

FIM_CaseKit::PropertyPreset2D BuildPreset2D(const PhysicalScenarioDef& scenario) {
    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset2D();

    preset.enable_rock_region = scenario.variable_properties;
    preset.rock_bg.phi_r = 0.10;
    preset.rock_bg.kxx = 1.0e-13;
    preset.rock_bg.kyy = 1.0e-13;
    preset.rock_bg.kzz = 1.0e-13;
    preset.rock_bg.compressibility = scenario.variable_properties ? 5.0e-9 : 1.0e-11;

    preset.frac.phi_f = 0.55;
    preset.frac.permeability = scenario.variable_properties ? 5.0e-11 : 1.0e-11;
    preset.frac.compressibility = scenario.variable_properties ? 8.0e-9 : 1.0e-11;
    preset.frac.aperture = 1.0e-3;

    if (scenario.variable_properties) {
        preset.rock_region = preset.rock_bg;
        preset.rock_region.kxx = 4.0e-13;
        preset.rock_region.kyy = 2.0e-13;
        preset.rock_region.kzz = 2.0e-13;
        preset.rock_region.compressibility = 1.0e-8;
        preset.rock_region_box = AABB(Vector(120.0, 8.0, -1.0), Vector(280.0, 32.0, 1.0));
    }

    if (!scenario.two_phase) {
        preset.single_phase_fluid = scenario.variable_properties
            ? FIM_Engine::SinglePhaseFluidModel::CO2
            : FIM_Engine::SinglePhaseFluidModel::Water;
    }

    return preset;
}

FIM_CaseKit::PropertyPreset3D BuildPreset3D(const PhysicalScenarioDef& scenario) {
    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset3D();

    if (!scenario.variable_properties) {
        preset.enable_rock_region = false;
        preset.rock_bg = RockPropertyParams(100.0, 100.0, 50.0, 0.12, 2600.0, 1000.0, 2.0, 1.0e-11);
        preset.rock_region = preset.rock_bg;
    }
    else {
        preset.enable_rock_region = true;
        preset.rock_bg = RockPropertyParams(100.0, 100.0, 50.0, 0.16, 2600.0, 1000.0, 2.0, 5.0e-9);
        preset.rock_region = RockPropertyParams(300.0, 200.0, 120.0, 0.18, 2600.0, 1000.0, 2.0, 1.0e-8);
        preset.rock_region_box = BoundingBox3D{30.0, 90.0, 8.0, 32.0, 8.0, 32.0};
    }

    if (!scenario.two_phase) {
        preset.single_phase_fluid = scenario.variable_properties
            ? FIM_Engine::SinglePhaseFluidModel::CO2
            : FIM_Engine::SinglePhaseFluidModel::Water;
    }

    return preset;
}

FIM_Engine::InitialConditions BuildInitialConditions(const PhysicalScenarioDef& scenario) {
    FIM_Engine::InitialConditions ic;
    ic.P_init = scenario.variable_properties ? 25.0e6 : 10.0e6;
    if (scenario.with_wells) {
        ic.P_init += 5.0e6;
    }
    ic.T_init = scenario.heat_enabled ? 380.0 : 360.0;
    ic.Sw_init = scenario.two_phase ? 0.80 : 1.0;
    return ic;
}

FIM_Engine::TransientSolverParams BuildScenarioParams(const PhysicalScenarioDef& scenario) {
    const double dtInit = scenario.two_phase ? 25.0 : 100.0;
    auto params = FIM_CaseKit::BuildSolverParams(scenario.two_phase, 60, dtInit);

    params.gravity_vector = Vector(0.0, 0.0, 0.0);
    params.enable_non_orthogonal_correction = true;
    params.lin_solver = scenario.two_phase
        ? FIM_Engine::LinearSolverType::AMGCL_CPR
        : FIM_Engine::LinearSolverType::AMGCL;

    params.max_newton_iter = scenario.two_phase ? 18 : 14;
    params.max_dT = scenario.heat_enabled ? 5.0 : 0.25;
    params.max_dSw = scenario.two_phase ? 0.08 : params.max_dSw;
    params.dt_min = 1.0e-4;

    return params;
}

std::vector<WellScheduleStep> BuildWells2D(
    const MeshManager& mgr,
    const PhysicalScenarioDef& scenario,
    const FIM_Engine::InitialConditions& ic,
    bool singlePhaseIsCO2) {

    if (!scenario.with_wells) {
        return {};
    }

    const int injCell = FindNearestCell2D(mgr, 80.0, 20.0);
    const int prodCell = FindNearestCell2D(mgr, 320.0, 20.0);

    WellScheduleStep inj;
    inj.well_name = "INJ_" + ToTxx(scenario.id);
    inj.domain = WellTargetDomain::Matrix;
    inj.control_mode = WellControlMode::BHP;
    inj.target_value = ic.P_init + 1.0e6;
    inj.completion_id = injCell;
    inj.injection_temperature = scenario.heat_enabled ? (ic.T_init - 20.0) : ic.T_init;

    WellScheduleStep prod;
    prod.well_name = "PROD_" + ToTxx(scenario.id);
    prod.domain = WellTargetDomain::Matrix;
    prod.control_mode = WellControlMode::BHP;
    prod.target_value = ic.P_init - 1.0e6;
    prod.completion_id = prodCell;

    if (scenario.two_phase) {
        inj.component_mode = scenario.variable_properties ? WellComponentMode::Gas : WellComponentMode::Water;
        if (inj.component_mode == WellComponentMode::Gas) {
            inj.frac_w = 0.0;
            inj.frac_g = 1.0;
            inj.injection_is_co2 = true;
        }
        else {
            inj.frac_w = 1.0;
            inj.frac_g = 0.0;
            inj.injection_is_co2 = false;
        }

        prod.component_mode = WellComponentMode::Total;
    }
    else {
        inj.component_mode = singlePhaseIsCO2 ? WellComponentMode::Gas : WellComponentMode::Water;
        inj.frac_w = singlePhaseIsCO2 ? 0.0 : 1.0;
        inj.frac_g = singlePhaseIsCO2 ? 1.0 : 0.0;
        inj.injection_is_co2 = singlePhaseIsCO2;

        prod.component_mode = inj.component_mode;
        prod.frac_w = inj.frac_w;
        prod.frac_g = inj.frac_g;
    }

    return { inj, prod };
}

std::vector<WellScheduleStep> BuildWells3D(
    const MeshManager_3D& mgr,
    const PhysicalScenarioDef& scenario,
    const FIM_Engine::InitialConditions& ic,
    bool singlePhaseIsCO2) {

    if (!scenario.with_wells) {
        return {};
    }

    const int injCell = FindNearestCell3D(mgr, 30.0, 20.0, 20.0);
    const int prodCell = FindNearestCell3D(mgr, 90.0, 20.0, 20.0);

    WellScheduleStep inj;
    inj.well_name = "INJ_" + ToTxx(scenario.id);
    inj.domain = WellTargetDomain::Matrix;
    inj.control_mode = WellControlMode::BHP;
    inj.target_value = ic.P_init + 1.0e6;
    inj.completion_id = injCell;
    inj.well_axis = WellAxis::Z;
    inj.injection_temperature = scenario.heat_enabled ? (ic.T_init - 20.0) : ic.T_init;

    WellScheduleStep prod;
    prod.well_name = "PROD_" + ToTxx(scenario.id);
    prod.domain = WellTargetDomain::Matrix;
    prod.control_mode = WellControlMode::BHP;
    prod.target_value = ic.P_init - 1.0e6;
    prod.completion_id = prodCell;
    prod.well_axis = WellAxis::Z;

    if (scenario.two_phase) {
        inj.component_mode = scenario.variable_properties ? WellComponentMode::Gas : WellComponentMode::Water;
        if (inj.component_mode == WellComponentMode::Gas) {
            inj.frac_w = 0.0;
            inj.frac_g = 1.0;
            inj.injection_is_co2 = true;
        }
        else {
            inj.frac_w = 1.0;
            inj.frac_g = 0.0;
            inj.injection_is_co2 = false;
        }

        prod.component_mode = WellComponentMode::Total;
    }
    else {
        inj.component_mode = singlePhaseIsCO2 ? WellComponentMode::Gas : WellComponentMode::Water;
        inj.frac_w = singlePhaseIsCO2 ? 0.0 : 1.0;
        inj.frac_g = singlePhaseIsCO2 ? 1.0 : 0.0;
        inj.injection_is_co2 = singlePhaseIsCO2;

        prod.component_mode = inj.component_mode;
        prod.frac_w = inj.frac_w;
        prod.frac_g = inj.frac_g;
    }

    return { inj, prod };
}

struct RunOutcome {
    bool success = false;
    bool final_vtk_ok = false;
    double elapsed_s = 0.0;
    std::string scenario_key;
    std::string case_folder;
    std::string message;
};

void WriteMetricsCsv(
    const PhysicalScenarioDef& scenario,
    TopologyAxis topology,
    int dim,
    const RunOutcome& outcome,
    const std::string& logPath,
    const std::string& finalVtkPath) {

    const std::string campaignDir = "Test/Transient/Day6/Campaign/" + DimLabel(dim) + "/" + outcome.scenario_key;
    EnsureDirRecursive(campaignDir);

    const std::string csvPath = campaignDir + "/metrics.csv";
    std::ofstream ofs(csvPath, std::ios::out | std::ios::trunc);
    ofs << "timestamp_utc,dimension,scenario_id,topology,case_folder,status,elapsed_s,final_vtk_ok,log_path,final_vtk_path,variable_properties,with_wells,two_phase,heat_enabled,message\n";
    ofs << CsvEscape(NowUtcIso8601()) << ","
        << CsvEscape(DimLabel(dim)) << ","
        << CsvEscape(ToTxx(scenario.id)) << ","
        << CsvEscape(ToFy(topology)) << ","
        << CsvEscape(outcome.case_folder) << ","
        << CsvEscape(outcome.success ? "PASS" : "FAIL") << ","
        << std::scientific << std::setprecision(6) << outcome.elapsed_s << ","
        << (outcome.final_vtk_ok ? "1" : "0") << ","
        << CsvEscape(logPath) << ","
        << CsvEscape(finalVtkPath) << ","
        << (scenario.variable_properties ? "1" : "0") << ","
        << (scenario.with_wells ? "1" : "0") << ","
        << (scenario.two_phase ? "1" : "0") << ","
        << (scenario.heat_enabled ? "1" : "0") << ","
        << CsvEscape(outcome.message)
        << "\n";
}

void ThrowIfOutcomeFailed(const RunOutcome& outcome, int dim) {
    if (outcome.success && outcome.final_vtk_ok) {
        return;
    }
    std::ostringstream oss;
    oss << "[DAY6-CAMPAIGN-FAIL] " << DimLabel(dim) << " " << outcome.scenario_key
        << " failed. status=" << (outcome.success ? "PASS" : "FAIL")
        << " final_vtk_ok=" << (outcome.final_vtk_ok ? "1" : "0")
        << " message=" << outcome.message;
    throw std::runtime_error(oss.str());
}

const PhysicalScenarioDef& GetScenarioOrThrow(int scenarioId) {
    for (const auto& scenario : kPhysicalScenarios) {
        if (scenario.id == scenarioId) {
            return scenario;
        }
    }
    throw std::runtime_error("Unknown Day6 physical scenario id: " + std::to_string(scenarioId));
}

RunOutcome RunScenario2D(const PhysicalScenarioDef& scenario, TopologyAxis topology) {
    RunOutcome outcome;
    outcome.scenario_key = BuildScenarioKey(scenario, topology);
    outcome.case_folder = BuildCaseFolderName(2, outcome.scenario_key);

    const std::string caseDir = "Test/Transient/Day6/" + outcome.case_folder;
    EnsureDirRecursive(caseDir);

    const std::string logPath = caseDir + "/convergence.log";
    const std::string finalVtkPath = caseDir + "/final.vtk";

    const auto begin = std::chrono::steady_clock::now();
    try {
        ScopedCoutLog scopedLog(logPath);

        std::cout << "\n[DAY6-CAMPAIGN] Start 2D " << outcome.scenario_key
                  << " | " << scenario.description << "\n";

        MeshManager mgr(400.0, 40.0, 0.0, 40, 4, 0, true, false);
        mgr.BuildSolidMatrixGrid_2D();
        BuildTopology2D(mgr, topology);
        mgr.setNumDOFs(scenario.two_phase ? 3 : 2);

        FieldManager_2D fm;
        FIM_CaseKit::InitFieldManager(mgr, fm);

        auto preset = BuildPreset2D(scenario);
        auto ic = BuildInitialConditions(scenario);

        BoundarySetting::BoundaryConditionManager bcP;
        BoundarySetting::BoundaryConditionManager bcT;
        BoundarySetting::BoundaryConditionManager bcS;

        ConfigurePressureBC2D(bcP, scenario.with_wells, ic.P_init);
        ConfigureTemperatureBC2D(bcT, scenario.heat_enabled, ic.T_init);

        BoundarySetting::BoundaryConditionManager* satBC = nullptr;
        if (scenario.two_phase) {
            ConfigureSaturationBC2D(bcS, scenario.with_wells);
            satBC = &bcS;
        }

        auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, satBC);
        const bool singlePhaseIsCO2 =
            (!scenario.two_phase) && (preset.single_phase_fluid == FIM_Engine::SinglePhaseFluidModel::CO2);
        auto wells = BuildWells2D(mgr, scenario, ic, singlePhaseIsCO2);
        auto params = BuildScenarioParams(scenario);

        if (scenario.two_phase) {
            FIM_Engine::RunGenericFIMTransient<3>(
                outcome.case_folder,
                mgr,
                fm,
                ic,
                wells,
                params,
                FIM_Engine::SolverRoute::FIM,
                modules);
        }
        else {
            FIM_Engine::RunGenericFIMTransient<2>(
                outcome.case_folder,
                mgr,
                fm,
                ic,
                wells,
                params,
                FIM_Engine::SolverRoute::FIM,
                modules);
        }

        outcome.success = true;
        outcome.message = "PASS";
    }
    catch (const std::exception& ex) {
        outcome.success = false;
        outcome.message = ex.what();
    }
    catch (...) {
        outcome.success = false;
        outcome.message = "unknown_exception";
    }

    const auto end = std::chrono::steady_clock::now();
    outcome.elapsed_s = std::chrono::duration<double>(end - begin).count();
    outcome.final_vtk_ok = FileExistsNonEmpty(finalVtkPath);

    WriteMetricsCsv(scenario, topology, 2, outcome, logPath, finalVtkPath);
    return outcome;
}

RunOutcome RunScenario3D(const PhysicalScenarioDef& scenario, TopologyAxis topology) {
    RunOutcome outcome;
    outcome.scenario_key = BuildScenarioKey(scenario, topology);
    outcome.case_folder = BuildCaseFolderName(3, outcome.scenario_key);

    const std::string caseDir = "Test/Transient/Day6/" + outcome.case_folder;
    EnsureDirRecursive(caseDir);

    const std::string logPath = caseDir + "/convergence.log";
    const std::string finalVtkPath = caseDir + "/final.vtk";

    const auto begin = std::chrono::steady_clock::now();
    try {
        ScopedCoutLog scopedLog(logPath);

        std::cout << "\n[DAY6-CAMPAIGN] Start 3D " << outcome.scenario_key
                  << " | " << scenario.description << "\n";

        MeshManager_3D mgr(120.0, 40.0, 40.0, 12, 4, 4, true, false);
        mgr.BuildSolidMatrixGrid_3D();
        BuildTopology3D(mgr, topology);
        mgr.setNumDOFs(scenario.two_phase ? 3 : 2);

        FieldManager_3D fm;
        FIM_CaseKit::InitFieldManager(mgr, fm);

        auto preset = BuildPreset3D(scenario);
        auto ic = BuildInitialConditions(scenario);

        BoundarySetting::BoundaryConditionManager bcP;
        BoundarySetting::BoundaryConditionManager bcT;
        BoundarySetting::BoundaryConditionManager bcS;

        ConfigurePressureBC3D(bcP, scenario.with_wells, ic.P_init);
        ConfigureTemperatureBC3D(bcT, scenario.heat_enabled, ic.T_init);

        BoundarySetting::BoundaryConditionManager* satBC = nullptr;
        if (scenario.two_phase) {
            ConfigureSaturationBC3D(bcS, scenario.with_wells);
            satBC = &bcS;
        }

        auto modules = FIM_CaseKit::BuildModules3D(preset, &bcP, &bcT, satBC);
        const bool singlePhaseIsCO2 =
            (!scenario.two_phase) && (preset.single_phase_fluid == FIM_Engine::SinglePhaseFluidModel::CO2);
        auto wells = BuildWells3D(mgr, scenario, ic, singlePhaseIsCO2);
        auto params = BuildScenarioParams(scenario);

        if (scenario.two_phase) {
            FIM_Engine::RunGenericFIMTransient<3>(
                outcome.case_folder,
                mgr,
                fm,
                ic,
                wells,
                params,
                FIM_Engine::SolverRoute::FIM,
                modules);
        }
        else {
            FIM_Engine::RunGenericFIMTransient<2>(
                outcome.case_folder,
                mgr,
                fm,
                ic,
                wells,
                params,
                FIM_Engine::SolverRoute::FIM,
                modules);
        }

        outcome.success = true;
        outcome.message = "PASS";
    }
    catch (const std::exception& ex) {
        outcome.success = false;
        outcome.message = ex.what();
    }
    catch (...) {
        outcome.success = false;
        outcome.message = "unknown_exception";
    }

    const auto end = std::chrono::steady_clock::now();
    outcome.elapsed_s = std::chrono::duration<double>(end - begin).count();
    outcome.final_vtk_ok = FileExistsNonEmpty(finalVtkPath);

    WriteMetricsCsv(scenario, topology, 3, outcome, logPath, finalVtkPath);
    return outcome;
}

void RunMilestoneRegression2D(int scenarioId) {
    std::cout << "\n[DAY6-CAMPAIGN] 2D milestone regression at " << ToTxx(scenarioId)
              << " (matrix_audit + issue11 + issue12)\n";
    Run_Day6_MatrixAudit_2D_EDFM();
    Test_Issue11::Run_All();
    Test_Issue12::Run_All();
    std::cout << "[DAY6-CAMPAIGN] 2D milestone regression pass at " << ToTxx(scenarioId) << "\n";
}

void RunMilestoneRegression3D(int scenarioId) {
    std::cout << "\n[DAY6-CAMPAIGN] 3D milestone regression at " << ToTxx(scenarioId)
              << " (matrix_audit + issue11 + issue12)\n";
    Run_Day6_MatrixAudit_3D_EDFM();
    Test_Issue11::Run_All();
    Test_Issue12::Run_All();
    std::cout << "[DAY6-CAMPAIGN] 3D milestone regression pass at " << ToTxx(scenarioId) << "\n";
}

void RunCampaign2DAllInternal() {
    std::cout << "\n================ Day6 Campaign 2D (Txx_Fy) ================\n";
    for (const auto& scenario : kPhysicalScenarios) {
        std::cout << "\n[DAY6-CAMPAIGN] Enter " << ToTxx(scenario.id)
                  << " -> require F0/F1/F2 all pass\n";

        for (const auto topology : kTopologyOrder) {
            const auto outcome = RunScenario2D(scenario, topology);
            ThrowIfOutcomeFailed(outcome, 2);
        }

        std::cout << "[DAY6-CAMPAIGN] PASS " << ToTxx(scenario.id)
                  << " (F0+F1+F2)\n";

        if (IsMilestoneScenario(scenario.id)) {
            RunMilestoneRegression2D(scenario.id);
        }
    }
    std::cout << "\n[DAY6-CAMPAIGN] PASS 2D full campaign T01..T16 with F0/F1/F2\n";
}

void RunCampaign3DAllInternal() {
    std::cout << "\n================ Day6 Campaign 3D (Txx_Fy) ================\n";
    for (const auto& scenario : kPhysicalScenarios) {
        std::cout << "\n[DAY6-CAMPAIGN] Enter " << ToTxx(scenario.id)
                  << " -> require F0/F1/F2 all pass\n";

        for (const auto topology : kTopologyOrder) {
            const auto outcome = RunScenario3D(scenario, topology);
            ThrowIfOutcomeFailed(outcome, 3);
        }

        std::cout << "[DAY6-CAMPAIGN] PASS " << ToTxx(scenario.id)
                  << " (F0+F1+F2)\n";

        if (IsMilestoneScenario(scenario.id)) {
            RunMilestoneRegression3D(scenario.id);
        }
    }
    std::cout << "\n[DAY6-CAMPAIGN] PASS 3D full campaign T01..T16 with F0/F1/F2\n";
}

void RunSingle2DOrThrow(int scenarioId, TopologyAxis topology) {
    const auto& scenario = GetScenarioOrThrow(scenarioId);
    const auto outcome = RunScenario2D(scenario, topology);
    ThrowIfOutcomeFailed(outcome, 2);
}

} // namespace


    void Run_Day6_T1_2D_SP_NoWell_Analytical() {
    std::cout << "\n=== Run_Day6_T1_2D_SP_NoWell_Analytical ===\n";

    constexpr double kLx = 400.0;
    constexpr double kLy = 40.0;
    constexpr int kNx = 40;
    constexpr int kNy = 4;

    constexpr double kPInit = 10.0e6;
    constexpr double kPLeft = 12.0e6;
    constexpr double kPRight = 8.0e6;
    constexpr double kTInit = 360.0;

    constexpr double kPhi = 0.10;
    constexpr double kPerm = 1.0e-13;
    constexpr double kCt = 5.0e-9;
    constexpr double kMu = 1.0e-3;

    constexpr double kDtInit = 100.0;
    constexpr double kTFinal = 2.0e5;
    constexpr int kFourierTerms = 200;
    constexpr double kErrorThreshold = 0.03;
    constexpr double kPi = 3.14159265358979323846;

    MeshManager mgr(kLx, kLy, 0.0, kNx, kNy, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset2D();
    preset.enable_rock_region = false;
    preset.rock_bg.phi_r = kPhi;
    preset.rock_bg.kxx = kPerm;
    preset.rock_bg.kyy = kPerm;
    preset.rock_bg.kzz = kPerm;
    preset.rock_bg.compressibility = kCt;
    preset.rock_region = preset.rock_bg;
    preset.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::ConstantWater;

    FIM_Engine::InitialConditions ic;
    ic.P_init = kPInit;
    ic.T_init = kTInit;
    ic.Sw_init = 1.0;

    BoundarySetting::BoundaryConditionManager bcP;
    BoundarySetting::BoundaryConditionManager bcT;
    bcP.Clear();
    bcT.Clear();

    bcP.SetDirichletBC(MeshTags::LEFT, kPLeft);
    bcP.SetDirichletBC(MeshTags::RIGHT, kPRight);
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcP.SetNeumannBC(MeshTags::TOP, 0.0);

    bcT.SetDirichletBC(MeshTags::LEFT, kTInit);
    bcT.SetDirichletBC(MeshTags::RIGHT, kTInit);
    bcT.SetDirichletBC(MeshTags::BOTTOM, kTInit);
    bcT.SetDirichletBC(MeshTags::TOP, kTInit);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, nullptr);

    auto params = FIM_CaseKit::BuildSolverParams(false, 2000, kDtInit);
    params.target_end_time_s = kTFinal;
    params.dt_min = 1.0;
    params.dt_max = 5.0e4;
    params.gravity_vector = Vector(0.0, 0.0, 0.0);
    params.enable_non_orthogonal_correction = false;
    params.enable_control_ramp = false;
    params.enable_ptc = false;
    params.max_newton_iter = 10;
    params.max_dP = 5.0e6;
    params.abs_res_tol = 1.0e-2;
    params.rel_res_tol = 5.0e-2;
    params.rel_update_tol = 1.0e-5;
    params.dt_relres_iter_grow_hi = 6;
    params.dt_relres_iter_neutral_hi = 8;
    params.dt_relres_iter_soft_shrink_hi = 10;
    params.dt_relres_grow_factor = 1.40;
    params.dt_relres_neutral_factor = 1.15;
    params.dt_relres_soft_shrink_factor = 0.95;
    params.dt_relres_hard_shrink_factor = 0.70;
    params.enable_ls_trace = false;
    params.diag_level = FIM_Engine::DiagLevel::Off;
    params.lin_solver = FIM_Engine::LinearSolverType::AMGCL;

    const std::string caseName = "day6_t1_2d_sp_nowell_analytical";
    FIM_Engine::RunGenericFIMTransient<2>(
        caseName,
        mgr,
        fm,
        ic,
        {},
        params,
        FIM_Engine::SolverRoute::FIM,
        modules);

    const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
    auto pField = fm.getMatrixScalar(pCfg.pressure_field);
    if (!pField) {
        throw std::runtime_error("[T1] pressure field not found after transient run.");
    }

    const auto& cells = mgr.mesh().getCells();
    if (pField->data.size() != cells.size()) {
        throw std::runtime_error("[T1] pressure field size mismatch with matrix cells.");
    }

    const double deltaP = kPLeft - kPRight;
    const double deltaPAbs = std::max(std::abs(deltaP), 1.0);
    const double D_h = kPerm / (kMu * kPhi * kCt);
    const double a0 = kPInit - kPLeft;

    auto pressureAnalytical = [&](double x) {
        const double p_ss = kPLeft - deltaP * (x / kLx);
        double phi_xt = 0.0;
        for (int n = 1; n <= kFourierTerms; ++n) {
            const double sign_n = (n % 2 == 0) ? 1.0 : -1.0;
            const double bn = (2.0 / (static_cast<double>(n) * kPi))
                * (a0 * (1.0 - sign_n) - deltaP * sign_n);
            const double lambda_n = static_cast<double>(n) * kPi / kLx;
            phi_xt += bn * std::sin(lambda_n * x)
                * std::exp(-D_h * lambda_n * lambda_n * kTFinal);
        }
        return p_ss + phi_xt;
    };

    const std::string caseDir = "Test/Transient/Day6/" + caseName;
    EnsureDirRecursive(caseDir);

    const std::string csvPath = caseDir + "/analytical_compare.csv";
    std::ofstream csv(csvPath, std::ios::out | std::ios::trunc);
    csv << "cell_id,x,y,p_num,p_ana,abs_err,abs_err_over_dP\n";

    double maxAbsErr = 0.0;
    int maxErrCell = -1;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        const double x = cells[i].center.m_x;
        const double y = cells[i].center.m_y;
        const double pNum = (*pField)[i];
        const double pAna = pressureAnalytical(x);
        const double absErr = std::abs(pNum - pAna);
        const double relErr = absErr / deltaPAbs;

        if (absErr > maxAbsErr) {
            maxAbsErr = absErr;
            maxErrCell = i;
        }

        csv << i << ","
            << std::setprecision(12) << x << ","
            << std::setprecision(12) << y << ","
            << std::setprecision(12) << pNum << ","
            << std::setprecision(12) << pAna << ","
            << std::setprecision(12) << absErr << ","
            << std::setprecision(12) << relErr << "\n";
    }

    const double normInfOverDP = maxAbsErr / deltaPAbs;
    const std::string summaryPath = caseDir + "/analytical_summary.txt";
    std::ofstream summary(summaryPath, std::ios::out | std::ios::trunc);
    summary << "case=" << caseName << "\n";
    summary << "D_h=" << std::setprecision(12) << D_h << "\n";
    summary << "t_final=" << std::setprecision(12) << kTFinal << "\n";
    summary << "max_abs_err=" << std::setprecision(12) << maxAbsErr << "\n";
    summary << "norm_inf_over_dP=" << std::setprecision(12) << normInfOverDP << "\n";
    summary << "threshold=" << std::setprecision(12) << kErrorThreshold << "\n";
    summary << "max_err_cell=" << maxErrCell << "\n";

    std::cout << "[T1] D_h = " << D_h << " m^2/s\n";
    std::cout << "[T1] ||e||_inf/DeltaP = " << normInfOverDP
        << " (threshold=" << kErrorThreshold << ")\n";
    std::cout << "[T1] analytical profile csv: " << csvPath << "\n";
    std::cout << "[T1] analytical summary: " << summaryPath << "\n";

    if (normInfOverDP > kErrorThreshold) {
        std::ostringstream oss;
        oss << "[T1] analytical check failed: ||e||_inf/DeltaP="
            << normInfOverDP << " > " << kErrorThreshold;
        throw std::runtime_error(oss.str());
    }
}
void Run_Day6_Transient_2D_SP_InjProd() {
        std::cout << "\n=== Run_Day6_Transient_2D_SP_InjProd ===\n";

        // -- Feature flags summary --
#ifdef USE_COOLPROP_EOS
        std::cout << "  [Step 1] EOS backend: CoolProp 6.x (Span-Wagner CO2 + IAPWS-95 Water)\n";
#else
        std::cout << "  [Step 1] EOS backend: Table-based (Catmull-Rom + 4-tier FD)\n";
#endif
        std::cout << "  [Step 2] Non-orthogonal correction: ENABLED\n";
        std::cout << "  [Step 3] Well independent DOF: ENABLED (WellDOFManager)\n";
        std::cout << "  [Step 4] Post-processing: VTU/PVD (ParaView XML)\n\n";

        // =====================================================================
        // 1. Mesh & Topology Setup
        // =====================================================================
        MeshManager mgr(100.0, 100.0, 0.0, 10, 10, 0, true, false);
        mgr.BuildSolidMatrixGrid_2D();
        mgr.addFracture(Vector(10.0, 10.0, 0.0), Vector(90.0, 90.0, 0.0));
        mgr.addFracture(Vector(10.0, 90.0, 0.0), Vector(90.0, 10.0, 0.0));
        mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
        mgr.BuildGlobalSystemIndexing();
        mgr.BuildFracturetoFractureTopology();
        mgr.setNumDOFs(2);

        // =====================================================================
        // 2. Field Manager Initialization
        // =====================================================================
        FieldManager_2D fm;
        FIM_CaseKit::InitFieldManager(mgr, fm);

        // =====================================================================
        // 3. Rock/Fracture Properties Setup
        // =====================================================================
        auto preset = FIM_CaseKit::MakeDefaultPropertyPreset2D();
		preset.enable_rock_region = false;  // 2D case with uniform rock properties for simplicity
        preset.rock_bg.phi_r = 0.12;        // 孔隙度
        preset.rock_bg.kxx = 1.0e-13; 
        preset.rock_bg.kyy = 1.0e-13;
		preset.rock_bg.kzz = 1.0e-13;       // 渗透率，单位 m^2，约等于 0.5 mD
		preset.rock_bg.compressibility = 5.0e-9;  // Pa^-1, 约等于 5e-6 MPa^-1，保持适度压缩性以观察压力传播和孔隙体积变化的影响
		preset.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::Water; // 使用 CO2 作为单相流体，具有较高的可压缩性和温度敏感性，能够更明显地展示压力和温度变化的耦合效应

        // =====================================================================
        // 4. Initial Conditions (IC)
        // =====================================================================
        FIM_Engine::InitialConditions ic;
        ic.P_init = 30.0e6;    // 30 MPa  (深储层，>Pc 裕量×4)
        ic.T_init = 380.0;     // 380 K   (≈107°C，深3km地热，>Tc 裕量+76K)
        ic.Sw_init = 1.0;

        // =====================================================================
        // 5. External BC (Dirichlet / Neumann / Robin)
        // =====================================================================
		BoundarySetting::BoundaryConditionManager bcP;  // 压力 BC 管理器
		BoundarySetting::BoundaryConditionManager bcT;  // 温度 BC 管理器
        bcP.Clear();
        bcT.Clear();

        // Pressure: closed (no-flow) on all sides
        bcP.SetNeumannBC(MeshTags::LEFT, 0.0);          // 压力左边界
		bcP.SetNeumannBC(MeshTags::RIGHT, 0.0);         // 压力右边界
		bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);        // 压力底边界
		bcP.SetNeumannBC(MeshTags::TOP, 0.0);           // 压力上边界

        // Temperature: adiabatic on all sides
        bcT.SetNeumannBC(MeshTags::LEFT, 0.0);          // 左边界    
		bcT.SetNeumannBC(MeshTags::RIGHT, 0.0);         // 右边界
		bcT.SetNeumannBC(MeshTags::BOTTOM, 0.0);        // 底边界
		bcT.SetNeumannBC(MeshTags::TOP, 0.0);           // 上边界

		auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, nullptr); // 创建模型

        // =====================================================================
        // 6. Wells Setup
        // =====================================================================
        const auto& cells = mgr.mesh().getCells();
        auto findNearestCell = [&](double x, double y) {
            int best = 0;
            double bestD2 = std::numeric_limits<double>::max();
            for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
                const double dx = cells[i].center.m_x - x;
                const double dy = cells[i].center.m_y - y;
                const double d2 = dx * dx + dy * dy;
                if (d2 < bestD2) {
                    bestD2 = d2;
                    best = i;
                }
            }
            return best;
            };

        const int injCell = findNearestCell(15.0, 15.0);   // 移至单元中心(i=1,j=1)，避开y=10单元边界，且位于对角裂缝路径上
        const int prodCell = findNearestCell(85.0, 95.0);  // 明确落于顶行单元中心

        WellScheduleStep w1, w2;
        w1.well_name = "INJ_CO2_2D";
        w1.domain = WellTargetDomain::Matrix;
        w1.control_mode = WellControlMode::BHP;
        w1.target_value = 31.0e6;  // 31 MPa
        w1.component_mode = WellComponentMode::Water;
        w1.completion_id = injCell;
        w1.frac_w = 1.0;
        w1.frac_g = 0.0;
        w1.injection_temperature = 340.0;   // 340 K，比储层低 40K，观察冷注降温效果
        w1.injection_is_co2 = false;


        w2.well_name = "PROD_CO2_2D";
        w2.domain = WellTargetDomain::Matrix;
        w2.control_mode = WellControlMode::BHP;
        w2.target_value = 29.0e6;  // 29 MPa
        w2.component_mode = WellComponentMode::Water;
        w2.completion_id = prodCell;
        w2.frac_w = 1.0;
        w2.frac_g = 0.0;
        

        // =====================================================================
        // 7. Solver Parameters Setup
        // =====================================================================
        const double horizon_years = 1.0;
        // Startup: 30天，让温度斜坡在startup内完整走完，避免切换长阶段时的28.8K硬跳变
        const double startup_days = 30.0;
        const double startup_vtk_interval_days = 0.5;   // 30天startup内每5天出一帧
        const double long_vtk_interval_days = 5.0;

        auto params = FIM_CaseKit::BuildLongRunTemplate(
            FIM_CaseKit::LongRunScenario::S2D_SP,
            horizon_years,
            startup_days);

        // Explicit local controls for long-run setup in this case.
        params.target_end_time_s = horizon_years * 365.0 * 86400.0;
        params.startup_end_time_s = startup_days * 86400.0;
        params.startup_vtk_output_interval_s = startup_vtk_interval_days * 86400.0;
        params.long_vtk_output_interval_s = long_vtk_interval_days * 86400.0;
        params.gravity_vector = Vector(0.0, 0.0, 0.0);

        // Keep startup short and allow dt growth when iter~14.
		params.startup_profile.max_newton_iter = 18; // Startup阶段允许更深入的非线性迭代，以更快地适应初始条件和边界条件，促进时间步长增长
		params.startup_profile.ptc_lambda_init = 1.0; // 初始PTC lambda较大，提供更强的正则化以稳定初始迭代，促进时间步长增长
        params.startup_profile.ptc_lambda_decay = 0.70;
        params.startup_profile.ptc_lambda_min = 0.02;
        params.startup_profile.dt_relres_iter_grow_hi = 12;
        params.startup_profile.dt_relres_iter_neutral_hi = 18;
        params.startup_profile.dt_relres_iter_soft_shrink_hi = 24;
        params.startup_profile.dt_relres_grow_factor = 1.10;
        params.startup_profile.dt_relres_neutral_factor = 1.02;
        params.startup_profile.dt_relres_soft_shrink_factor = 0.995;
        params.startup_profile.dt_relres_hard_shrink_factor = 0.92;
        // 100步完成斜坡（~4天 @ dt_max=3600s），远早于30天startup结束，消除切换时硬跳变
        params.startup_profile.control_ramp_steps = 10;
		params.startup_profile.rel_res_tol = 1e-2; // [Fix-A] 放宽至1e-2：dt_init=1e-3时growth~3.77e-3>1e-3导致NL-DIVERGED；1e-2允许rel_res_update触发，不影响dt=3600路径


        params.long_profile.max_newton_iter = 36; // Allow deeper nonlinear progress in late stiff periods.
        params.long_profile.rel_update_tol = 1e-3;
        params.long_profile.ptc_lambda_init = 0.5;
        params.long_profile.ptc_lambda_decay = 0.7;
        params.long_profile.ptc_lambda_min = 5.0e-2; // 加强long stage对角正则化，防止能量方程在高刚性时失稳
        params.long_profile.dt_relres_iter_grow_hi = 8;
        params.long_profile.dt_relres_iter_neutral_hi = 14;
        params.long_profile.dt_relres_iter_soft_shrink_hi = 18;
        params.long_profile.dt_relres_grow_factor = 1.05;
        params.long_profile.dt_relres_neutral_factor = 1.02;
        params.long_profile.dt_relres_soft_shrink_factor = 0.90;
        params.long_profile.dt_relres_hard_shrink_factor = 0.90;
        params.rollback_shrink_factor = 0.85;
        params.clamp_state_to_eos_bounds = true; // Keep single-phase CO2 state inside EOS table bounds.
        params.long_profile.rel_res_tol = 1e-2;
        
        // 降低Newton溫度更新上限：预设2K/iter过大，在偽临界區域（?h/?T極大）會過衝
        params.max_dT = 5;   // 从 0.5 放宽到 5K/iter，远离拟临界不再需要极小步

        // Issue#12: 切换至 CPR-AMG 求解器，验证全流程稳定性
        params.lin_solver = FIM_Engine::LinearSolverType::AMGCL_CPR;

        // Step 2: Enable deferred non-orthogonal correction (T-vector gradient term).
        // On this 10x10 orthogonal grid |vectorT|≈0, so the correction has near-zero effect
        // but exercises the full code path. Use skewed meshes for accuracy improvement.
        params.enable_non_orthogonal_correction = true;

        // =====================================================================
        // 8. Execution
        // =====================================================================
        // Step 3: Wells w1, w2 are passed to RunGenericFIMTransient, which internally
        // invokes WellDOFManager to add independent P_wbh DOF blocks per well.
        // The sparse matrix is extended to (nMatrix + nFrac + nWells) blocks;
        // well equations (Peaceman model) are assembled each Newton iteration.
        FIM_Engine::RunGenericFIMTransient<2>(
            "day6_transient_2d_sp_injprod", mgr, fm, ic, { w1, w2 }, params, FIM_Engine::SolverRoute::FIM, modules);

        // =====================================================================
        // 9. VTU/PVD Post-Processing (Step 4)
        // =====================================================================
        // Export final state as ParaView XML Unstructured Grid (.vtu).
        // RunGenericFIMTransient already synced field data at the last VTK snapshot,
        // so fm contains the final-time solution.
        {
            PostProcess_2D pp(mgr, fm);
            const std::string baseDir = "Test/Transient/Day6/day6_transient_2d_sp_injprod/";
            const std::string vtu_file = baseDir + "final_state.vtu";
            pp.ExportVTU(vtu_file, params.target_end_time_s);
            std::cout << "    [VTU] Exported: " << vtu_file << "\n";

            // PVD time-series collection (single-snapshot demo).
            // In production, RunGenericFIMTransient would collect all intermediate
            // VTU snapshots and produce a complete .pvd for ParaView animation.
            const std::string pvd_file = baseDir + "timeseries.pvd";
            PostProcess_2D::ExportPVD(pvd_file,
                                      { vtu_file },
                                      { params.target_end_time_s });
            std::cout << "    [PVD] Exported: " << pvd_file << "\n";
        }
    }

void Run_Day6_Transient_2D_TP_InjProd() {
    // =====================================================================
    // 1. Mesh & Topology Setup
    // =====================================================================
    MeshManager mgr(10.0, 10.0, 0.0, 10, 10, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    mgr.addFracture(Vector(1, 1, 0), Vector(9, 9, 0));
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(3);

    // =====================================================================
    // 2. Field Manager Initialization
    // =====================================================================
    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    // =====================================================================
    // 3. Rock/Fracture Properties Setup
    // =====================================================================
    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset2D();

    // =====================================================================
    // 4. Initial Conditions (IC)
    // =====================================================================
    FIM_Engine::InitialConditions ic;
    ic.P_init = 25.0e6;
    ic.T_init = 360.0;
    ic.Sw_init = 0.25;

    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    BoundarySetting::BoundaryConditionManager bcP;
    BoundarySetting::BoundaryConditionManager bcT;
    bcP.Clear();
    bcT.Clear();

    // Pressure: closed (no-flow) on all sides
    bcP.SetNeumannBC(MeshTags::LEFT, 0.0);
    bcP.SetNeumannBC(MeshTags::RIGHT, 0.0);
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcP.SetNeumannBC(MeshTags::TOP, 0.0);

    // Temperature: adiabatic on all sides
    bcT.SetNeumannBC(MeshTags::LEFT, 0.0);
    bcT.SetNeumannBC(MeshTags::RIGHT, 0.0);
    bcT.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcT.SetNeumannBC(MeshTags::TOP, 0.0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, nullptr);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    const int nMat = FIM_CaseKit::MatrixBlockCount(mgr);
    const int injCell = nMat / 5;
    const int prodCell = (4 * nMat) / 5;

    WellScheduleStep w1, w2;
    w1.well_name = "INJ_2D_TP";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::BHP;
    w1.target_value = 25.15e6;
    w1.component_mode = WellComponentMode::Water;
    w1.injection_temperature = 330.0;
    w1.completion_id = injCell;
    w1.frac_w = 1.0;
    w1.frac_g = 0.0;

    w2.well_name = "PROD_2D_TP";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 24.90e6;
    w2.component_mode = WellComponentMode::Total;
    w2.completion_id = prodCell;
    // Keep frac_w/frac_g unset for Total producer:
    // BoundaryAssembler now auto-splits by local mobility.

    // =====================================================================
    // 7. Solver Parameters Setup
    // =====================================================================
    auto params = FIM_CaseKit::BuildSolverParams(true);

    // =====================================================================
    // 8. Execution
    // =====================================================================
    FIM_Engine::RunGenericFIMTransient<3>(
        "day6_transient_2d_tp_injprod", mgr, fm, ic, { w1, w2 }, params, FIM_Engine::SolverRoute::FIM, modules);
}

void Run_Day6_Transient_2D_TP_Multiwell() {
    // =====================================================================
    // 1. Mesh & Topology Setup
    // =====================================================================
    MeshManager mgr(10.0, 10.0, 0.0, 10, 10, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    mgr.addFracture(Vector(1, 1, 0), Vector(9, 9, 0));
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(3);

    // =====================================================================
    // 2. Field Manager Initialization
    // =====================================================================
    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    // =====================================================================
    // 3. Rock/Fracture Properties Setup
    // =====================================================================
    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset2D();

    // =====================================================================
    // 4. Initial Conditions (IC)
    // =====================================================================
    FIM_Engine::InitialConditions ic;
    ic.P_init = 25.0e6;
    ic.T_init = 360.0;
    ic.Sw_init = 0.25;

    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    // Closed boundary for Day6 robustness: no external BC injection into solver.
    auto modules = FIM_CaseKit::BuildModules2D(preset, nullptr, nullptr, nullptr);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    const int nMat = FIM_CaseKit::MatrixBlockCount(mgr);
    const int injCell = nMat / 5;
    const int prodCell = (4 * nMat) / 5;
    const int midCell = nMat / 2;

    WellScheduleStep w1, w2, w3;
    w1.well_name = "INJ_W";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::Rate;
    w1.target_value = -0.04;
    w1.component_mode = WellComponentMode::Water;
    w1.injection_temperature = 330.0;
    w1.completion_id = injCell;
    w1.frac_w = 1.0;
    w1.frac_g = 0.0;

    w2.well_name = "PROD_1";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 24.92e6;
    w2.component_mode = WellComponentMode::Total;
    w2.completion_id = prodCell;

    w3.well_name = "PROD_2";
    w3.domain = WellTargetDomain::Matrix;
    w3.control_mode = WellControlMode::BHP;
    w3.target_value = 24.92e6;
    w3.component_mode = WellComponentMode::Total;
    w3.completion_id = midCell;

    // =====================================================================
    // 7. Solver Parameters Setup
    // =====================================================================
    auto params = FIM_CaseKit::BuildSolverParams(true);

    // =====================================================================
    // 8. Execution
    // =====================================================================
    FIM_Engine::RunGenericFIMTransient<3>(
        "day6_transient_2d_tp_multiwell", mgr, fm, ic, { w1, w2, w3 }, params, FIM_Engine::SolverRoute::FIM, modules);
}

void Run_Day6_Transient_3D_SP_InjProd() {
    // =====================================================================
    // 1. Mesh & Topology Setup
    // =====================================================================
    MeshManager_3D mgr(10.0, 10.0, 10.0, 20, 20, 20, true, false);
    mgr.BuildSolidMatrixGrid_3D();
    std::vector<Vector> pts = { Vector(1, 1, 1), Vector(9, 1, 1), Vector(9, 9, 1), Vector(1, 9, 1) };
    mgr.addFracturetoFractureNetwork(Fracture_2D(0, pts));
    mgr.meshAllFracturesinNetwork(2, 2);

    mgr.setupGlobalIndices();
    mgr.DetectFractureFractureIntersectionsInNetwork();
    mgr.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);
    mgr.fracture_network().rebuildEdgeProperties();
    mgr.removeDuplicateInteractions();
    mgr.resolveCoplanarInteractions();
    mgr.buildTopologyMaps();
    mgr.setNumDOFs(2);

    // =====================================================================
    // 2. Field Manager Initialization
    // =====================================================================
    FieldManager_3D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    // =====================================================================
    // 3. Rock/Fracture Properties Setup
    // =====================================================================
    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset3D();
    preset.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::CO2;

    // =====================================================================
    // 4. Initial Conditions (IC)
    // =====================================================================
    FIM_Engine::InitialConditions ic;
    ic.P_init = 25.0e6;
    ic.T_init = 360.0;
    ic.Sw_init = 1.0;

    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    // Closed boundary for Day6 robustness: no external BC injection into solver.
    auto modules = FIM_CaseKit::BuildModules3D(preset, nullptr, nullptr, nullptr);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    const int nMat = FIM_CaseKit::MatrixBlockCount(mgr);
    const int injCell = nMat / 5;
    const int prodCell = (4 * nMat) / 5;

    WellScheduleStep w1, w2;
    w1.well_name = "INJ_3D";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::Rate;
    w1.target_value = -0.08;
    w1.component_mode = WellComponentMode::Gas;
    w1.injection_temperature = 330.0;
    w1.injection_is_co2 = true;
    w1.completion_id = injCell;
    w1.frac_w = 0.0;
    w1.frac_g = 1.0;
    w1.well_axis = WellAxis::Z;

    w2.well_name = "PROD_3D";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 24.90e6;
    w2.component_mode = WellComponentMode::Gas;
    w2.completion_id = prodCell;
    w2.frac_w = 0.0;
    w2.frac_g = 1.0;
    w2.well_axis = WellAxis::Z;

    // =====================================================================
    // 7. Solver Parameters Setup
    // =====================================================================
    auto params = FIM_CaseKit::BuildSolverParams(false);

    // =====================================================================
    // 8. Execution
    // =====================================================================
    FIM_Engine::RunGenericFIMTransient<2>(
        "day6_transient_3d_sp_injprod", mgr, fm, ic, { w1, w2 }, params, FIM_Engine::SolverRoute::FIM, modules);
}

void Run_Day6_Transient_3D_TP_InjProd() {
    // =====================================================================
    // 1. Mesh & Topology Setup
    // =====================================================================
    MeshManager_3D mgr(10.0, 10.0, 10.0, 10, 10, 10, true, false);
    mgr.BuildSolidMatrixGrid_3D();
    std::vector<Vector> pts = { Vector(1, 1, 1), Vector(9, 1, 1), Vector(9, 9, 1), Vector(1, 9, 1) };
    mgr.addFracturetoFractureNetwork(Fracture_2D(0, pts));
    mgr.meshAllFracturesinNetwork(2, 2);

    mgr.setupGlobalIndices();
    mgr.DetectFractureFractureIntersectionsInNetwork();
    mgr.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);
    mgr.fracture_network().rebuildEdgeProperties();
    mgr.removeDuplicateInteractions();
    mgr.resolveCoplanarInteractions();
    mgr.buildTopologyMaps();
    mgr.setNumDOFs(3);

    // =====================================================================
    // 2. Field Manager Initialization
    // =====================================================================
    FieldManager_3D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    // =====================================================================
    // 3. Rock/Fracture Properties Setup
    // =====================================================================
    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset3D();
    preset.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::CO2;

    // =====================================================================
    // 4. Initial Conditions (IC)
    // =====================================================================
    FIM_Engine::InitialConditions ic;
    ic.P_init = 25.0e6;
    ic.T_init = 360.0;
    ic.Sw_init = 0.25;


    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    // Closed boundary for Day6 robustness: no external BC injection into solver.
    auto modules = FIM_CaseKit::BuildModules3D(preset, nullptr, nullptr, nullptr);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    const int nMat = FIM_CaseKit::MatrixBlockCount(mgr);
    const int injCell = nMat / 5;
    const int prodCell = (4 * nMat) / 5;

    WellScheduleStep w1, w2;
    w1.well_name = "INJ_3D_TP";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::BHP;
    w1.target_value = 25.15e6;
    w1.component_mode = WellComponentMode::Water;
    w1.injection_temperature = 330.0;
    w1.completion_id = injCell;
    w1.frac_w = 1.0;
    w1.frac_g = 0.0;
    w1.well_axis = WellAxis::Z;

    w2.well_name = "PROD_3D_TP";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 24.90e6;
    w2.component_mode = WellComponentMode::Total;
    w2.completion_id = prodCell;
    w2.well_axis = WellAxis::Z;


    // =====================================================================
    // 7. Solver Parameters Setup
    // =====================================================================
    auto params = FIM_CaseKit::BuildSolverParams(true);

    // =====================================================================
    // 8. Execution
    // =====================================================================
    FIM_Engine::RunGenericFIMTransient<3>(
        "day6_transient_3d_tp_injprod", mgr, fm, ic, { w1, w2 }, params, FIM_Engine::SolverRoute::FIM, modules);
}

void Run_Day6_Transient_3D_TP_Multiwell() {
    // =====================================================================
    // 1. Mesh & Topology Setup
    // =====================================================================
    MeshManager_3D mgr(10.0, 10.0, 10.0, 10, 10, 10, true, false);
    mgr.BuildSolidMatrixGrid_3D();
    std::vector<Vector> pts = { Vector(1, 1, 1), Vector(9, 1, 1), Vector(9, 9, 1), Vector(1, 9, 1) };
    mgr.addFracturetoFractureNetwork(Fracture_2D(0, pts));
    mgr.meshAllFracturesinNetwork(2, 2);

    mgr.setupGlobalIndices();
    mgr.DetectFractureFractureIntersectionsInNetwork();
    mgr.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);
    mgr.fracture_network().rebuildEdgeProperties();
    mgr.removeDuplicateInteractions();
    mgr.resolveCoplanarInteractions();
    mgr.buildTopologyMaps();
    mgr.setNumDOFs(3);

    // =====================================================================
    // 2. Field Manager Initialization
    // =====================================================================
    FieldManager_3D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    // =====================================================================
    // 3. Rock/Fracture Properties Setup
    // =====================================================================
    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset3D();
    preset.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::CO2;

    // =====================================================================
    // 4. Initial Conditions (IC)
    // =====================================================================
    FIM_Engine::InitialConditions ic;
    ic.P_init = 25.0e6;
    ic.T_init = 360.0;
    ic.Sw_init = 0.25;

    // =====================================================================
    // 5. External BC (Dirichlet / Neumann / Robin)
    // =====================================================================
    // Closed boundary for Day6 robustness: no external BC injection into solver.
    auto modules = FIM_CaseKit::BuildModules3D(preset, nullptr, nullptr, nullptr);

    // =====================================================================
    // 6. Wells Setup
    // =====================================================================
    const int nMat = FIM_CaseKit::MatrixBlockCount(mgr);
    const int injCell = nMat / 5;
    const int prodCell = (4 * nMat) / 5;
    const int midCell = nMat / 2;

    WellScheduleStep w1, w2, w3;
    w1.well_name = "INJ_W_3D";
    w1.domain = WellTargetDomain::Matrix;
    w1.control_mode = WellControlMode::Rate;
    w1.target_value = -0.08;
    w1.component_mode = WellComponentMode::Water;
    w1.injection_temperature = 330.0;
    w1.completion_id = injCell;
    w1.frac_w = 1.0;
    w1.frac_g = 0.0;
    w1.well_axis = WellAxis::Z;

    w2.well_name = "PROD_1_3D";
    w2.domain = WellTargetDomain::Matrix;
    w2.control_mode = WellControlMode::BHP;
    w2.target_value = 24.92e6;
    w2.component_mode = WellComponentMode::Total;
    w2.completion_id = prodCell;
    w2.well_axis = WellAxis::Z;

    w3.well_name = "PROD_2_3D";
    w3.domain = WellTargetDomain::Matrix;
    w3.control_mode = WellControlMode::BHP;
    w3.target_value = 24.92e6;
    w3.component_mode = WellComponentMode::Total;
    w3.completion_id = midCell;
    w3.well_axis = WellAxis::Z;

    // =====================================================================
    // 7. Solver Parameters Setup
    // =====================================================================
    auto params = FIM_CaseKit::BuildSolverParams(true);

    // =====================================================================
    // 8. Execution
    // =====================================================================
    FIM_Engine::RunGenericFIMTransient<3>(
        "day6_transient_3d_tp_multiwell", mgr, fm, ic, { w1, w2, w3 }, params, FIM_Engine::SolverRoute::FIM, modules);
}

// ---------- 新增 矩阵审计 测试代码 ----------
void Run_Day6_MatrixAudit_2D_EDFM() {
    MeshManager mgr(100.0, 100.0, 0.0, 10, 10, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    mgr.addFracture(Vector(10.0, 10.0, 0.0), Vector(90.0, 90.0, 0.0));
    mgr.addFracture(Vector(10.0, 90.0, 0.0), Vector(90.0, 10.0, 0.0));
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset2D();
    FIM_Engine::InitialConditions ic;
    auto modules = FIM_CaseKit::BuildModules2D(preset, nullptr, nullptr, nullptr);

    auto params = FIM_CaseKit::BuildSolverParams(false, 1, 1.0);
    params.enable_matrix_audit = true;
    params.matrix_audit_strict = true;
    params.matrix_audit_step = 1;
    params.matrix_audit_iter = 1;
    params.max_steps = 1;
    params.max_newton_iter = 1;
    params.abs_res_tol = 1e30; // 仅过一次残差及装配，不管牛顿收敛
    params.gravity_vector = Vector(0.0, 0.0, 0.0);

    FIM_Engine::RunGenericFIMTransient<2>(
        "day6_matrix_audit_2d_edfm", mgr, fm, ic, {}, params, FIM_Engine::SolverRoute::FIM, modules);
}

void Run_Day6_MatrixAudit_3D_EDFM() {
    MeshManager_3D mgr(10.0, 10.0, 10.0, 10, 10, 10, true, false);
    mgr.BuildSolidMatrixGrid_3D();
    std::vector<Vector> pts1 = { Vector(1, 1, 1), Vector(9, 1, 1), Vector(9, 9, 9), Vector(1, 9, 9) };
    std::vector<Vector> pts2 = { Vector(1, 9, 1), Vector(9, 9, 1), Vector(9, 1, 9), Vector(1, 1, 9) };
    mgr.addFracturetoFractureNetwork(Fracture_2D(0, pts1));
    mgr.addFracturetoFractureNetwork(Fracture_2D(1, pts2));
    mgr.meshAllFracturesinNetwork(2, 2);

    mgr.setupGlobalIndices();
    mgr.DetectFractureFractureIntersectionsInNetwork();
    mgr.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);
    mgr.fracture_network().rebuildEdgeProperties();
    mgr.removeDuplicateInteractions();
    mgr.resolveCoplanarInteractions();
    mgr.buildTopologyMaps();
    mgr.setNumDOFs(2);

    FieldManager_3D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = FIM_CaseKit::MakeDefaultPropertyPreset3D();
    preset.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::CO2;
    FIM_Engine::InitialConditions ic;
    auto modules = FIM_CaseKit::BuildModules3D(preset, nullptr, nullptr, nullptr);

    auto params = FIM_CaseKit::BuildSolverParams(false, 1, 1.0);
    params.enable_matrix_audit = true;
    params.matrix_audit_strict = true;
    params.matrix_audit_step = 1;
    params.matrix_audit_iter = 1;
    params.max_steps = 1;
    params.max_newton_iter = 1;
    params.abs_res_tol = 1e30;

    FIM_Engine::RunGenericFIMTransient<2>(
        "day6_matrix_audit_3d_edfm", mgr, fm, ic, {}, params, FIM_Engine::SolverRoute::FIM, modules);
}

void Run_Day6_Campaign_2D_All_TxxFy() {
    RunCampaign2DAllInternal();
}

void Run_Day6_Campaign_3D_All_TxxFy() {
    RunCampaign3DAllInternal();
}

void Run_Day6_Campaign_2D_T01_F0() {
    RunSingle2DOrThrow(1, TopologyAxis::F0);
}

void Run_Day6_Campaign_2D_T01_F1() {
    RunSingle2DOrThrow(1, TopologyAxis::F1);
}

void Run_Day6_Campaign_2D_T01_F2() {
    RunSingle2DOrThrow(1, TopologyAxis::F2);
}

} // namespace Test_Day6



#pragma once
#include <string>

namespace FullCaseTest {

enum class TopologyVariant {
    NoFrac = 0,
    SingleFrac = 1,
    CrossFrac = 2
};

// 2D topology dispatcher (mapped to migrated validation const variants).
void Run2DCase(int scenarioId, TopologyVariant topology);

// 3D dispatcher entrypoint (migrated through FullCase dispatcher).
void Run3DCase(int scenarioId, TopologyVariant topology);

void Run2DAll();
void Run3DAll();
void RunAll();

// N=1 pressure-only registry-backed presets.
void RunN1L1ConstNoWellAnalytical();
void RunN1L2ConstNoWellGrid();
void RunN1L3ConstNoWellSolver();
void RunN1L4VarPropNoWell();
void RunN1L4VarPropNoWellSingleFrac();
void RunN1LadderAll();

void RunN1L1ConstNoWellAnalyticalLegacy();
void RunN1L2ConstNoWellGridLegacy();
void RunN1L3ConstNoWellSolverLegacy();
void RunN1L4VarPropNoWellLegacy();
void RunN1LadderAllLegacy();

// Reusable N=1 template bodies (copy these for new N=1 scenarios).
void RunN1TemplateConstNoWellNoFrac();
void RunN1TemplateConstNoWellSingleFrac();
void RunN1TemplateConstNoWellCrossFrac();

// Migrated N=2/N=3 validation entrypoints.
void RunValidationCaseByKey(const std::string& case_key);
void RunValidationGroup(int t_id);   // t_id in [1,8]
void RunValidationAllMigrated();

} // namespace FullCaseTest

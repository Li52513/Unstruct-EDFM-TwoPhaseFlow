#include "Test_2D_EDFM_H_T_CO2_VaryPP_NoFrac_NoWell.h"

#include "Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_NoWell.h"

namespace Test_H_T_CO2_VaryPP_NoFrac {
namespace {

constexpr const char* kPlanKey = "h_t_co2_varypp_nofrac_nowell";
constexpr const char* kSmokePlanKey = "h_t_co2_varypp_nofrac_nowell_smoke";

} // namespace

void RunTestCase() {
    RunFullWorkflow();
}

void ExecutePlanByKey(const std::string& key) {
    Test_H_T_CO2_ConstPP_NoFrac::ExecutePlanByKey(key);
}

void RunSolveOnly() {
    Test_H_T_CO2_ConstPP_NoFrac::RunStageByKey(kSmokePlanKey, CaseCommon::CaseStage::SolveOnly);
}

void RunPrepareReference() {
    Test_H_T_CO2_ConstPP_NoFrac::RunStageByKey(kPlanKey, CaseCommon::CaseStage::PrepareReference);
}

void RunValidateOnly() {
    Test_H_T_CO2_ConstPP_NoFrac::RunStageByKey(kPlanKey, CaseCommon::CaseStage::ValidateOnly);
}

void RunFullWorkflow() {
    Test_H_T_CO2_ConstPP_NoFrac::RunStageByKey(kPlanKey, CaseCommon::CaseStage::FullWorkflow);
}

} // namespace Test_H_T_CO2_VaryPP_NoFrac

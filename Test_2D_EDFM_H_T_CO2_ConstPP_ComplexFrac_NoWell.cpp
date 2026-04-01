#include "Test_2D_EDFM_H_T_CO2_ConstPP_ComplexFrac_NoWell.h"

#include "Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.h"

namespace Test_H_T_CO2_ConstPP_ComplexFrac {
namespace {

constexpr const char* kPlanKey = "h_t_co2_constpp_complexfrac_nowell";
constexpr const char* kSmokeKey = "h_t_co2_constpp_complexfrac_nowell_smoke";
constexpr const char* kPrepareKey = "h_t_co2_constpp_complexfrac_nowell_prepare_reference";
constexpr const char* kValidateKey = "h_t_co2_constpp_complexfrac_nowell_validate_reference";
constexpr const char* kFullKey = "h_t_co2_constpp_complexfrac_nowell_full_comsol";

} // namespace

void RunTestCase() {
    RunFullWorkflow();
}

void ExecutePlanByKey(const std::string& key) {
    Test_H_T_CO2_ConstPP_SingleFrac::ExecutePlanByKey(key);
}

void RunSolveOnly() {
    Test_H_T_CO2_ConstPP_SingleFrac::ExecutePlanByKey(kSmokeKey);
}

void RunPrepareReference() {
    Test_H_T_CO2_ConstPP_SingleFrac::ExecutePlanByKey(kPrepareKey);
}

void RunValidateOnly() {
    Test_H_T_CO2_ConstPP_SingleFrac::ExecutePlanByKey(kValidateKey);
}

void RunFullWorkflow() {
    Test_H_T_CO2_ConstPP_SingleFrac::ExecutePlanByKey(kFullKey);
}

} // namespace Test_H_T_CO2_ConstPP_ComplexFrac

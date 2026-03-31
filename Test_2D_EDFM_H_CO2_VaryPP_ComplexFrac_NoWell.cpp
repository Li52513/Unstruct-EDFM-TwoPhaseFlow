/**
 * @file Test_2D_EDFM_H_CO2_VaryPP_ComplexFrac_NoWell.cpp
 * @brief Standalone test: 2D single-phase CO2 variable-property (EOS), complex-fracture, no-well.
 */

#include "Test_2D_EDFM_H_CO2_VaryPP_ComplexFrac_NoWell.h"

#include "Test_2D_EDFM_H_CO2_VaryPP_SingleFrac_NoWell.h"

namespace Test_H_CO2_VaryPP_ComplexFrac {
namespace {

constexpr const char* kPlanKey = "h_co2_varypp_complexfrac_nowell";

} // namespace

void RunTestCase() {
    Test_H_CO2_VaryPP_SingleFrac::ExecutePlanByKey(kPlanKey);
}

void ExecutePlanByKey(const std::string& key) {
    Test_H_CO2_VaryPP_SingleFrac::ExecutePlanByKey(key);
}

void RunSolveOnly() {
    Test_H_CO2_VaryPP_SingleFrac::ExecutePlanByKey(kPlanKey);
}

void RunPrepareReference() {
    Test_H_CO2_VaryPP_SingleFrac::ExecutePlanByKey(kPlanKey);
}

void RunValidateOnly() {
    Test_H_CO2_VaryPP_SingleFrac::ExecutePlanByKey(kPlanKey);
}

void RunFullWorkflow() {
    Test_H_CO2_VaryPP_SingleFrac::ExecutePlanByKey(kPlanKey);
}

} // namespace Test_H_CO2_VaryPP_ComplexFrac

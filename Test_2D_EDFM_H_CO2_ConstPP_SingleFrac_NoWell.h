#pragma once

#include <string>

/**
 * @file Test_2D_EDFM_H_CO2_ConstPP_SingleFrac_NoWell.h
 * @brief Standalone test: 2D single-phase CO2 constant-property, single-fracture, no-well.
 */

namespace Test_H_CO2_ConstPP_SingleFrac {

void RunTestCase();
void ExecutePlanByKey(const std::string& key);
void RunSolveOnly();
void RunPrepareReference();
void RunValidateOnly();
void RunFullWorkflow();

} // namespace Test_H_CO2_ConstPP_SingleFrac

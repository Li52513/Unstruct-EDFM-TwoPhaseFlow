#pragma once

/**
 * @file  Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.h
 * @brief Standalone analytical-validation test: 2D single-phase CO2, constant properties,
 *        no-fracture, no-well pressure diffusion.
 */

#include <string>

namespace Test_H_CO2_ConstPP {

/// Single public entry routed by the internal registry.
void RunTestCase();
void ExecutePlanByKey(const std::string& key);
void RunSolveOnly();
void RunPrepareReference();
void RunValidateOnly();
void RunFullWorkflow();

} // namespace Test_H_CO2_ConstPP

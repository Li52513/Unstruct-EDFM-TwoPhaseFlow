#pragma once

#include <string>

/**
 * @file Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.h
 * @brief Standalone test: 2D single-phase CO2 constant-property P-T coupled, single-fracture, no-well.
 */

namespace Test_H_T_CO2_ConstPP_SingleFrac {

void RunTestCase();
void ExecutePlanByKey(const std::string& key);

} // namespace Test_H_T_CO2_ConstPP_SingleFrac

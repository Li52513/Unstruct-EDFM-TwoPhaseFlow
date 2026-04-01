#pragma once

#include "CaseCommon_Catalog.h"

#include <string>

/**
 * @file Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_NoWell.h
 * @brief Standalone test: 2D single-phase CO2 constant-property P-T coupled, no-fracture, no-well.
 */

namespace Test_H_T_CO2_ConstPP_NoFrac {

void RunTestCase();
void ExecutePlanByKey(const std::string& key);
void RunStageByKey(const std::string& key, CaseCommon::CaseStage stage);
void RunSolveOnly();
void RunPrepareReference();
void RunValidateOnly();
void RunFullWorkflow();

} // namespace Test_H_T_CO2_ConstPP_NoFrac

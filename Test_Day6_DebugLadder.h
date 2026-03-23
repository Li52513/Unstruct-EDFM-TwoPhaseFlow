/**
 * @file Test_Day6_DebugLadder.h
 * @brief Day6 debug ladder: pressure-only (N=1), single-phase CO2 (constant/variable properties), no wells.
 */

#pragma once

namespace Test_Day6 {

void Run_Day6L1_2D_SP_CO2Const_NoWell_Analytical();
void Run_Day6L2_2D_SP_CO2Const_NoWell_GridConvergence();
void Run_Day6L3_2D_SP_CO2Const_NoWell_SolverRobustness();
void Run_Day6L4_2D_SP_CO2VarProp_NoWell_Nonlinear();
void Run_Day6Ladder_2D_SP_CO2Const_NoWell_All();
void Run_Day6L1_2D_SP_CO2Const_NoWell_Analytical_Legacy();
void Run_Day6L2_2D_SP_CO2Const_NoWell_GridConvergence_Legacy();
void Run_Day6L3_2D_SP_CO2Const_NoWell_SolverRobustness_Legacy();
void Run_Day6L4_2D_SP_CO2VarProp_NoWell_Nonlinear_Legacy();
void Run_Day6Ladder_2D_SP_CO2Const_NoWell_All_Legacy();

} // namespace Test_Day6

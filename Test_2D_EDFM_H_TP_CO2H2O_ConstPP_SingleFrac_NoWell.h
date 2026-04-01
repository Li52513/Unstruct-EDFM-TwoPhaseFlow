#pragma once

#include <string>

namespace Test_H_TP_CO2H2O_ConstPP_SingleFrac {

void RunTestCase();
void ExecutePlanByKey(const std::string& key);
void RunSolveOnly();
void RunPrepareReference();
void RunValidateOnly();
void RunFullWorkflow();

} // namespace Test_H_TP_CO2H2O_ConstPP_SingleFrac

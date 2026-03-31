#pragma once

#include <string>
#include <vector>

namespace Case2DMatlab {

struct ValidationPlotScriptSpec {
    std::string script_path;
    std::vector<std::string> profile_families;
    std::string final_tag = "t100pct";
};

void WriteNoFracPTValidationPlotScript(const ValidationPlotScriptSpec& spec);
void WriteSingleFracPTValidationPlotScript(const ValidationPlotScriptSpec& spec);

} // namespace Case2DMatlab

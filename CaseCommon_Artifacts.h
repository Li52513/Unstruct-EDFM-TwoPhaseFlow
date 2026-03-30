#pragma once

#include <string>

namespace CaseCommon {

enum class CaseStage {
    SolveOnly,
    PrepareReference,
    ValidateOnly,
    FullWorkflow
};

struct CaseArtifactPaths {
    std::string root_dir;
    std::string case_dir;
    std::string studies_dir;
    std::string figures_dir;
    std::string engineering_dir;
    std::string reference_dir;
    std::string report_dir;
    std::string report_scripts_dir;
};

void EnsureDirRecursive(const std::string& rawPath);
CaseArtifactPaths BuildArtifactPaths(const std::string& rootDir,
                                     const std::string& caseCode,
                                     const std::string& caseSlug);
const char* ToString(CaseStage stage);
bool ParseCaseStage(const std::string& text, CaseStage* stage);

} // namespace CaseCommon

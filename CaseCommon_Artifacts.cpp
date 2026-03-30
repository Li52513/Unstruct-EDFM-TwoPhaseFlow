#include "CaseCommon_Artifacts.h"

#include <algorithm>
#include <sstream>

#ifdef _WIN32
#include <direct.h>
#define CASECOMMON_MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define CASECOMMON_MKDIR(path) mkdir(path, 0777)
#endif

namespace CaseCommon {

void EnsureDirRecursive(const std::string& rawPath) {
    if (rawPath.empty()) return;
    std::string path = rawPath;
    std::replace(path.begin(), path.end(), '\\', '/');

    std::stringstream ss(path);
    std::string token;
    std::string current;
    while (std::getline(ss, token, '/')) {
        if (token.empty() || token == ".") continue;
        if (!current.empty()) current += "/";
        current += token;
        CASECOMMON_MKDIR(current.c_str());
    }
}

CaseArtifactPaths BuildArtifactPaths(const std::string& rootDir,
                                     const std::string& caseCode,
                                     const std::string& caseSlug) {
    CaseArtifactPaths out;
    out.root_dir = rootDir;
    out.case_dir = rootDir + "/" + caseCode + "_" + caseSlug;
    out.studies_dir = out.case_dir + "/studies";
    out.figures_dir = out.case_dir + "/figures";
    out.engineering_dir = out.case_dir + "/engineering";
    out.reference_dir = out.case_dir + "/reference";
    out.report_dir = out.case_dir + "/report";
    out.report_scripts_dir = out.report_dir + "/scripts";
    return out;
}

const char* ToString(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: return "solve_only";
    case CaseStage::PrepareReference: return "prepare_reference";
    case CaseStage::ValidateOnly: return "validate_only";
    case CaseStage::FullWorkflow: return "full_workflow";
    default: return "unknown";
    }
}

bool ParseCaseStage(const std::string& text, CaseStage* stage) {
    if (!stage) return false;
    if (text == "solve_only" || text == "solve") {
        *stage = CaseStage::SolveOnly;
        return true;
    }
    if (text == "prepare_reference" || text == "prepare") {
        *stage = CaseStage::PrepareReference;
        return true;
    }
    if (text == "validate_only" || text == "validate") {
        *stage = CaseStage::ValidateOnly;
        return true;
    }
    if (text == "full_workflow" || text == "full") {
        *stage = CaseStage::FullWorkflow;
        return true;
    }
    return false;
}

} // namespace CaseCommon

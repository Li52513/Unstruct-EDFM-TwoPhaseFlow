#pragma once

#include "CaseCommon_Artifacts.h"

#include <string>
#include <vector>

namespace CaseCommon {

struct SkeletonTemplateContext {
    std::string case_code;
    std::string case_slug;
    std::string title;
    std::string output_root = "Test/Transient/A1_F12_TemplateSystem";
    std::vector<std::string> notes;
};

CaseArtifactPaths PrepareSkeletonArtifacts(const SkeletonTemplateContext& ctx, CaseStage stage);
[[noreturn]] void ThrowStageNotImplemented(const SkeletonTemplateContext& ctx,
                                          CaseStage stage,
                                          const std::string& detail);

} // namespace CaseCommon

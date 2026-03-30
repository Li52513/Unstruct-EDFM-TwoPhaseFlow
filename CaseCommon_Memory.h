#pragma once

#include <string>
#include <vector>

namespace CaseCommon {

struct TemplateSystemMemorySnapshot {
    std::string title;
    std::string current_phase;
    std::vector<std::string> frozen_decisions;
    std::vector<std::string> completed_items;
    std::vector<std::string> replaced_legacy_entries;
    std::vector<std::string> blockers;
    std::vector<std::string> next_steps;
    std::vector<std::string> key_files;
};

bool WriteTemplateSystemMemory(const std::string& path,
                               const TemplateSystemMemorySnapshot& snapshot,
                               std::string* errorMessage = nullptr);

} // namespace CaseCommon

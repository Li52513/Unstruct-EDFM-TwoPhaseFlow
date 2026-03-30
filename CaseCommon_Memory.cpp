#include "CaseCommon_Memory.h"

#include <fstream>

namespace CaseCommon {

namespace {

void WriteSection(std::ofstream& out,
                  const char* title,
                  const std::vector<std::string>& items) {
    out << "## " << title << "\n";
    if (items.empty()) {
        out << "- (none)\n\n";
        return;
    }
    for (const std::string& item : items) {
        out << "- " << item << "\n";
    }
    out << "\n";
}

} // namespace

bool WriteTemplateSystemMemory(const std::string& path,
                               const TemplateSystemMemorySnapshot& snapshot,
                               std::string* errorMessage) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.good()) {
        if (errorMessage) *errorMessage = "failed to open memory file: " + path;
        return false;
    }

    out << "# " << snapshot.title << "\n\n";
    out << "## Current Phase\n";
    out << "- " << snapshot.current_phase << "\n\n";
    WriteSection(out, "Frozen Decisions", snapshot.frozen_decisions);
    WriteSection(out, "Completed Items", snapshot.completed_items);
    WriteSection(out, "Replaced Legacy Entries", snapshot.replaced_legacy_entries);
    WriteSection(out, "Blockers", snapshot.blockers);
    WriteSection(out, "Next Steps", snapshot.next_steps);
    WriteSection(out, "Key Files", snapshot.key_files);
    return true;
}

} // namespace CaseCommon

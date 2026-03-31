#pragma once

#include <string>
#include <vector>

namespace Case2DStudies {

struct SweepStudyRow {
    std::string label;
    std::string case_dir;
    int nx = 0;
    int ny = 0;
    int steps = 0;
    double dt_init = 0.0;
    double h_char = 0.0;
    double t_end = 0.0;
    double l1_norm = 0.0;
    double l2_norm = 0.0;
    double linf_norm = 0.0;
};

void WriteStudyCSV(const std::vector<SweepStudyRow>& rows,
                   const std::string& csvPath,
                   const std::string& study);

bool IsMonotonicNonIncreasingL2(const std::vector<SweepStudyRow>& rows);

} // namespace Case2DStudies

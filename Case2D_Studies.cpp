#include "Case2D_Studies.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace Case2DStudies {

void WriteStudyCSV(const std::vector<SweepStudyRow>& rows,
                   const std::string& csvPath,
                   const std::string& study) {
    std::ofstream csv(csvPath.c_str(), std::ios::out | std::ios::trunc);
    if (!csv.good()) {
        throw std::runtime_error("[Case2D_Studies] failed to open study csv: " + csvPath);
    }

    csv << "study,label,case_dir,nx,ny,dt_init,h_char,t_end,steps,l1_norm,l2_norm,linf_norm\n";
    for (const auto& row : rows) {
        csv << study << "," << row.label << "," << row.case_dir << ","
            << row.nx << "," << row.ny << ","
            << std::setprecision(12) << row.dt_init << ","
            << row.h_char << "," << row.t_end << "," << row.steps << ","
            << row.l1_norm << "," << row.l2_norm << "," << row.linf_norm << "\n";
    }
}

bool IsMonotonicNonIncreasingL2(const std::vector<SweepStudyRow>& rows) {
    const double absTol = 1.0e-10;
    const double relTol = 1.0e-6;
    for (std::size_t i = 1; i < rows.size(); ++i) {
        const double prev = rows[i - 1].l2_norm;
        const double allowed = std::max(absTol, relTol * std::max(prev, 1.0));
        if (rows[i].l2_norm > prev + allowed) {
            return false;
        }
    }
    return true;
}

} // namespace Case2DStudies

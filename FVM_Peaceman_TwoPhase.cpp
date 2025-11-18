#include "FVM_Peaceman_TwoPhase.h"
#include "FieldAcessForDiscre.h"   // 假设的场访问工具头文件
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cassert>
#include <iostream> // For std::cout and std::cerr

// An anonymous namespace to keep helper functions local to this file.
namespace {

    constexpr double kPI = 3.14159265358979323846;

    // --- Helper functions for mesh data access ---

    static inline int cell_id_by_index(const Mesh& mesh, int cidx) {
        // Ensures we get the cell's public ID, which might be different from its index in the vector.
        return mesh.getCells()[cidx].id;
    }

    static inline double cell_measure_2D(const Mesh& mesh, int cidx) {
        // In your framework, cell.volume stores the area for 2D cells.
        const int cid = cell_id_by_index(mesh, cidx);
        return mesh.getCellArea(cid);
    }

    static inline Vector cell_center(const Mesh& mesh, int cidx) {
        return mesh.getCells()[cidx].center;
    }

    // --- Helper to determine effective horizontal permeability from fields ---

    static inline double kh_from_fields(const FieldRegistry& reg, const Mesh& mesh,
        int cid, double fallbackKh)
    {
        // Prioritize anisotropic permeabilities if available
        const double kxx = cellScalar(reg, mesh, "kxx", cid, std::numeric_limits<double>::quiet_NaN());
        const double kyy = cellScalar(reg, mesh, "kyy", cid, std::numeric_limits<double>::quiet_NaN());
        if (std::isfinite(kxx) && std::isfinite(kyy) && kxx > 0.0 && kyy > 0.0) {
            return std::sqrt(kxx * kyy); // Geometric mean for effective permeability
        }

        // Fallback to isotropic permeability
        const double k = cellScalar(reg, mesh, "k", cid, std::numeric_limits<double>::quiet_NaN());
        if (std::isfinite(k) && k > 0.0) {
            return k;
        }

        // If no permeability field is found, use the default fallback value
        return fallbackKh;
    }

    // --- Helper to select perforated cells based on the well specification ---

    static std::vector<int> pick_perforated_cells(const Mesh& mesh, const PeacemanTwoPhase::WellSpec& spec)
    {
        const int Nc = (int)mesh.getCells().size();
        std::vector<int> hit;

        // Strategy 1: Select all cells within the perforation radius
        if (spec.perfRadius > 0.0) {
            const double R2 = spec.perfRadius * spec.perfRadius;
            for (int i = 0; i < Nc; ++i) {
                const Vector cc = cell_center(mesh, i);
                const double d2 = (cc.m_x - spec.pos.m_x) * (cc.m_x - spec.pos.m_x) + (cc.m_y - spec.pos.m_y) * (cc.m_y - spec.pos.m_y);
                if (d2 <= R2) {
                    hit.push_back(i);
                }
            }
            if (!hit.empty()) {
                return hit;
            }
            // If no cells are found within the radius, fall through to the "nearest N" strategy
            std::cerr << "[Peaceman_TwoPhase] Warning: No cells found within perfRadius for well '" << spec.name << "'. Falling back to nearest cell(s).\n";
        }

        // Strategy 2: Find the N nearest cells to the well position
        struct Item { int idx; double d2; };
        std::vector<Item> all_distances(Nc);
        for (int i = 0; i < Nc; ++i) {
            const Vector cc = cell_center(mesh, i);
            const double d2 = (cc.m_x - spec.pos.m_x) * (cc.m_x - spec.pos.m_x) + (cc.m_y - spec.pos.m_y) * (cc.m_y - spec.pos.m_y);
            all_distances[i] = { i, d2 };
        }

        const int take = std::max(1, std::min(spec.maxHitCells, Nc));
        std::nth_element(all_distances.begin(), all_distances.begin() + (take - 1), all_distances.end(),
            [](const Item& a, const Item& b) { return a.d2 < b.d2; });

        hit.resize(take);
        for (int j = 0; j < take; ++j) {
            hit[j] = all_distances[j].idx;
        }
        return hit;
    }
} // end anonymous namespace


// --- Implementation of the public function in the PeacemanTwoPhase namespace ---

namespace PeacemanTwoPhase {

    void build_well_mask_and_WI(
        MeshManager& mgr, FieldRegistry& reg,
        const WellSpec& spec, // This is PeacemanTwoPhase::WellSpec
        const PeacemanParams_TwoPhase& prm,
        const std::string& mask_name, const std::string& WI_name)
    {
        Mesh& mesh = mgr.mesh();
        const int Nc = (int)mesh.getCells().size();

        // 1) Find indices of perforated cells
        std::vector<int> perf_indices = pick_perforated_cells(mesh, spec);

        // 2) Create or retrieve and then clear the output fields
        auto maskFld = reg.getOrCreate<volScalarField>(mask_name, Nc, 0.0);
        maskFld->data.assign(Nc, 0.0);

        auto WIFld = reg.getOrCreate<volScalarField>(WI_name, Nc, 0.0);
        WIFld->data.assign(Nc, 0.0);

        // 3) Calculate the geometric Well Index (WI) for each perforated cell.
        //    WI_i = 2 * pi * k_h_i * H / ( ln(r_e_i / r_w) + skin )
        double WI_total = 0.0;
        for (int cidx : perf_indices) {
            const int cid = cell_id_by_index(mesh, cidx);
            const double Ai = cell_measure_2D(mesh, cidx);
            const double kh = kh_from_fields(reg, mesh, cid, prm.fallbackKh);

            const double re = prm.reFactor * std::sqrt(std::max(Ai, 1e-30) / kPI);
            const double lnre = std::log(std::max(re / spec.rw, 1.0 + 1e-12));
            const double denom = lnre + spec.skin; // Denominator can be zero or negative with large negative skin

            // Ensure denominator is positive for a physically meaningful WI
            if (denom > 1e-12) {
                const double WI_i = (2.0 * kPI) * kh * spec.H / denom;

                (*maskFld)[cidx] = 1.0;
                (*WIFld)[cidx] = WI_i;
                WI_total += WI_i;
            }
            else {
                // Handle non-physical denominator by setting WI to a very large number or zero
                (*maskFld)[cidx] = 1.0; // Still mark as a perforated cell
                (*WIFld)[cidx] = 0.0;   // Set WI to 0 to avoid division by zero or negative flow
                std::cerr << "[Peaceman_TwoPhase] Warning: Non-positive denominator (ln(re/rw)+skin) for well '" << spec.name
                    << "' in cell " << cid << ". Setting WI to 0 for this cell.\n";
            }
        }

        // 4) Provide a summary of the operation
        std::cout << "[Peaceman_TwoPhase] Well '" << spec.name << "' pre-processed: \n"
            << "  - Mask field '" << mask_name << "' created/updated for " << perf_indices.size() << " cells.\n"
            << "  - WI field '" << WI_name << "' created/updated with total Σ(WI) = " << WI_total << ".\n";
    }

} // namespace PeacemanTwoPhase
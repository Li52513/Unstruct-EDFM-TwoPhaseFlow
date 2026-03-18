/**
 * @file WellDOFManager.h
 * @brief Manages independent well bottom-hole-pressure (BHP) DOFs in the FIM global sparse system.
 *
 * @details
 * For each unique well name in the active schedule, one extra block is appended
 * to FIM_BlockSparseMatrix beyond the reservoir (matrix + fracture) blocks.
 * The well block stores the well BHP as DOF index 0.  DOFs 1..N-1 are
 * identity-constrained (residual=0, J_diag=1) to prevent singular rows.
 *
 * Well equations assembled per Newton iteration:
 *   BHP mode:  R_well[0] = P_wbh - P_bhp_target         (Dirichlet)
 *   Rate mode: R_well[0] = WI_mob*(P_res - P_wbh) - q_spec
 *
 * Reservoir equation augmentation (both modes):
 *   J_off(comp, well_block)[0,0] += -WI_mob  (coupling dR_res/dP_wbh)
 *
 * Usage order:
 *   1. Setup(wells, totalReservoirBlocks)          -- before matrix construction
 *   2. RegisterPatternConnections(global_mat)       -- before FreezePattern()
 *   3. InitWellState(state, ic, active_wells)       -- after IC assignment
 *   4. [Newton iter] AssembleWellEquations(...)     -- after BoundaryAssembler
 *   5. [post-solve]  UpdateWellState(state, dx, a)  -- after ApplyTrialUpdate
 */

#pragma once

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <Eigen/Core>

#include "FIM_BlockSparseMatrix.h"
#include "FIM_StateMap.h"
#include "Well_WellControlTypes.h"

namespace FIM_Engine {

template <int N>
class WellDOFManager {
public:
    struct WellEntry {
        std::string name;
        int block_idx;   ///< global block index (= reservoir_block_count_ + w)
        int comp_cell;   ///< completion reservoir block index (from WellScheduleStep::completion_id)
    };

    WellDOFManager() : reservoir_block_count_(0) {}

    // =========================================================
    // Phase 1: Setup (before FIM_BlockSparseMatrix construction)
    // =========================================================

    /**
     * @brief Register unique well names from active schedule and assign block indices.
     * @param active_wells    Active WellScheduleStep list.
     * @param reservoir_blocks  Number of reservoir (matrix+fracture) blocks.
     */
    void Setup(const std::vector<WellScheduleStep>& active_wells, int reservoir_blocks) {
        reservoir_block_count_ = reservoir_blocks;
        wells_.clear();
        name_to_idx_.clear();
        for (const auto& w : active_wells) {
            if (name_to_idx_.count(w.well_name)) continue; // already registered
            WellEntry e;
            e.name      = w.well_name;
            e.block_idx = reservoir_blocks + static_cast<int>(wells_.size());
            e.comp_cell = w.completion_id;
            name_to_idx_[w.well_name] = static_cast<int>(wells_.size());
            wells_.push_back(e);
        }
    }

    int NumWells()              const { return static_cast<int>(wells_.size()); }
    int WellBlockOffset()       const { return reservoir_block_count_; }
    int TotalBlocksWithWells()  const { return reservoir_block_count_ + NumWells(); }
    const WellEntry& GetEntry(int w_idx) const { return wells_[w_idx]; }
    bool Empty()                const { return wells_.empty(); }

    // =========================================================
    // Phase 2: Pattern registration (before FreezePattern())
    // =========================================================

    /**
     * @brief Register off-diagonal zero-blocks between each completion cell and its
     *        well block. Call this before global_mat.FreezePattern().
     */
    void RegisterPatternConnections(FIM_BlockSparseMatrix<N>& mat) const {
        for (const auto& e : wells_) {
            if (e.comp_cell < 0 || e.comp_cell >= reservoir_block_count_) {
                std::cerr << "[WellDOF] Well '" << e.name
                    << "' comp_cell=" << e.comp_cell << " out of range [0,"
                    << reservoir_block_count_ << "). Skipping pattern.\n";
                continue;
            }
            // reservoir block → well block
            mat.AddOffDiagBlock(e.comp_cell, e.block_idx,
                                Eigen::Matrix<double, N, N>::Zero());
            // well block → reservoir block
            mat.AddOffDiagBlock(e.block_idx, e.comp_cell,
                                Eigen::Matrix<double, N, N>::Zero());
        }
    }

    // =========================================================
    // Phase 3: State initialization
    // =========================================================

    /**
     * @brief Initialize well BHP state.
     * BHP-mode wells start at the target pressure; Rate-mode wells start at p_init.
     * @param state         FIM state (must already be sized to TotalBlocksWithWells()).
     * @param p_init        Initial reservoir pressure [Pa].
     * @param t_init        Initial temperature [K].
     * @param sw_init       Initial water saturation (ignored when N != 3).
     * @param active_wells  Active schedule steps.
     */
    void InitWellState(FIM_StateMap<N>& state,
                       double p_init, double t_init, double sw_init,
                       const std::vector<WellScheduleStep>& active_wells) const {
        for (const auto& w : active_wells) {
            auto it = name_to_idx_.find(w.well_name);
            if (it == name_to_idx_.end()) continue;
            const int wb = wells_[it->second].block_idx;
            if (wb >= static_cast<int>(state.P.size())) continue;
            state.P[wb] = (w.control_mode == WellControlMode::BHP)
                          ? w.target_value : p_init;
            state.T[wb] = t_init;
            if constexpr (N == 3) state.Sw[wb] = sw_init;
        }
    }

    // =========================================================
    // Phase 4: Per-Newton-iteration assembly
    // =========================================================

    /**
     * @brief Assemble well block equations and reservoir/well cross-coupling.
     *
     * Call AFTER BoundaryAssembler::Assemble_Wells_*_FullJac has been applied to
     * the reservoir blocks. This method adds:
     *   1) Off-diagonal coupling dR_res[0]/dP_wbh = -WI_mob  (reservoir → well)
     *   2) Well block equation row (R_well, J_well_diag, J_well_vs_Pres for Rate)
     *   3) Identity constraints for DOFs 1..N-1 of well block (prevent singular rows)
     *
     * WI_mob is extracted from kSign * w_jac3[g_eq_p][0] which the BoundaryAssembler
     * already computed as the effective coupling coefficient dR/dP_res.
     *
     * @tparam MeshMgrType  MeshManager or MeshManager_3D
     * @param mat           FIM_BlockSparseMatrix (pattern includes well blocks)
     * @param state         Current FIM state (includes P_wbh in state.P[well_block])
     * @param active_wells  Active WellScheduleStep (same list passed to Setup())
     * @param w_res         Residual from BoundaryAssembler [dim = totalReservoirEq]
     * @param w_jac3        Jacobian [dim = totalReservoirEq][3] from BoundaryAssembler
     * @param mgr           MeshManager (for getEquationIndex)
     * @param kSign         Well source sign convention (params.well_source_sign)
     */
    template<typename MeshMgrType>
    void AssembleWellEquations(
        FIM_BlockSparseMatrix<N>& mat,
        const FIM_StateMap<N>& state,
        const std::vector<WellScheduleStep>& active_wells,
        const std::vector<double>& w_res,
        const std::vector<std::array<double, 3>>& w_jac3,
        const MeshMgrType& mgr,
        double kSign) const
    {
        for (const auto& w : active_wells) {
            auto it = name_to_idx_.find(w.well_name);
            if (it == name_to_idx_.end()) continue;

            const auto& e   = wells_[it->second];
            const int comp  = e.comp_cell;
            const int wb    = e.block_idx;

            if (comp < 0 || comp >= reservoir_block_count_) continue;

            // Extract WI_mob = kSign * dR_res/dP_res from BoundaryAssembler output.
            // This is the effective pressure-coupling coefficient (WI × mobility).
            const int g_eq_p = mgr.getEquationIndex(comp, 0);
            if (g_eq_p < 0 || g_eq_p >= static_cast<int>(w_jac3.size())) continue;
            const double WI_mob = kSign * w_jac3[g_eq_p][0];

            // Guard against degenerate wells (zero WI or small coupling).
            if (std::abs(WI_mob) < 1.0e-30) {
                // Still add identity diagonal to prevent singular well block row.
                mat.AddDiagJacobian(wb, 0, 0, 1.0);
                for (int dof = 1; dof < N; ++dof)
                    mat.AddDiagJacobian(wb, dof, dof, 1.0);
                continue;
            }

            // ── 1. Off-diagonal coupling: dR_res[0]/dP_wbh = -WI_mob ──────────
            // The existing BoundaryAssembler already added:
            //   R_res  += kSign * WI_mob*(P_res - P_bhp_fixed)  [source term]
            //   dR_res/dP_res += WI_mob                         [diagonal]
            // With P_wbh now a DOF, we additionally need:
            //   dR_res/dP_wbh = -WI_mob                         [off-diagonal]
            mat.AddOffDiagJacobian(comp, wb, 0, 0, -WI_mob);

            // ── 2. Well block equation ───────────────────────────────────────────
            const double P_res = state.P[comp];
            const double P_wbh = state.P[wb];

            if (w.control_mode == WellControlMode::BHP) {
                // Dirichlet constraint: R_well = P_wbh - P_target = 0
                // dR_well/dP_wbh = 1.0
                const double R_well = P_wbh - w.target_value;
                mat.AddResidual(wb, 0, R_well);
                mat.AddDiagJacobian(wb, 0, 0, 1.0);
                // No coupling from reservoir into well row for BHP mode.
            }
            else {
                // Rate mode: Peaceman model   R_well = WI_mob*(P_res - P_wbh) - q_spec = 0
                //   dR_well/dP_wbh = -WI_mob   (diagonal of well block)
                //   dR_well/dP_res = +WI_mob   (off-diagonal: well → reservoir)
                //
                // q_spec = kSign * w.target_value  (outflow-positive convention,
                //          same as kSign * w_res[g_eq_p] when converged for rate-mode wells).
                const double q_spec = kSign * w.target_value;
                const double R_well = WI_mob * (P_res - P_wbh) - q_spec;
                mat.AddResidual(wb, 0, R_well);
                mat.AddDiagJacobian(wb, 0, 0, -WI_mob);                // dR_well/dP_wbh
                mat.AddOffDiagJacobian(wb, comp, 0, 0, WI_mob);        // dR_well/dP_res
            }

            // ── 3. Identity constraints for inactive DOFs 1..N-1 ────────────────
            // Prevents singular rows in the well block for T (and Sw in N=3).
            // Residual for these DOFs remains zero (already SetZero'd each Newton iter).
            for (int dof = 1; dof < N; ++dof) {
                mat.AddDiagJacobian(wb, dof, dof, 1.0);
            }
        }
    }

    // =========================================================
    // Phase 5: Post-solve state update
    // =========================================================

    /**
     * @brief Update well BHP state from the Newton update vector.
     * Uses interleaved storage convention: equation index of DOF 0 of well block wb
     * is (wb * N + 0).  T and Sw of well blocks are not physical; they are not updated.
     * @param state  FIM state to update in-place.
     * @param dx     Newton update vector (size = TotalBlocksWithWells() * N).
     * @param alpha  Accepted line-search step size.
     */
    void UpdateWellState(FIM_StateMap<N>& state,
                         const Eigen::VectorXd& dx,
                         double alpha) const {
        for (const auto& e : wells_) {
            const int wb   = e.block_idx;
            const int eq_p = wb * N;  // DOF 0 index in global solution vector
            if (eq_p >= static_cast<int>(dx.size())) continue;
            state.P[wb] += alpha * dx[eq_p];
            // Clamp to physical range (100 Pa .. 200 MPa) for robustness.
            state.P[wb] = std::max(1.0e2, std::min(2.0e8, state.P[wb]));
        }
    }

    // =========================================================
    // Utility
    // =========================================================

    void PrintSummary() const {
        if (wells_.empty()) return;
        std::cout << "    [WellDOF] " << wells_.size() << " independent well block(s) added:\n";
        for (const auto& e : wells_) {
            std::cout << "      well='" << e.name
                      << "' block_idx=" << e.block_idx
                      << " comp_cell=" << e.comp_cell << "\n";
        }
    }

private:
    int reservoir_block_count_;
    std::vector<WellEntry> wells_;
    std::map<std::string, int> name_to_idx_;
};

} // namespace FIM_Engine

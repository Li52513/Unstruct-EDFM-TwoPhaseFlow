/**
 * @file WellDOFManager.h
 * @brief Manages explicit well BHP DOFs with multi-completion coupling in FIM.
 */

#pragma once

#include <vector>
#include <map>
#include <set>
#include <string>
#include <array>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <Eigen/Core>

#include "FIM_BlockSparseMatrix.h"
#include "FIM_StateMap.h"
#include "Well_WellControlTypes.h"
#include "BoundaryAssembler.h"

namespace FIM_Engine {

template <int N>
class WellDOFManager {
public:
    struct WellEntry {
        std::string name;
        int block_idx = -1;                    ///< global block index (= reservoir_block_count_ + well_id)
        std::vector<int> completion_cells;     ///< normalized completion solver block indices
    };

    WellDOFManager() : reservoir_block_count_(0) {}

    /**
     * @brief Register unique wells from normalized schedule and assign block indices.
     * @param normalized_wells schedule list where completion_solver_index is already normalized.
     * @param reservoir_blocks number of reservoir blocks (matrix + fracture).
     */
    void Setup(const std::vector<WellScheduleStep>& normalized_wells, int reservoir_blocks) {
        reservoir_block_count_ = reservoir_blocks;
        wells_.clear();
        name_to_idx_.clear();

        for (const auto& w : normalized_wells) {
            if (w.well_name.empty()) {
                throw std::runtime_error("[WellDOF] Empty well_name is not allowed.");
            }
            const int comp = (w.completion_solver_index >= 0) ? w.completion_solver_index : w.completion_id;
            if (comp < 0 || comp >= reservoir_block_count_) {
                throw std::runtime_error("[WellDOF] completion solver index out of range for well '" + w.well_name + "'.");
            }

            int idx = -1;
            auto it = name_to_idx_.find(w.well_name);
            if (it == name_to_idx_.end()) {
                idx = static_cast<int>(wells_.size());
                WellEntry e;
                e.name = w.well_name;
                e.block_idx = reservoir_blocks + idx;
                wells_.push_back(e);
                name_to_idx_[w.well_name] = idx;
            }
            else {
                idx = it->second;
            }

            auto& comps = wells_[idx].completion_cells;
            if (std::find(comps.begin(), comps.end(), comp) == comps.end()) {
                comps.push_back(comp);
            }
        }

        for (auto& e : wells_) {
            std::sort(e.completion_cells.begin(), e.completion_cells.end());
            if (e.completion_cells.empty()) {
                throw std::runtime_error("[WellDOF] Well '" + e.name + "' has no valid completion.");
            }
        }
    }

    int NumWells() const { return static_cast<int>(wells_.size()); }
    int WellBlockOffset() const { return reservoir_block_count_; }
    int TotalBlocksWithWells() const { return reservoir_block_count_ + NumWells(); }
    const WellEntry& GetEntry(int w_idx) const { return wells_.at(w_idx); }
    bool Empty() const { return wells_.empty(); }

    int FindWellIndex(const std::string& well_name) const {
        auto it = name_to_idx_.find(well_name);
        if (it == name_to_idx_.end()) return -1;
        return it->second;
    }

    int FindWellBlockIndex(const std::string& well_name) const {
        const int idx = FindWellIndex(well_name);
        return (idx >= 0) ? wells_[idx].block_idx : -1;
    }

    int GetPrimaryCompletionCell(int w_idx) const {
        const auto& e = wells_.at(w_idx);
        return e.completion_cells.empty() ? -1 : e.completion_cells.front();
    }

    std::unordered_map<std::string, double> BuildWellBhpMap(const FIM_StateMap<N>& state) const {
        std::unordered_map<std::string, double> out;
        out.reserve(wells_.size());
        for (const auto& e : wells_) {
            if (e.block_idx >= 0 && e.block_idx < static_cast<int>(state.P.size())) {
                out[e.name] = state.P[e.block_idx];
            }
        }
        return out;
    }

    /**
     * @brief Register all completion<->well pattern blocks before FreezePattern().
     */
    void RegisterPatternConnections(FIM_BlockSparseMatrix<N>& mat) const {
        for (const auto& e : wells_) {
            for (int comp : e.completion_cells) {
                if (comp < 0 || comp >= reservoir_block_count_) {
                    throw std::runtime_error("[WellDOF] Pattern registration got out-of-range completion.");
                }
                mat.AddOffDiagBlock(comp, e.block_idx, Eigen::Matrix<double, N, N>::Zero());
                mat.AddOffDiagBlock(e.block_idx, comp, Eigen::Matrix<double, N, N>::Zero());
            }
        }
    }

    /**
     * @brief Initialize well BHP state.
     *        Active BHP wells start from target; others start from p_init.
     */
    void InitWellState(FIM_StateMap<N>& state,
                       double p_init, double t_init, double sw_init,
                       const std::vector<WellScheduleStep>& active_wells) const {
        for (const auto& e : wells_) {
            if (e.block_idx >= static_cast<int>(state.P.size())) continue;
            state.P[e.block_idx] = p_init;
            state.T[e.block_idx] = t_init;
            if constexpr (N == 3) state.Sw[e.block_idx] = sw_init;
        }

        for (const auto& w : active_wells) {
            auto it = name_to_idx_.find(w.well_name);
            if (it == name_to_idx_.end()) continue;
            const int wb = wells_[it->second].block_idx;
            if (wb < 0 || wb >= static_cast<int>(state.P.size())) continue;
            if (w.control_mode == WellControlMode::BHP) {
                state.P[wb] = w.target_value;
            }
        }
    }

    /**
     * @brief Assemble reservoir<->well couplings and well equations from per-completion linearization.
     *
     * Reservoir rows receive: dR_res/dPbh = kSign * dq/dPbh for all affected equations.
     * Well row:
     *   BHP mode : R = Pbh - Ptarget
     *   Rate mode: R = kSign*(sum(q_i) - Qtarget)
     */
    void AssembleWellEquations(
        FIM_BlockSparseMatrix<N>& mat,
        const FIM_StateMap<N>& state,
        const std::vector<WellScheduleStep>& active_wells,
        const std::vector<WellCompletionLinearization>& completion_lins,
        double kSign,
        bool strict = true) const
    {
        std::unordered_map<std::string, std::vector<size_t>> step_indices_by_well;
        step_indices_by_well.reserve(wells_.size());
        for (size_t i = 0; i < active_wells.size(); ++i) {
            const auto& s = active_wells[i];
            if (name_to_idx_.count(s.well_name)) {
                step_indices_by_well[s.well_name].push_back(i);
            }
        }

        std::vector<const WellCompletionLinearization*> lin_by_step(active_wells.size(), nullptr);
        for (const auto& lin : completion_lins) {
            if (lin.step_index >= active_wells.size()) continue;
            if (lin_by_step[lin.step_index] != nullptr && strict) {
                throw std::runtime_error("[WellDOF] Duplicate completion linearization for one schedule step.");
            }
            lin_by_step[lin.step_index] = &lin;
        }

        for (const auto& e : wells_) {
            const int wb = e.block_idx;
            for (int dof = 1; dof < N; ++dof) {
                mat.AddDiagJacobian(wb, dof, dof, 1.0);
            }

            auto it_steps = step_indices_by_well.find(e.name);
            if (it_steps == step_indices_by_well.end() || it_steps->second.empty()) {
                // Inactive this step: decouple as identity row.
                mat.AddDiagJacobian(wb, 0, 0, 1.0);
                continue;
            }

            const auto& step_ids = it_steps->second;
            const WellScheduleStep* ctrl = &active_wells[step_ids.front()];
            for (size_t sid : step_ids) {
                const auto& ws = active_wells[sid];
                if (ws.control_mode != ctrl->control_mode ||
                    ws.component_mode != ctrl->component_mode ||
                    ws.rate_target_type != ctrl->rate_target_type ||
                    ws.injection_is_co2 != ctrl->injection_is_co2 ||
                    !NearlyEqual(ws.injection_temperature, ctrl->injection_temperature) ||
                    !NearlyEqual(ws.target_value, ctrl->target_value)) {
                    throw std::runtime_error("[WellDOF] Inconsistent control config across completions of well '" + e.name + "'.");
                }
            }

            const double pbh = (wb >= 0 && wb < static_cast<int>(state.P.size())) ? state.P[wb] : 0.0;

            double sum_q = 0.0;
            double sum_dq_dPbh = 0.0;
            std::unordered_map<int, std::array<double, 3>> sum_dq_dcell_by_comp;

            for (size_t sid : step_ids) {
                const auto& ws = active_wells[sid];
                const int comp = ws.completion_solver_index;
                if (comp < 0 || comp >= reservoir_block_count_) {
                    if (strict) throw std::runtime_error("[WellDOF] Active completion index out of range.");
                    continue;
                }
                if (std::find(e.completion_cells.begin(), e.completion_cells.end(), comp) == e.completion_cells.end()) {
                    if (strict) throw std::runtime_error("[WellDOF] Active completion is missing from well pattern superset.");
                    continue;
                }

                const WellCompletionLinearization* lin = lin_by_step[sid];
                if (!lin) {
                    if (strict) throw std::runtime_error("[WellDOF] Missing per-completion linearization.");
                    continue;
                }

                auto add_res_to_well_coupling = [&](const WellEquationLinearization& eq_lin) {
                    if (!eq_lin.valid || eq_lin.eq_dof < 0 || eq_lin.eq_dof >= N) return;
                    mat.AddOffDiagJacobian(comp, wb, eq_lin.eq_dof, 0, kSign * eq_lin.dq_dPbh);
                };
                add_res_to_well_coupling(lin->water);
                add_res_to_well_coupling(lin->gas);
                add_res_to_well_coupling(lin->energy);

                double q_sel = 0.0;
                double dq_sel_dPbh = 0.0;
                std::array<double, 3> dq_sel_dcell{ {0.0, 0.0, 0.0} };
                if (ctrl->component_mode == WellComponentMode::Water) {
                    q_sel = lin->water.q;
                    dq_sel_dPbh = lin->water.dq_dPbh;
                    dq_sel_dcell = lin->water.dq_dcell;
                }
                else if (ctrl->component_mode == WellComponentMode::Gas) {
                    q_sel = lin->gas.q;
                    dq_sel_dPbh = lin->gas.dq_dPbh;
                    dq_sel_dcell = lin->gas.dq_dcell;
                }
                else {
                    q_sel = lin->q_total;
                    dq_sel_dPbh = lin->dqtotal_dPbh;
                    dq_sel_dcell = lin->dqtotal_dcell;
                }

                sum_q += q_sel;
                sum_dq_dPbh += dq_sel_dPbh;
                auto& acc = sum_dq_dcell_by_comp[comp];
                for (int g = 0; g < 3; ++g) acc[g] += dq_sel_dcell[g];
            }

            if (ctrl->control_mode == WellControlMode::BHP) {
                const double R_well = pbh - ctrl->target_value;
                mat.AddResidual(wb, 0, R_well);
                mat.AddDiagJacobian(wb, 0, 0, 1.0);
            }
            else {
                const double q_target = ConvertRateTargetToMass(*ctrl);
                const double R_well = kSign * (sum_q - q_target);
                mat.AddResidual(wb, 0, R_well);

                if (std::abs(sum_dq_dPbh) < 1.0e-30) {
                    throw std::runtime_error("[WellDOF] Degenerate Rate well equation (dR/dPbh ~ 0). Check WI/completion configuration.");
                }
                mat.AddDiagJacobian(wb, 0, 0, kSign * sum_dq_dPbh);

                for (const auto& kv : sum_dq_dcell_by_comp) {
                    const int comp = kv.first;
                    const auto& grad = kv.second;
                    for (int g = 0; g < 3; ++g) {
                        const int col_dof = GradIndexToReservoirDOF(g);
                        if (col_dof < 0) continue;
                        if (std::abs(grad[g]) <= 0.0) continue;
                        mat.AddOffDiagJacobian(wb, comp, 0, col_dof, kSign * grad[g]);
                    }
                }
            }
        }
    }

    /**
     * @brief Update well BHP state from Newton increment.
     */
    void UpdateWellState(FIM_StateMap<N>& state,
                         const Eigen::VectorXd& dx,
                         double alpha) const {
        for (const auto& e : wells_) {
            const int wb = e.block_idx;
            const int eq_p = wb * N;
            if (eq_p >= static_cast<int>(dx.size())) continue;
            state.P[wb] += alpha * dx[eq_p];
            state.P[wb] = std::max(1.0e2, std::min(2.0e8, state.P[wb]));
        }
    }

    void PrintSummary() const {
        if (wells_.empty()) return;
        std::cout << "    [WellDOF] " << wells_.size() << " independent well block(s) added:\n";
        for (const auto& e : wells_) {
            std::cout << "      well='" << e.name
                      << "' block_idx=" << e.block_idx
                      << " completions=" << e.completion_cells.size();
            if (!e.completion_cells.empty()) {
                std::cout << " [";
                for (size_t k = 0; k < e.completion_cells.size(); ++k) {
                    if (k) std::cout << ",";
                    std::cout << e.completion_cells[k];
                }
                std::cout << "]";
            }
            std::cout << "\n";
        }
    }

private:
    static bool NearlyEqual(double a, double b, double rel = 1e-10, double abs = 1e-12) {
        return std::abs(a - b) <= std::max(abs, rel * std::max(std::abs(a), std::abs(b)));
    }

    static int GradIndexToReservoirDOF(int grad_idx) {
        if constexpr (N == 3) {
            if (grad_idx == 0) return 0; // P
            if (grad_idx == 1) return 1; // Sw
            if (grad_idx == 2) return 2; // T
            return -1;
        }
        else {
            if (grad_idx == 0) return 0; // P
            if (grad_idx == 2) return 1; // T
            return -1;
        }
    }

    static double ConvertRateTargetToMass(const WellScheduleStep& step) {
        if (step.rate_target_type == WellRateTargetType::MassRate) {
            return step.target_value;
        }

        // Backward-compatible conversion for StdVolumeRate when explicit stream density is absent.
        constexpr double rho_w_std = 1000.0;
        constexpr double rho_g_std = 700.0;
        if (step.component_mode == WellComponentMode::Water) {
            return step.target_value * rho_w_std;
        }
        if (step.component_mode == WellComponentMode::Gas) {
            return step.target_value * rho_g_std;
        }

        double fw = std::max(0.0, step.frac_w);
        double fg = std::max(0.0, step.frac_g);
        if ((fw + fg) <= 1.0e-12) {
            fw = step.injection_is_co2 ? 0.0 : 1.0;
            fg = step.injection_is_co2 ? 1.0 : 0.0;
        }
        else {
            const double s = fw + fg;
            fw /= s;
            fg /= s;
        }
        return step.target_value * (fw * rho_w_std + fg * rho_g_std);
    }

private:
    int reservoir_block_count_;
    std::vector<WellEntry> wells_;
    std::map<std::string, int> name_to_idx_;
};

} // namespace FIM_Engine

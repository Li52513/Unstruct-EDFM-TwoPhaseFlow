#pragma once

#include "Types.hpp"
#include "../FIM_ConnectionManager.h"

#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace FIM_Engine {

    struct EqContrib {
        double R_acc = 0.0, R_flux = 0.0, R_well = 0.0, R_bc = 0.0;
        double D_acc = 0.0, D_flux = 0.0, D_well = 0.0, D_bc = 0.0;

        void reset() {
            R_acc = R_flux = R_well = R_bc = 0.0;
            D_acc = D_flux = D_well = D_bc = 0.0;
        }
        double R_total() const { return R_acc + R_flux + R_well + R_bc; }
        double D_total() const { return D_acc + D_flux + D_well + D_bc; }
    };

    inline const char* ConnectionTypeLabel(ConnectionType type) {
        switch (type) {
        case ConnectionType::Matrix_Matrix: return "MM";
        case ConnectionType::Matrix_Fracture: return "MF";
        case ConnectionType::Fracture_Fracture: return "FF";
        case ConnectionType::Fracture_Internal: return "FI";
        default: return "Unknown";
        }
    }

    template<typename MeshMgrType>
    inline void Run3DDiagnosticPrecheck(
        MeshMgrType& mgr,
        const std::vector<Connection>& conns,
        const TransientSolverParams& params) {
        if constexpr (std::is_same_v<MeshMgrType, MeshManager_3D>) {
            if (params.diag_level == DiagLevel::Off) return;

            std::cout << "\n=========================================================\n"
                << "[PRE3D-DIAG] V3 Diagnostic Pre-check Starting...\n";

            int invalid_idx = 0, invalid_area = 0, invalid_dist = 0;
            const auto& pairs = mgr.getInteractionPairs();
            for (const auto& p : pairs) {
                if (p.matrixSolverIndex < 0 || p.fracCellSolverIndex < 0) invalid_idx++;
                if (p.intersectionArea <= 1e-12) invalid_area++;
                if (p.distMatrixToFracPlane <= 1e-12) invalid_dist++;
            }

            std::cout << "  [PRE3D-PAIR] Total Pairs: " << pairs.size() << "\n"
                << "               Invalid Index: " << invalid_idx << "\n"
                << "               Area <= eps  : " << invalid_area << "\n"
                << "               Dist <= eps  : " << invalid_dist << "\n";

            int mm = 0, mf = 0, ff = 0, fi = 0;
            int neg_t_flow = 0, zero_t_flow = 0;
            int neg_t_heat = 0, zero_t_heat = 0;

            for (const auto& c : conns) {
                if (c.type == ConnectionType::Matrix_Matrix) mm++;
                else if (c.type == ConnectionType::Matrix_Fracture) mf++;
                else if (c.type == ConnectionType::Fracture_Fracture) ff++;
                else if (c.type == ConnectionType::Fracture_Internal) fi++;

                if (c.T_Flow < 0.0) neg_t_flow++;
                if (std::abs(c.T_Flow) < 1e-16) zero_t_flow++;

                if (c.T_Heat < 0.0) neg_t_heat++;
                if (std::abs(c.T_Heat) < 1e-16) zero_t_heat++;
            }

            std::cout << "  [PRE3D-CONN] Connections Count -> MM: " << mm
                << ", MF(NNC): " << mf
                << ", FF: " << ff
                << ", FI: " << fi << "\n"
                << "               T_Flow Anomaly  -> Negative: " << neg_t_flow
                << ", Zero: " << zero_t_flow << "\n"
                << "               T_Heat Anomaly  -> Negative: " << neg_t_heat
                << ", Zero: " << zero_t_heat << "\n"
                << "=========================================================\n\n";
        }
    }

    struct MatrixAuditSummary {
        int total_checked = 0;
        int nnc_checked = 0;
        int ff_checked = 0;
        int nnc_flow_coupled = 0;
        int ff_flow_coupled = 0;
        int nnc_heat_coupled = 0;
        int ff_heat_coupled = 0;
        int nnc_zero_tflow = 0;
        int ff_zero_tflow = 0;
        int nnc_zero_theat = 0;
        int ff_zero_theat = 0;
        int nnc_missing_flow = 0;
        int ff_missing_flow = 0;
        int nnc_missing_heat = 0;
        int ff_missing_heat = 0;

        bool Passed(bool strict) const {
            if (!strict) return true;
            if (total_checked <= 0) return false;

            const int nnc_active_flow = nnc_checked - nnc_zero_tflow;
            const int ff_active_flow = ff_checked - ff_zero_tflow;
            const int nnc_active_heat = nnc_checked - nnc_zero_theat;
            const int ff_active_heat = ff_checked - ff_zero_theat;

            const bool nnc_flow_ok = (nnc_active_flow <= 0) || (nnc_flow_coupled > 0);
            const bool ff_flow_ok = (ff_active_flow <= 0) || (ff_flow_coupled > 0);
            const bool nnc_heat_ok = (nnc_active_heat <= 0) || (nnc_heat_coupled > 0);
            const bool ff_heat_ok = (ff_active_heat <= 0) || (ff_heat_coupled > 0);
            return nnc_flow_ok && ff_flow_ok && nnc_heat_ok && ff_heat_ok;
        }
    };

    template <int N>
    inline MatrixAuditSummary RunMatrixAssemblyAudit(
        const Eigen::SparseMatrix<double>& A,
        const std::vector<Connection>& conns,
        double tol,
        int max_detail) {

        MatrixAuditSummary s;
        const int flow_eq = 0;
        const int heat_eq = (N == 3) ? 2 : ((N == 2) ? 1 : -1);
        int detail_left = std::max(0, max_detail);

        auto abs_coeff = [&](int row, int col) -> double {
            if (row < 0 || col < 0 || row >= A.rows() || col >= A.cols()) return 0.0;
            return std::abs(A.coeff(row, col));
            };

        auto print_missing = [&](const Connection& c, const char* kind) {
            if (detail_left <= 0) return;
            std::cout << "      [MATRIX-AUDIT-MISS] type=" << ConnectionTypeLabel(c.type)
                << " i=" << c.nodeI << " j=" << c.nodeJ
                << " kind=" << kind << "\n";
            --detail_left;
            };

        const double eps = std::max(0.0, tol);
        const double trans_eps = std::max(1.0e-12, eps);

        for (const auto& c : conns) {
            if (c.type != ConnectionType::Matrix_Fracture && c.type != ConnectionType::Fracture_Fracture) {
                continue;
            }

            const bool is_nnc = (c.type == ConnectionType::Matrix_Fracture);
            ++s.total_checked;
            if (is_nnc) ++s.nnc_checked;
            else ++s.ff_checked;

            const bool active_flow = std::abs(c.T_Flow) > trans_eps;
            if (!active_flow) {
                if (is_nnc) ++s.nnc_zero_tflow;
                else ++s.ff_zero_tflow;
            }

            bool active_heat = false;
            if (heat_eq < 0) {
                // N=1 has no heat equation row; mark as inactive to skip strict heat coupling checks.
                if (is_nnc) ++s.nnc_zero_theat;
                else ++s.ff_zero_theat;
            }
            else {
                active_heat = std::abs(c.T_Heat) > trans_eps;
                if (!active_heat) {
                    if (is_nnc) ++s.nnc_zero_theat;
                    else ++s.ff_zero_theat;
                }
            }

            if (active_flow) {
                const int i_flow = c.nodeI * N + flow_eq;
                const int j_flow = c.nodeJ * N + flow_eq;
                const bool has_flow = (abs_coeff(i_flow, j_flow) > eps) || (abs_coeff(j_flow, i_flow) > eps);
                if (has_flow) {
                    if (is_nnc) ++s.nnc_flow_coupled;
                    else ++s.ff_flow_coupled;
                }
                else {
                    if (is_nnc) ++s.nnc_missing_flow;
                    else ++s.ff_missing_flow;
                    print_missing(c, "flow");
                }
            }

            if (heat_eq >= 0 && active_heat) {
                const int i_heat = c.nodeI * N + heat_eq;
                const int j_heat = c.nodeJ * N + heat_eq;
                const bool has_heat = (abs_coeff(i_heat, j_heat) > eps) || (abs_coeff(j_heat, i_heat) > eps);
                if (has_heat) {
                    if (is_nnc) ++s.nnc_heat_coupled;
                    else ++s.ff_heat_coupled;
                }
                else {
                    if (is_nnc) ++s.nnc_missing_heat;
                    else ++s.ff_missing_heat;
                    print_missing(c, "heat");
                }
            }
        }

        return s;
    }

} // namespace FIM_Engine

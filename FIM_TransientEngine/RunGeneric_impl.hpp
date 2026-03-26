#pragma once

#include "StepKernels.hpp"
#include "../FIM_BlockSparseMatrix.h"
#include "../FIM_GlobalAssembler.h"
#include "../FIM_ConnectionManager.h"
#include "../FVM_Ops_AD.h"
#include "../FIM_TopologyBuilder2D.h"
#include "../FIM_TopologyBuilder3D.h"
#include "../TransmissibilitySolver_2D.h"
#include "../TransmissibilitySolver_3D.h"
#include "../FVM_Grad.h"
#include "../2D_PostProcess.h"
#include "../3D_PostProcess.h"
#include "../WellDOFManager.h"

#include <algorithm>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <nlohmann/json.hpp>


namespace FIM_Engine {

    namespace detail {

        template <int N, typename MeshMgrType, typename FieldMgrType>
        inline void EmitStepAccepted(
            const TransientOptionalModules<MeshMgrType, FieldMgrType>& modules,
            int step,
            double time_s,
            double dt_used_s,
            int newton_iters,
            double residual_inf,
            int total_rollbacks,
            const std::string& converge_mode,
            const FIM_StateMap<N>& state)
        {
            if (!modules.on_step_accepted) return;
            const std::vector<double>* sw_ptr = nullptr;
            if constexpr (N == 3) sw_ptr = &state.Sw;
            modules.on_step_accepted(
                step, time_s, dt_used_s, newton_iters, residual_inf, total_rollbacks, converge_mode, state.P, state.T, sw_ptr);
        }

        template <typename MeshMgrType, typename FieldMgrType>
        inline void EmitSnapshot(
            const TransientOptionalModules<MeshMgrType, FieldMgrType>& modules,
            const std::string& tag,
            int step,
            double time_s,
            const std::string& vtk_path)
        {
            if (!modules.on_snapshot_written) return;
            modules.on_snapshot_written(tag, step, time_s, vtk_path);
        }

        inline bool NearlyEqualRelAbs(double a, double b, double rel = 1e-10, double abs = 1e-12) {
            return std::abs(a - b) <= std::max(abs, rel * std::max(std::abs(a), std::abs(b)));
        }

        inline const char* CompletionIdSpaceName(CompletionIdSpace s) {
            switch (s) {
            case CompletionIdSpace::SolverIndex: return "SolverIndex";
            case CompletionIdSpace::FractureLocalIndex: return "FractureLocalIndex";
            case CompletionIdSpace::AutoLegacy: return "AutoLegacy";
            default: return "Unknown";
            }
        }

        struct FluxCellPropsCache {
            double rho_w = 0.0, mu_w = 1e-3, h_w = 0.0;
            double rho_g = 0.0, mu_g = 1e-5, h_g = 0.0;
            double Pc = 0.0, dPc_dSw = 0.0, krw = 1.0, krg = 0.0;
        };

        template <typename MeshMgrType>
        inline bool ComputeMatrixGreenGaussGradient(
            const MeshMgrType& mgr,
            int nMat,
            const std::vector<double>& values,
            const BoundarySetting::BoundaryConditionManager* bc_mgr,
            const std::string& tmp_name,
            std::vector<Vector>& grad_out)
        {
            try {
                volScalarField tmp(tmp_name, static_cast<size_t>(nMat), 0.0);
                for (int ci = 0; ci < nMat; ++ci) tmp[ci] = values[ci];
                FVM_Grad grad_solver(mgr.mesh(), nullptr, nullptr, bc_mgr);
                auto grad_field = grad_solver.compute(tmp, FVM_Grad::Method::GreenGauss);
                if (!grad_field || grad_field->data.size() < static_cast<size_t>(nMat)) return false;
                for (int ci = 0; ci < nMat; ++ci) grad_out[ci] = (*grad_field)[ci];
                return true;
            }
            catch (...) {
                return false;
            }
        }

        template <typename MeshMgrType>
        inline void AccumulateFallbackMatrixFaceGradients(
            const MeshMgrType& mgr,
            int nMat,
            const std::vector<double>& vols,
            const std::vector<double>& pressure_values,
            const std::vector<double>& temperature_values,
            const std::vector<double>* saturation_values,
            std::vector<Vector>& grad_p,
            std::vector<Vector>& grad_t,
            std::vector<Vector>& grad_sw)
        {
            const auto& mesh_faces = mgr.mesh().getFaces();
            for (const auto& mf : mesh_faces) {
                if (mf.isBoundary()) continue;
                const int nI = mf.ownerCell_index;
                const int nJ = mf.neighborCell_index;
                if (nI < 0 || nJ < 0 || nI >= nMat || nJ >= nMat) continue;

                const double vol_I = std::max(vols[nI], 1.0e-30);
                const double vol_J = std::max(vols[nJ], 1.0e-30);
                const double Pf = 0.5 * (pressure_values[nI] + pressure_values[nJ]);
                const double Tf = 0.5 * (temperature_values[nI] + temperature_values[nJ]);
                const double Ax = mf.normal.m_x * mf.length;
                const double Ay = mf.normal.m_y * mf.length;
                const double Az = mf.normal.m_z * mf.length;

                grad_p[nI].m_x += Pf * Ax / vol_I;
                grad_p[nI].m_y += Pf * Ay / vol_I;
                grad_p[nI].m_z += Pf * Az / vol_I;
                grad_p[nJ].m_x -= Pf * Ax / vol_J;
                grad_p[nJ].m_y -= Pf * Ay / vol_J;
                grad_p[nJ].m_z -= Pf * Az / vol_J;

                grad_t[nI].m_x += Tf * Ax / vol_I;
                grad_t[nI].m_y += Tf * Ay / vol_I;
                grad_t[nI].m_z += Tf * Az / vol_I;
                grad_t[nJ].m_x -= Tf * Ax / vol_J;
                grad_t[nJ].m_y -= Tf * Ay / vol_J;
                grad_t[nJ].m_z -= Tf * Az / vol_J;

                if (saturation_values) {
                    const double Swf = 0.5 * ((*saturation_values)[nI] + (*saturation_values)[nJ]);
                    grad_sw[nI].m_x += Swf * Ax / vol_I;
                    grad_sw[nI].m_y += Swf * Ay / vol_I;
                    grad_sw[nI].m_z += Swf * Az / vol_I;
                    grad_sw[nJ].m_x -= Swf * Ax / vol_J;
                    grad_sw[nJ].m_y -= Swf * Ay / vol_J;
                    grad_sw[nJ].m_z -= Swf * Az / vol_J;
                }
            }
        }

        template <int N, typename MeshMgrType>
        inline void PrepareNonOrthogonalGradients(
            const MeshMgrType& mgr,
            const FIM_StateMap<N>& state,
            int nMat,
            const std::vector<double>& vols,
            const BoundarySetting::BoundaryConditionManager* pressure_bc,
            const BoundarySetting::BoundaryConditionManager* temperature_bc,
            const BoundarySetting::BoundaryConditionManager* saturation_bc,
            const std::string& pressure_tmp_name,
            const std::string& temperature_tmp_name,
            const std::string& saturation_tmp_name,
            std::vector<Vector>& grad_p,
            std::vector<Vector>& grad_t,
            std::vector<Vector>& grad_sw)
        {
            std::fill(grad_p.begin(), grad_p.end(), Vector(0.0, 0.0, 0.0));
            std::fill(grad_t.begin(), grad_t.end(), Vector(0.0, 0.0, 0.0));
            std::fill(grad_sw.begin(), grad_sw.end(), Vector(0.0, 0.0, 0.0));

            const bool ok_p = ComputeMatrixGreenGaussGradient(
                mgr, nMat, state.P, pressure_bc, pressure_tmp_name, grad_p);
            const bool ok_t = ComputeMatrixGreenGaussGradient(
                mgr, nMat, state.T, temperature_bc, temperature_tmp_name, grad_t);
            bool ok_sw = true;
            if constexpr (N == 3) {
                ok_sw = ComputeMatrixGreenGaussGradient(
                    mgr, nMat, state.Sw, saturation_bc, saturation_tmp_name, grad_sw);
            }

            if (!(ok_p && ok_t && ok_sw)) {
                const std::vector<double>* sw_values = nullptr;
                if constexpr (N == 3) sw_values = &state.Sw;
                AccumulateFallbackMatrixFaceGradients(
                    mgr, nMat, vols, state.P, state.T, sw_values, grad_p, grad_t, grad_sw);
            }
        }

        inline std::string BuildFailureSnapshotPath(
            const std::string& output_root,
            const std::string& caseName,
            int step,
            int iter_used)
        {
            return BuildCaseOutputFile(
                output_root,
                caseName,
                "fail_snapshot_step" + std::to_string(step) + "_iter" + std::to_string(iter_used) + ".json");
        }

        inline std::string BuildVtkTagPath(
            const std::string& output_root,
            const std::string& caseName,
            const std::string& tag)
        {
            return BuildCaseOutputFile(output_root, caseName, tag + ".vtk");
        }

        template <typename MeshMgrType, typename FieldMgrType, typename SyncFn>
        inline std::string ExportVtkSnapshotFile(
            MeshMgrType& mgr,
            FieldMgrType& fm,
            const VTKBoundaryVisualizationContext* vtk_bc_ctx_ptr,
            const std::string& output_root,
            const std::string& caseName,
            const std::string& tag,
            double time_s,
            bool check_sw,
            SyncFn&& sync_fn)
        {
            sync_fn();
            const std::string fname = BuildVtkTagPath(output_root, caseName, tag);
            if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                PostProcess_2D(mgr, fm, vtk_bc_ctx_ptr).ExportVTK(fname, time_s);
            }
            else {
                PostProcess_3D(mgr, fm, vtk_bc_ctx_ptr).ExportVTK(fname, time_s);
            }
            VerifyVtkExport(fname, check_sw);
            return fname;
        }

        template <int N, typename EvalWPhaseFn, typename EvalGPhaseFn>
        inline std::vector<ADVar<N>> EvaluateConnectionFlux(
            const Connection& conn,
            int i,
            int j,
            int nMat,
            const FIM_StateMap<N>& state,
            const std::vector<Vector>& blockCenters,
            const Vector& gravityVec,
            bool wrt_i,
            const std::vector<FluxCellPropsCache>* passive_cache,
            EvalWPhaseFn&& eval_w_phase,
            EvalGPhaseFn&& eval_g_phase,
            const CapRelPerm::VGParams& vg_cfg,
            const CapRelPerm::RelPermParams& rp_cfg,
            double sw_constitutive_eps)
        {
            std::vector<ADVar<N>> F(N);
            ADVar<N> P_i(state.P[i]), T_i(state.T[i]), P_j(state.P[j]), T_j(state.T[j]);
            const Vector& x_i = blockCenters[i];
            const Vector& x_j = blockCenters[j];
            const bool use_passive_cache =
                (passive_cache != nullptr) &&
                (static_cast<int>(passive_cache->size()) > std::max(i, j));

            if constexpr (N == 2) {
                if (wrt_i) { P_i.grad(0) = 1.0; T_i.grad(1) = 1.0; }
                else { P_j.grad(0) = 1.0; T_j.grad(1) = 1.0; }

                AD_Fluid::ADFluidProperties<N> pW_i, pW_j;
                if (use_passive_cache) {
                    if (wrt_i) {
                        pW_i = eval_w_phase(P_i, T_i);
                        pW_j.rho = ADVar<N>((*passive_cache)[j].rho_w);
                        pW_j.mu = ADVar<N>((*passive_cache)[j].mu_w);
                        pW_j.h = ADVar<N>((*passive_cache)[j].h_w);
                    }
                    else {
                        pW_i.rho = ADVar<N>((*passive_cache)[i].rho_w);
                        pW_i.mu = ADVar<N>((*passive_cache)[i].mu_w);
                        pW_i.h = ADVar<N>((*passive_cache)[i].h_w);
                        pW_j = eval_w_phase(P_j, T_j);
                    }
                }
                else {
                    pW_i = eval_w_phase(P_i, T_i);
                    pW_j = eval_w_phase(P_j, T_j);
                }

                ADVar<N> rho_avg_w = ADVar<N>(0.5) * (pW_i.rho + pW_j.rho);
                ADVar<N> dPhi = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, rho_avg_w, x_i, x_j, gravityVec);
                ADVar<N> mob_i = pW_i.rho / pW_i.mu;
                ADVar<N> mob_j = pW_j.rho / pW_j.mu;
                ADVar<N> up_mob = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi, mob_i, mob_j);
                F[0] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mob, dPhi);
                ADVar<N> up_h = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi, pW_i.h, pW_j.h);
                F[1] = FVM_Ops::Compute_Heat_Flux<N, ADVar<N>>(conn.T_Heat, T_i, T_j, F[0], up_h);
                return F;
            }

            ADVar<N> Sw_i(state.Sw[i]), Sw_j(state.Sw[j]);
            if (wrt_i) { P_i.grad(0) = 1.0; Sw_i.grad(1) = 1.0; T_i.grad(2) = 1.0; }
            else { P_j.grad(0) = 1.0; Sw_j.grad(1) = 1.0; T_j.grad(2) = 1.0; }

            ADVar<N> Sw_i_const, Pc_i, krw_i, krg_i;
            ADVar<N> Sw_j_const, Pc_j, krw_j, krg_j;
            if (use_passive_cache) {
                if (wrt_i) {
                    if (i >= nMat) {
                        Sw_i_const = Sw_i;
                        Pc_i = ADVar<N>(0.0); krw_i = Sw_i_const; krg_i = ADVar<N>(1.0) - Sw_i_const;
                    }
                    else {
                        Sw_i_const = ClampSwForConstitutive<N>(Sw_i, vg_cfg, sw_constitutive_eps);
                        Pc_i = CapRelPerm::pc_vG<N>(Sw_i_const, vg_cfg);
                        CapRelPerm::kr_Mualem_vG<N>(Sw_i_const, vg_cfg, rp_cfg, krw_i, krg_i);
                    }
                    Pc_j = ADVar<N>((*passive_cache)[j].Pc);
                    krw_j = ADVar<N>((*passive_cache)[j].krw);
                    krg_j = ADVar<N>((*passive_cache)[j].krg);
                }
                else {
                    Pc_i = ADVar<N>((*passive_cache)[i].Pc);
                    krw_i = ADVar<N>((*passive_cache)[i].krw);
                    krg_i = ADVar<N>((*passive_cache)[i].krg);
                    if (j >= nMat) {
                        Sw_j_const = Sw_j;
                        Pc_j = ADVar<N>(0.0); krw_j = Sw_j_const; krg_j = ADVar<N>(1.0) - Sw_j_const;
                    }
                    else {
                        Sw_j_const = ClampSwForConstitutive<N>(Sw_j, vg_cfg, sw_constitutive_eps);
                        Pc_j = CapRelPerm::pc_vG<N>(Sw_j_const, vg_cfg);
                        CapRelPerm::kr_Mualem_vG<N>(Sw_j_const, vg_cfg, rp_cfg, krw_j, krg_j);
                    }
                }
            }
            else {
                if (i >= nMat) {
                    Sw_i_const = Sw_i;
                    Pc_i = ADVar<N>(0.0); krw_i = Sw_i_const; krg_i = ADVar<N>(1.0) - Sw_i_const;
                }
                else {
                    Sw_i_const = ClampSwForConstitutive<N>(Sw_i, vg_cfg, sw_constitutive_eps);
                    Pc_i = CapRelPerm::pc_vG<N>(Sw_i_const, vg_cfg);
                    CapRelPerm::kr_Mualem_vG<N>(Sw_i_const, vg_cfg, rp_cfg, krw_i, krg_i);
                }
                if (j >= nMat) {
                    Sw_j_const = Sw_j;
                    Pc_j = ADVar<N>(0.0); krw_j = Sw_j_const; krg_j = ADVar<N>(1.0) - Sw_j_const;
                }
                else {
                    Sw_j_const = ClampSwForConstitutive<N>(Sw_j, vg_cfg, sw_constitutive_eps);
                    Pc_j = CapRelPerm::pc_vG<N>(Sw_j_const, vg_cfg);
                    CapRelPerm::kr_Mualem_vG<N>(Sw_j_const, vg_cfg, rp_cfg, krw_j, krg_j);
                }
            }

            ADVar<N> Pg_i = P_i + Pc_i;
            ADVar<N> Pg_j = P_j + Pc_j;
            AD_Fluid::ADFluidProperties<N> pW_i, pW_j, pG_i, pG_j;
            if (use_passive_cache) {
                if (wrt_i) {
                    pW_i = eval_w_phase(P_i, T_i);
                    pG_i = eval_g_phase(Pg_i, T_i);
                    pW_j.rho = ADVar<N>((*passive_cache)[j].rho_w);
                    pW_j.mu = ADVar<N>((*passive_cache)[j].mu_w);
                    pW_j.h = ADVar<N>((*passive_cache)[j].h_w);
                    pG_j.rho = ADVar<N>((*passive_cache)[j].rho_g);
                    pG_j.mu = ADVar<N>((*passive_cache)[j].mu_g);
                    pG_j.h = ADVar<N>((*passive_cache)[j].h_g);
                }
                else {
                    pW_i.rho = ADVar<N>((*passive_cache)[i].rho_w);
                    pW_i.mu = ADVar<N>((*passive_cache)[i].mu_w);
                    pW_i.h = ADVar<N>((*passive_cache)[i].h_w);
                    pG_i.rho = ADVar<N>((*passive_cache)[i].rho_g);
                    pG_i.mu = ADVar<N>((*passive_cache)[i].mu_g);
                    pG_i.h = ADVar<N>((*passive_cache)[i].h_g);
                    pW_j = eval_w_phase(P_j, T_j);
                    pG_j = eval_g_phase(Pg_j, T_j);
                }
            }
            else {
                pW_i = eval_w_phase(P_i, T_i);
                pW_j = eval_w_phase(P_j, T_j);
                pG_i = eval_g_phase(Pg_i, T_i);
                pG_j = eval_g_phase(Pg_j, T_j);
            }

            ADVar<N> rho_avg_w = ADVar<N>(0.5) * (pW_i.rho + pW_j.rho);
            ADVar<N> rho_avg_g = ADVar<N>(0.5) * (pG_i.rho + pG_j.rho);
            ADVar<N> dPhi_w = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, rho_avg_w, x_i, x_j, gravityVec);
            ADVar<N> dPhi_g = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, Pc_i, Pc_j, rho_avg_g, x_i, x_j, gravityVec);

            ADVar<N> mobW_i = krw_i * pW_i.rho / pW_i.mu;
            ADVar<N> mobW_j = krw_j * pW_j.rho / pW_j.mu;
            ADVar<N> mobG_i = krg_i * pG_i.rho / pG_i.mu;
            ADVar<N> mobG_j = krg_j * pG_j.rho / pG_j.mu;
            ADVar<N> up_mobW = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_w, mobW_i, mobW_j);
            ADVar<N> up_mobG = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_g, mobG_i, mobG_j);

            F[0] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mobW, dPhi_w);
            F[1] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mobG, dPhi_g);
            ADVar<N> up_h_w = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_w, pW_i.h, pW_j.h);
            ADVar<N> up_h_g = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_g, pG_i.h, pG_j.h);
            F[2] = FVM_Ops::Compute_Heat_Flux<N, ADVar<N>>(conn.T_Heat, T_i, T_j, F[0], F[1], up_h_w, up_h_g);
            return F;
        }

        template<typename MeshMgrType>
        inline int NormalizeCompletionSolverIndex(
            const WellScheduleStep& step,
            const WellCompletionSpec& comp,
            const MeshMgrType& mgr,
            bool& used_legacy_mapping)
        {
            const int nMat = MatrixBlockCount(mgr);
            const int totalBlocks = mgr.getTotalDOFCount();
            if (comp.completion_id < 0) {
                throw std::runtime_error("[WellNormalize] completion_id < 0 for well '" + step.well_name + "'.");
            }

            int solverIdx = -1;
            if (step.domain == WellTargetDomain::Matrix) {
                solverIdx = comp.completion_id;
            }
            else {
                if (comp.completion_id_space == CompletionIdSpace::SolverIndex) {
                    solverIdx = comp.completion_id;
                }
                else if (comp.completion_id_space == CompletionIdSpace::FractureLocalIndex) {
                    solverIdx = nMat + comp.completion_id;
                }
                else {
                    used_legacy_mapping = true;
                    solverIdx = (comp.completion_id >= nMat)
                        ? comp.completion_id
                        : (nMat + comp.completion_id);
                }
            }

            if (solverIdx < 0 || solverIdx >= totalBlocks) {
                throw std::runtime_error("[WellNormalize] normalized solver index out of range for well '" + step.well_name + "'.");
            }
            if (step.domain == WellTargetDomain::Matrix && solverIdx >= nMat) {
                throw std::runtime_error("[WellNormalize] Matrix completion mapped into fracture range for well '" + step.well_name + "'.");
            }
            if (step.domain == WellTargetDomain::Fracture && solverIdx < nMat) {
                throw std::runtime_error("[WellNormalize] Fracture completion mapped into matrix range for well '" + step.well_name + "'.");
            }
            return solverIdx;
        }

        template<typename MeshMgrType>
        inline std::vector<WellScheduleStep> SelectActiveAndNormalizeWells(
            const std::vector<WellScheduleStep>& wells,
            const MeshMgrType& mgr,
            double t_now,
            bool filter_by_time = true)
        {
            struct CtrlSig {
                WellControlMode control_mode = WellControlMode::BHP;
                WellComponentMode component_mode = WellComponentMode::Total;
                WellRateTargetType rate_target_type = WellRateTargetType::MassRate;
                double target_value = 0.0;
                double injection_temperature = -1.0;
                bool injection_is_co2 = false;
                bool initialized = false;
            };

            std::unordered_map<std::string, CtrlSig> ctrl_by_well;
            std::unordered_set<std::string> warned_legacy;
            std::unordered_set<std::string> seen_well_completion;
            std::vector<WellScheduleStep> out;
            out.reserve(wells.size());
            const int nMat = MatrixBlockCount(mgr);
            const int totalBlocks = mgr.getTotalDOFCount();

            for (const auto& w : wells) {
                const bool active = (t_now >= w.t_start && t_now < w.t_end);
                if (filter_by_time && !active) continue;

                std::vector<WellCompletionSpec> completion_specs = w.completions;
                if (completion_specs.empty()) {
                    WellCompletionSpec c;
                    c.completion_id = w.completion_id;
                    c.completion_id_space = w.completion_id_space;
                    c.completion_solver_index = w.completion_solver_index;
                    completion_specs.push_back(c);
                }
                if (completion_specs.empty()) {
                    throw std::runtime_error("[WellNormalize] No completion found for well '" + w.well_name + "'.");
                }

                if (filter_by_time) {
                    auto& sig = ctrl_by_well[w.well_name];
                    if (!sig.initialized) {
                        sig.control_mode = w.control_mode;
                        sig.component_mode = w.component_mode;
                        sig.rate_target_type = w.rate_target_type;
                        sig.target_value = w.target_value;
                        sig.injection_temperature = w.injection_temperature;
                        sig.injection_is_co2 = w.injection_is_co2;
                        sig.initialized = true;
                    }
                    else {
                        if (sig.control_mode != w.control_mode ||
                            sig.component_mode != w.component_mode ||
                            sig.rate_target_type != w.rate_target_type ||
                            sig.injection_is_co2 != w.injection_is_co2 ||
                            !NearlyEqualRelAbs(sig.injection_temperature, w.injection_temperature) ||
                            !NearlyEqualRelAbs(sig.target_value, w.target_value)) {
                            throw std::runtime_error("[WellNormalize] Active same-name well has conflicting controls: '" + w.well_name + "'.");
                        }
                    }
                }

                for (const auto& comp : completion_specs) {
                    bool used_legacy_mapping = false;
                    int normalized_solver = (comp.completion_solver_index >= 0)
                        ? comp.completion_solver_index
                        : NormalizeCompletionSolverIndex(w, comp, mgr, used_legacy_mapping);

                    if (normalized_solver < 0 || normalized_solver >= totalBlocks) {
                        throw std::runtime_error("[WellNormalize] completion_solver_index out of range for well '" + w.well_name + "'.");
                    }
                    if (w.domain == WellTargetDomain::Matrix && normalized_solver >= nMat) {
                        throw std::runtime_error("[WellNormalize] Matrix completion_solver_index out of matrix range for well '" + w.well_name + "'.");
                    }
                    if (w.domain == WellTargetDomain::Fracture && normalized_solver < nMat) {
                        throw std::runtime_error("[WellNormalize] Fracture completion_solver_index out of fracture range for well '" + w.well_name + "'.");
                    }

                    if (used_legacy_mapping) {
                        if (warned_legacy.insert(w.well_name).second) {
                            std::cout << "[WellNormalize] well='" << w.well_name
                                << "' completion_id_space=" << CompletionIdSpaceName(comp.completion_id_space)
                                << "' uses AutoLegacy completion mapping; prefer explicit completion_id_space.\n";
                        }
                    }

                    std::ostringstream key;
                    key << w.well_name << "#" << normalized_solver;
                    if (!seen_well_completion.insert(key.str()).second) {
                        continue;
                    }

                    WellScheduleStep n = w;
                    n.completions.clear();
                    n.completion_id = comp.completion_id;
                    n.completion_id_space = comp.completion_id_space;
                    n.completion_solver_index = normalized_solver;
                    WellCompletionSpec normalized_comp;
                    normalized_comp.completion_id = comp.completion_id;
                    normalized_comp.completion_id_space = comp.completion_id_space;
                    normalized_comp.completion_solver_index = normalized_solver;
                    n.completions.push_back(normalized_comp);
                    out.push_back(n);
                }
            }

            return out;
        }

        template <typename MeshMgrType, typename FieldMgrType>
        inline void RunGenericFIMTransientPressureOnlyN1(
            const std::string& caseName,
            MeshMgrType& mgr,
            FieldMgrType& fm,
            const InitialConditions& ic,
            const std::vector<WellScheduleStep>& wells,
            const TransientSolverParams& params,
            const TransientOptionalModules<MeshMgrType, FieldMgrType>& modules)
        {
            if (!wells.empty()) {
                throw std::runtime_error("[N=1] unsupported: wells are not enabled in pressure-only AD route.");
            }

            std::cout << "\n========== Starting Transient Scenario (N=1 AD): " << caseName << " ==========\n";
            MakePath(params.output_root_dir, caseName);
            InjectStaticProperties(fm);
            if (modules.property_initializer) {
                modules.property_initializer(mgr, fm);
                std::cout << "[Init] External property module injected.\n";
            }

            const int totalBlocks = mgr.getTotalDOFCount();
            const int totalEq = totalBlocks;
            if (totalBlocks <= 0 || totalEq <= 0) {
                throw std::runtime_error("[N=1] invalid system size.");
            }

            const int nMat = MatrixBlockCount(mgr);
            FIM_StateMap<1> state;
            state.InitSizes(totalBlocks);
            for (int i = 0; i < totalBlocks; ++i) {
                state.P[i] = ic.P_init;
                state.T[i] = ic.T_init;
            }

            std::vector<double> vols(totalBlocks, 1.0);
            for (size_t i = 0; i < mgr.mesh().getCells().size(); ++i) vols[i] = mgr.mesh().getCells()[i].volume;
            for (size_t i = 0; i < mgr.fracture_network().getOrderedFractureElements().size(); ++i) {
                auto* elem = mgr.fracture_network().getOrderedFractureElements()[i];
                if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                    vols[nMat + static_cast<int>(i)] = elem->length * std::max(elem->aperture, 1e-6);
                }
                else {
                    const auto* frac2d = mgr.findFractureByID(mgr.fracture_network(), elem->parentFractureID);
                    const double ap = frac2d ? frac2d->aperture : 1e-3;
                    vols[nMat + static_cast<int>(i)] = std::max(elem->area, 1e-12) * std::max(ap, 1e-8);
                }
            }

            std::vector<Vector> blockCenters(totalBlocks, Vector(0.0, 0.0, 0.0));
            for (size_t i = 0; i < mgr.mesh().getCells().size(); ++i) blockCenters[i] = mgr.mesh().getCells()[i].center;
            if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                const auto& orderedFrac = mgr.fracture_network().getOrderedFractureElements();
                const auto& macros = mgr.fracture_network().fractures;
                for (size_t i = 0; i < orderedFrac.size(); ++i) {
                    const auto* elem = orderedFrac[i];
                    Vector c(0.0, 0.0, 0.0);
                    if (elem && elem->parentFractureID >= 0 && elem->parentFractureID < static_cast<int>(macros.size())) {
                        const auto& frac = macros[elem->parentFractureID];
                        const double u = std::max(0.0, std::min(1.0, 0.5 * (elem->param0 + elem->param1)));
                        c = frac.start + (frac.end - frac.start) * u;
                    }
                    blockCenters[nMat + static_cast<int>(i)] = c;
                }
            }
            else {
                const auto& orderedFrac = mgr.fracture_network().getOrderedFractureElements();
                for (size_t i = 0; i < orderedFrac.size(); ++i) {
                    blockCenters[nMat + static_cast<int>(i)] = orderedFrac[i] ? orderedFrac[i]->centroid : Vector(0.0, 0.0, 0.0);
                }
            }

            if (modules.state_initializer) {
                modules.state_initializer(mgr, blockCenters, nMat, state.P, state.T, nullptr);
                std::cout << "[Init] External state initializer injected.\n";
            }

            FIM_ConnectionManager connMgr;
            if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                TransmissibilitySolver_2D::Calculate_Transmissibility_Matrix(mgr, fm);
                TransmissibilitySolver_2D::Calculate_Transmissibility_FractureInternal(mgr, fm);
                TransmissibilitySolver_2D::Calculate_Transmissibility_NNC(mgr, fm);
                TransmissibilitySolver_2D::Calculate_Transmissibility_FF(mgr, fm);
                FIM_TopologyBuilder2D::LoadAllConnections(connMgr, mgr, fm);
            }
            else {
                TransmissibilitySolver_3D::Calculate_Transmissibility_Matrix(mgr, fm);
                TransmissibilitySolver_3D::Calculate_Transmissibility_FractureInternal(mgr, fm);
                TransmissibilitySolver_3D::Calculate_Transmissibility_NNC(mgr, fm);
                TransmissibilitySolver_3D::Calculate_Transmissibility_FF(mgr, fm);
                FIM_TopologyBuilder3D::LoadAllConnections(connMgr, mgr, fm);
            }
            connMgr.FinalizeAndAggregate();

            FIM_BlockSparseMatrix<1> global_mat(totalBlocks);
            for (const auto& conn : connMgr.GetConnections()) {
                global_mat.AddOffDiagBlock(conn.nodeI, conn.nodeJ, Eigen::Matrix<double, 1, 1>::Zero());
                global_mat.AddOffDiagBlock(conn.nodeJ, conn.nodeI, Eigen::Matrix<double, 1, 1>::Zero());
            }
            global_mat.FreezePattern();

            const auto str_rock = PhysicalProperties_string_op::Rock();
            const auto str_frac = PhysicalProperties_string_op::Fracture_string();
            auto phi_mat = fm.getMatrixScalar(str_rock.phi_tag);
            auto cr_mat = fm.getMatrixScalar(str_rock.c_r_tag);
            auto kxx_mat = fm.getMatrixScalar(str_rock.k_xx_tag);
            auto kyy_mat = fm.getMatrixScalar(str_rock.k_yy_tag);
            auto kzz_mat = fm.getMatrixScalar(str_rock.k_zz_tag);
            auto phi_frac = fm.getFractureScalar(str_frac.phi_tag);
            auto cr_frac = fm.getFractureScalar(str_frac.c_r_tag);
            auto kt_frac = fm.getFractureScalar(str_frac.k_t_tag);
            if (!cr_mat) std::cout << "    [Warning] Matrix c_r field not found, defaulting to 0.0.\n";
            if (nMat < totalBlocks && !cr_frac) std::cout << "    [Warning] Fracture c_r field not found, defaulting to 0.0.\n";

            const double t_eval = (modules.pressure_only_temperature_k > 0.0)
                ? modules.pressure_only_temperature_k
                : ic.T_init;
            const bool use_eos = (modules.pressure_only_property_mode == PressureOnlyPropertyMode::CO2_EOS);
            double baseline_rho = modules.pressure_only_baseline_rho;
            double baseline_mu = modules.pressure_only_baseline_mu;
            if (!(baseline_rho > 0.0) || !(baseline_mu > 0.0)) {
                ADVar<1> p0(ic.P_init), t0(t_eval);
                AD_Fluid::ADFluidProperties<1> base_props;
                if (modules.single_phase_fluid == SinglePhaseFluidModel::CO2 || use_eos) {
                    base_props = AD_Fluid::Evaluator::evaluateCO2<1>(p0, t0);
                }
                else {
                    base_props = EvalPrimaryFluid<1>(modules.single_phase_fluid, modules.fluid_property_eval, p0, t0);
                }
                if (!(baseline_rho > 0.0)) baseline_rho = std::max(base_props.rho.val, 1.0e-8);
                if (!(baseline_mu > 0.0)) baseline_mu = std::max(base_props.mu.val, 1.0e-12);
            }

            auto eval_pressure_props = [&](const ADVar<1>& P) {
                if (use_eos) {
                    return AD_Fluid::Evaluator::evaluateCO2<1>(P, ADVar<1>(t_eval));
                }
                AD_Fluid::ADFluidProperties<1> props;
                props.rho = ADVar<1>(baseline_rho);
                props.mu = ADVar<1>(baseline_mu);
                props.cp = ADVar<1>(0.0);
                props.cv = ADVar<1>(0.0);
                props.h = ADVar<1>(0.0);
                props.k = ADVar<1>(0.0);
                props.isFallback = false;
                props.near_bound = false;
                return props;
            };

            const auto p_cfg_n1 = PhysicalProperties_string_op::PressureEquation_String::FIM();
            const auto t_cfg_n1 = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
            const PhysicalProperties_string_op::Water w_cfg_n1;
            VTKBoundaryVisualizationContext vtk_bc_ctx;
            vtk_bc_ctx.water_family_policy = VTKBCWaterFamilyDerivePolicy::FollowPrimaryFluid;
            vtk_bc_ctx.primary_fluid_model =
                (use_eos || modules.single_phase_fluid == SinglePhaseFluidModel::CO2)
                ? VTKBCPrimaryFluidModel::CO2
                : VTKBCPrimaryFluidModel::Water;
            if (modules.pressure_bc) {
                vtk_bc_ctx.bindings.push_back(
                    VTKBCVariableBinding{ p_cfg_n1.pressure_field, modules.pressure_bc, VTKBCTransportKind::Pressure });
            }
            if (modules.temperature_bc) {
                vtk_bc_ctx.bindings.push_back(
                    VTKBCVariableBinding{ t_cfg_n1.temperatue_field, modules.temperature_bc, VTKBCTransportKind::Temperature });
            }
            const VTKBoundaryVisualizationContext* vtk_bc_ctx_ptr =
                vtk_bc_ctx.bindings.empty() ? nullptr : &vtk_bc_ctx;
            const PhysicalProperties_string_op::TransientInternalFieldNames internalNamesN1;
            auto sync_pressure_only_state_to_fields = [&]() {
                if (use_eos) {
                    SyncStateToFieldManager(state, fm, mgr, modules.single_phase_fluid, modules.fluid_property_eval, modules.vg_params, modules.rp_params);
                    return;
                }

                auto f_pw = fm.getOrCreateMatrixScalar(p_cfg_n1.pressure_field, 0.0);
                auto f_T = fm.getOrCreateMatrixScalar(t_cfg_n1.temperatue_field, t_eval);
                auto f_rhow = fm.getOrCreateMatrixScalar(w_cfg_n1.rho_tag, baseline_rho);
                auto f_muw = fm.getOrCreateMatrixScalar(w_cfg_n1.mu_tag, baseline_mu);
                auto f_hw = fm.getOrCreateMatrixScalar(w_cfg_n1.h_tag, 0.0);
                auto f_lamw_mob = fm.getOrCreateMatrixScalar(w_cfg_n1.lambda_w_tag, 1.0 / std::max(baseline_mu, 1.0e-20));
                auto f_kw = fm.getOrCreateMatrixScalar(w_cfg_n1.k_tag, 0.6);
                auto f_P_viz = fm.getOrCreateMatrixScalar(internalNamesN1.pressure_viz, 0.0);

                auto frac_pw = fm.getOrCreateFractureScalar(p_cfg_n1.pressure_field, 0.0);
                auto frac_T = fm.getOrCreateFractureScalar(t_cfg_n1.temperatue_field, t_eval);
                auto frac_rhow = fm.getOrCreateFractureScalar(w_cfg_n1.rho_tag, baseline_rho);
                auto frac_muw = fm.getOrCreateFractureScalar(w_cfg_n1.mu_tag, baseline_mu);
                auto frac_hw = fm.getOrCreateFractureScalar(w_cfg_n1.h_tag, 0.0);
                auto frac_lamw_mob = fm.getOrCreateFractureScalar(w_cfg_n1.lambda_w_tag, 1.0 / std::max(baseline_mu, 1.0e-20));
                auto frac_kw = fm.getOrCreateFractureScalar(w_cfg_n1.k_tag, 0.6);
                auto frac_P_viz = fm.getOrCreateFractureScalar(internalNamesN1.pressure_viz, 0.0);

                const double lam_w_mob = 1.0 / std::max(baseline_mu, 1.0e-20);
                for (int i = 0; i < totalBlocks; ++i) {
                    const double p = state.P[i];
                    if (i < nMat) {
                        if (f_pw && i < static_cast<int>(f_pw->data.size())) f_pw->data[i] = p;
                        if (f_T && i < static_cast<int>(f_T->data.size())) f_T->data[i] = t_eval;
                        if (f_rhow && i < static_cast<int>(f_rhow->data.size())) f_rhow->data[i] = baseline_rho;
                        if (f_muw && i < static_cast<int>(f_muw->data.size())) f_muw->data[i] = baseline_mu;
                        if (f_hw && i < static_cast<int>(f_hw->data.size())) f_hw->data[i] = 0.0;
                        if (f_lamw_mob && i < static_cast<int>(f_lamw_mob->data.size())) f_lamw_mob->data[i] = lam_w_mob;
                        if (f_kw && i < static_cast<int>(f_kw->data.size())) f_kw->data[i] = 0.6;
                        if (f_P_viz && i < static_cast<int>(f_P_viz->data.size())) f_P_viz->data[i] = p;
                    }
                    else {
                        const int fi = i - nMat;
                        if (frac_pw && fi < static_cast<int>(frac_pw->data.size())) frac_pw->data[fi] = p;
                        if (frac_T && fi < static_cast<int>(frac_T->data.size())) frac_T->data[fi] = t_eval;
                        if (frac_rhow && fi < static_cast<int>(frac_rhow->data.size())) frac_rhow->data[fi] = baseline_rho;
                        if (frac_muw && fi < static_cast<int>(frac_muw->data.size())) frac_muw->data[fi] = baseline_mu;
                        if (frac_hw && fi < static_cast<int>(frac_hw->data.size())) frac_hw->data[fi] = 0.0;
                        if (frac_lamw_mob && fi < static_cast<int>(frac_lamw_mob->data.size())) frac_lamw_mob->data[fi] = lam_w_mob;
                        if (frac_kw && fi < static_cast<int>(frac_kw->data.size())) frac_kw->data[fi] = 0.6;
                        if (frac_P_viz && fi < static_cast<int>(frac_P_viz->data.size())) frac_P_viz->data[fi] = p;
                    }
                }
            };

            std::vector<double> P_ref(totalBlocks, ic.P_init);
            for (int i = 0; i < totalBlocks; ++i) P_ref[i] = state.P[i];

            double p_floor = 1.0e4;
            double p_ceil = std::numeric_limits<double>::max();
            double t_floor = 273.15;
            double t_ceil = std::numeric_limits<double>::max();
            if (params.clamp_state_to_eos_bounds && use_eos) {
                const auto& table = CO2PropertyTable::instance();
                p_floor = std::max(p_floor, table.minPressure() + 1.0);
                p_ceil = std::min(p_ceil, table.maxPressure() - 1.0);
                t_floor = std::max(t_floor, table.minTemperature() + 1.0e-6);
                t_ceil = std::min(t_ceil, table.maxTemperature() - 1.0e-6);
                if (!(p_floor < p_ceil)) { p_floor = 1.0e4; p_ceil = std::numeric_limits<double>::max(); }
                if (!(t_floor < t_ceil)) { t_floor = 273.15; t_ceil = std::numeric_limits<double>::max(); }
            }

            double t = 0.0;
            double dt = params.dt_init;
            int completed_steps = 0;
            int total_rollbacks = 0;
            int total_limiters = 0;
            int vtk_export_count = 0;
            bool final_vtk_exported = false;
            double final_residual = std::numeric_limits<double>::quiet_NaN();
            const double target_end_time_s = params.target_end_time_s;
            const Vector gravityVec = params.gravity_vector;

            EmitStepAccepted<1>(modules, 0, t, 0.0, 0, 0.0, total_rollbacks, "initial", state);

            auto export_vtk_snapshot = [&](const std::string& tag, int step_idx) {
                if (modules.disable_default_vtk_output) return;
                const std::string fname = ExportVtkSnapshotFile(
                    mgr, fm, vtk_bc_ctx_ptr, params.output_root_dir, caseName, tag, t, false,
                    [&]() { sync_pressure_only_state_to_fields(); });
                if (tag == "final") final_vtk_exported = true;
                ++vtk_export_count;
                EmitSnapshot(modules, tag, step_idx, t, fname);
                std::cout << "    [VTK Export PASS] " << fname << "\n";
            };

            for (int step = 1; step <= params.max_steps && (target_end_time_s <= 0.0 || t < target_end_time_s - 1.0e-12); ++step) {
                if (target_end_time_s > 0.0) {
                    const double t_remaining = target_end_time_s - t;
                    if (t_remaining <= 1.0e-14) break;
                    dt = std::min(dt, t_remaining);
                }
                dt = std::max(params.dt_min, std::min(dt, params.dt_max));

                FIM_StateMap<1> old_state = state;
                LinearSolverCache<1> linear_solver_cache;
                linear_solver_cache.Configure(params);

                bool converged = false;
                std::string converge_mode = "none";
                std::string fail_reason;
                int iter_used = 0;
                double res_iter1 = -1.0;
                double last_rel_update = std::numeric_limits<double>::infinity();
                double step_final_residual = std::numeric_limits<double>::quiet_NaN();

                for (int iter = 1; iter <= std::max(1, params.max_newton_iter); ++iter) {
                    iter_used = iter;
                    global_mat.SetZero();
                    std::vector<EqContrib> eq_contribs(totalEq);

                    struct CellCache { double rho = 0.0; double mu = 0.0; double mob = 0.0; };
                    std::vector<CellCache> cprop(totalBlocks);
                    for (int bi = 0; bi < totalBlocks; ++bi) {
                        const auto props = eval_pressure_props(ADVar<1>(state.P[bi]));
                        const double rho = std::max(props.rho.val, 1.0e-12);
                        const double mu = std::max(props.mu.val, 1.0e-20);
                        cprop[bi].rho = rho;
                        cprop[bi].mu = mu;
                        cprop[bi].mob = rho / mu;
                    }

                    std::vector<Vector> nonorth_grad_P(nMat, Vector(0.0, 0.0, 0.0));
                    if (params.enable_non_orthogonal_correction && nMat > 0) {
                        bool gg_ok = false;
                        try {
                            volScalarField p_tmp("n1_p_tmp", static_cast<size_t>(nMat), 0.0);
                            for (int ci = 0; ci < nMat; ++ci) p_tmp[ci] = state.P[ci];
                            const BoundarySetting::BoundaryConditionManager* bc_grad = modules.pressure_bc;
                            FVM_Grad grad_solver(mgr.mesh(), nullptr, nullptr, bc_grad);
                            auto grad_field = grad_solver.compute(p_tmp, FVM_Grad::Method::GreenGauss);
                            if (grad_field && grad_field->data.size() >= static_cast<size_t>(nMat)) {
                                for (int ci = 0; ci < nMat; ++ci) {
                                    nonorth_grad_P[ci] = (*grad_field)[ci];
                                }
                                gg_ok = true;
                            }
                        }
                        catch (...) {
                            gg_ok = false;
                        }

                        if (!gg_ok) {
                            const auto& mesh_faces = mgr.mesh().getFaces();
                            for (const auto& mf : mesh_faces) {
                                if (mf.isBoundary()) continue;
                                const int nI = mf.ownerCell_index;
                                const int nJ = mf.neighborCell_index;
                                if (nI < 0 || nJ < 0 || nI >= nMat || nJ >= nMat) continue;
                                const double volI = std::max(vols[nI], 1.0e-30);
                                const double volJ = std::max(vols[nJ], 1.0e-30);
                                const double pFace = 0.5 * (state.P[nI] + state.P[nJ]);
                                double area = std::max(mf.vectorE.Mag(), 0.0);
                                if (area <= 1.0e-20) area = std::max(mf.length, 0.0);
                                const double Ax = mf.normal.m_x * area;
                                const double Ay = mf.normal.m_y * area;
                                const double Az = mf.normal.m_z * area;
                                nonorth_grad_P[nI].m_x += pFace * Ax / volI;
                                nonorth_grad_P[nI].m_y += pFace * Ay / volI;
                                nonorth_grad_P[nI].m_z += pFace * Az / volI;
                                nonorth_grad_P[nJ].m_x -= pFace * Ax / volJ;
                                nonorth_grad_P[nJ].m_y -= pFace * Ay / volJ;
                                nonorth_grad_P[nJ].m_z -= pFace * Az / volJ;
                            }
                        }
                    }

                    for (int bi = 0; bi < totalBlocks; ++bi) {
                        double phi_ref = 0.2;
                        double c_r = 0.0;
                        if (bi < nMat) {
                            if (phi_mat) phi_ref = (*phi_mat)[bi];
                            if (cr_mat) c_r = (*cr_mat)[bi];
                        }
                        else {
                            const int fi = bi - nMat;
                            if (phi_frac) phi_ref = (*phi_frac)[fi];
                            if (cr_frac) c_r = (*cr_frac)[fi];
                        }

                        ADVar<1> P(state.P[bi]); P.grad(0) = 1.0;
                        ADVar<1> P_old(old_state.P[bi]);
                        ADVar<1> acc(0.0);

                        if (!use_eos) {
                            // Constant-baseline pressure-only route:
                            // keep linear transient coefficient consistent with legacy Day6 ladder
                            // (phi*ct) and avoid introducing extra phi(P) nonlinearity.
                            const double ct_eff = std::max(c_r, 0.0);
                            const double coeff = baseline_rho * phi_ref * ct_eff
                                * (vols[bi] / std::max(dt, 1.0e-20));
                            acc = ADVar<1>(coeff) * (P - P_old);
                        }
                        else {
                            auto props = eval_pressure_props(P);
                            auto props_old = eval_pressure_props(P_old);
                            ADVar<1> phi = ADVar<1>(phi_ref) * (ADVar<1>(1.0) + ADVar<1>(c_r) * (P - ADVar<1>(P_ref[bi])));
                            ADVar<1> phi_old = ADVar<1>(phi_ref) * (ADVar<1>(1.0) + ADVar<1>(c_r) * (P_old - ADVar<1>(P_ref[bi])));
                            acc = (props.rho * phi - props_old.rho * phi_old) * (vols[bi] / std::max(dt, 1.0e-20));
                        }

                        std::vector<ADVar<1>> acc_eqs(1);
                        acc_eqs[0] = acc;
                        FIM_GlobalAssembler<1, ADVar<1>>::AssembleAccumulation(bi, acc_eqs, global_mat);

                        const int eq = mgr.getEquationIndex(bi, 0);
                        if (eq >= 0 && eq < totalEq) {
                            eq_contribs[eq].R_acc += acc.val;
                            eq_contribs[eq].D_acc += acc.grad(0);
                        }
                    }

                    for (const auto& conn : connMgr.GetConnections()) {
                        const int i = conn.nodeI;
                        const int j = conn.nodeJ;
                        const Vector& x_i = blockCenters[i];
                        const Vector& x_j = blockCenters[j];

                        auto eval_flux = [&](bool wrt_i) {
                            std::vector<ADVar<1>> F(1);
                            ADVar<1> P_i(state.P[i]);
                            ADVar<1> P_j(state.P[j]);
                            if (wrt_i) P_i.grad(0) = 1.0;
                            else P_j.grad(0) = 1.0;

                            const auto p_i = eval_pressure_props(P_i);
                            const auto p_j = eval_pressure_props(P_j);
                            ADVar<1> rho_avg = ADVar<1>(0.5) * (p_i.rho + p_j.rho);
                            ADVar<1> dPhi = FVM_Ops::Compute_Potential_Diff<1, ADVar<1>, Vector>(P_i, P_j, rho_avg, x_i, x_j, gravityVec);
                            ADVar<1> mob_i = p_i.rho / p_i.mu;
                            ADVar<1> mob_j = p_j.rho / p_j.mu;
                            ADVar<1> up_mob = FVM_Ops::Op_Upwind_AD<1, ADVar<1>>(dPhi, mob_i, mob_j);
                            F[0] = FVM_Ops::Compute_Mass_Flux<1, ADVar<1>>(conn.T_Flow, up_mob, dPhi);
                            return F;
                        };

                        auto f_i = eval_flux(true);
                        auto f_j = eval_flux(false);
                        FIM_GlobalAssembler<1, ADVar<1>>::AssembleFlux(i, j, f_i, f_j, global_mat);

                        const int eq_i = mgr.getEquationIndex(i, 0);
                        const int eq_j = mgr.getEquationIndex(j, 0);
                        if (eq_i >= 0 && eq_i < totalEq) {
                            eq_contribs[eq_i].R_flux += -f_i[0].val;
                            eq_contribs[eq_i].D_flux += -f_i[0].grad(0);
                        }
                        if (eq_j >= 0 && eq_j < totalEq) {
                            eq_contribs[eq_j].R_flux += f_i[0].val;
                            eq_contribs[eq_j].D_flux += f_j[0].grad(0);
                        }

                        if (params.enable_non_orthogonal_correction &&
                            conn.type == ConnectionType::Matrix_Matrix &&
                            i < nMat && j < nMat) {
                            if constexpr (!std::is_same_v<MeshMgrType, MeshManager>) {
                                const double vT_mag = conn.vectorT.Mag();
                                if (vT_mag > 1.0e-15) {
                                    const double gPx = 0.5 * (nonorth_grad_P[i].m_x + nonorth_grad_P[j].m_x);
                                    const double gPy = 0.5 * (nonorth_grad_P[i].m_y + nonorth_grad_P[j].m_y);
                                    const double gPz = 0.5 * (nonorth_grad_P[i].m_z + nonorth_grad_P[j].m_z);
                                    const double corr_P = gPx * conn.vectorT.m_x + gPy * conn.vectorT.m_y + gPz * conn.vectorT.m_z;
                                    const double K_eff = (conn.aux_area > 1.0e-12)
                                        ? conn.T_Flow * conn.aux_dist / conn.aux_area
                                        : 0.0;
                                    const double gdx = gravityVec * (blockCenters[j] - blockCenters[i]);
                                    const double dPhi_sign = (state.P[j] - state.P[i]) - 0.5 * (cprop[i].rho + cprop[j].rho) * gdx;
                                    const double mob = (dPhi_sign < 0.0) ? cprop[i].mob : cprop[j].mob;
                                    const double mass_corr = K_eff * mob * corr_P;
                                    global_mat.AddResidual(i, 0, -mass_corr);
                                    global_mat.AddResidual(j, 0, mass_corr);
                                    if (eq_i >= 0 && eq_i < totalEq) eq_contribs[eq_i].R_flux += -mass_corr;
                                    if (eq_j >= 0 && eq_j < totalEq) eq_contribs[eq_j].R_flux += mass_corr;
                                }
                            }
                        }
                    }

                    if (params.enable_non_orthogonal_correction && nMat > 0) {
                        if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                            const auto& mesh_faces = mgr.mesh().getFaces();
                            for (const auto& face : mesh_faces) {
                                if (face.isBoundary()) continue;
                                const int owner = face.ownerCell_index;
                                const int nei = face.neighborCell_index;
                                if (owner < 0 || nei < 0 || owner >= nMat || nei >= nMat) continue;
                                if (face.vectorT.Mag() <= 1.0e-15) continue;

                                const double w = std::min(1.0, std::max(0.0, face.f_linearInterpolationCoef));
                                const Vector grad_f = nonorth_grad_P[owner] * w + nonorth_grad_P[nei] * (1.0 - w);
                                const double corr_p = grad_f * face.vectorT;

                                double kxx_o = 1.0e-13, kyy_o = 1.0e-13, kzz_o = 1.0e-13;
                                double kxx_n = 1.0e-13, kyy_n = 1.0e-13, kzz_n = 1.0e-13;
                                if (kxx_mat && owner < static_cast<int>(kxx_mat->data.size())) kxx_o = std::max(kxx_mat->data[owner], 0.0);
                                if (kyy_mat && owner < static_cast<int>(kyy_mat->data.size())) kyy_o = std::max(kyy_mat->data[owner], 0.0);
                                if (kzz_mat && owner < static_cast<int>(kzz_mat->data.size())) kzz_o = std::max(kzz_mat->data[owner], 0.0);
                                if (kxx_mat && nei < static_cast<int>(kxx_mat->data.size())) kxx_n = std::max(kxx_mat->data[nei], 0.0);
                                if (kyy_mat && nei < static_cast<int>(kyy_mat->data.size())) kyy_n = std::max(kyy_mat->data[nei], 0.0);
                                if (kzz_mat && nei < static_cast<int>(kzz_mat->data.size())) kzz_n = std::max(kzz_mat->data[nei], 0.0);

                                const double nx = face.normal.m_x;
                                const double ny = face.normal.m_y;
                                const double nz = face.normal.m_z;
                                double kn_o = nx * nx * kxx_o + ny * ny * kyy_o + nz * nz * kzz_o;
                                double kn_n = nx * nx * kxx_n + ny * ny * kyy_n + nz * nz * kzz_n;
                                if (kn_o <= 0.0) kn_o = std::max({ kxx_o, kyy_o, kzz_o, 1.0e-20 });
                                if (kn_n <= 0.0) kn_n = std::max({ kxx_n, kyy_n, kzz_n, 1.0e-20 });

                                const double m_o = kn_o * cprop[owner].rho / std::max(cprop[owner].mu, 1.0e-20);
                                const double m_n = kn_n * cprop[nei].rho / std::max(cprop[nei].mu, 1.0e-20);
                                const double m1 = std::max(m_o, 1.0e-30);
                                const double m2 = std::max(m_n, 1.0e-30);
                                const double mob_face = 2.0 * m1 * m2 / (m1 + m2);

                                const double corr_mass = mob_face * corr_p;
                                // Legacy-equivalent deferred correction in residual form:
                                // rhs_owner -= corr, rhs_nei += corr  <=>  R_owner += corr, R_nei -= corr.
                                global_mat.AddResidual(owner, 0, corr_mass);
                                global_mat.AddResidual(nei, 0, -corr_mass);

                                const int eq_owner = mgr.getEquationIndex(owner, 0);
                                const int eq_nei = mgr.getEquationIndex(nei, 0);
                                if (eq_owner >= 0 && eq_owner < totalEq) eq_contribs[eq_owner].R_flux += corr_mass;
                                if (eq_nei >= 0 && eq_nei < totalEq) eq_contribs[eq_nei].R_flux += -corr_mass;
                            }
                        }
                    }

                    if (modules.pressure_bc) {
                        if (use_eos) {
                            sync_pressure_only_state_to_fields();
                            std::vector<double> bc_res(totalEq, 0.0);
                            std::vector<std::array<double, 3>> bc_jac3(totalEq, std::array<double, 3>{ 0.0, 0.0, 0.0 });
                            FluidPropertyEvalConfig bc_fluid_cfg = modules.fluid_property_eval;
                            bc_fluid_cfg.enable_single_phase_constant = false;
                            bc_fluid_cfg.enable_two_phase_constant = false;
                            bc_fluid_cfg.single_phase_is_co2 = true;
                            if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                                BoundaryAssembler::Assemble_2D_FullJac(
                                    mgr, *modules.pressure_bc, 0, fm,
                                    PhysicalProperties_string_op::PressureEquation_String::FIM().pressure_field,
                                    bc_res, bc_jac3, modules.pressure_bc, nullptr, bc_fluid_cfg, modules.vg_params, modules.rp_params);
                            }
                            else {
                                BoundaryAssembler::Assemble_3D_FullJac(
                                    mgr, *modules.pressure_bc, 0, fm,
                                    PhysicalProperties_string_op::PressureEquation_String::FIM().pressure_field,
                                    bc_res, bc_jac3, modules.pressure_bc, nullptr, bc_fluid_cfg, modules.vg_params, modules.rp_params);
                            }

                            for (int bi = 0; bi < totalBlocks; ++bi) {
                                const int eq = mgr.getEquationIndex(bi, 0);
                                if (eq < 0 || eq >= totalEq) continue;
                                const double r_bc = bc_res[eq];
                                const double dRdP = bc_jac3[eq][0];
                                if (std::abs(r_bc) <= 1.0e-16 && std::abs(dRdP) <= 1.0e-16) continue;
                                global_mat.AddResidual(bi, 0, r_bc);
                                global_mat.AddDiagJacobian(bi, 0, 0, dRdP);
                                eq_contribs[eq].R_bc += r_bc;
                                eq_contribs[eq].D_bc += dRdP;
                            }
                        }
                        else {
                            constexpr double kBcEps = 1.0e-14;
                            const auto& mesh = mgr.mesh();
                            const auto& cells = mesh.getCells();
                            const auto& faces = mesh.getFaces();
                            for (const auto& face : faces) {
                                if (!face.isBoundary()) continue;
                                if (!modules.pressure_bc->HasBC(face.physicalGroupId)) continue;
                                const int owner = face.ownerCell_index;
                                if (owner < 0 || owner >= totalBlocks || owner >= static_cast<int>(cells.size())) continue;

                                const auto bc = modules.pressure_bc->GetBCCoefficients(face.physicalGroupId, face.midpoint);
                                ADVar<1> P_cell(state.P[owner]);
                                P_cell.grad(0) = 1.0;
                                const auto props = eval_pressure_props(P_cell);
                                const ADVar<1> mob = props.rho / props.mu;

                                double k_n = 1.0e-13;
                                if (owner < nMat) {
                                    double kxx = 1.0e-13, kyy = 1.0e-13, kzz = 1.0e-13;
                                    if (kxx_mat && owner < static_cast<int>(kxx_mat->data.size())) kxx = std::max(kxx_mat->data[owner], 0.0);
                                    if (kyy_mat && owner < static_cast<int>(kyy_mat->data.size())) kyy = std::max(kyy_mat->data[owner], 0.0);
                                    if (kzz_mat && owner < static_cast<int>(kzz_mat->data.size())) kzz = std::max(kzz_mat->data[owner], 0.0);
                                    const double nx = face.normal.m_x;
                                    const double ny = face.normal.m_y;
                                    const double nz = face.normal.m_z;
                                    k_n = nx * nx * kxx + ny * ny * kyy + nz * nz * kzz;
                                    if (k_n <= 0.0) k_n = std::max({ kxx, kyy, kzz, 1.0e-20 });
                                }
                                else {
                                    const int fi = owner - nMat;
                                    if (kt_frac && fi >= 0 && fi < static_cast<int>(kt_frac->data.size())) {
                                        k_n = std::max(kt_frac->data[fi], 1.0e-20);
                                    }
                                }

                                const Vector d_owner = face.midpoint - cells[owner].center;
                                const double dist = std::max(std::abs(d_owner * face.normal), 1.0e-12);
                                const double area = std::max(face.vectorE.Mag(), 0.0);
                                const double T_geom = area / dist;
                                const ADVar<1> phys_coeff = ADVar<1>(std::max(k_n, 1.0e-20)) * mob;

                                ADVar<1> q_bc(0.0);
                                bool applied = false;
                                if (bc.type == BoundarySetting::BoundaryType::Dirichlet) {
                                    if (std::abs(bc.a) > kBcEps) {
                                        q_bc = FVM_Ops::Op_Boundary_Dirichlet_AD<1, ADVar<1>>(T_geom, P_cell, bc.c / bc.a);
                                        q_bc = FVM_Ops::Op_Boundary_ScaleFlux_AD<1, ADVar<1>>(q_bc, phys_coeff);
                                        applied = true;
                                    }
                                }
                                else if (bc.type == BoundarySetting::BoundaryType::Neumann) {
                                    q_bc = FVM_Ops::Op_Boundary_Neumann_AD<1, ADVar<1>>(area, bc.c);
                                    applied = true;
                                }
                                else if (bc.type == BoundarySetting::BoundaryType::Robin) {
                                    if (std::abs(bc.a) > kBcEps) {
                                        const double c_l = -bc.a;
                                        const double far_field_value = bc.c / bc.a;
                                        q_bc = ADVar<1>(area) * FVM_Ops::Op_Leakoff_Source_AD<1, ADVar<1>>(true, c_l, P_cell, far_field_value);
                                        q_bc = FVM_Ops::Op_Boundary_ScaleFlux_AD<1, ADVar<1>>(q_bc, phys_coeff);
                                        applied = true;
                                    }
                                }
                                if (!applied) continue;

                                global_mat.AddResidual(owner, 0, q_bc.val);
                                global_mat.AddDiagJacobian(owner, 0, 0, q_bc.grad(0));
                                const int eq = mgr.getEquationIndex(owner, 0);
                                if (eq >= 0 && eq < totalEq) {
                                    eq_contribs[eq].R_bc += q_bc.val;
                                    eq_contribs[eq].D_bc += q_bc.grad(0);
                                }
                            }
                        }
                    }

                    auto A = global_mat.GetFrozenMatrix();
                    auto b = global_mat.ExportEigenResidual();
                    const double conv_res = b.lpNorm<Eigen::Infinity>();
                    if (iter == 1) res_iter1 = conv_res;
                    const double rel_res = conv_res / std::max(res_iter1, 1.0e-30);
                    step_final_residual = conv_res;

                    if (conv_res <= params.abs_res_tol) {
                        converged = true;
                        converge_mode = "abs_res";
                        break;
                    }
                    if (rel_res <= params.rel_res_tol && last_rel_update <= params.rel_update_tol) {
                        converged = true;
                        converge_mode = "rel_res_update";
                        break;
                    }

                    auto linear_out = SolveLinearSystem<1>(A, b, eq_contribs, totalEq, params, linear_solver_cache);
                    if (!linear_out.compute_ok || !linear_out.solve_ok || linear_out.dx.size() != totalEq) {
                        fail_reason = "linear_solve_fail";
                        break;
                    }
                    const Eigen::VectorXd& dx = linear_out.dx;
                    if (!dx.allFinite()) {
                        fail_reason = "dx_nan_inf";
                        break;
                    }

                    // N=1 pressure-only route:
                    // do NOT couple max_dP to dt shrinking, otherwise rollback to small dt
                    // can stall Newton updates and trigger perpetual nonlinear_max_iter.
                    const double max_dP_eff = std::max(1.0e-3, params.max_dP);
                    double alpha = 1.0;
                    for (int bi = 0; bi < totalBlocks; ++bi) {
                        const int eqP = mgr.getEquationIndex(bi, 0);
                        if (eqP < 0 || eqP >= dx.size()) { fail_reason = "invalid_eq_index"; alpha = 0.0; break; }
                        alpha = std::min(alpha, max_dP_eff / (std::abs(dx[eqP]) + 1.0e-14));
                    }
                    alpha = std::min(1.0, alpha);
                    if (!std::isfinite(alpha) || alpha < params.min_alpha) {
                        fail_reason = "alpha_too_small";
                        break;
                    }

                    FIM_StateMap<1> trial_state = state;
                    int trial_limiter = 0;
                    double trial_rel_update = std::numeric_limits<double>::infinity();
                    if (!ApplyTrialUpdate<1>(
                        state, alpha, dx, mgr, totalBlocks, params,
                        p_floor, p_ceil, t_floor, t_ceil,
                        trial_state, trial_limiter, trial_rel_update, 0)) {
                        fail_reason = "state_nan_inf";
                        break;
                    }

                    state = std::move(trial_state);
                    total_limiters += trial_limiter;
                    last_rel_update = trial_rel_update;
                }

                if (!converged) {
                    if (fail_reason.empty()) fail_reason = "nonlinear_max_iter";
                    if (fail_reason == "nonlinear_max_iter") {
                        std::cout << "    [N1-DIAG] step=" << step
                            << " iter_used=" << iter_used
                            << " residual_inf=" << step_final_residual
                            << " last_rel_update=" << last_rel_update
                            << " max_newton_iter=" << params.max_newton_iter
                            << "\n";
                    }
                    ++total_rollbacks;
                    const double rollback_fac = std::min(1.0, std::max(0.05, params.rollback_shrink_factor));
                    dt = std::max(dt * rollback_fac, params.dt_min);
                    state = old_state;
                    --step;
                    std::cout << "    [Rollback] step=" << (step + 1) << " new_dt=" << dt << " reason=" << fail_reason << "\n";
                    if (dt <= params.dt_min && total_rollbacks > 20) {
                        throw std::runtime_error("[FAIL] dt reached lower bound with repeated rollback.");
                    }
                    if (total_rollbacks > 80) {
                        throw std::runtime_error("[FAIL] Max rollbacks exceeded.");
                    }
                    continue;
                }

                t += dt;
                completed_steps = step;
                final_residual = step_final_residual;
                std::cout << "  [Step Success] step=" << std::setw(3) << step
                    << " dt=" << std::scientific << std::setprecision(3) << dt
                    << " iter_used=" << iter_used
                    << " conv_mode=" << converge_mode
                    << " rollback_count=" << total_rollbacks
                    << " limiter_count=" << total_limiters << "\n";

                EmitStepAccepted<1>(modules, step, t, dt, iter_used, step_final_residual, total_rollbacks, converge_mode, state);

                // N=1 route should recover quickly from dt_min after successful steps.
                // Previous rule shrank dt for iter_used >= 5, which can lock dt at dt_min.
                if (iter_used <= 2) dt = std::min(dt * 1.20, params.dt_max);
                else if (iter_used <= 5) dt = std::min(dt * 1.10, params.dt_max);
                else if (iter_used <= 8) dt = std::min(dt * 1.02, params.dt_max);
                else dt = std::max(dt * 0.85, params.dt_min);

                const bool reached_target_end = (target_end_time_s > 0.0) && (t >= target_end_time_s - 1.0e-12);
                if (!modules.disable_default_vtk_output) {
                    if (vtk_export_count == 0) {
                        export_vtk_snapshot("step_" + std::to_string(step), step);
                    }
                    if (step % 10 == 0 || step == params.max_steps || reached_target_end) {
                        export_vtk_snapshot(reached_target_end ? "final" : ("step_" + std::to_string(step)), step);
                    }
                }
            }

            if (!final_vtk_exported && !modules.disable_default_vtk_output) {
                export_vtk_snapshot("final", completed_steps);
            }

            std::cout << "[PASS] Day6 case completed: " << caseName
                << " | steps=" << completed_steps
                << " | rollbacks=" << total_rollbacks
                << " | limiters=" << total_limiters
                << " | final_residual=" << std::scientific << final_residual
                << " | vtk_exports=" << vtk_export_count
                << " | t_end=" << std::scientific << t << " s\n";
        }

    } // namespace detail

    template <int N, typename MeshMgrType, typename FieldMgrType>
    inline void RunGenericFIMTransient(
        const std::string& caseName,
        MeshMgrType& mgr,
        FieldMgrType& fm,
        const InitialConditions& ic,
        const std::vector<WellScheduleStep>& wells,
        const TransientSolverParams& params,
        SolverRoute route,
        const TransientOptionalModules<MeshMgrType, FieldMgrType>& modules)
    {

        if (route == SolverRoute::IMPES) {
            throw std::runtime_error("[TODO] IMPES explicit route is reserved but currently bypassed.");
        }

        if constexpr (N == 1) {
            detail::RunGenericFIMTransientPressureOnlyN1(
                caseName, mgr, fm, ic, wells, params, modules);
            return;
        }

        std::cout << "\n========== Starting Transient Scenario: " << caseName << " ==========\n";
        MakePath(params.output_root_dir, caseName);
        InjectStaticProperties(fm);
        if (modules.property_initializer) {
            modules.property_initializer(mgr, fm);
            std::cout << "[Init] External property module injected.\n";
        }

        const SinglePhaseFluidModel sp_model = modules.single_phase_fluid;
        const FluidPropertyEvalContext fluid_ctx =
            BuildFluidPropertyEvalContext(sp_model, modules.fluid_property_eval);
        const FluidPropertyEvalConfig& fluid_cfg = fluid_ctx.config;
        const PhysicalProperties_string_op::TransientInternalFieldNames internalNames;
        const bool sp_use_co2 = (N == 2) && fluid_cfg.single_phase_is_co2;
        auto eval_w_phase = [&](const ADVar<N>& P, const ADVar<N>& T) {
            if constexpr (N == 3) {
                return EvalTwoPhaseWaterFluid<N>(fluid_ctx, P, T);
            }
            else {
                return EvalSinglePhaseFluid<N>(fluid_ctx, P, T);
            }
        };
        auto eval_g_phase = [&](const ADVar<N>& P, const ADVar<N>& T) {
            return EvalTwoPhaseGasFluid<N>(fluid_ctx, P, T);
        };
        if constexpr (N == 2) {
            for (const auto& w : wells) {
                if (sp_use_co2 && w.component_mode == WellComponentMode::Water) {
                    std::cout << "[WellModeWarn] Single-phase CO2 case uses Water component mode for well '" << w.well_name
                        << "'. Recommend WellComponentMode::Gas for explicit phase semantics.\n";
                }
                if (!sp_use_co2 && w.component_mode == WellComponentMode::Gas) {
                    std::cout << "[WellModeWarn] Single-phase Water case uses Gas component mode for well '" << w.well_name
                        << "'. Recommend WellComponentMode::Water for explicit phase semantics.\n";
                }
            }
        }
        const auto& vg_cfg = modules.vg_params;
        const auto& rp_cfg = modules.rp_params;

        const int totalBlocks = mgr.getTotalDOFCount();
        const int totalEq = mgr.getTotalEquationDOFs();

        const std::vector<WellScheduleStep> normalized_wells_superset =
            detail::SelectActiveAndNormalizeWells(wells, mgr, 0.0, false);

        WellDOFManager<N> well_mgr;
        well_mgr.Setup(normalized_wells_superset, totalBlocks);
        const int totalBlocksWithWells = well_mgr.TotalBlocksWithWells();
        const int totalEqWithWells     = totalBlocksWithWells * N;

        FIM_StateMap<N> state;
        state.InitSizes(totalBlocksWithWells);

        const auto pEqCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const auto tEqCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        const auto sEqCfg = PhysicalProperties_string_op::SaturationEquation_String::FIM();
        const int pressureDof = 0;
        const int saturationDof = (N == 3) ? 1 : -1;
        const int temperatureDof = (N == 3) ? 2 : 1;
        VTKBoundaryVisualizationContext vtk_bc_ctx;
        vtk_bc_ctx.water_family_policy = VTKBCWaterFamilyDerivePolicy::FollowPrimaryFluid;
        if constexpr (N == 3) {
            vtk_bc_ctx.primary_fluid_model = VTKBCPrimaryFluidModel::Water;
        }
        else {
            vtk_bc_ctx.primary_fluid_model = sp_use_co2
                ? VTKBCPrimaryFluidModel::CO2
                : VTKBCPrimaryFluidModel::Water;
        }
        if (modules.pressure_bc) {
            vtk_bc_ctx.bindings.push_back(
                VTKBCVariableBinding{ pEqCfg.pressure_field, modules.pressure_bc, VTKBCTransportKind::Pressure });
        }
        if constexpr (N == 3) {
            if (modules.saturation_bc) {
                vtk_bc_ctx.bindings.push_back(
                    VTKBCVariableBinding{ sEqCfg.saturation, modules.saturation_bc, VTKBCTransportKind::Saturation });
            }
        }
        if (modules.temperature_bc) {
            vtk_bc_ctx.bindings.push_back(
                VTKBCVariableBinding{ tEqCfg.temperatue_field, modules.temperature_bc, VTKBCTransportKind::Temperature });
        }
        const VTKBoundaryVisualizationContext* vtk_bc_ctx_ptr =
            vtk_bc_ctx.bindings.empty() ? nullptr : &vtk_bc_ctx;

        for (int i = 0; i < totalBlocksWithWells; ++i) {
            state.P[i] = ic.P_init;
            state.T[i] = ic.T_init;
            if constexpr (N == 3) state.Sw[i] = ic.Sw_init;
        }

        std::vector<double> vols(totalBlocksWithWells, 1.0);
        const int nMat = MatrixBlockCount(mgr);
        for (size_t i = 0; i < mgr.mesh().getCells().size(); ++i) vols[i] = mgr.mesh().getCells()[i].volume;
        for (size_t i = 0; i < mgr.fracture_network().getOrderedFractureElements().size(); ++i) {
            auto* elem = mgr.fracture_network().getOrderedFractureElements()[i];
            if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                vols[nMat + static_cast<int>(i)] = elem->length * std::max(elem->aperture, 1e-6);
            }
            else {
                const auto* frac2d = mgr.findFractureByID(mgr.fracture_network(), elem->parentFractureID);
                const double ap = frac2d ? frac2d->aperture : 1e-3;
                vols[nMat + static_cast<int>(i)] = std::max(elem->area, 1e-12) * std::max(ap, 1e-8);
            }
        }

        std::vector<Vector> blockCenters(totalBlocksWithWells, Vector(0.0, 0.0, 0.0));
        for (size_t i = 0; i < mgr.mesh().getCells().size(); ++i) blockCenters[i] = mgr.mesh().getCells()[i].center;
        if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
            const auto& orderedFrac = mgr.fracture_network().getOrderedFractureElements();
            const auto& macros = mgr.fracture_network().fractures;
            for (size_t i = 0; i < orderedFrac.size(); ++i) {
                const auto* elem = orderedFrac[i];
                Vector c(0.0, 0.0, 0.0);
                if (elem && elem->parentFractureID >= 0 && elem->parentFractureID < static_cast<int>(macros.size())) {
                    const auto& frac = macros[elem->parentFractureID];
                    double u = std::max(0.0, std::min(1.0, 0.5 * (elem->param0 + elem->param1)));
                    c = frac.start + (frac.end - frac.start) * u;
                }
                blockCenters[nMat + static_cast<int>(i)] = c;
            }
        }
        else {
            const auto& orderedFrac = mgr.fracture_network().getOrderedFractureElements();
            for (size_t i = 0; i < orderedFrac.size(); ++i) {
                blockCenters[nMat + static_cast<int>(i)] = orderedFrac[i] ? orderedFrac[i]->centroid : Vector(0.0, 0.0, 0.0);
            }
        }
        if (modules.state_initializer) {
            modules.state_initializer(mgr, blockCenters, nMat, state.P, state.T, (N == 3) ? &state.Sw : nullptr);
            std::cout << "[Init] External state initializer injected.\n";
        }

        // 閳光偓閳光偓 Step 3: init well BHP state; set well block centers to comp_cell center 閳光偓閳光偓
        const std::vector<WellScheduleStep> initial_active_wells =
            detail::SelectActiveAndNormalizeWells(wells, mgr, 0.0, true);
        well_mgr.InitWellState(state, ic.P_init, ic.T_init, ic.Sw_init, initial_active_wells);
        for (int w = 0; w < well_mgr.NumWells(); ++w) {
            const auto& e = well_mgr.GetEntry(w);
            const int comp_cell = well_mgr.GetPrimaryCompletionCell(w);
            if (comp_cell >= 0 && comp_cell < static_cast<int>(blockCenters.size())) {
                blockCenters[e.block_idx] = blockCenters[comp_cell];
            }
        }
        well_mgr.PrintSummary();
        const Vector gravityVec = params.gravity_vector;

        FIM_ConnectionManager connMgr;
        if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
            TransmissibilitySolver_2D::Calculate_Transmissibility_Matrix(mgr, fm);
            TransmissibilitySolver_2D::Calculate_Transmissibility_FractureInternal(mgr, fm);
            TransmissibilitySolver_2D::Calculate_Transmissibility_NNC(mgr, fm);
            TransmissibilitySolver_2D::Calculate_Transmissibility_FF(mgr, fm);
            FIM_TopologyBuilder2D::LoadAllConnections(connMgr, mgr, fm);
        }
        else {
            TransmissibilitySolver_3D::Calculate_Transmissibility_Matrix(mgr, fm);
            TransmissibilitySolver_3D::Calculate_Transmissibility_FractureInternal(mgr, fm);
            TransmissibilitySolver_3D::Calculate_Transmissibility_NNC(mgr, fm);
            TransmissibilitySolver_3D::Calculate_Transmissibility_FF(mgr, fm);
            FIM_TopologyBuilder3D::LoadAllConnections(connMgr, mgr, fm);
        }
        connMgr.FinalizeAndAggregate();

        FIM_BlockSparseMatrix<N> global_mat(totalBlocksWithWells);
        for (const auto& conn : connMgr.GetConnections()) {
            global_mat.AddOffDiagBlock(conn.nodeI, conn.nodeJ, Eigen::Matrix<double, N, N>::Zero());
            global_mat.AddOffDiagBlock(conn.nodeJ, conn.nodeI, Eigen::Matrix<double, N, N>::Zero());
        }
        well_mgr.RegisterPatternConnections(global_mat);
        
        global_mat.FreezePattern();

        const auto str_rock = PhysicalProperties_string_op::Rock();
        const auto str_frac = PhysicalProperties_string_op::Fracture_string();

        auto phi_mat = fm.getMatrixScalar(str_rock.phi_tag);
        auto rho_mat = fm.getMatrixScalar(str_rock.rho_tag);
        auto cp_mat = fm.getMatrixScalar(str_rock.cp_tag);
        auto cr_mat = fm.getMatrixScalar(str_rock.c_r_tag);
        if (!cr_mat) std::cout << "    [Warning] Matrix c_r field not found, defaulting to 0.0.\n";

        auto phi_frac = fm.getFractureScalar(str_frac.phi_tag);
        auto rho_frac = fm.getFractureScalar(str_frac.rho_tag);
        auto cp_frac = fm.getFractureScalar(str_frac.cp_tag);
        auto cr_frac = fm.getFractureScalar(str_frac.c_r_tag);
        if (nMat < totalBlocks && !cr_frac) std::cout << "    [Warning] Fracture c_r field not found, defaulting to 0.0.\n";

        std::vector<double> P_ref(totalBlocksWithWells, ic.P_init);
        for (int i = 0; i < totalBlocksWithWells; ++i) P_ref[i] = state.P[i];

        double t = 0.0;
        double dt = params.dt_init;
        int total_rollbacks = 0;
        int total_limiters = 0;
        int completed_steps = 0;
        int vtk_export_count = 0;
        bool final_vtk_exported = false;
        double final_residual = std::numeric_limits<double>::quiet_NaN();

        detail::EmitStepAccepted<N>(
            modules, 0, t, 0.0, 0, 0.0, total_rollbacks, "initial", state);

        // State guard rails:
        // 1) always keep conservative physical lower bounds;
        // 2) optionally enforce EOS-table hard clipping (off by default).
        double p_floor = 1.0e4;
        double p_ceil = std::numeric_limits<double>::max();
        double t_floor = 273.15;
        double t_ceil = std::numeric_limits<double>::max();
        if (params.clamp_state_to_eos_bounds) {
            auto merge_table_bounds = [&](const auto& table) {
                p_floor = std::max(p_floor, table.minPressure() + 1.0);
                p_ceil = std::min(p_ceil, table.maxPressure() - 1.0);
                t_floor = std::max(t_floor, table.minTemperature() + 1.0e-6);
                t_ceil = std::min(t_ceil, table.maxTemperature() - 1.0e-6);
                };

            if constexpr (N == 3) {
                // Two-phase: enforce the intersection of water and CO2 EOS domains.
                merge_table_bounds(WaterPropertyTable::instance());
                merge_table_bounds(CO2PropertyTable::instance());
            }
            else if (sp_use_co2) {
                // Single-phase CO2 (N=2): clamp to CO2 EOS bounds.
                merge_table_bounds(CO2PropertyTable::instance());
            }
            else {
                // Single-phase water (N=2): clamp to water EOS bounds.
                merge_table_bounds(WaterPropertyTable::instance());
            }
            if (!(p_floor < p_ceil)) { p_floor = 1.0e4; p_ceil = std::numeric_limits<double>::max(); }
            if (!(t_floor < t_ceil)) { t_floor = 273.15; t_ceil = std::numeric_limits<double>::max(); }
        }
        const double sw_constitutive_eps = std::max(1.0e-12, std::min(1.0e-2, params.sw_safe_eps));

        auto build_profile_from_params = [&]() {
            TransientStageProfile prof;
            prof.dt_max = params.dt_max;
            prof.max_newton_iter = params.max_newton_iter;
            prof.rel_res_tol = params.rel_res_tol;
            prof.rel_update_tol = params.rel_update_tol;
            prof.enable_ptc = params.enable_ptc;
            prof.ptc_lambda_init = params.ptc_lambda_init;
            prof.ptc_lambda_decay = params.ptc_lambda_decay;
            prof.ptc_lambda_min = params.ptc_lambda_min;
            prof.enable_control_ramp = params.enable_control_ramp;
            prof.control_ramp_steps = params.control_ramp_steps;
            prof.control_ramp_min = params.control_ramp_min;
            prof.control_ramp_apply_rate = params.control_ramp_apply_rate;
            prof.control_ramp_apply_bhp = params.control_ramp_apply_bhp;
            prof.dt_relres_iter_grow_hi = params.dt_relres_iter_grow_hi;
            prof.dt_relres_iter_neutral_hi = params.dt_relres_iter_neutral_hi;
            prof.dt_relres_iter_soft_shrink_hi = params.dt_relres_iter_soft_shrink_hi;
            prof.dt_relres_grow_factor = params.dt_relres_grow_factor;
            prof.dt_relres_neutral_factor = params.dt_relres_neutral_factor;
            prof.dt_relres_soft_shrink_factor = params.dt_relres_soft_shrink_factor;
            prof.dt_relres_hard_shrink_factor = params.dt_relres_hard_shrink_factor;
            return prof;
            };

        auto sanitize_profile = [&](TransientStageProfile prof, const TransientStageProfile& fallback) {
            if (!(prof.dt_max > 0.0)) prof.dt_max = fallback.dt_max;
            prof.dt_max = std::max(params.dt_min, prof.dt_max);
            if (prof.max_newton_iter <= 0) prof.max_newton_iter = fallback.max_newton_iter;
            prof.max_newton_iter = std::max(1, prof.max_newton_iter);
            if (!(prof.rel_res_tol > 0.0)) prof.rel_res_tol = fallback.rel_res_tol;
            if (!(prof.rel_update_tol > 0.0)) prof.rel_update_tol = fallback.rel_update_tol;
            if (!(prof.ptc_lambda_init >= 0.0)) prof.ptc_lambda_init = fallback.ptc_lambda_init;
            if (!(prof.ptc_lambda_decay >= 0.0 && prof.ptc_lambda_decay <= 1.0)) prof.ptc_lambda_decay = fallback.ptc_lambda_decay;
            if (!(prof.ptc_lambda_min >= 0.0)) prof.ptc_lambda_min = fallback.ptc_lambda_min;
            if (prof.control_ramp_steps <= 0) prof.control_ramp_steps = fallback.control_ramp_steps;
            prof.control_ramp_steps = std::max(1, prof.control_ramp_steps);
            if (!(prof.control_ramp_min > 0.0 && prof.control_ramp_min <= 1.0)) prof.control_ramp_min = fallback.control_ramp_min;
            prof.dt_relres_iter_grow_hi = std::max(1, prof.dt_relres_iter_grow_hi);
            prof.dt_relres_iter_neutral_hi = std::max(prof.dt_relres_iter_grow_hi, prof.dt_relres_iter_neutral_hi);
            prof.dt_relres_iter_soft_shrink_hi = std::max(prof.dt_relres_iter_neutral_hi, prof.dt_relres_iter_soft_shrink_hi);
            prof.dt_relres_grow_factor = std::max(1.0, prof.dt_relres_grow_factor);
            prof.dt_relres_neutral_factor = std::max(0.1, prof.dt_relres_neutral_factor);
            prof.dt_relres_soft_shrink_factor = std::min(1.0, std::max(0.1, prof.dt_relres_soft_shrink_factor));
            prof.dt_relres_hard_shrink_factor = std::min(1.0, std::max(0.1, prof.dt_relres_hard_shrink_factor));
            return prof;
            };

        const TransientStageProfile base_profile = sanitize_profile(build_profile_from_params(), build_profile_from_params());
        TransientStageProfile startup_profile = base_profile;
        TransientStageProfile long_profile = base_profile;
        if (params.enable_two_stage_profile) {
            startup_profile = sanitize_profile(params.startup_profile, base_profile);
            long_profile = sanitize_profile(params.long_profile, base_profile);
        }

        const double target_end_time_s = params.target_end_time_s;
        bool stage_initialized = false;
        bool last_stage_startup = true;
        double next_vtk_output_time_s = -1.0;

        int cache_step = -1;
        std::vector<double> line_search_hist_cache;
        int ls_fail_iter1_cache = 0;
        int ptc_rescue_used_cache = 0;
        double ptc_rescue_boost_cache = 1.0;
        int controlled_accept_used_cache = 0;

        for (int step = 1; step <= params.max_steps && (target_end_time_s <= 0.0 || t < target_end_time_s - 1.0e-12); ++step) {
            bool converged = false;
            FIM_StateMap<N> old_state = state;
            FIM_StateMap<N> best_state = state;
            std::string fail_reason;
            std::string converge_mode = "none";
            int iter_used = 0;
            double res_iter1 = -1.0;
            double step_final_residual = std::numeric_limits<double>::quiet_NaN();
            double best_res = std::numeric_limits<double>::infinity();
            int best_iter = -1;

            bool incident_dumped_this_step = false;
            int prev_hot_idx = -1;
            double prev_hot_res = -1.0;
            int hot_repeat_count = 0;
            double prev_max_mass_flux = 0.0;
            double prev_max_heat_flux = 0.0;

            if (target_end_time_s > 0.0) {
                const double t_remaining = target_end_time_s - t;
                if (t_remaining <= 1.0e-14) {
                    break;
                }
                dt = std::min(dt, t_remaining);
                if (dt <= 0.0) {
                    break;
                }
            }

            const bool use_startup_stage = params.enable_two_stage_profile && (t < params.startup_end_time_s);
            const TransientStageProfile& active_profile = use_startup_stage ? startup_profile : long_profile;
            const int active_max_newton_iter = std::max(1, active_profile.max_newton_iter);
            const double active_rel_res_tol = std::max(1.0e-16, active_profile.rel_res_tol);
            const double active_rel_update_tol = std::max(1.0e-16, active_profile.rel_update_tol);
            const bool active_enable_ptc = active_profile.enable_ptc;
            const double active_dt_max = std::max(params.dt_min, active_profile.dt_max);
            const bool active_enable_control_ramp = active_profile.enable_control_ramp;
            const int active_control_ramp_steps = std::max(1, active_profile.control_ramp_steps);
            const double active_control_ramp_min = std::max(1.0e-6, std::min(1.0, active_profile.control_ramp_min));
            const bool active_control_ramp_apply_rate = active_profile.control_ramp_apply_rate;
            const bool active_control_ramp_apply_bhp = active_profile.control_ramp_apply_bhp;

            if (params.enable_two_stage_profile && (!stage_initialized || use_startup_stage != last_stage_startup)) {
                std::cout << "    [PROFILE] stage=" << (use_startup_stage ? "startup" : "long")
                    << " t=" << std::scientific << t
                    << " dt=" << dt
                    << " dt_max=" << active_dt_max
                    << " max_iter=" << active_max_newton_iter << "\n";
                next_vtk_output_time_s = -1.0;
            }
            stage_initialized = true;
            last_stage_startup = use_startup_stage;


            if (step == 1) {
                Run3DDiagnosticPrecheck(mgr, connMgr.GetConnections(), params);
            }

            // Dual convergence uses previous accepted relative update.
            double last_rel_update = std::numeric_limits<double>::infinity();

            if (step != cache_step) {
                cache_step = step;
                line_search_hist_cache.clear();
                ls_fail_iter1_cache = 0;
                ptc_rescue_used_cache = 0;
                ptc_rescue_boost_cache = 1.0;
                controlled_accept_used_cache = 0;
            }

            std::vector<double>& line_search_hist = line_search_hist_cache;
            if (line_search_hist.capacity() < static_cast<size_t>(std::max(2, active_max_newton_iter + 1))) {
                line_search_hist.reserve(std::max(2, active_max_newton_iter + 1));
            }

            std::vector<WellScheduleStep> active_wells =
                detail::SelectActiveAndNormalizeWells(wells, mgr, t, true);
            double control_ramp = 1.0;
            double linear_solve_ms_sum = 0.0;
            int linear_solve_calls = 0;
            detail::LinearSolverCache<N> linear_solver_cache;
            linear_solver_cache.Configure(params);
            if (active_enable_control_ramp)
            {
                const int ramp_steps = active_control_ramp_steps;
                const double r0 = active_control_ramp_min;
                if (step <= ramp_steps && ramp_steps > 1) {
                    const double s = static_cast<double>(step - 1) / static_cast<double>(ramp_steps - 1);
                    control_ramp = r0 + (1.0 - r0) * s;
                }
                else if (step <= ramp_steps) {
                    control_ramp = r0;
                }
                control_ramp = std::max(1.0e-6, std::min(1.0, control_ramp));

                if (control_ramp < 1.0 - 1.0e-12) {
                    for (auto& w : active_wells) {
                        const int wbidx = well_mgr.FindWellBlockIndex(w.well_name);
                        const int cidx = w.completion_solver_index;
                        const double p_anchor =
                            (wbidx >= 0 && wbidx < static_cast<int>(state.P.size())) ? state.P[wbidx] :
                            ((cidx >= 0 && cidx < static_cast<int>(state.P.size())) ? state.P[cidx] : ic.P_init);
                        const double t_anchor =
                            (wbidx >= 0 && wbidx < static_cast<int>(state.T.size())) ? state.T[wbidx] :
                            ((cidx >= 0 && cidx < static_cast<int>(state.T.size())) ? state.T[cidx] : ic.T_init);

                        if (w.control_mode == WellControlMode::Rate && active_control_ramp_apply_rate) {
                            w.target_value *= control_ramp;
                        }
                        else if (w.control_mode == WellControlMode::BHP && active_control_ramp_apply_bhp) {
                            w.target_value = p_anchor + control_ramp * (w.target_value - p_anchor);
                        }

                        if (w.injection_temperature > 0.0) {
                            w.injection_temperature = t_anchor + control_ramp * (w.injection_temperature - t_anchor);
                        }
                    }
                    std::cout << "    [CTRL-RAMP] step=" << step << " factor=" << control_ramp << "\n";
                }
            }

            // Residual probe used by Armijo line search.
            auto compute_residual_inf_for_state = [&](const FIM_StateMap<N>& eval_state) -> double {
                FIM_BlockSparseMatrix<N> probe_mat(totalBlocksWithWells);
                probe_mat.SetZero();

                // 1) accumulation
                for (int bi = 0; bi < totalBlocks; ++bi) {
                    double phi_ref = 0.2, c_pr = 1000.0, rho_r = 2600.0, c_r = 0.0;
                    if (bi < nMat) {
                        if (phi_mat) phi_ref = (*phi_mat)[bi];
                        if (cp_mat) c_pr = (*cp_mat)[bi];
                        if (rho_mat) rho_r = (*rho_mat)[bi];
                        if (cr_mat) c_r = (*cr_mat)[bi];
                    }
                    else {
                        int fi = bi - nMat;
                        if (phi_frac) phi_ref = (*phi_frac)[fi];
                        if (cp_frac) c_pr = (*cp_frac)[fi];
                        if (rho_frac) rho_r = (*rho_frac)[fi];
                        if (cr_frac) c_r = (*cr_frac)[fi];
                    }
                    ADVar<N> P(eval_state.P[bi]); P.grad(0) = 1.0;
                    ADVar<N> T(eval_state.T[bi]); T.grad((N == 2) ? 1 : 2) = 1.0;
                    auto pW = eval_w_phase(P, T);
                    ADVar<N> P_old(old_state.P[bi]), T_old(old_state.T[bi]);
                    auto pW_old = eval_w_phase(P_old, T_old);

                    ADVar<N> phi = ADVar<N>(phi_ref) * (ADVar<N>(1.0) + ADVar<N>(c_r) * (P - ADVar<N>(P_ref[bi])));
                    ADVar<N> phi_old = ADVar<N>(phi_ref) * (ADVar<N>(1.0) + ADVar<N>(c_r) * (P_old - ADVar<N>(P_ref[bi])));

                    std::vector<ADVar<N>> acc_eqs(N);
                    if constexpr (N == 2) {
                        ADVar<N> m_w = pW.rho * phi, m_w_old = pW_old.rho * phi_old;
                        acc_eqs[0] = (m_w - m_w_old) * (vols[bi] / dt);
                        ADVar<N> e_w = m_w * (pW.h - P / pW.rho) + ADVar<N>((1.0 - phi_ref) * rho_r * c_pr) * T;
                        ADVar<N> e_w_old = m_w_old * (pW_old.h - P_old / pW_old.rho) + ADVar<N>((1.0 - phi_ref) * rho_r * c_pr) * T_old;
                        acc_eqs[1] = (e_w - e_w_old) * (vols[bi] / dt);
                    }
                    else {
                        ADVar<N> Sw(eval_state.Sw[bi]); Sw.grad(1) = 1.0;
                        ADVar<N> Sw_old(old_state.Sw[bi]);
                        ADVar<N> Sw_const = ClampSwForConstitutive<N>(Sw, vg_cfg, sw_constitutive_eps);
                        ADVar<N> Sw_old_const = ClampSwForConstitutive<N>(Sw_old, vg_cfg, sw_constitutive_eps);
                        ADVar<N> Sg = ADVar<N>(1.0) - Sw;
                        ADVar<N> Sg_old = ADVar<N>(1.0) - Sw_old;
                        ADVar<N> Pc     = (bi >= nMat) ? ADVar<N>(0.0) : CapRelPerm::pc_vG<N>(Sw_const,     vg_cfg);
                        ADVar<N> Pc_old = (bi >= nMat) ? ADVar<N>(0.0) : CapRelPerm::pc_vG<N>(Sw_old_const, vg_cfg);
                        ADVar<N> Pg = P + Pc;
                        ADVar<N> Pg_old = P_old + Pc_old;
                        auto pG = eval_g_phase(Pg, T);
                        auto pG_old = eval_g_phase(Pg_old, T_old);
                        acc_eqs[0] = (pW.rho * phi * Sw - pW_old.rho * phi_old * Sw_old) * (vols[bi] / dt);
                        acc_eqs[1] = (pG.rho * phi * Sg - pG_old.rho * phi_old * Sg_old) * (vols[bi] / dt);
                        ADVar<N> e_fluid = pW.rho * Sw * (pW.h - P / pW.rho) + pG.rho * Sg * (pG.h - Pg / pG.rho);
                        ADVar<N> e_fluid_old = pW_old.rho * Sw_old * (pW_old.h - P_old / pW_old.rho) + pG_old.rho * Sg_old * (pG_old.h - Pg_old / pG_old.rho);
                        ADVar<N> e_rock = ADVar<N>((1.0 - phi_ref) * rho_r * c_pr) * T;
                        ADVar<N> e_rock_old = ADVar<N>((1.0 - phi_ref) * rho_r * c_pr) * T_old;
                        acc_eqs[2] = ((e_fluid * phi + e_rock) - (e_fluid_old * phi_old + e_rock_old)) * (vols[bi] / dt);
                    }
                    FIM_GlobalAssembler<N, ADVar<N>>::AssembleAccumulation(bi, acc_eqs, probe_mat);
                }

                std::vector<detail::FluxCellPropsCache> cprop(totalBlocks);
                for (int ci = 0; ci < totalBlocks; ++ci) {
                    const double p_ci = eval_state.P[ci], t_ci = eval_state.T[ci];
                    {
                        ADVar<N> P0(p_ci), T0(t_ci);
                        auto pw = eval_w_phase(P0, T0);
                        cprop[ci].rho_w = pw.rho.val;
                        cprop[ci].mu_w  = pw.mu.val;
                        cprop[ci].h_w   = pw.h.val;
                    }
                    if constexpr (N == 3) {
                        const double sw_ci = eval_state.Sw[ci];
                        double pc_s = 0.0, dpc_dsw_s = 0.0, krw_s = 1.0, krg_s = 0.0;
                        if (ci >= nMat) {
                            const double sw_c = std::max(0.0, std::min(1.0, sw_ci));
                            krw_s = sw_c; krg_s = 1.0 - sw_c;
                        }
                        else {
                            ADVar<N> SwC(sw_ci);
                            SwC.grad(1) = 1.0;
                            ADVar<N> SwCC = ClampSwForConstitutive<N>(SwC, vg_cfg, sw_constitutive_eps);
                            ADVar<N> pc_ad = CapRelPerm::pc_vG<N>(SwCC, vg_cfg);
                            pc_s = pc_ad.val;
                            dpc_dsw_s = pc_ad.grad(1);
                            ADVar<N> krw0, krg0;
                            CapRelPerm::kr_Mualem_vG<N>(SwCC, vg_cfg, rp_cfg, krw0, krg0);
                            krw_s = krw0.val; krg_s = krg0.val;
                        }
                        cprop[ci].Pc = pc_s;
                        cprop[ci].dPc_dSw = dpc_dsw_s;
                        cprop[ci].krw = krw_s;
                        cprop[ci].krg = krg_s;
                        ADVar<N> Pg0(p_ci + pc_s), T0(t_ci);
                        auto pg = eval_g_phase(Pg0, T0);
                        cprop[ci].rho_g = pg.rho.val;
                        cprop[ci].mu_g = pg.mu.val;
                        cprop[ci].h_g = pg.h.val;
                    }
                }

                std::vector<Vector> nonorth_grad_P(nMat, Vector(0.0, 0.0, 0.0));
                std::vector<Vector> nonorth_grad_T_cell(nMat, Vector(0.0, 0.0, 0.0));
                std::vector<Vector> nonorth_grad_Sw(nMat, Vector(0.0, 0.0, 0.0));
                if (params.enable_non_orthogonal_correction && nMat > 0) {
                    detail::PrepareNonOrthogonalGradients(
                        mgr, eval_state, nMat, vols,
                        modules.pressure_bc, modules.temperature_bc, modules.saturation_bc,
                        internalNames.nonorth_probe_p_tmp,
                        internalNames.nonorth_probe_t_tmp,
                        internalNames.nonorth_probe_sw_tmp,
                        nonorth_grad_P, nonorth_grad_T_cell, nonorth_grad_Sw);
                }

                // 2) flux
                for (const auto& conn : connMgr.GetConnections()) {
                    int i = conn.nodeI, j = conn.nodeJ;
                    auto f_wrt_i = detail::EvaluateConnectionFlux(
                        conn, i, j, nMat, eval_state, blockCenters, gravityVec, true, nullptr,
                        eval_w_phase, eval_g_phase, vg_cfg, rp_cfg, sw_constitutive_eps);
                    auto f_wrt_j = detail::EvaluateConnectionFlux(
                        conn, i, j, nMat, eval_state, blockCenters, gravityVec, false, nullptr,
                        eval_w_phase, eval_g_phase, vg_cfg, rp_cfg, sw_constitutive_eps);
                    FIM_GlobalAssembler<N, ADVar<N>>::AssembleFlux(i, j, f_wrt_i, f_wrt_j, probe_mat);

                    if (params.enable_non_orthogonal_correction &&
                        conn.type == ConnectionType::Matrix_Matrix &&
                        i < nMat && j < nMat)
                    {
                        const double vT_mag = conn.vectorT.Mag();
                        if (vT_mag > 1.0e-15) {
                            const double gP_x = 0.5 * (nonorth_grad_P[i].m_x + nonorth_grad_P[j].m_x);
                            const double gP_y = 0.5 * (nonorth_grad_P[i].m_y + nonorth_grad_P[j].m_y);
                            const double gP_z = 0.5 * (nonorth_grad_P[i].m_z + nonorth_grad_P[j].m_z);
                            const double gT_x = 0.5 * (nonorth_grad_T_cell[i].m_x + nonorth_grad_T_cell[j].m_x);
                            const double gT_y = 0.5 * (nonorth_grad_T_cell[i].m_y + nonorth_grad_T_cell[j].m_y);
                            const double gT_z = 0.5 * (nonorth_grad_T_cell[i].m_z + nonorth_grad_T_cell[j].m_z);
                            const double corr_P = gP_x * conn.vectorT.m_x + gP_y * conn.vectorT.m_y + gP_z * conn.vectorT.m_z;
                            const double corr_T = gT_x * conn.vectorT.m_x + gT_y * conn.vectorT.m_y + gT_z * conn.vectorT.m_z;
                            const double K_eff = (conn.aux_area > 1e-12)
                                ? conn.T_Flow * conn.aux_dist / conn.aux_area : 0.0;
                            const double g_dx_w = gravityVec * (blockCenters[j] - blockCenters[i]);
                            const double dPhi_w_sign = (eval_state.P[j] - eval_state.P[i])
                                - 0.5 * (cprop[i].rho_w + cprop[j].rho_w) * g_dx_w;
                            const double mob_w = (dPhi_w_sign < 0.0)
                                ? cprop[i].krw * cprop[i].rho_w / std::max(cprop[i].mu_w, 1e-20)
                                : cprop[j].krw * cprop[j].rho_w / std::max(cprop[j].mu_w, 1e-20);
                            const double mass_w_corr = K_eff * mob_w * corr_P;
                            const double heat_corr = conn.T_Heat * corr_T;
                            const double h_w_upwind = (dPhi_w_sign < 0.0) ? cprop[i].h_w : cprop[j].h_w;
                            double conv_heat_corr = mass_w_corr * h_w_upwind;
                            probe_mat.AddResidual(i, 0, -mass_w_corr);
                            probe_mat.AddResidual(j, 0, mass_w_corr);
                            if constexpr (N == 3) {
                                const double dPhi_g_sign = (eval_state.P[j] - eval_state.P[i])
                                    + (cprop[j].Pc - cprop[i].Pc)
                                    - 0.5 * (cprop[i].rho_g + cprop[j].rho_g) * g_dx_w;
                                const double gSw_x = 0.5 * (nonorth_grad_Sw[i].m_x + nonorth_grad_Sw[j].m_x);
                                const double gSw_y = 0.5 * (nonorth_grad_Sw[i].m_y + nonorth_grad_Sw[j].m_y);
                                const double gSw_z = 0.5 * (nonorth_grad_Sw[i].m_z + nonorth_grad_Sw[j].m_z);
                                const double corr_Sw = gSw_x * conn.vectorT.m_x + gSw_y * conn.vectorT.m_y + gSw_z * conn.vectorT.m_z;
                                const double dPc_dSw = 0.5 * (cprop[i].dPc_dSw + cprop[j].dPc_dSw);
                                const double corr_Pc = corr_Sw * dPc_dSw;
                                const double mob_g = (dPhi_g_sign < 0.0)
                                    ? cprop[i].krg * cprop[i].rho_g / std::max(cprop[i].mu_g, 1e-20)
                                    : cprop[j].krg * cprop[j].rho_g / std::max(cprop[j].mu_g, 1e-20);
                                const double mass_g_corr = K_eff * mob_g * (corr_P + corr_Pc);
                                const double h_g_upwind = (dPhi_g_sign < 0.0) ? cprop[i].h_g : cprop[j].h_g;
                                conv_heat_corr += mass_g_corr * h_g_upwind;
                                probe_mat.AddResidual(i, 1, -mass_g_corr);
                                probe_mat.AddResidual(j, 1, mass_g_corr);
                            }
                            const double total_heat_corr = heat_corr + conv_heat_corr;
                            probe_mat.AddResidual(i, N - 1, -total_heat_corr);
                            probe_mat.AddResidual(j, N - 1, total_heat_corr);
                        }
                    }
                }

                // 3) well source
                SyncStateToFieldManager(eval_state, fm, mgr, sp_model, fluid_cfg, vg_cfg, rp_cfg);
                std::vector<double> w_res(totalEq, 0.0);
                std::vector<std::array<double, 3>> w_jac3(totalEq, std::array<double, 3>{ 0.0, 0.0, 0.0 });
                std::vector<WellCompletionLinearization> w_lin;
                const auto well_bhp_probe = well_mgr.BuildWellBhpMap(eval_state);
                const int well_dof_w = (N == 3) ? 0 : (sp_use_co2 ? -1 : 0);
                const int well_dof_g = (N == 3) ? 1 : (sp_use_co2 ? 0 : -1);
                const int well_dof_e = (N == 3) ? 2 : 1;
                if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                    BoundaryAssembler::Assemble_Wells_2D_FullJac(
                        mgr, fm, active_wells, 0, well_dof_w, well_dof_g, well_dof_e,
                        w_res, w_jac3, fluid_cfg, vg_cfg, rp_cfg, &w_lin, &well_bhp_probe);
                }
                else {
                    BoundaryAssembler::Assemble_Wells_3D_FullJac(
                        mgr, fm, active_wells, 0, well_dof_w, well_dof_g, well_dof_e,
                        w_res, w_jac3, fluid_cfg, vg_cfg, rp_cfg, &w_lin, &well_bhp_probe);
                }
                const double kWellSourceSignProbe = params.well_source_sign;
                for (int bi = 0; bi < totalBlocks; ++bi) {
                    for (int eq = 0; eq < N; ++eq) {
                        int g_eq = mgr.getEquationIndex(bi, eq);
                        if (g_eq < 0 || g_eq >= totalEq) continue;
                        const double rWell = kWellSourceSignProbe * w_res[g_eq];
                        if (std::abs(rWell) <= 1e-16) continue;
                        probe_mat.AddResidual(bi, eq, rWell);
                    }
                }
                if (!well_mgr.Empty()) {
                    well_mgr.AssembleWellEquations(
                        probe_mat, eval_state, active_wells,
                        w_lin, kWellSourceSignProbe);
                }

                // 4) boundary source
                if constexpr (N == 3) {
                    std::vector<double> bc_res(totalEq, 0.0);
                    std::vector<std::array<double, 3>> bc_jac3(totalEq, std::array<double, 3>{ 0.0, 0.0, 0.0 });
                    BoundaryAssemblyStats bc_stats;
                    if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                        bc_stats = BoundaryAssembler::Assemble_2D_CoupledN3_FullJac(
                            mgr, fm, modules.pressure_bc, modules.saturation_bc, modules.temperature_bc,
                            pressureDof, pressureDof, saturationDof, temperatureDof,
                            bc_res, bc_jac3, fluid_cfg, vg_cfg, rp_cfg);
                    }
                    else {
                        bc_stats = BoundaryAssembler::Assemble_3D_CoupledN3_FullJac(
                            mgr, fm, modules.pressure_bc, modules.saturation_bc, modules.temperature_bc,
                            pressureDof, pressureDof, saturationDof, temperatureDof,
                            bc_res, bc_jac3, fluid_cfg, vg_cfg, rp_cfg);
                    }
                    (void)bc_stats;
                    for (int bi = 0; bi < totalBlocks; ++bi) {
                        for (int eq = 0; eq < 3; ++eq) {
                            const int eqIdx = mgr.getEquationIndex(bi, eq);
                            if (eqIdx < 0 || eqIdx >= totalEq) continue;
                            const double r_bc = bc_res[eqIdx];
                            if (std::abs(r_bc) <= 1e-16) continue;
                            probe_mat.AddResidual(bi, eq, r_bc);
                        }
                    }
                }
                else {
                    auto assembleBoundaryFieldProbe = [&](const BoundarySetting::BoundaryConditionManager* bcMgr, int dofOffset, const std::string& fieldName) {
                        if (!bcMgr || dofOffset < 0) return;
                        std::vector<double> bc_res(totalEq, 0.0);
                        std::vector<std::array<double, 3>> bc_jac3(totalEq, std::array<double, 3>{ 0.0, 0.0, 0.0 });
                        BoundaryAssemblyStats bc_stats;
                        if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                            bc_stats = BoundaryAssembler::Assemble_2D_FullJac(
                                mgr, *bcMgr, dofOffset, fm, fieldName, bc_res, bc_jac3,
                                modules.pressure_bc, modules.saturation_bc, fluid_cfg, vg_cfg, rp_cfg);
                        }
                        else {
                            bc_stats = BoundaryAssembler::Assemble_3D_FullJac(
                                mgr, *bcMgr, dofOffset, fm, fieldName, bc_res, bc_jac3,
                                modules.pressure_bc, modules.saturation_bc, fluid_cfg, vg_cfg, rp_cfg);
                        }
                        (void)bc_stats;
                        for (int bi = 0; bi < totalBlocks; ++bi) {
                            const int eqIdx = mgr.getEquationIndex(bi, dofOffset);
                            if (eqIdx < 0 || eqIdx >= totalEq) continue;
                            const double r_bc = bc_res[eqIdx];
                            if (std::abs(r_bc) <= 1e-16) continue;
                            probe_mat.AddResidual(bi, dofOffset, r_bc);
                        }
                        };
                    assembleBoundaryFieldProbe(modules.pressure_bc, pressureDof, pEqCfg.pressure_field);
                    assembleBoundaryFieldProbe(modules.temperature_bc, temperatureDof, tEqCfg.temperatue_field);
                }

                Eigen::VectorXd b_probe = probe_mat.ExportEigenResidual();
                // 閳光偓閳光偓 ConstantWaterNoConvection: freeze T probe rows to match main Newton system 閳光偓閳光偓
                // Without this, probe_res >> conv_res 閳?LS-BASE-CHECK corrupts use_ref 閳?Armijo fails.
                if constexpr (N == 2) {
                    if (IsSinglePhaseNoConvectionActive(fluid_ctx)) {
                        for (int bi = 0; bi < totalBlocks; ++bi) {
                            const int t_row = mgr.getEquationIndex(bi, 1);
                            if (t_row >= 0 && t_row < static_cast<int>(b_probe.size()))
                                b_probe(t_row) = eval_state.T[bi] - ic.T_init;
                        }
                    }
                }
                int max_probe_idx = 0;
                const double max_probe = b_probe.cwiseAbs().maxCoeff(&max_probe_idx);
                (void)max_probe_idx;
                return max_probe; // [閿熸枻鎷?2] 閿熸枻鎷烽敓鏂ゆ嫹閿熸枻鎷锋帰閿熸枻鎷烽敓鏂ゆ嫹閿熸枻鎷烽敓鏂ゆ嫹鍘熷閿熷彨璇ф嫹 (raw_res) 閿熺殕鎲嬫嫹璇侀敓鏂ゆ嫹閿熸枻鎷烽敓鏂ゆ嫹閿熸枻鎷?
                };

            for (int iter = 0; iter < active_max_newton_iter; ++iter) {
                iter_used++;
                global_mat.SetZero();

                std::vector<EqContrib> eq_contribs(totalEqWithWells);

                // [V3] Diagnostics counters for summary and incident triggers
                const bool track_eos_domain = params.enforce_eos_domain || (params.diag_level != DiagLevel::Off);
                int eos_fallback_water = 0;
                int eos_fallback_co2 = 0;
                int eos_near_bound_count = 0;
                int eos_total_samples = 0;

                double max_mass_flux = 0.0;
                double max_heat_flux = 0.0;
                int hot_mass_i = -1, hot_mass_j = -1;
                int hot_heat_i = -1, hot_heat_j = -1;
                ConnectionType hot_mass_type = ConnectionType::Matrix_Matrix;
                ConnectionType hot_heat_type = ConnectionType::Matrix_Matrix;


                for (int bi = 0; bi < totalBlocks; ++bi) {
                    double phi_ref = 0.2, c_pr = 1000.0, rho_r = 2600.0, c_r = 0.0;
                    if (bi < nMat) {
                        if (phi_mat) phi_ref = (*phi_mat)[bi];
                        if (cp_mat) c_pr = (*cp_mat)[bi];
                        if (rho_mat) rho_r = (*rho_mat)[bi];
                        if (cr_mat) c_r = (*cr_mat)[bi];
                    }
                    else {
                        int fi = bi - nMat;
                        if (phi_frac) phi_ref = (*phi_frac)[fi];
                        if (cp_frac) c_pr = (*cp_frac)[fi];
                        if (rho_frac) rho_r = (*rho_frac)[fi];
                        if (cr_frac) c_r = (*cr_frac)[fi];
                    }
                    ADVar<N> P(state.P[bi]); P.grad(0) = 1.0;
                    ADVar<N> T(state.T[bi]); T.grad((N == 2) ? 1 : 2) = 1.0;
                    auto pW = eval_w_phase(P, T);
                    ADVar<N> P_old(old_state.P[bi]), T_old(old_state.T[bi]);
                    auto pW_old = eval_w_phase(P_old, T_old);
                    if (track_eos_domain) {
                        ++eos_total_samples;
                        if (pW.isFallback) ++eos_fallback_water;
                        if (pW.near_bound) ++eos_near_bound_count;
                    }

                    ADVar<N> phi = ADVar<N>(phi_ref) * (ADVar<N>(1.0) + ADVar<N>(c_r) * (P - ADVar<N>(P_ref[bi])));
                    ADVar<N> phi_old = ADVar<N>(phi_ref) * (ADVar<N>(1.0) + ADVar<N>(c_r) * (P_old - ADVar<N>(P_ref[bi])));

                    std::vector<ADVar<N>> acc_eqs(N);
                    if constexpr (N == 2) {
                        ADVar<N> m_w = pW.rho * phi, m_w_old = pW_old.rho * phi_old;
                        acc_eqs[0] = (m_w - m_w_old) * (vols[bi] / dt);
                        ADVar<N> e_w = m_w * (pW.h - P / pW.rho) + ADVar<N>((1.0 - phi_ref) * rho_r * c_pr) * T;
                        ADVar<N> e_w_old = m_w_old * (pW_old.h - P_old / pW_old.rho) + ADVar<N>((1.0 - phi_ref) * rho_r * c_pr) * T_old;
                        acc_eqs[1] = (e_w - e_w_old) * (vols[bi] / dt);
                    }
                    else {
                        ADVar<N> Sw(state.Sw[bi]); Sw.grad(1) = 1.0;
                        ADVar<N> Sw_old(old_state.Sw[bi]);
                        ADVar<N> Sw_const = ClampSwForConstitutive<N>(Sw, vg_cfg, sw_constitutive_eps);
                        ADVar<N> Sw_old_const = ClampSwForConstitutive<N>(Sw_old, vg_cfg, sw_constitutive_eps);
                        ADVar<N> Sg = ADVar<N>(1.0) - Sw;
                        ADVar<N> Sg_old = ADVar<N>(1.0) - Sw_old;

                        // [DAY6-08] 鐟佸倻绱?Pc=0閿涘苯鐔€鐠愩劋绻氶幐?vG
                        ADVar<N> Pc     = (bi >= nMat) ? ADVar<N>(0.0) : CapRelPerm::pc_vG<N>(Sw_const,     vg_cfg);
                        ADVar<N> Pc_old = (bi >= nMat) ? ADVar<N>(0.0) : CapRelPerm::pc_vG<N>(Sw_old_const, vg_cfg);
                        ADVar<N> Pg = P + Pc;
                        ADVar<N> Pg_old = P_old + Pc_old;

                        auto pG = eval_g_phase(Pg, T);
                        auto pG_old = eval_g_phase(Pg_old, T_old);
                        if (track_eos_domain) {
                            ++eos_total_samples;
                            if (pG.isFallback) ++eos_fallback_co2;
                            if (pG.near_bound) ++eos_near_bound_count;
                        }

                        acc_eqs[0] = (pW.rho * phi * Sw - pW_old.rho * phi_old * Sw_old) * (vols[bi] / dt);
                        acc_eqs[1] = (pG.rho * phi * Sg - pG_old.rho * phi_old * Sg_old) * (vols[bi] / dt);

                        ADVar<N> e_fluid = pW.rho * Sw * (pW.h - P / pW.rho) + pG.rho * Sg * (pG.h - Pg / pG.rho);
                        ADVar<N> e_fluid_old = pW_old.rho * Sw_old * (pW_old.h - P_old / pW_old.rho) + pG_old.rho * Sg_old * (pG_old.h - Pg_old / pG_old.rho);
                        ADVar<N> e_rock = ADVar<N>((1.0 - phi_ref) * rho_r * c_pr) * T;
                        ADVar<N> e_rock_old = ADVar<N>((1.0 - phi_ref) * rho_r * c_pr) * T_old;
                        acc_eqs[2] = ((e_fluid * phi + e_rock) - (e_fluid_old * phi_old + e_rock_old)) * (vols[bi] / dt);
                    }
                    FIM_GlobalAssembler<N, ADVar<N>>::AssembleAccumulation(bi, acc_eqs, global_mat);
                    // Keep accumulation contributions always available for row scaling / diagnostics.
                    for (int eq = 0; eq < N; ++eq) {
                        int g_eq = mgr.getEquationIndex(bi, eq);
                        if (g_eq >= 0 && g_eq < eq_contribs.size()) {
                            eq_contribs[g_eq].R_acc += acc_eqs[eq].val;
                            eq_contribs[g_eq].D_acc += acc_eqs[eq].grad[eq];
                        }
                    }
                }

                struct ConnAuditStat {
                    int conn_count = 0;
                    double sum_abs_tflow = 0.0;
                    double sum_abs_theat = 0.0;
                    int offdiag_nnz = 0;
                    double sum_abs_offdiag = 0.0;
                };
                std::map<ConnectionType, ConnAuditStat> audit_stats;

                // [DAY6-09] 閸楁洝鍑禒锝囧⒖閹呯处鐎涙﹫绱版０鍕吀缁犳鎮嘽ell閺嶅洭鍣洪悧鈺傗偓褝绱濋柆鍨帳passive娓氀冨晳娴ｆOS鐠嬪啰鏁?
                std::vector<detail::FluxCellPropsCache> cprop(totalBlocks);
                for (int ci = 0; ci < totalBlocks; ++ci) {
                    const double p_ci = state.P[ci], t_ci = state.T[ci];
                    {
                        ADVar<N> P0(p_ci), T0(t_ci);
                        auto pw = eval_w_phase(P0, T0);
                        cprop[ci].rho_w = pw.rho.val;
                        cprop[ci].mu_w  = pw.mu.val;
                        cprop[ci].h_w   = pw.h.val;
                    }
                    if constexpr (N == 3) {
                        const double sw_ci = state.Sw[ci];
                        double pc_s = 0.0, dpc_dsw_s = 0.0, krw_s = 1.0, krg_s = 0.0;
                        if (ci >= nMat) {
                            const double sw_c = (sw_ci < 0.0) ? 0.0 : (sw_ci > 1.0) ? 1.0 : sw_ci;
                            krw_s = sw_c; krg_s = 1.0 - sw_c;
                        } else {
                            ADVar<N> SwC(sw_ci);
                            SwC.grad(1) = 1.0;
                            ADVar<N> SwCC = ClampSwForConstitutive<N>(SwC, vg_cfg, sw_constitutive_eps);
                            ADVar<N> pc_ad = CapRelPerm::pc_vG<N>(SwCC, vg_cfg);
                            pc_s  = pc_ad.val;
                            dpc_dsw_s = pc_ad.grad(1);
                            ADVar<N> krw0, krg0;
                            CapRelPerm::kr_Mualem_vG<N>(SwCC, vg_cfg, rp_cfg, krw0, krg0);
                            krw_s = krw0.val; krg_s = krg0.val;
                        }
                        cprop[ci].Pc  = pc_s;
                        cprop[ci].dPc_dSw = dpc_dsw_s;
                        cprop[ci].krw = krw_s;
                        cprop[ci].krg = krg_s;
                        {
                            ADVar<N> Pg0(p_ci + pc_s), T0(t_ci);
                            auto pg = eval_g_phase(Pg0, T0);
                            cprop[ci].rho_g = pg.rho.val;
                            cprop[ci].mu_g  = pg.mu.val;
                            cprop[ci].h_g   = pg.h.val;
                        }
                    }
                }

                // Pre-compute cell-centre gradients for deferred non-orth correction.
                // Align N=2/N=3 with the corrected N=1 style:
                // prefer Green-Gauss + BC handling, fallback to geometric accumulation.
                std::vector<Vector> nonorth_grad_P(nMat, Vector(0.0, 0.0, 0.0));
                std::vector<Vector> nonorth_grad_T_cell(nMat, Vector(0.0, 0.0, 0.0));
                std::vector<Vector> nonorth_grad_Sw(nMat, Vector(0.0, 0.0, 0.0));
                if (params.enable_non_orthogonal_correction && nMat > 0) {
                    detail::PrepareNonOrthogonalGradients(
                        mgr, state, nMat, vols,
                        modules.pressure_bc, modules.temperature_bc, modules.saturation_bc,
                        internalNames.nonorth_main_p_tmp,
                        internalNames.nonorth_main_t_tmp,
                        internalNames.nonorth_main_sw_tmp,
                        nonorth_grad_P, nonorth_grad_T_cell, nonorth_grad_Sw);
                }
                for (const auto& conn : connMgr.GetConnections()) {
                    audit_stats[conn.type].conn_count++;
                    audit_stats[conn.type].sum_abs_tflow += std::abs(conn.T_Flow);
                    audit_stats[conn.type].sum_abs_theat += std::abs(conn.T_Heat);

                    int i = conn.nodeI, j = conn.nodeJ;
                    auto f_wrt_i = detail::EvaluateConnectionFlux(
                        conn, i, j, nMat, state, blockCenters, gravityVec, true, &cprop,
                        eval_w_phase, eval_g_phase, vg_cfg, rp_cfg, sw_constitutive_eps);
                    auto f_wrt_j = detail::EvaluateConnectionFlux(
                        conn, i, j, nMat, state, blockCenters, gravityVec, false, &cprop,
                        eval_w_phase, eval_g_phase, vg_cfg, rp_cfg, sw_constitutive_eps);
                    FIM_GlobalAssembler<N, ADVar<N>>::AssembleFlux(i, j, f_wrt_i, f_wrt_j, global_mat);

                    // 閳光偓閳光偓 Step 2: non-orthogonal deferred correction 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
                    if (params.enable_non_orthogonal_correction &&
                        conn.type == ConnectionType::Matrix_Matrix &&
                        i < nMat && j < nMat)
                    {
                        const double vT_mag = conn.vectorT.Mag();
                        if (vT_mag > 1e-15) {
                            // Face-interpolated cell-centre gradients (arithmetic mean)
                            const double gP_x = 0.5*(nonorth_grad_P[i].m_x + nonorth_grad_P[j].m_x);
                            const double gP_y = 0.5*(nonorth_grad_P[i].m_y + nonorth_grad_P[j].m_y);
                            const double gP_z = 0.5*(nonorth_grad_P[i].m_z + nonorth_grad_P[j].m_z);
                            const double gT_x = 0.5*(nonorth_grad_T_cell[i].m_x + nonorth_grad_T_cell[j].m_x);
                            const double gT_y = 0.5*(nonorth_grad_T_cell[i].m_y + nonorth_grad_T_cell[j].m_y);
                            const double gT_z = 0.5*(nonorth_grad_T_cell[i].m_z + nonorth_grad_T_cell[j].m_z);
                            // Tangential correction scalars: grad 璺?vectorT  [Pa], [K]
                            const double corr_P = gP_x * conn.vectorT.m_x + gP_y * conn.vectorT.m_y + gP_z * conn.vectorT.m_z;
                            const double corr_T = gT_x * conn.vectorT.m_x + gT_y * conn.vectorT.m_y + gT_z * conn.vectorT.m_z;
                            // Effective permeability: K_eff = T_Flow * dist / |vectorE|
                            const double K_eff = (conn.aux_area > 1e-12)
                                ? conn.T_Flow * conn.aux_dist / conn.aux_area : 0.0;
                            // Upwind water mobility (from per-cell cache; potential-based)
                            const double g_dx_w =
                                gravityVec * (blockCenters[j] - blockCenters[i]);
                            const double dPhi_w_sign = (state.P[j] - state.P[i])
                                - 0.5*(cprop[i].rho_w + cprop[j].rho_w) * g_dx_w;
                            const double mob_w = (dPhi_w_sign < 0.0)
                                ? cprop[i].krw * cprop[i].rho_w / std::max(cprop[i].mu_w, 1e-20)
                                : cprop[j].krw * cprop[j].rho_w / std::max(cprop[j].mu_w, 1e-20);
                            // Correction contributions (explicit 閳?residual only, no Jacobian)
                            const double mass_w_corr = K_eff * mob_w * corr_P;
                            const double heat_corr   = conn.T_Heat * corr_T;
                            const double h_w_upwind = (dPhi_w_sign < 0.0) ? cprop[i].h_w : cprop[j].h_w;
                            double conv_heat_corr = mass_w_corr * h_w_upwind;
                            global_mat.AddResidual(i, 0,  -mass_w_corr);
                            global_mat.AddResidual(j, 0,   mass_w_corr);
                            if constexpr (N == 3) {
                                // CO2 mass correction
                                const double dPhi_g_sign = (state.P[j] - state.P[i])
                                    + (cprop[j].Pc - cprop[i].Pc)
                                    - 0.5*(cprop[i].rho_g + cprop[j].rho_g) * g_dx_w;
                                const double gSw_x = 0.5 * (nonorth_grad_Sw[i].m_x + nonorth_grad_Sw[j].m_x);
                                const double gSw_y = 0.5 * (nonorth_grad_Sw[i].m_y + nonorth_grad_Sw[j].m_y);
                                const double gSw_z = 0.5 * (nonorth_grad_Sw[i].m_z + nonorth_grad_Sw[j].m_z);
                                const double corr_Sw = gSw_x * conn.vectorT.m_x + gSw_y * conn.vectorT.m_y + gSw_z * conn.vectorT.m_z;
                                const double dPc_dSw = 0.5 * (cprop[i].dPc_dSw + cprop[j].dPc_dSw);
                                const double corr_Pc = corr_Sw * dPc_dSw;
                                const double mob_g = (dPhi_g_sign < 0.0)
                                    ? cprop[i].krg * cprop[i].rho_g / std::max(cprop[i].mu_g, 1e-20)
                                    : cprop[j].krg * cprop[j].rho_g / std::max(cprop[j].mu_g, 1e-20);
                                const double mass_g_corr = K_eff * mob_g * (corr_P + corr_Pc);
                                const double h_g_upwind = (dPhi_g_sign < 0.0) ? cprop[i].h_g : cprop[j].h_g;
                                conv_heat_corr += mass_g_corr * h_g_upwind;
                                global_mat.AddResidual(i, 1, -mass_g_corr);
                                global_mat.AddResidual(j, 1,  mass_g_corr);
                            }
                            const double total_heat_corr = heat_corr + conv_heat_corr;
                            global_mat.AddResidual(i, N - 1, -total_heat_corr);
                            global_mat.AddResidual(j, N - 1,  total_heat_corr);
                        }
                    }
                    // 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

                    const double abs_mass_flux = [&]() {
                        if constexpr (N == 2) return std::abs(f_wrt_i[0].val);
                        else return std::max(std::abs(f_wrt_i[0].val), std::abs(f_wrt_i[1].val));
                        }();
                    const double abs_heat_flux = [&]() {
                        if constexpr (N == 2) return std::abs(f_wrt_i[1].val);
                        else return std::abs(f_wrt_i[2].val);
                        }();

                    if (abs_mass_flux > max_mass_flux) {
                        max_mass_flux = abs_mass_flux;
                        hot_mass_i = i;
                        hot_mass_j = j;
                        hot_mass_type = conn.type;
                    }
                    if (abs_heat_flux > max_heat_flux) {
                        max_heat_flux = abs_heat_flux;
                        hot_heat_i = i;
                        hot_heat_j = j;
                        hot_heat_type = conn.type;
                    }

                    if (params.diag_level != DiagLevel::Off) {
                        for (int eq = 0; eq < N; ++eq) {
                            int g_eq_i = mgr.getEquationIndex(i, eq);
                            int g_eq_j = mgr.getEquationIndex(j, eq);
                            if (g_eq_i >= 0 && g_eq_i < eq_contribs.size()) {
                                eq_contribs[g_eq_i].R_flux -= f_wrt_i[eq].val;
                                eq_contribs[g_eq_i].D_flux -= f_wrt_i[eq].grad[eq];
                            }
                            if (g_eq_j >= 0 && g_eq_j < eq_contribs.size()) {
                                eq_contribs[g_eq_j].R_flux += f_wrt_j[eq].val;
                                eq_contribs[g_eq_j].D_flux += f_wrt_j[eq].grad[eq];
                            }
                        }
                    }
                }


                SyncStateToFieldManager(state, fm, mgr, sp_model, fluid_cfg, vg_cfg, rp_cfg);
                std::vector<double> w_res(totalEq, 0.0);
                std::vector<std::array<double, 3>> w_jac3(totalEq, std::array<double, 3>{ 0.0, 0.0, 0.0 });
                std::vector<WellCompletionLinearization> w_lin;
                const auto well_bhp_now = well_mgr.BuildWellBhpMap(state);
                const int well_dof_w = (N == 3) ? 0 : (sp_use_co2 ? -1 : 0);
                const int well_dof_g = (N == 3) ? 1 : (sp_use_co2 ? 0 : -1);
                const int well_dof_e = (N == 3) ? 2 : 1;

                if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                    BoundaryAssembler::Assemble_Wells_2D_FullJac(
                        mgr, fm, active_wells, 0, well_dof_w, well_dof_g, well_dof_e,
                        w_res, w_jac3, fluid_cfg, vg_cfg, rp_cfg, &w_lin, &well_bhp_now);
                }
                else {
                    BoundaryAssembler::Assemble_Wells_3D_FullJac(
                        mgr, fm, active_wells, 0, well_dof_w, well_dof_g, well_dof_e,
                        w_res, w_jac3, fluid_cfg, vg_cfg, rp_cfg, &w_lin, &well_bhp_now);
                }

                double max_abs_well_dsw = 0.0;
                double max_abs_well_dt = 0.0;

                // BoundaryAssembler returns well terms in outflow-positive convention.
                // Residual here is assembled as Acc - Flux_inflow + Q_out = 0.
                // Keep sign configurable for compatibility with custom conventions.
                const double kWellSourceSign = params.well_source_sign;

                for (int bi = 0; bi < totalBlocks; ++bi) {
                    for (int eq = 0; eq < N; ++eq) {
                        int g_eq = mgr.getEquationIndex(bi, eq);
                        if (g_eq < 0 || g_eq >= totalEq) continue;

                        const double dRdP = kWellSourceSign * w_jac3[g_eq][0];
                        const double dRdSw = kWellSourceSign * w_jac3[g_eq][1];
                        const double dRdT = kWellSourceSign * w_jac3[g_eq][2];
                        const double rWell = kWellSourceSign * w_res[g_eq];

                        max_abs_well_dsw = std::max(max_abs_well_dsw, std::abs(dRdSw));
                        max_abs_well_dt = std::max(max_abs_well_dt, std::abs(dRdT));

                        if (std::abs(rWell) <= 1e-16 &&
                            std::abs(dRdP) <= 1e-16 &&
                            std::abs(dRdSw) <= 1e-16 &&
                            std::abs(dRdT) <= 1e-16) {
                            continue;
                        }

                        global_mat.AddResidual(bi, eq, rWell);
                        global_mat.AddDiagJacobian(bi, eq, 0, dRdP);

                        if constexpr (N == 2) {
                            global_mat.AddDiagJacobian(bi, eq, 1, dRdT);
                        }
                        else {
                            global_mat.AddDiagJacobian(bi, eq, 1, dRdSw);
                            global_mat.AddDiagJacobian(bi, eq, 2, dRdT);
                        }

                        if (params.diag_level != DiagLevel::Off && g_eq >= 0 && g_eq < eq_contribs.size()) {
                            eq_contribs[g_eq].R_well += rWell;
                            if (eq == 0) eq_contribs[g_eq].D_well += dRdP;
                            else if constexpr (N == 2) { if (eq == 1) eq_contribs[g_eq].D_well += dRdT; }
                            else if constexpr (N == 3) {
                                if (eq == 1) eq_contribs[g_eq].D_well += dRdSw;
                                else if (eq == 2) eq_contribs[g_eq].D_well += dRdT;
                            }
                        }

                    }
                }

                if constexpr (N == 3) {
                    std::cout << "    [WellJac] max|dR/dSw|=" << std::scientific << max_abs_well_dsw
                        << " max|dR/dT|=" << max_abs_well_dt << "\n";
                }
                else {
                    std::cout << "    [WellJac] max|dR/dT|=" << std::scientific << max_abs_well_dt << "\n";
                }

                // 閳光偓閳光偓 Step 3: assemble independent well DOF equations 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
                if (!well_mgr.Empty()) {
                    well_mgr.AssembleWellEquations(
                        global_mat, state, active_wells,
                        w_lin, kWellSourceSign);
                }
                // 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

                if constexpr (N == 3) {
                    std::vector<double> bc_res(totalEq, 0.0);
                    std::vector<std::array<double, 3>> bc_jac3(totalEq, std::array<double, 3>{ 0.0, 0.0, 0.0 });
                    BoundaryAssemblyStats bc_stats;
                    if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                        bc_stats = BoundaryAssembler::Assemble_2D_CoupledN3_FullJac(
                            mgr, fm, modules.pressure_bc, modules.saturation_bc, modules.temperature_bc,
                            pressureDof, pressureDof, saturationDof, temperatureDof,
                            bc_res, bc_jac3, fluid_cfg, vg_cfg, rp_cfg);
                    }
                    else {
                        bc_stats = BoundaryAssembler::Assemble_3D_CoupledN3_FullJac(
                            mgr, fm, modules.pressure_bc, modules.saturation_bc, modules.temperature_bc,
                            pressureDof, pressureDof, saturationDof, temperatureDof,
                            bc_res, bc_jac3, fluid_cfg, vg_cfg, rp_cfg);
                    }

                    int appliedEq = 0;
                    for (int bi = 0; bi < totalBlocks; ++bi) {
                        for (int eq = 0; eq < 3; ++eq) {
                            const int eqIdx = mgr.getEquationIndex(bi, eq);
                            if (eqIdx < 0 || eqIdx >= totalEq) continue;
                            const double r_bc = bc_res[eqIdx];
                            const double dRdP = bc_jac3[eqIdx][0];
                            const double dRdSw = bc_jac3[eqIdx][1];
                            const double dRdT = bc_jac3[eqIdx][2];
                            if (std::abs(r_bc) <= 1e-16 &&
                                std::abs(dRdP) <= 1e-16 &&
                                std::abs(dRdSw) <= 1e-16 &&
                                std::abs(dRdT) <= 1e-16) {
                                continue;
                            }

                            global_mat.AddResidual(bi, eq, r_bc);
                            global_mat.AddDiagJacobian(bi, eq, 0, dRdP);
                            global_mat.AddDiagJacobian(bi, eq, 1, dRdSw);
                            global_mat.AddDiagJacobian(bi, eq, 2, dRdT);
                            ++appliedEq;

                            if (params.diag_level != DiagLevel::Off && eqIdx >= 0 && eqIdx < eq_contribs.size()) {
                                eq_contribs[eqIdx].R_bc += r_bc;
                                if (eq == 0) eq_contribs[eqIdx].D_bc += dRdP;
                                else if (eq == 1) eq_contribs[eqIdx].D_bc += dRdSw;
                                else eq_contribs[eqIdx].D_bc += dRdT;
                            }
                        }
                    }

                    if (params.diag_level != DiagLevel::Off) {
                        std::cout << "    [BC-SUM] field=N3-Coupled"
                            << " faces=" << bc_stats.matrixBCCount + bc_stats.fractureBCCount
                            << " applied_eq=" << appliedEq
                            << " (visited=" << bc_stats.visitedEqRows
                            << ", nonzero=" << bc_stats.nonzeroEqRows
                            << ", zero_row=" << bc_stats.zeroEqRows
                            << ", invalid=" << bc_stats.invalidEqRows
                            << ", infer_tag_fail=" << bc_stats.inferredTagFail
                            << ", sat_comp_used=" << bc_stats.satCompositionUsed
                            << ", sat_drive_ignored=" << bc_stats.satDriveIgnored << ")\n";
                    }
                    if (bc_stats.satDriveIgnored > 0) {
                        std::cout << "    [BC-WARN] N=3 saturation_bc Neumann/Robin is ignored as independent drive under pressure-driven-two-phase semantics; count="
                            << bc_stats.satDriveIgnored << "\n";
                    }
                }
                else {
                    auto assembleBoundaryField = [&](const BoundarySetting::BoundaryConditionManager* bcMgr,
                        int dofOffset,
                        const std::string& fieldName,
                        const char* fieldLabel) {
                            if (!bcMgr || dofOffset < 0) return;

                            std::vector<double> bc_res(totalEq, 0.0);
                            std::vector<std::array<double, 3>> bc_jac3(totalEq, std::array<double, 3>{ 0.0, 0.0, 0.0 });
                            BoundaryAssemblyStats bc_stats;

                            if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                                bc_stats = BoundaryAssembler::Assemble_2D_FullJac(
                                    mgr, *bcMgr, dofOffset, fm, fieldName, bc_res, bc_jac3,
                                    modules.pressure_bc, modules.saturation_bc, fluid_cfg, vg_cfg, rp_cfg);
                            }
                            else {
                                bc_stats = BoundaryAssembler::Assemble_3D_FullJac(
                                    mgr, *bcMgr, dofOffset, fm, fieldName, bc_res, bc_jac3,
                                    modules.pressure_bc, modules.saturation_bc, fluid_cfg, vg_cfg, rp_cfg);
                            }

                            int appliedEq = 0;
                            for (int bi = 0; bi < totalBlocks; ++bi) {
                                const int eqIdx = mgr.getEquationIndex(bi, dofOffset);
                                if (eqIdx < 0 || eqIdx >= totalEq) continue;

                                const double r_bc = bc_res[eqIdx];
                                const double dRdP = bc_jac3[eqIdx][0];
                                const double dRdT = bc_jac3[eqIdx][2];
                                if (std::abs(r_bc) <= 1e-16 &&
                                    std::abs(dRdP) <= 1e-16 &&
                                    std::abs(dRdT) <= 1e-16) {
                                    continue;
                                }

                                global_mat.AddResidual(bi, dofOffset, r_bc);
                                global_mat.AddDiagJacobian(bi, dofOffset, 0, dRdP);
                                global_mat.AddDiagJacobian(bi, dofOffset, 1, dRdT);
                                ++appliedEq;

                                if (params.diag_level != DiagLevel::Off && eqIdx >= 0 && eqIdx < eq_contribs.size()) {
                                    eq_contribs[eqIdx].R_bc += r_bc;
                                    eq_contribs[eqIdx].D_bc += (dofOffset == 0) ? dRdP : dRdT;
                                }
                            }
                            if (params.diag_level != DiagLevel::Off) {
                                std::cout << "    [BC-SUM] field=" << fieldLabel
                                    << " faces=" << bc_stats.matrixBCCount + bc_stats.fractureBCCount
                                    << " applied_eq=" << appliedEq
                                    << " (visited=" << bc_stats.visitedEqRows
                                    << ", nonzero=" << bc_stats.nonzeroEqRows
                                    << ", zero_row=" << bc_stats.zeroEqRows
                                    << ", invalid=" << bc_stats.invalidEqRows
                                    << ", infer_tag_fail=" << bc_stats.inferredTagFail << ")\n";
                            }
                        };

                    assembleBoundaryField(modules.pressure_bc, pressureDof, pEqCfg.pressure_field, pEqCfg.pressure_field.c_str());
                    assembleBoundaryField(modules.temperature_bc, temperatureDof, tEqCfg.temperatue_field, tEqCfg.temperatue_field.c_str());
                }
                auto A = global_mat.GetFrozenMatrix(); // Issue#11: O(nnz) CSR value-update, no Triplet sort
                auto b = global_mat.ExportEigenResidual();

                // 閳光偓閳光偓 ConstantWaterNoConvection: freeze T DOF 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
                // Replace the energy equation row with a Dirichlet constraint T = T_init.
                // R_T = T - T_init, J_TT = 1, all off-diagonals = 0.
                // This makes the energy equation trivially satisfied every Newton step
                // (T stays at T_init), isolating the pressure equation for clean convergence.
                if constexpr (N == 2) {
                    if (IsSinglePhaseNoConvectionActive(fluid_ctx)) {
                        using SpMat = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;
                        for (int bi = 0; bi < totalBlocks; ++bi) {
                            const int t_row = mgr.getEquationIndex(bi, 1);
                            if (t_row < 0 || t_row >= static_cast<int>(b.size())) continue;
                            for (SpMat::InnerIterator it(A, t_row); it; ++it)
                                it.valueRef() = 0.0;
                            A.coeffRef(t_row, t_row) = 1.0;
                            b(t_row) = state.T[bi] - ic.T_init;
                            // Reset D_acc so row scaling uses denominator=1 for T rows
                            // (original D_acc_energy is large 閳?would scale T row to ~0 閳?no precision gain)
                            if (t_row < static_cast<int>(eq_contribs.size()))
                                eq_contribs[t_row].D_acc = 1.0;
                        }
                    }
                }

                if (params.enable_matrix_audit && step == params.matrix_audit_step && iter_used == params.matrix_audit_iter) {
                    for (const auto& conn : connMgr.GetConnections()) {
                        int i = conn.nodeI, j = conn.nodeJ;
                        for (int eq = 0; eq < N; ++eq) {
                            int g_eq_i = mgr.getEquationIndex(i, eq);
                            int g_eq_j = mgr.getEquationIndex(j, eq);
                            for (int var = 0; var < N; ++var) {
                                int g_var_i = mgr.getEquationIndex(i, var);
                                int g_var_j = mgr.getEquationIndex(j, var);

                                if (g_eq_i >= 0 && g_var_j >= 0) {
                                    double val = A.coeff(g_eq_i, g_var_j);
                                    if (std::abs(val) > 0.0) {
                                        audit_stats[conn.type].offdiag_nnz++;
                                        audit_stats[conn.type].sum_abs_offdiag += std::abs(val);
                                    }
                                }
                                if (g_eq_j >= 0 && g_var_i >= 0) {
                                    double val = A.coeff(g_eq_j, g_var_i);
                                    if (std::abs(val) > 0.0) {
                                        audit_stats[conn.type].offdiag_nnz++;
                                        audit_stats[conn.type].sum_abs_offdiag += std::abs(val);
                                    }
                                }
                            }
                        }
                    }

                    for (const auto& kv : audit_stats) {
                        std::cout << "    [MATRIX-AUDIT] type=" << ConnectionTypeLabel(kv.first)
                            << " conn=" << kv.second.conn_count
                            << " sum|Tflow|=" << std::scientific << kv.second.sum_abs_tflow
                            << " sum|Theat|=" << kv.second.sum_abs_theat
                            << " offdiag_nnz=" << kv.second.offdiag_nnz
                            << " sum|Aoff|=" << kv.second.sum_abs_offdiag << "\n";
                    }

                    if (params.matrix_audit_strict) {
                        auto check_conn = [&](ConnectionType type, const std::string& name, bool required) {
                            if (!required) return;
                            auto it = audit_stats.find(type);
                            if (it == audit_stats.end() || it->second.conn_count == 0) {
                                throw std::runtime_error("[MATRIX-AUDIT-FAIL] Strict mode requires " + name + " connections, but none found.");
                            }
                            if ((it->second.sum_abs_tflow + it->second.sum_abs_theat) <= params.matrix_audit_eps) {
                                throw std::runtime_error("[MATRIX-AUDIT-FAIL] " + name + " TFlow/THeat sum is zero or extremely small.");
                            }
                            if (it->second.offdiag_nnz == 0 || it->second.sum_abs_offdiag <= params.matrix_audit_eps) {
                                throw std::runtime_error("[MATRIX-AUDIT-FAIL] " + name + " coefficients were not successfully assembled into off-diagonal Jacobian.");
                            }
                            };
                        check_conn(ConnectionType::Matrix_Fracture, "MF(NNC)", params.matrix_audit_require_nnc);
                        check_conn(ConnectionType::Fracture_Fracture, "FF", params.matrix_audit_require_ff);
                        std::cout << "    [MATRIX-AUDIT-PASS] Strict audit completed successfully.\n";
                    }
                }

                if (params.enable_matrix_audit &&
                    step == std::max(1, params.matrix_audit_step) &&
                    iter_used == std::max(1, params.matrix_audit_iter)) {
                    const auto audit = RunMatrixAssemblyAudit<N>(
                        A,
                        connMgr.GetConnections(),
                        params.matrix_audit_coeff_tol,
                        params.matrix_audit_max_detail);

                    std::cout << "    [MATRIX-AUDIT] checked=" << audit.total_checked
                        << " nnc=" << audit.nnc_checked
                        << " ff=" << audit.ff_checked
                        << " coupled_flow(nnc/ff)=" << audit.nnc_flow_coupled << "/" << audit.ff_flow_coupled
                        << " coupled_heat(nnc/ff)=" << audit.nnc_heat_coupled << "/" << audit.ff_heat_coupled
                        << " miss_flow(nnc/ff)=" << audit.nnc_missing_flow << "/" << audit.ff_missing_flow
                        << " miss_heat(nnc/ff)=" << audit.nnc_missing_heat << "/" << audit.ff_missing_heat
                        << " zeroTflow(nnc/ff)=" << audit.nnc_zero_tflow << "/" << audit.ff_zero_tflow
                        << " zeroTheat(nnc/ff)=" << audit.nnc_zero_theat << "/" << audit.ff_zero_theat
                        << " strict=" << (params.matrix_audit_strict ? "true" : "false") << "\n";

                    if (!audit.Passed(params.matrix_audit_strict)) {
                        throw std::runtime_error("[FAIL] matrix_audit_strict: missing NNC/FF coupling detected or no NNC/FF connection found.");
                    }
                }

                double ptc_lambda_iter = 0.0;
                if (active_enable_ptc) {
                    const double decay = std::max(0.0, std::min(1.0, active_profile.ptc_lambda_decay));
                    const double base_ptc = std::max(active_profile.ptc_lambda_min,
                        active_profile.ptc_lambda_init * std::pow(decay, static_cast<double>(iter)));
                    const double rescue_boost = std::max(1.0, ptc_rescue_boost_cache);
                    ptc_lambda_iter = base_ptc * rescue_boost;
                    if (ptc_lambda_iter > 0.0) {
                        auto select_ptc_sign = [](double d_acc, double diag_now) -> double {
                            if (d_acc > 0.0) return 1.0;
                            if (d_acc < 0.0) return -1.0;
                            if (diag_now > 0.0) return 1.0;
                            if (diag_now < 0.0) return -1.0;
                            return 1.0;
                            };

                        const double row_floor = std::max(std::abs(params.row_scale_floor), 1.0e-30);
                        for (int r = 0; r < totalEq; ++r) {
                            const double d_acc = (r >= 0 && r < static_cast<int>(eq_contribs.size())) ? eq_contribs[r].D_acc : 0.0;
                            const double m_ptc = std::max(std::abs(d_acc), row_floor);
                            const double sign_ptc = select_ptc_sign(d_acc, A.coeff(r, r));
                            A.coeffRef(r, r) += sign_ptc * ptc_lambda_iter * m_ptc;
                        }
                    }
                }

                int max_idx = 0;
                double max_res = b.cwiseAbs().maxCoeff(&max_idx);
                double max_res_scaled = max_res;
                if (params.enable_row_scaling) {
                    max_res_scaled = 0.0;
                    for (int r = 0; r < totalEq; ++r) {
                        const double diag_acc_abs = (r >= 0 && r < static_cast<int>(eq_contribs.size())) ? std::abs(eq_contribs[r].D_acc) : 0.0;
                        const double diag_abs = std::abs(A.coeff(r, r));
                        const double denom = std::max({ diag_acc_abs, diag_abs, params.row_scale_floor });
                        const double scaled = std::abs(b[r]) / std::max(denom, 1.0e-30);
                        if (scaled > max_res_scaled) {
                            max_res_scaled = scaled;
                        }
                    }
                }
                const double conv_res = max_res;

                if (iter == 0) res_iter1 = conv_res;
                step_final_residual = conv_res;
                if (conv_res < best_res) {
                    best_res = conv_res;
                    best_state = state;
                    best_iter = iter_used;
                }

                std::cout << "    [NL] step=" << step << " iter=" << iter_used
                    << " dt=" << std::scientific << dt << " ramp=" << control_ramp << " ptc=" << ptc_lambda_iter << " res_inf=" << max_res
                    << " res_scaled=" << max_res_scaled
                    << " (at DOF=" << max_idx << ")\n";

                if (!std::isfinite(conv_res)) { fail_reason = "residual_nan_inf"; break; }
                if (params.enforce_eos_domain &&
                    (eos_fallback_water > 0 || eos_fallback_co2 > 0 || eos_near_bound_count > 0)) {
                    fail_reason = "eos_domain_violation";
                    std::cout << "    [EOS-STRICT] fallback_w=" << eos_fallback_water
                        << " fallback_g=" << eos_fallback_co2
                        << " near_bound=" << eos_near_bound_count << "/" << eos_total_samples
                        << " -> reject step\n";
                    break;
                }
                if (conv_res < params.abs_res_tol) {
                    converged = true; converge_mode = "abs_res";
                    std::cout << "    [NL-CONVERGED] step=" << step << " mode=" << converge_mode << " res_scaled=" << conv_res << "\n";
                    break;
                }
                const double rel_res = (res_iter1 > 0.0) ? (conv_res / std::max(res_iter1, 1e-30)) : 1.0;
                if (iter_used > 1 && rel_res <= active_rel_res_tol && last_rel_update <= active_rel_update_tol) {
                    converged = true; converge_mode = "rel_res_update";
                    std::cout << "    [NL-CONVERGED] step=" << step << " mode=" << converge_mode
                        << " rel_res=" << rel_res << " rel_update=" << last_rel_update << "\n";
                    break;
                }
                // General best-iterate guard:
                // accept the in-step best state when later Newton iterations
                // clearly move away from it after an already acceptable drop.
                if (params.enable_best_iter_guard && iter_used >= params.best_iter_guard_min_iter && best_iter > 1 && res_iter1 > 0.0) {
                    const double grow_from_best = conv_res / std::max(best_res, 1e-30);
                    const double best_growth = best_res / std::max(res_iter1, 1e-30);
                    const double best_drop = 1.0 - best_growth;
                    const bool best_quality_ok = (best_res <= params.stagnation_abs_res_tol) &&
                        (best_res < res_iter1) &&
                        (best_drop >= params.stagnation_min_drop);
                    if (best_quality_ok && grow_from_best > params.best_iter_growth_trigger) {
                        state = best_state;
                        step_final_residual = best_res;
                        converged = true;
                        converge_mode = "best_iter_guard";
                        std::cout << "    [NL-CONVERGED] step=" << step
                            << " mode=" << converge_mode
                            << " best_iter=" << best_iter
                            << " best_res=" << best_res
                            << " grow_from_best=" << grow_from_best << "\n";
                        break;
                    }
                }
                // At dt floor, recover best iterate in-step when later Newton iterations
                // start blowing up after an already acceptable decrease.
                const bool at_dt_floor = (dt <= std::max(params.dt_min, 1.0e-18) * (1.0 + 1.0e-12));
                if (params.enable_dt_floor_best && at_dt_floor && iter_used >= 3 && best_iter > 1) {
                    const double grow_from_best = conv_res / std::max(best_res, 1e-30);
                    const double best_growth = (res_iter1 > 0.0) ? (best_res / std::max(res_iter1, 1e-30)) : 1.0;
                    const double best_drop = 1.0 - best_growth;
                    const bool best_quality_ok = (best_res <= params.stagnation_abs_res_tol) &&
                        (best_res < res_iter1) &&
                        (best_drop >= params.stagnation_min_drop);
                    const bool accept_best = (grow_from_best > 2.0) && best_quality_ok;
                    if (accept_best) {
                        state = best_state;
                        step_final_residual = best_res;
                        converged = true;
                        converge_mode = "dt_floor_best";
                        std::cout << "    [NL-CONVERGED] step=" << step
                            << " mode=" << converge_mode
                            << " best_iter=" << best_iter
                            << " best_res=" << best_res
                            << " grow_from_best=" << grow_from_best << "\n";
                        break;
                    }
                }
                // Guarded dt-floor hold:
                // when residual at iter-1 is already within stagnation tolerance but Newton
                // corrections keep increasing residual, accept iter-1 state and move on.
                // This avoids repeated rollback loops at dt_min while still rejecting high residuals.
                if (params.enable_dt_floor_hold && at_dt_floor && iter_used >= 2 && best_iter == 1) {
                    const double grow_from_best = conv_res / std::max(best_res, 1e-30);
                    const bool best_quality_ok = (best_res <= params.stagnation_abs_res_tol);
                    const bool accept_hold = best_quality_ok && (grow_from_best > params.best_iter_growth_trigger);
                    if (accept_hold) {
                        state = best_state;
                        step_final_residual = best_res;
                        converged = true;
                        converge_mode = "dt_floor_hold";
                        std::cout << "    [NL-CONVERGED] step=" << step
                            << " mode=" << converge_mode
                            << " best_iter=" << best_iter
                            << " best_res=" << best_res
                            << " grow_from_best=" << grow_from_best << "\n";
                        break;
                    }
                }
                if (iter == active_max_newton_iter - 1) {
                    double growth = (res_iter1 > 0.0) ? (conv_res / std::max(res_iter1, 1e-30)) : 1.0;
                    const double rel_drop = 1.0 - growth;
                    const bool stagnation_ok = (growth <= params.stagnation_growth_tol) &&
                        (conv_res <= params.stagnation_abs_res_tol) &&
                        (conv_res < res_iter1) &&
                        (rel_drop >= params.stagnation_min_drop);
                    if (params.enable_stagnation_accept && stagnation_ok) {
                        converged = true; converge_mode = "stagnation";
                        std::cout << "    [NL-CONVERGED] step=" << step << " mode=" << converge_mode << " growth=" << growth << " rel_drop=" << rel_drop << " res_scaled=" << conv_res << "\n";
                    }
                    else {
                        fail_reason = "nonlinear_diverged";
                        std::cout << "    [NL-DIVERGED] step=" << step << " growth=" << growth << " rel_drop=" << rel_drop << " res_scaled=" << conv_res << "\n";
                    }
                    break;
                }

                if (params.enable_nonmonotone_line_search) {
                    line_search_hist.push_back(conv_res);
                    const size_t keep = static_cast<size_t>(std::max(1, params.nonmonotone_window));
                    if (line_search_hist.size() > keep) {
                        line_search_hist.erase(line_search_hist.begin(), line_search_hist.begin() + (line_search_hist.size() - keep));
                    }
                }

                const auto linear_t0 = std::chrono::steady_clock::now();
                const auto linear_result = detail::SolveLinearSystem<N>(
                    A, b, eq_contribs, totalEqWithWells, params, linear_solver_cache);
                const auto linear_t1 = std::chrono::steady_clock::now();

                Eigen::VectorXd dx = linear_result.dx;
                const bool compute_ok = linear_result.compute_ok;
                const bool solve_ok = linear_result.solve_ok;
                const bool row_scaling_applied = linear_result.row_scaling_applied;
                const std::string solver_log = linear_result.solver_log;
                linear_solve_ms_sum += std::chrono::duration<double, std::milli>(linear_t1 - linear_t0).count();
                linear_solve_calls++;

                if (params.diag_level != DiagLevel::Off) {
                    const int printEvery = std::max(1, params.diag_print_every_iter);
                    const bool printSummary = ((iter % printEvery) == 0) || !solve_ok;
                    const double near_ratio = (eos_total_samples > 0) ? (static_cast<double>(eos_near_bound_count) / static_cast<double>(eos_total_samples)) : 0.0;
                    const bool eos_warn = (eos_fallback_water > 0) || (eos_fallback_co2 > 0) || (near_ratio >= params.diag_eos_near_bound_ratio);

                    const double mass_flux_spike = (prev_max_mass_flux > 1e-30) ? (max_mass_flux / prev_max_mass_flux) : 1.0;
                    const double heat_flux_spike = (prev_max_heat_flux > 1e-30) ? (max_heat_flux / prev_max_heat_flux) : 1.0;
                    const bool flux_spike = (mass_flux_spike > params.diag_flux_spike_factor) || (heat_flux_spike > params.diag_flux_spike_factor);

                    if (printSummary) {
                        std::cout << "    [SOLVER] " << solver_log << "\n";
                        std::cout << "    [FLUX-GLOBAL] max_mass_flux=" << max_mass_flux
                            << " (" << ConnectionTypeLabel(hot_mass_type) << ", i=" << hot_mass_i << ", j=" << hot_mass_j << ")"
                            << " max_heat_flux=" << max_heat_flux
                            << " (" << ConnectionTypeLabel(hot_heat_type) << ", i=" << hot_heat_i << ", j=" << hot_heat_j << ")\n";
                        if (eos_warn) {
                            std::cout << "    [EOS-WARN] water_fallback=" << eos_fallback_water
                                << " co2_fallback=" << eos_fallback_co2
                                << " near_bound=" << eos_near_bound_count << "/" << eos_total_samples
                                << " ratio=" << near_ratio << "\n";
                        }
                    }

                    bool trigger_hotspot = !solve_ok;
                    const double growth_factor = (res_iter1 > 0.0) ? (conv_res / std::max(res_iter1, 1e-30)) : 1.0;
                    if (iter_used > 1 && growth_factor > params.diag_blowup_factor) {
                        trigger_hotspot = true;
                    }
                    if (flux_spike) {
                        trigger_hotspot = true;
                    }

                    if (prev_hot_idx == max_idx) {
                        const double rel_change = std::abs(max_res - prev_hot_res) / std::max(std::abs(prev_hot_res), 1e-30);
                        if (rel_change <= params.diag_hot_res_change_tol) {
                            ++hot_repeat_count;
                        }
                        else {
                            hot_repeat_count = 1;
                        }
                    }
                    else {
                        hot_repeat_count = 1;
                    }
                    if (hot_repeat_count >= params.diag_hot_repeat_iters) {
                        trigger_hotspot = true;
                    }

                    if (trigger_hotspot) {
                        const int hot_block = max_idx / N;
                        const int hot_var = max_idx % N;
                        const std::string domain_type = (hot_block < MatrixBlockCount(mgr)) ? "Matrix" : "Fracture";
                        std::cout << "    >> [HOT] trigger=1 eq=" << max_idx
                            << " (block=" << hot_block << ", domain=" << domain_type << ", var=" << hot_var << ")\n";

                        const bool hot_block_valid = (hot_block >= 0 && hot_block < totalBlocks);
                        if (!hot_block_valid) {
                            std::cout << "    >> [HOT] skipped: invalid block index " << hot_block << "\n";
                        }

                        if (hot_block_valid && max_idx >= 0 && max_idx < static_cast<int>(eq_contribs.size())) {
                            std::cout << "    >> [HOT-R&J] R_total=" << eq_contribs[max_idx].R_total()
                                << " = acc(" << eq_contribs[max_idx].R_acc << ") + flux(" << eq_contribs[max_idx].R_flux
                                << ") + well(" << eq_contribs[max_idx].R_well << ") + bc(" << eq_contribs[max_idx].R_bc << ")\n";
                            std::cout << "                 D_total=" << eq_contribs[max_idx].D_total()
                                << " = acc(" << eq_contribs[max_idx].D_acc << ") + flux(" << eq_contribs[max_idx].D_flux
                                << ") + well(" << eq_contribs[max_idx].D_well << ") + bc(" << eq_contribs[max_idx].D_bc << ")\n";
                        }
                        else {
                            std::cout << "    >> [HOT-R&J] skipped: invalid eq index " << max_idx << "\n";
                        }

                        if (hot_block_valid && (!incident_dumped_this_step || !params.diag_incident_once_per_step)) {
                            nlohmann::json snap;
                            snap["DiagnosticLevel"] = "Forensic_V3";
                            snap["TimeControl"]["step"] = step;
                            snap["TimeControl"]["iter"] = iter_used;
                            snap["TimeControl"]["dt"] = dt;

                            snap["LinearSolver"] = solver_log;
                            snap["Hotspot"]["block_id"] = hot_block;
                            snap["Hotspot"]["domain"] = domain_type;
                            snap["Hotspot"]["P"] = state.P[hot_block];
                            snap["Hotspot"]["T"] = state.T[hot_block];
                            if constexpr (N == 3) {
                                snap["Hotspot"]["Sw"] = state.Sw[hot_block];
                            }
                            if (solve_ok && max_idx >= 0 && max_idx < dx.size()) {
                                snap["Hotspot"]["dx_raw"] = dx[max_idx];
                            }

                            if constexpr (N == 3) {
                                try {
                                    const double sw_val = state.Sw[hot_block];
                                    ADVar<3> Sw_ad = ADVar<3>(sw_val);
                                    Sw_ad.grad[0] = 0.0;
                                    Sw_ad.grad[1] = 1.0;
                                    Sw_ad.grad[2] = 0.0;
                                    const auto& vg = vg_cfg;
                                    const auto& rp = rp_cfg;

                                    ADVar<3> Sw_const_ad = ClampSwForConstitutive<3>(Sw_ad, vg, sw_constitutive_eps);
                                    ADVar<3> krw, krg;
                                    CapRelPerm::kr_Mualem_vG<3>(Sw_const_ad, vg, rp, krw, krg);
                                    ADVar<3> pc = CapRelPerm::pc_vG<3>(Sw_const_ad, vg);

                                    nlohmann::json const_snap;
                                    const_snap["Sw"] = Sw_ad.val;
                                    const_snap["Sw_const"] = Sw_const_ad.val;
                                    const_snap["krw"] = krw.val;
                                    const_snap["krg"] = krg.val;
                                    const_snap["Pc"] = pc.val;
                                    const_snap["dkrw_dSw"] = krw.grad[1];
                                    const_snap["dkrg_dSw"] = krg.grad[1];
                                    const_snap["dPc_dSw"] = pc.grad[1];
                                    if (std::isnan(krw.val) || std::isnan(pc.grad[1]) || std::isinf(pc.grad[1])) {
                                        const_snap["flag"] = "invalid_constitutive";
                                    }
                                    snap["IncidentConst"] = const_snap;
                                    std::cout << "    >> [INCIDENT-CONST] Sw=" << const_snap["Sw"]
                                        << " krw=" << const_snap["krw"]
                                        << " krg=" << const_snap["krg"]
                                        << " Pc=" << const_snap["Pc"] << "\n";
                                }
                                catch (...) {
                                    snap["IncidentConst"]["flag"] = "invalid_constitutive";
                                }
                            }

                            const std::string snap_path = detail::BuildFailureSnapshotPath(
                                params.output_root_dir, caseName, step, iter_used);
                            std::ofstream ofs(snap_path);
                            if (ofs.is_open()) {
                                ofs << std::setw(4) << snap << std::endl;
                                incident_dumped_this_step = true;
                                std::cout << "    >> [INCIDENT] snapshot=" << snap_path << "\n";
                            }
                            else {
                                std::cout << "    >> [INCIDENT] failed to write snapshot file\n";
                            }
                        }
                    }

                    prev_hot_idx = max_idx;
                    prev_hot_res = max_res;
                    prev_max_mass_flux = max_mass_flux;
                    prev_max_heat_flux = max_heat_flux;
                }

                if (!solve_ok) { fail_reason = "linear_solve_fail"; break; }
                if (!dx.allFinite()) { fail_reason = "dx_nan_inf"; break; }

                bool state_valid = true;
                double alpha = 1.0;
                // Time-step-aware damping:
                // accumulation scales with 1/dt, so update caps should shrink at least linearly with dt.
                const double dt_ref = std::max(params.dt_init, params.dt_min);
                const double dt_eff = std::max(dt, params.dt_min);
                const double damp_scale = std::min(1.0, std::max(1.0e-12, dt_eff / dt_ref));
                // Keep tiny but non-zero floors for stiff late-time rollback loops.
                const double max_dP_eff = std::max(1.0e-3, params.max_dP * damp_scale);
                const double max_dT_eff = std::max(1.0e-5, params.max_dT * damp_scale);
                const double max_dSw_eff = std::max(1.0e-6, params.max_dSw * damp_scale);
                for (int bi = 0; bi < totalBlocks; ++bi) {
                    int eqP = mgr.getEquationIndex(bi, 0), eqT = mgr.getEquationIndex(bi, (N == 3) ? 2 : 1);
                    if (eqP < 0 || eqP >= dx.size() || eqT < 0 || eqT >= dx.size()) { state_valid = false; fail_reason = "invalid_eq_index"; break; }
                    alpha = std::min(alpha, max_dP_eff / (std::abs(dx[eqP]) + 1e-14));
                    alpha = std::min(alpha, max_dT_eff / (std::abs(dx[eqT]) + 1e-14));
                    if constexpr (N == 3) {
                        int eqSw = mgr.getEquationIndex(bi, 1);
                        if (eqSw < 0 || eqSw >= dx.size()) { state_valid = false; fail_reason = "invalid_eq_index"; break; }
                        alpha = std::min(alpha, max_dSw_eff / (std::abs(dx[eqSw]) + 1e-14));
                    }
                }
                if (!state_valid) break;

                if constexpr (N == 3) {
                    if (params.enable_alpha_safe_two_phase) {
                        const double sw_eps = std::max(1.0e-12, std::min(1.0e-2, params.sw_safe_eps));
                        const double sw_upper = 1.0 - sw_eps;
                        const double shrink = std::max(0.5, std::min(1.0, params.sw_alpha_shrink));
                        for (int bi = 0; bi < totalBlocks; ++bi) {
                            const int eqSw = mgr.getEquationIndex(bi, 1);
                            if (eqSw < 0 || eqSw >= dx.size()) { state_valid = false; fail_reason = "invalid_eq_index"; break; }
                            const double sw0 = state.Sw[bi];
                            const double dsw = dx[eqSw];
                            if (dsw > 0.0) {
                                const double amax = (sw_upper - sw0) / (dsw + 1.0e-30);
                                if (std::isfinite(amax)) alpha = std::min(alpha, shrink * amax);
                            }
                            else if (dsw < 0.0) {
                                const double amax = (sw0 - sw_eps) / (-dsw + 1.0e-30);
                                if (std::isfinite(amax)) alpha = std::min(alpha, shrink * amax);
                            }
                        }
                    }
                }

                if (!state_valid) break;
                alpha = std::min(1.0, alpha);
                if (!std::isfinite(alpha) || alpha < params.min_alpha) {
                    fail_reason = "alpha_too_small";
                    break;
                }

                const auto apply_trial_update = [&](const FIM_StateMap<N>& base_state,
                    double alpha_try,
                    FIM_StateMap<N>& out_state,
                    int& limiter_added_local,
                    double& rel_update_inf) -> bool {
                        return detail::ApplyTrialUpdate<N>(
                            base_state,
                            alpha_try,
                            dx,
                            mgr,
                            totalBlocks,
                            params,
                            p_floor,
                            p_ceil,
                            t_floor,
                            t_ceil,
                            out_state,
                            limiter_added_local,
                            rel_update_inf,
                            well_mgr.NumWells()); // Step 3: update well BHP DOFs
                    };
                const FIM_StateMap<N> state_before_update = state;
                double accepted_alpha = alpha;
                int accepted_limiter_added = 0;
                double accepted_rel_update = std::numeric_limits<double>::infinity();
                bool update_accepted = false;

                double conv_res_for_line_search = conv_res;
                double conv_res_probe = std::numeric_limits<double>::quiet_NaN();
                if (params.enable_armijo_line_search && params.enable_ls_base_check) {
                    conv_res_probe = compute_residual_inf_for_state(state_before_update);
                    if (std::isfinite(conv_res_probe) && conv_res_probe > 0.0) {
                        conv_res_for_line_search = conv_res_probe;
                        const double mismatch = std::abs(conv_res_probe - conv_res) / std::max(std::abs(conv_res), 1.0e-30);
                        if (mismatch > params.ls_base_mismatch_tol) {
                            std::cout << "    [LS-BASE-CHECK] step=" << step
                                << " iter=" << iter_used
                                << " conv_res=" << conv_res
                                << " probe_res=" << conv_res_probe
                                << " mismatch=" << mismatch
                                << " use_ref=" << conv_res_for_line_search << "\n";
                        }
                    }
                }

                if (params.enable_armijo_line_search) {
                    const int max_bt = std::max(1, params.armijo_max_backtracks);
                    const double bt_beta = std::max(0.1, std::min(0.95, params.armijo_beta));
                    const double c1 = std::max(1.0e-8, std::min(1.0e-2, params.armijo_c1));

                    int ls_trials = 0;
                    int ls_reject_state_invalid = 0;
                    int ls_reject_nan_inf = 0;
                    int ls_reject_armijo = 0;
                    int ls_reject_alpha_floor = 0;
                    std::string ls_last_reason = "none";
                    std::vector<std::string> ls_trace_lines;
                    ls_trace_lines.reserve(static_cast<size_t>(max_bt + 1));

                    auto line_search_ref_value = [&]() {
                        double ref = conv_res_for_line_search;
                        if (params.enable_nonmonotone_line_search && !line_search_hist.empty()) {
                            ref = std::max(ref, *std::max_element(line_search_hist.begin(), line_search_hist.end()));
                        }
                        return ref;
                        };

                    const bool ls_trace_live = params.enable_ls_trace && (step == 1);

                    double alpha_try = alpha;
                    double best_trial_res = std::numeric_limits<double>::infinity();
                    double best_trial_alpha = alpha;
                    int best_trial_limiter = 0;
                    double best_trial_rel_update = std::numeric_limits<double>::infinity();
                    FIM_StateMap<N> best_trial_state = state_before_update;
                    bool best_trial_valid = false;

                    for (int bt = 0; bt < max_bt; ++bt) {
                        ++ls_trials;
                        const double line_search_ref = line_search_ref_value();
                        const double armijo_rhs = std::max(0.0, line_search_ref - c1 * alpha_try * conv_res_for_line_search);

                        FIM_StateMap<N> trial_state = state_before_update;
                        int trial_limiter = 0;
                        double trial_rel_update = std::numeric_limits<double>::infinity();
                        double trial_res = std::numeric_limits<double>::quiet_NaN();

                        if (!apply_trial_update(state_before_update, alpha_try, trial_state, trial_limiter, trial_rel_update)) {
                            ++ls_reject_state_invalid;
                            ls_last_reason = "state_invalid";
                        }
                        else {
                            trial_res = compute_residual_inf_for_state(trial_state);
                            if (!std::isfinite(trial_res)) {
                                ++ls_reject_nan_inf;
                                ls_last_reason = "trial_res_nan_inf";
                            }
                            else {
                                if (trial_res < best_trial_res) {
                                    best_trial_res = trial_res;
                                    best_trial_alpha = alpha_try;
                                    best_trial_limiter = trial_limiter;
                                    best_trial_rel_update = trial_rel_update;
                                    best_trial_state = trial_state;
                                    best_trial_valid = true;
                                }

                                if (trial_res <= armijo_rhs) {
                                    state = std::move(trial_state);
                                    accepted_alpha = alpha_try;
                                    accepted_limiter_added = trial_limiter;
                                    accepted_rel_update = trial_rel_update;
                                    update_accepted = true;
                                    ls_last_reason = "accepted";
                                }
                                else {
                                    ++ls_reject_armijo;
                                    ls_last_reason = "armijo_not_satisfied";
                                }
                            }
                        }

                        if (params.enable_ls_trace) {
                            std::ostringstream oss;
                            oss << "    [LS-TRACE] step=" << step
                                << " iter=" << iter_used
                                << " bt=" << (bt + 1)
                                << " alpha_try=" << alpha_try
                                << " trial_res=" << trial_res
                                << " armijo_rhs=" << armijo_rhs
                                << " ref_res=" << line_search_ref
                                << " limiter_count=" << trial_limiter
                                << " reject_reason=" << ls_last_reason;
                            ls_trace_lines.push_back(oss.str());
                            if (ls_trace_live) {
                                std::cout << oss.str() << "\n";
                            }
                        }

                        if (update_accepted) {
                            break;
                        }

                        alpha_try *= bt_beta;
                        if (alpha_try < params.min_alpha) {
                            ++ls_reject_alpha_floor;
                            ls_last_reason = "alpha_too_small";
                            if (params.enable_ls_trace) {
                                std::ostringstream oss;
                                oss << "    [LS-TRACE] step=" << step
                                    << " iter=" << iter_used
                                    << " bt=" << (bt + 1)
                                    << " alpha_try=" << alpha_try
                                    << " trial_res=" << std::numeric_limits<double>::quiet_NaN()
                                    << " armijo_rhs=" << armijo_rhs
                                    << " ref_res=" << line_search_ref
                                    << " limiter_count=" << trial_limiter
                                    << " reject_reason=" << ls_last_reason;
                                ls_trace_lines.push_back(oss.str());
                                if (ls_trace_live) {
                                    std::cout << oss.str() << "\n";
                                }
                            }
                            break;
                        }
                    }

                    if (!update_accepted) {
                        const double line_search_ref = line_search_ref_value();
                        const double fallback_relax = 1.10; // allow mild non-monotone acceptance to escape hard rejects
                        if (best_trial_valid && std::isfinite(best_trial_res) && best_trial_res <= line_search_ref * fallback_relax) {
                            state = std::move(best_trial_state);
                            accepted_alpha = best_trial_alpha;
                            accepted_limiter_added = best_trial_limiter;
                            accepted_rel_update = best_trial_rel_update;
                            update_accepted = true;
                            if (params.diag_level != DiagLevel::Off) {
                                std::cout << "    [LS-FALLBACK] step=" << step
                                    << " iter=" << iter_used
                                    << " best_trial_res=" << best_trial_res
                                    << " ref=" << line_search_ref
                                    << " alpha=" << accepted_alpha << "\n";
                            }
                        }

                        if (!update_accepted) {
                            const bool allow_controlled_accept =
                                params.enable_controlled_accept_iter1 &&
                                (iter_used == 1) &&
                                (ls_reject_armijo >= std::max(1, params.controlled_accept_min_armijo_reject)) &&
                                (ls_reject_nan_inf == 0) &&
                                (ls_reject_state_invalid == 0) &&
                                best_trial_valid && std::isfinite(best_trial_res) &&
                                (controlled_accept_used_cache < std::max(0, params.controlled_accept_max_per_step));

                            if (allow_controlled_accept) {
                                const double relax = std::max(1.0, params.controlled_accept_relax);
                                const double rel_update_cap = params.controlled_accept_max_rel_update;
                                const bool res_ok = (best_trial_res <= line_search_ref * relax);
                                const bool upd_ok = (rel_update_cap <= 0.0) ||
                                    (std::isfinite(best_trial_rel_update) && best_trial_rel_update <= rel_update_cap);
                                if (res_ok && upd_ok) {
                                    state = std::move(best_trial_state);
                                    accepted_alpha = best_trial_alpha;
                                    accepted_limiter_added = best_trial_limiter;
                                    accepted_rel_update = best_trial_rel_update;
                                    update_accepted = true;
                                    controlled_accept_used_cache += 1;
                                    std::cout << "    [LS-CONTROLLED-ACCEPT] step=" << step
                                        << " iter=" << iter_used
                                        << " best_trial_res=" << best_trial_res
                                        << " ref=" << line_search_ref
                                        << " relax=" << relax
                                        << " alpha=" << accepted_alpha
                                        << " rel_update=" << accepted_rel_update
                                        << " used=" << controlled_accept_used_cache << "\n";
                                }
                            }
                        }

                        if (!update_accepted) {
                            if (params.enable_ls_trace && !ls_trace_live) {
                                for (const auto& ln : ls_trace_lines) {
                                    std::cout << ln << "\n";
                                }
                            }
                            std::cout << "    [LS-FAIL] step=" << step
                                << " iter=" << iter_used
                                << " trials=" << ls_trials
                                << " reject_state_invalid=" << ls_reject_state_invalid
                                << " reject_nan_inf=" << ls_reject_nan_inf
                                << " reject_armijo=" << ls_reject_armijo
                                << " reject_alpha_floor=" << ls_reject_alpha_floor
                                << " best_trial_res=" << best_trial_res
                                << " ref=" << line_search_ref
                                << " last_reason=" << ls_last_reason << "\n";

                            state = state_before_update;
                            if (ls_reject_nan_inf > 0) {
                                fail_reason = "line_search_fail_nan_inf";
                            }
                            else if (ls_reject_state_invalid > 0 && ls_reject_armijo == 0) {
                                fail_reason = "line_search_fail_state_invalid";
                            }
                            else if (ls_reject_alpha_floor > 0) {
                                fail_reason = "line_search_fail_alpha_floor";
                            }
                            else {
                                fail_reason = "line_search_fail_armijo";
                            }
                            break;
                        }
                    }
                }
                else {
                    FIM_StateMap<N> trial_state = state_before_update;
                    int trial_limiter = 0;
                    double trial_rel_update = std::numeric_limits<double>::infinity();
                    if (!apply_trial_update(state_before_update, alpha, trial_state, trial_limiter, trial_rel_update)) {
                        state = state_before_update;
                        fail_reason = "state_nan_inf";
                        break;
                    }
                    state = std::move(trial_state);
                    accepted_alpha = alpha;
                    accepted_limiter_added = trial_limiter;
                    accepted_rel_update = trial_rel_update;
                    update_accepted = true;
                }

                if (!update_accepted) { fail_reason = "state_nan_inf"; break; }
                total_limiters += accepted_limiter_added;
                last_rel_update = accepted_rel_update;

                if (params.diag_level != DiagLevel::Off && accepted_limiter_added >= params.diag_clamp_trigger) {
                    std::cout << "    [Limiter] added=" << accepted_limiter_added << " alpha=" << accepted_alpha << "\n";
                }
            }

            if (!converged && fail_reason.empty()) fail_reason = "nonlinear_max_iter";

            if (!converged) {
                const bool ls_fail = (fail_reason.rfind("line_search_fail", 0) == 0);
                if (ls_fail && iter_used <= 1) {
                    ls_fail_iter1_cache += 1;
                }
                else {
                    ls_fail_iter1_cache = 0;
                }

                bool rescue_applied = false;
                const int rescue_threshold = std::max(1, params.ls_fail_rescue_threshold);
                const int rescue_max = std::max(0, params.ls_fail_rescue_max);
                const double rescue_boost_factor = std::max(1.0, params.ptc_rescue_boost);
                if (active_enable_ptc && ls_fail && iter_used <= 1 &&
                    ls_fail_iter1_cache >= rescue_threshold &&
                    ptc_rescue_used_cache < rescue_max &&
                    rescue_boost_factor > 1.0) {
                    ptc_rescue_used_cache += 1;
                    ptc_rescue_boost_cache = std::max(1.0, ptc_rescue_boost_cache) * rescue_boost_factor;
                    state = old_state;
                    step--;
                    rescue_applied = true;
                    std::cout << "    [PTC-RESCUE] step=" << (step + 1)
                        << " iter=" << iter_used
                        << " ls_fail_iter1_count=" << ls_fail_iter1_cache
                        << " boost=" << ptc_rescue_boost_cache
                        << " reason=" << fail_reason
                        << " action=retry_without_dt_cut\n";
                }

                if (!rescue_applied) {
                    total_rollbacks++;
                    const double rollback_fac = std::min(1.0, std::max(0.05, params.rollback_shrink_factor));
                    dt = std::max(dt * rollback_fac, params.dt_min);
                    state = old_state;
                    step--;
                    std::cout << "    [Rollback] step=" << (step + 1) << " new_dt=" << dt << " reason=" << fail_reason << "\n";
                    auto export_crash_vtk = [&]() {
                        try {
                            const std::string cfname = detail::ExportVtkSnapshotFile(
                                mgr,
                                fm,
                                vtk_bc_ctx_ptr,
                                params.output_root_dir,
                                caseName,
                                "crash",
                                t,
                                N == 3,
                                [&]() {
                                    SyncStateToFieldManager(state, fm, mgr, sp_model, fluid_cfg, vg_cfg, rp_cfg);
                                });
                            std::cout << "    [VTK Export PASS] " << cfname << " (crash snapshot t=" << std::scientific << t << "s)\n";
                        }
                        catch (const std::exception& ex) {
                            std::cout << "    [VTK Export FAIL] crash snapshot reason=" << ex.what() << "\n";
                        }
                        catch (...) {
                            std::cout << "    [VTK Export FAIL] crash snapshot reason=unknown\n";
                        }
                        };
                    if (dt <= params.dt_min && total_rollbacks > 20) {
                        export_crash_vtk();
                        throw std::runtime_error("[FAIL] dt reached lower bound with repeated rollback.");
                    }
                    if (total_rollbacks > 80) {
                        export_crash_vtk();
                        throw std::runtime_error("[FAIL] Max rollbacks exceeded.");
                    }
                }
            }
            else {
                ls_fail_iter1_cache = 0;
                ptc_rescue_used_cache = 0;
                ptc_rescue_boost_cache = 1.0;
                controlled_accept_used_cache = 0;
                line_search_hist.clear();
                t += dt;
                completed_steps = step;
                final_residual = step_final_residual;
                std::cout << "  [Step Success] step=" << std::setw(3) << step << " dt=" << std::scientific << std::setprecision(3) << dt
                    << " iter_used=" << iter_used << " conv_mode=" << converge_mode
                    << " rollback_count=" << total_rollbacks << " limiter_count=" << total_limiters << "\n";

                const double dt_used = dt;
                detail::EmitStepAccepted<N>(
                    modules, step, t, dt_used, iter_used, step_final_residual, total_rollbacks, converge_mode, state);

                if (converge_mode == "abs_res") {
                    if (iter_used <= 3) dt = std::min(dt * 1.2, active_dt_max);
                    else if (iter_used <= 5) dt = std::min(dt * 1.05, active_dt_max);
                    else dt = std::max(dt * 0.8, params.dt_min);
                }
                else if (converge_mode == "rel_res_update") {
                    const int grow_hi = std::max(1, active_profile.dt_relres_iter_grow_hi);
                    const int neutral_hi = std::max(grow_hi, active_profile.dt_relres_iter_neutral_hi);
                    const int soft_hi = std::max(neutral_hi, active_profile.dt_relres_iter_soft_shrink_hi);
                    const double fac_grow = std::max(1.0, active_profile.dt_relres_grow_factor);
                    const double fac_neutral = std::max(0.1, active_profile.dt_relres_neutral_factor);
                    const double fac_soft = std::min(1.0, std::max(0.1, active_profile.dt_relres_soft_shrink_factor));
                    const double fac_hard = std::min(1.0, std::max(0.1, active_profile.dt_relres_hard_shrink_factor));

                    double fac = fac_hard;
                    if (iter_used <= grow_hi) fac = fac_grow;
                    else if (iter_used <= neutral_hi) fac = fac_neutral;
                    else if (iter_used <= soft_hi) fac = fac_soft;

                    dt = std::max(std::min(dt * fac, active_dt_max), params.dt_min);
                }
                else if (converge_mode == "stagnation" || converge_mode == "dt_floor_best" || converge_mode == "best_iter_guard" || converge_mode == "dt_floor_hold") {
                    dt = std::max(dt * 0.85, params.dt_min);
                }
                else {
                    const double rollback_fac = std::min(1.0, std::max(0.05, params.rollback_shrink_factor));
                    dt = std::max(dt * rollback_fac, params.dt_min);
                }

                if (params.enable_ramp_dt_protection && active_enable_control_ramp) {
                    const int ramp_steps = active_control_ramp_steps;
                    if (step <= ramp_steps) {
                        const double floor_scale = std::min(1.0, std::max(0.0, params.ramp_dt_min_scale));
                        const double dt_floor = std::max(dt_used * floor_scale, params.dt_min);
                        if (dt < dt_floor) dt = dt_floor;
                    }
                }

                const double dt_scale_factor = (dt_used > 0.0) ? (dt / dt_used) : 1.0;
                const double linear_solve_ms_avg = (linear_solve_calls > 0) ? (linear_solve_ms_sum / static_cast<double>(linear_solve_calls)) : 0.0;
                std::ostringstream speed_ss;
                speed_ss << std::scientific << std::setprecision(3)
                    << "    [SPEED] dt_next=" << dt
                    << " dt_scale_factor=" << dt_scale_factor
                    << " linear_solve_ms_avg_this_step=" << linear_solve_ms_avg << " ms";
                std::cout << speed_ss.str() << "\n";

                if (target_end_time_s > 0.0) {
                    const double t_remaining = target_end_time_s - t;
                    if (t_remaining > 0.0) {
                        dt = std::min(dt, t_remaining);
                    }
                }

                auto export_vtk_snapshot = [&](const std::string& tag) {
                    const std::string fname = detail::ExportVtkSnapshotFile(
                        mgr,
                        fm,
                        vtk_bc_ctx_ptr,
                        params.output_root_dir,
                        caseName,
                        tag,
                        t,
                        N == 3,
                        [&]() {
                            SyncStateToFieldManager(state, fm, mgr, sp_model, fluid_cfg, vg_cfg, rp_cfg);
                        });
                    if (tag == "final") final_vtk_exported = true;
                    vtk_export_count++;
                    detail::EmitSnapshot(modules, tag, step, t, fname);
                    std::cout << "    [VTK Export PASS] " << fname << "\n";
                    };

                const bool reached_target_end = (target_end_time_s > 0.0) && (t >= target_end_time_s - 1.0e-12);
                const double active_vtk_interval_s =
                    params.enable_two_stage_profile
                    ? (use_startup_stage ? params.startup_vtk_output_interval_s : params.long_vtk_output_interval_s)
                    : ((params.long_vtk_output_interval_s > 0.0) ? params.long_vtk_output_interval_s : params.startup_vtk_output_interval_s);

                // Ensure at least one recoverable snapshot exists even if the run fails
                // before the first time-based VTK interval is reached.
                if (!modules.disable_default_vtk_output) {
                    if (vtk_export_count == 0) {
                        export_vtk_snapshot("step_" + std::to_string(step));
                    }

                    if (active_vtk_interval_s > 0.0) {
                        if (!(next_vtk_output_time_s > 0.0)) {
                            next_vtk_output_time_s = t + active_vtk_interval_s;
                        }
                        if (t + 1.0e-12 >= next_vtk_output_time_s || reached_target_end) {
                            export_vtk_snapshot("step_" + std::to_string(step));
                            do {
                                next_vtk_output_time_s += active_vtk_interval_s;
                            } while (next_vtk_output_time_s <= t + 1.0e-12);
                        }
                    }
                    else if (step % 10 == 0 || step == params.max_steps || reached_target_end) {
                        export_vtk_snapshot((step == params.max_steps || reached_target_end)
                            ? "final"
                            : ("step_" + std::to_string(step)));
                    }
                }
            }
        }
        if (!final_vtk_exported && !modules.disable_default_vtk_output) {
            const std::string fname = detail::ExportVtkSnapshotFile(
                mgr,
                fm,
                vtk_bc_ctx_ptr,
                params.output_root_dir,
                caseName,
                "final",
                t,
                N == 3,
                [&]() {
                    SyncStateToFieldManager(state, fm, mgr, sp_model, fluid_cfg, vg_cfg, rp_cfg);
                });
            vtk_export_count++;
            detail::EmitSnapshot(modules, "final", completed_steps, t, fname);
            std::cout << "    [VTK Export PASS] " << fname << "\n";
        }

        std::cout << "[PASS] Day6 case completed: " << caseName
            << " | steps=" << completed_steps
            << " | rollbacks=" << total_rollbacks
            << " | limiters=" << total_limiters
            << " | final_residual=" << std::scientific << final_residual
            << " | vtk_exports=" << vtk_export_count
            << " | t_end=" << std::scientific << t << " s\n";
    }

} // namespace FIM_Engine





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
#include "../2D_PostProcess.h"
#include "../3D_PostProcess.h"

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
#include <memory>
#include <sstream>
#include <nlohmann/json.hpp>


namespace FIM_Engine {

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

        std::cout << "\n========== Starting Transient Scenario: " << caseName << " ==========\n";
        MakePath(caseName);
        InjectStaticProperties(fm);
        if (modules.property_initializer) {
            modules.property_initializer(mgr, fm);
            std::cout << "[Init] External property module injected.\n";
        }

        const SinglePhaseFluidModel sp_model = (N == 2) ? modules.single_phase_fluid : SinglePhaseFluidModel::Water;
        const bool sp_use_co2 = (N == 2) && (sp_model == SinglePhaseFluidModel::CO2);
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
        FIM_StateMap<N> state;
        state.InitSizes(totalBlocks);

        const auto pEqCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const auto tEqCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        const auto sEqCfg = PhysicalProperties_string_op::SaturationEquation_String::FIM();
        const int pressureDof = 0;
        const int saturationDof = (N == 3) ? 1 : -1;
        const int temperatureDof = (N == 3) ? 2 : 1;

        for (int i = 0; i < totalBlocks; ++i) {
            state.P[i] = ic.P_init;
            state.T[i] = ic.T_init;
            if constexpr (N == 3) state.Sw[i] = ic.Sw_init;
        }

        std::vector<double> vols(totalBlocks, 1.0);
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

        FIM_BlockSparseMatrix<N> global_mat(totalBlocks);
        for (const auto& conn : connMgr.GetConnections()) {
            global_mat.AddOffDiagBlock(conn.nodeI, conn.nodeJ, Eigen::Matrix<double, N, N>::Zero());
            global_mat.AddOffDiagBlock(conn.nodeJ, conn.nodeI, Eigen::Matrix<double, N, N>::Zero());
        }
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

        std::vector<double> P_ref(totalBlocks, ic.P_init);
        for (int i = 0; i < totalBlocks; ++i) P_ref[i] = state.P[i]; // 使用初始场作为压缩系数计算的参考压力

        double t = 0.0;
        double dt = params.dt_init;
        int total_rollbacks = 0;
        int total_limiters = 0;
        int completed_steps = 0;
        int vtk_export_count = 0;
        bool final_vtk_exported = false;
        double final_residual = std::numeric_limits<double>::quiet_NaN();
        // State guard rails:
        // 1) always keep conservative physical lower bounds;
        // 2) optionally enforce EOS-table hard clipping (off by default).
        double p_floor = 1.0e4;
        double p_ceil = std::numeric_limits<double>::max();
        double t_floor = 273.15;
        double t_ceil = std::numeric_limits<double>::max();
        if (params.clamp_state_to_eos_bounds) {
            const auto& wt_table = WaterPropertyTable::instance();
            p_floor = std::max(p_floor, wt_table.minPressure() + 1.0);
            p_ceil = std::min(p_ceil, wt_table.maxPressure() - 1.0);
            t_floor = std::max(t_floor, wt_table.minTemperature() + 1.0e-6);
            t_ceil = std::min(t_ceil, wt_table.maxTemperature() - 1.0e-6);
            if constexpr (N == 3) {
                const auto& gt_table = CO2PropertyTable::instance();
                p_floor = std::max(p_floor, gt_table.minPressure() + 1.0);
                p_ceil = std::min(p_ceil, gt_table.maxPressure() - 1.0);
                t_floor = std::max(t_floor, gt_table.minTemperature() + 1.0e-6);
                t_ceil = std::min(t_ceil, gt_table.maxTemperature() - 1.0e-6);
            }
            if (!(p_floor < p_ceil)) { p_floor = 1.0e4; p_ceil = std::numeric_limits<double>::max(); }
            if (!(t_floor < t_ceil)) { t_floor = 273.15; t_ceil = std::numeric_limits<double>::max(); }
        }

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

            std::vector<WellScheduleStep> active_wells = wells;
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
                    int nMat = 0;
                    if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                        nMat = mgr.getMatrixDOFCount();
                    }
                    else {
                        nMat = mgr.fracture_network().getSolverIndexOffset();
                    }

                    for (auto& w : active_wells) {
                        int bidx = -1;
                        if (w.domain == WellTargetDomain::Matrix) {
                            bidx = w.completion_id;
                        }
                        else if (w.domain == WellTargetDomain::Fracture) {
                            bidx = nMat + w.completion_id;
                        }

                        const double p_anchor = (bidx >= 0 && bidx < static_cast<int>(state.P.size())) ? state.P[bidx] : ic.P_init;
                        const double t_anchor = (bidx >= 0 && bidx < static_cast<int>(state.T.size())) ? state.T[bidx] : ic.T_init;

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
                FIM_BlockSparseMatrix<N> probe_mat(totalBlocks);
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
                    auto pW = EvalPrimaryFluid<N>(sp_model, P, T);
                    ADVar<N> P_old(old_state.P[bi]), T_old(old_state.T[bi]);
                    auto pW_old = EvalPrimaryFluid<N>(sp_model, P_old, T_old);

                    // 引入岩石压缩性，由于 P 含有梯度 1.0，phi 自动继承对压力的导数 c_r * phi_ref
                    ADVar<N> phi = ADVar<N>(phi_ref) * (ADVar<N>(1.0) + ADVar<N>(c_r) * (P - ADVar<N>(P_ref[bi])));
                    ADVar<N> phi_old = ADVar<N>(phi_ref) * (ADVar<N>(1.0) + ADVar<N>(c_r) * (P_old - ADVar<N>(P_ref[bi])));

                    std::vector<ADVar<N>> acc_eqs(N);
                    if constexpr (N == 2) {
                        ADVar<N> m_w = pW.rho * phi, m_w_old = pW_old.rho * phi_old;
                        acc_eqs[0] = (m_w - m_w_old) * (vols[bi] / dt);
                        ADVar<N> e_w = m_w * (pW.h - P / pW.rho) + ADVar<N>((1.0 - phi) * rho_r * c_pr) * T;
                        ADVar<N> e_w_old = m_w_old * (pW_old.h - P_old / pW_old.rho) + ADVar<N>((1.0 - phi_old) * rho_r * c_pr) * T_old;
                        acc_eqs[1] = (e_w - e_w_old) * (vols[bi] / dt);
                    }
                    else {
                        ADVar<N> Sw(eval_state.Sw[bi]); Sw.grad(1) = 1.0;
                        ADVar<N> Sg = ADVar<N>(1.0) - Sw;
                        ADVar<N> Sw_old(old_state.Sw[bi]), Sg_old = ADVar<N>(1.0) - Sw_old;
                        auto pG = AD_Fluid::Evaluator::evaluateCO2<N>(P, T);
                        auto pG_old = AD_Fluid::Evaluator::evaluateCO2<N>(P_old, T_old);
                        acc_eqs[0] = (pW.rho * phi * Sw - pW_old.rho * phi_old * Sw_old) * (vols[bi] / dt);
                        acc_eqs[1] = (pG.rho * phi * Sg - pG_old.rho * phi_old * Sg_old) * (vols[bi] / dt);
                        ADVar<N> e_fluid = pW.rho * Sw * (pW.h - P / pW.rho) + pG.rho * Sg * (pG.h - P / pG.rho);
                        ADVar<N> e_fluid_old = pW_old.rho * Sw_old * (pW_old.h - P_old / pW_old.rho) + pG_old.rho * Sg_old * (pG_old.h - P_old / pG_old.rho);
                        ADVar<N> e_rock = ADVar<N>((1.0 - phi) * rho_r * c_pr) * T;
                        ADVar<N> e_rock_old = ADVar<N>((1.0 - phi_old) * rho_r * c_pr) * T_old;
                        acc_eqs[2] = ((e_fluid * phi + e_rock) - (e_fluid_old * phi_old + e_rock_old)) * (vols[bi] / dt);
                    }
                    FIM_GlobalAssembler<N, ADVar<N>>::AssembleAccumulation(bi, acc_eqs, probe_mat);
                }

                // 2) flux
                for (const auto& conn : connMgr.GetConnections()) {
                    int i = conn.nodeI, j = conn.nodeJ;
                    auto evalFlux = [&](bool wrt_i) -> std::vector<ADVar<N>> {
                        std::vector<ADVar<N>> F(N);
                        ADVar<N> P_i(eval_state.P[i]), T_i(eval_state.T[i]), P_j(eval_state.P[j]), T_j(eval_state.T[j]);
                        const Vector& x_i = blockCenters[i]; const Vector& x_j = blockCenters[j];

                        if constexpr (N == 2) {
                            if (wrt_i) { P_i.grad(0) = 1.0; T_i.grad(1) = 1.0; }
                            else { P_j.grad(0) = 1.0; T_j.grad(1) = 1.0; }
                            auto pW_i = EvalPrimaryFluid<N>(sp_model, P_i, T_i);
                            auto pW_j = EvalPrimaryFluid<N>(sp_model, P_j, T_j);
                            ADVar<N> rho_avg_w = ADVar<N>(0.5) * (pW_i.rho + pW_j.rho);
                            ADVar<N> dPhi = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, rho_avg_w, x_i, x_j, gravityVec);
                            ADVar<N> mob_i = pW_i.rho / pW_i.mu, mob_j = pW_j.rho / pW_j.mu;
                            ADVar<N> up_mob = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi, mob_i, mob_j);
                            F[0] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mob, dPhi);
                            ADVar<N> up_h = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi, pW_i.h, pW_j.h);
                            F[1] = FVM_Ops::Compute_Heat_Flux<N, ADVar<N>>(conn.T_Heat, T_i, T_j, F[0], up_h);
                        }
                        else {
                            ADVar<N> Sw_i(eval_state.Sw[i]), Sw_j(eval_state.Sw[j]);
                            if (wrt_i) { P_i.grad(0) = 1.0; Sw_i.grad(1) = 1.0; T_i.grad(2) = 1.0; }
                            else { P_j.grad(0) = 1.0; Sw_j.grad(1) = 1.0; T_j.grad(2) = 1.0; }
                            auto pW_i = AD_Fluid::Evaluator::evaluateWater<N>(P_i, T_i), pW_j = AD_Fluid::Evaluator::evaluateWater<N>(P_j, T_j);
                            auto pG_i = AD_Fluid::Evaluator::evaluateCO2<N>(P_i, T_i), pG_j = AD_Fluid::Evaluator::evaluateCO2<N>(P_j, T_j);
                            const auto& vg = vg_cfg;
                            const auto& rp = rp_cfg;
                            ADVar<N> krw_i, krg_i, krw_j, krg_j;
                            CapRelPerm::kr_Mualem_vG<N>(Sw_i, vg, rp, krw_i, krg_i);
                            CapRelPerm::kr_Mualem_vG<N>(Sw_j, vg, rp, krw_j, krg_j);
                            ADVar<N> rho_avg_w = ADVar<N>(0.5) * (pW_i.rho + pW_j.rho), rho_avg_g = ADVar<N>(0.5) * (pG_i.rho + pG_j.rho);
                            ADVar<N> Pc_i = CapRelPerm::pc_vG<N>(Sw_i, vg), Pc_j = CapRelPerm::pc_vG<N>(Sw_j, vg);
                            ADVar<N> dPhi_w = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, rho_avg_w, x_i, x_j, gravityVec);
                            ADVar<N> dPhi_g = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, Pc_i, Pc_j, rho_avg_g, x_i, x_j, gravityVec);
                            ADVar<N> mobW_i = krw_i * pW_i.rho / pW_i.mu, mobW_j = krw_j * pW_j.rho / pW_j.mu;
                            ADVar<N> mobG_i = krg_i * pG_i.rho / pG_i.mu, mobG_j = krg_j * pG_j.rho / pG_j.mu;
                            ADVar<N> up_mobW = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_w, mobW_i, mobW_j);
                            ADVar<N> up_mobG = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_g, mobG_i, mobG_j);
                            F[0] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mobW, dPhi_w);
                            F[1] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mobG, dPhi_g);
                            ADVar<N> up_h_w = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_w, pW_i.h, pW_j.h);
                            ADVar<N> up_h_g = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_g, pG_i.h, pG_j.h);
                            F[2] = FVM_Ops::Compute_Heat_Flux<N, ADVar<N>>(conn.T_Heat, T_i, T_j, F[0], F[1], up_h_w, up_h_g);
                        }
                        return F;
                        };
                    auto f_wrt_i = evalFlux(true);
                    auto f_wrt_j = evalFlux(false);
                    FIM_GlobalAssembler<N, ADVar<N>>::AssembleFlux(i, j, f_wrt_i, f_wrt_j, probe_mat);
                }

                // 3) well source
                SyncStateToFieldManager(eval_state, fm, mgr, sp_model, vg_cfg, rp_cfg);
                std::vector<double> w_res(totalEq, 0.0);
                std::vector<std::array<double, 3>> w_jac3(totalEq, std::array<double, 3>{ 0.0, 0.0, 0.0 });
                const int well_dof_w = (N == 3) ? 0 : (sp_use_co2 ? -1 : 0);
                const int well_dof_g = (N == 3) ? 1 : (sp_use_co2 ? 0 : -1);
                const int well_dof_e = (N == 3) ? 2 : 1;
                if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                    BoundaryAssembler::Assemble_Wells_2D_FullJac(mgr, fm, active_wells, 0, well_dof_w, well_dof_g, well_dof_e, w_res, w_jac3, sp_use_co2, vg_cfg, rp_cfg);
                }
                else {
                    BoundaryAssembler::Assemble_Wells_3D_FullJac(mgr, fm, active_wells, 0, well_dof_w, well_dof_g, well_dof_e, w_res, w_jac3, sp_use_co2, vg_cfg, rp_cfg);
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

                // 4) boundary source
                auto assembleBoundaryFieldProbe = [&](const BoundarySetting::BoundaryConditionManager* bcMgr, int dofOffset, const std::string& fieldName) {
                    if (!bcMgr || dofOffset < 0) return;
                    std::vector<double> bc_res(totalEq, 0.0);
                    std::vector<double> bc_diag(totalEq, 0.0);
                    BoundaryAssemblyStats bc_stats;
                    if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                        bc_stats = BoundaryAssembler::Assemble_2D(mgr, *bcMgr, dofOffset, fm, fieldName, bc_res, bc_diag);
                    }
                    else {
                        bc_stats = BoundaryAssembler::Assemble_3D(mgr, *bcMgr, dofOffset, fm, fieldName, bc_res, bc_diag);
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
                if constexpr (N == 3) {
                    assembleBoundaryFieldProbe(modules.saturation_bc, saturationDof, sEqCfg.saturation);
                }
                assembleBoundaryFieldProbe(modules.temperature_bc, temperatureDof, tEqCfg.temperatue_field);

                const Eigen::VectorXd b_probe = probe_mat.ExportEigenResidual();
                int max_probe_idx = 0;
                const double max_probe = b_probe.cwiseAbs().maxCoeff(&max_probe_idx);
                (void)max_probe_idx;
                return max_probe; // [包02] 线搜索探测必须基于原始残差 (raw_res) 以保证物理意义
                };

            for (int iter = 0; iter < active_max_newton_iter; ++iter) {
                iter_used++;
                global_mat.SetZero();

                std::vector<EqContrib> eq_contribs(totalEq);

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
                    auto pW = EvalPrimaryFluid<N>(sp_model, P, T);
                    ADVar<N> P_old(old_state.P[bi]), T_old(old_state.T[bi]);
                    auto pW_old = EvalPrimaryFluid<N>(sp_model, P_old, T_old);
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
                        ADVar<N> e_w = m_w * (pW.h - P / pW.rho) + ADVar<N>((1.0 - phi) * rho_r * c_pr) * T;
                        ADVar<N> e_w_old = m_w_old * (pW_old.h - P_old / pW_old.rho) + ADVar<N>((1.0 - phi_old) * rho_r * c_pr) * T_old;
                        acc_eqs[1] = (e_w - e_w_old) * (vols[bi] / dt);
                    }
                    else {
                        ADVar<N> Sw(state.Sw[bi]); Sw.grad(1) = 1.0;
                        ADVar<N> Sg = ADVar<N>(1.0) - Sw;
                        ADVar<N> Sw_old(old_state.Sw[bi]), Sg_old = ADVar<N>(1.0) - Sw_old;

                        auto pG = AD_Fluid::Evaluator::evaluateCO2<N>(P, T);
                        auto pG_old = AD_Fluid::Evaluator::evaluateCO2<N>(P_old, T_old);
                        if (track_eos_domain) {
                            ++eos_total_samples;
                            if (pG.isFallback) ++eos_fallback_co2;
                            if (pG.near_bound) ++eos_near_bound_count;
                        }

                        acc_eqs[0] = (pW.rho * phi * Sw - pW_old.rho * phi_old * Sw_old) * (vols[bi] / dt);
                        acc_eqs[1] = (pG.rho * phi * Sg - pG_old.rho * phi_old * Sg_old) * (vols[bi] / dt);

                        ADVar<N> e_fluid = pW.rho * Sw * (pW.h - P / pW.rho) + pG.rho * Sg * (pG.h - P / pG.rho);
                        ADVar<N> e_fluid_old = pW_old.rho * Sw_old * (pW_old.h - P_old / pW_old.rho) + pG_old.rho * Sg_old * (pG_old.h - P_old / pG_old.rho);
                        ADVar<N> e_rock = ADVar<N>((1.0 - phi) * rho_r * c_pr) * T;
                        ADVar<N> e_rock_old = ADVar<N>((1.0 - phi_old) * rho_r * c_pr) * T_old;
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


                for (const auto& conn : connMgr.GetConnections()) {
                    audit_stats[conn.type].conn_count++;
                    audit_stats[conn.type].sum_abs_tflow += std::abs(conn.T_Flow);
                    audit_stats[conn.type].sum_abs_theat += std::abs(conn.T_Heat);

                    int i = conn.nodeI, j = conn.nodeJ;
                    auto evalFlux = [&](bool wrt_i) -> std::vector<ADVar<N>> {
                        std::vector<ADVar<N>> F(N);
                        ADVar<N> P_i(state.P[i]), T_i(state.T[i]), P_j(state.P[j]), T_j(state.T[j]);
                        const Vector& x_i = blockCenters[i]; const Vector& x_j = blockCenters[j];

                        if constexpr (N == 2) {
                            if (wrt_i) { P_i.grad(0) = 1.0; T_i.grad(1) = 1.0; }
                            else { P_j.grad(0) = 1.0; T_j.grad(1) = 1.0; }
                            auto pW_i = EvalPrimaryFluid<N>(sp_model, P_i, T_i);
                            auto pW_j = EvalPrimaryFluid<N>(sp_model, P_j, T_j);
                            ADVar<N> rho_avg_w = ADVar<N>(0.5) * (pW_i.rho + pW_j.rho);
                            ADVar<N> dPhi = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, rho_avg_w, x_i, x_j, gravityVec);
                            ADVar<N> mob_i = pW_i.rho / pW_i.mu, mob_j = pW_j.rho / pW_j.mu;
                            ADVar<N> up_mob = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi, mob_i, mob_j);
                            F[0] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mob, dPhi);
                            ADVar<N> up_h = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi, pW_i.h, pW_j.h);
                            F[1] = FVM_Ops::Compute_Heat_Flux<N, ADVar<N>>(conn.T_Heat, T_i, T_j, F[0], up_h);
                        }
                        else {
                            ADVar<N> Sw_i(state.Sw[i]), Sw_j(state.Sw[j]);
                            if (wrt_i) { P_i.grad(0) = 1.0; Sw_i.grad(1) = 1.0; T_i.grad(2) = 1.0; }
                            else { P_j.grad(0) = 1.0; Sw_j.grad(1) = 1.0; T_j.grad(2) = 1.0; }
                            auto pW_i = AD_Fluid::Evaluator::evaluateWater<N>(P_i, T_i), pW_j = AD_Fluid::Evaluator::evaluateWater<N>(P_j, T_j);
                            auto pG_i = AD_Fluid::Evaluator::evaluateCO2<N>(P_i, T_i), pG_j = AD_Fluid::Evaluator::evaluateCO2<N>(P_j, T_j);
                            const auto& vg = vg_cfg;
                            const auto& rp = rp_cfg;
                            ADVar<N> krw_i, krg_i, krw_j, krg_j;
                            CapRelPerm::kr_Mualem_vG<N>(Sw_i, vg, rp, krw_i, krg_i);
                            CapRelPerm::kr_Mualem_vG<N>(Sw_j, vg, rp, krw_j, krg_j);

                            ADVar<N> rho_avg_w = ADVar<N>(0.5) * (pW_i.rho + pW_j.rho), rho_avg_g = ADVar<N>(0.5) * (pG_i.rho + pG_j.rho);
                            ADVar<N> Pc_i = CapRelPerm::pc_vG<N>(Sw_i, vg), Pc_j = CapRelPerm::pc_vG<N>(Sw_j, vg);
                            ADVar<N> dPhi_w = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, rho_avg_w, x_i, x_j, gravityVec);
                            ADVar<N> dPhi_g = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, Pc_i, Pc_j, rho_avg_g, x_i, x_j, gravityVec);

                            ADVar<N> mobW_i = krw_i * pW_i.rho / pW_i.mu, mobW_j = krw_j * pW_j.rho / pW_j.mu;
                            ADVar<N> mobG_i = krg_i * pG_i.rho / pG_i.mu, mobG_j = krg_j * pG_j.rho / pG_j.mu;
                            ADVar<N> up_mobW = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_w, mobW_i, mobW_j);
                            ADVar<N> up_mobG = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_g, mobG_i, mobG_j);

                            F[0] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mobW, dPhi_w);
                            F[1] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mobG, dPhi_g);
                            ADVar<N> up_h_w = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_w, pW_i.h, pW_j.h);
                            ADVar<N> up_h_g = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_g, pG_i.h, pG_j.h);
                            F[2] = FVM_Ops::Compute_Heat_Flux<N, ADVar<N>>(conn.T_Heat, T_i, T_j, F[0], F[1], up_h_w, up_h_g);
                        }
                        return F;
                        };
                    auto f_wrt_i = evalFlux(true);
                    auto f_wrt_j = evalFlux(false);
                    FIM_GlobalAssembler<N, ADVar<N>>::AssembleFlux(i, j, f_wrt_i, f_wrt_j, global_mat);

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
                                eq_contribs[g_eq_i].R_flux += f_wrt_i[eq].val;
                                eq_contribs[g_eq_i].D_flux += f_wrt_i[eq].grad[eq];
                            }
                            if (g_eq_j >= 0 && g_eq_j < eq_contribs.size()) {
                                eq_contribs[g_eq_j].R_flux -= f_wrt_j[eq].val;      // ??? j ???? -F
                                eq_contribs[g_eq_j].D_flux -= f_wrt_j[eq].grad[eq]; // ??????????????
                            }
                        }
                    }
                }


                SyncStateToFieldManager(state, fm, mgr, sp_model, vg_cfg, rp_cfg);
                std::vector<double> w_res(totalEq, 0.0);
                std::vector<std::array<double, 3>> w_jac3(totalEq, std::array<double, 3>{ 0.0, 0.0, 0.0 });
                const int well_dof_w = (N == 3) ? 0 : (sp_use_co2 ? -1 : 0);
                const int well_dof_g = (N == 3) ? 1 : (sp_use_co2 ? 0 : -1);
                const int well_dof_e = (N == 3) ? 2 : 1;

                if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                    BoundaryAssembler::Assemble_Wells_2D_FullJac(mgr, fm, active_wells, 0, well_dof_w, well_dof_g, well_dof_e, w_res, w_jac3, sp_use_co2, vg_cfg, rp_cfg);
                }
                else {
                    BoundaryAssembler::Assemble_Wells_3D_FullJac(mgr, fm, active_wells, 0, well_dof_w, well_dof_g, well_dof_e, w_res, w_jac3, sp_use_co2, vg_cfg, rp_cfg);
                }

                double max_abs_well_dsw = 0.0;
                double max_abs_well_dt = 0.0;

                // BoundaryAssembler returns well terms in outflow-positive convention.
                // Residual here is assembled as Acc + Flux + Q_out = 0 by default.
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
                auto assembleBoundaryField = [&](const BoundarySetting::BoundaryConditionManager* bcMgr,
                    int dofOffset,
                    const std::string& fieldName,
                    const char* fieldLabel) {
                        if (!bcMgr || dofOffset < 0) return;

                        std::vector<double> bc_res(totalEq, 0.0);
                        std::vector<double> bc_diag(totalEq, 0.0);
                        BoundaryAssemblyStats bc_stats;

                        if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                            bc_stats = BoundaryAssembler::Assemble_2D(mgr, *bcMgr, dofOffset, fm, fieldName, bc_res, bc_diag);
                        }
                        else {
                            bc_stats = BoundaryAssembler::Assemble_3D(mgr, *bcMgr, dofOffset, fm, fieldName, bc_res, bc_diag);
                        }

                        int appliedEq = 0;
                        for (int bi = 0; bi < totalBlocks; ++bi) {
                            const int eqIdx = mgr.getEquationIndex(bi, dofOffset);
                            if (eqIdx < 0 || eqIdx >= totalEq) continue;

                            const double r_bc = bc_res[eqIdx];
                            const double d_bc = bc_diag[eqIdx];
                            if (std::abs(r_bc) <= 1e-16 && std::abs(d_bc) <= 1e-16) continue;

                            global_mat.AddResidual(bi, dofOffset, r_bc);
                            global_mat.AddDiagJacobian(bi, dofOffset, dofOffset, d_bc);
                            ++appliedEq;

                            if (params.diag_level != DiagLevel::Off && eqIdx >= 0 && eqIdx < eq_contribs.size()) {
                                eq_contribs[eqIdx].R_bc += r_bc;
                                eq_contribs[eqIdx].D_bc += d_bc;
                            }
                        }
                        if (params.diag_level != DiagLevel::Off) {
                            std::cout << "    [BC-SUM] field=" << fieldLabel
                                << " faces=" << bc_stats.matrixBCCount + bc_stats.fractureBCCount
                                << " applied_eq=" << appliedEq
                                << " (visited=" << bc_stats.visitedEqRows
                                << ", nonzero=" << bc_stats.nonzeroEqRows
                                << ", zero_row=" << bc_stats.zeroEqRows
                                << ", invalid=" << bc_stats.invalidEqRows << ")\n";
                            for (const auto& kv : bc_stats.perTagType) {
                                const auto& v = kv.second;
                                std::cout << "      [BC-TAG] field=" << fieldLabel
                                    << " key=" << kv.first
                                    << " faces=" << v.faces
                                    << " applied=" << v.applied
                                    << " skipped=" << v.skipped
                                    << " |sumR|=" << v.sumR
                                    << " |sumDiag|=" << v.sumDiag << "\n";
                            }
                        }
                    };

                assembleBoundaryField(modules.pressure_bc, pressureDof, pEqCfg.pressure_field, "P");
                if constexpr (N == 3) {
                    assembleBoundaryField(modules.saturation_bc, saturationDof, sEqCfg.saturation, "Sw");
                }
                assembleBoundaryField(modules.temperature_bc, temperatureDof, tEqCfg.temperatue_field, "T");
                auto A = global_mat.ExportEigenSparseMatrix();
                auto b = global_mat.ExportEigenResidual();
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
                        for (int r = 0; r < totalEq; ++r) {
                            const double diag_acc_abs = (r >= 0 && r < static_cast<int>(eq_contribs.size())) ? std::abs(eq_contribs[r].D_acc) : 0.0;
                            const double m_ptc = std::max(diag_acc_abs, params.row_scale_floor);
                            A.coeffRef(r, r) += ptc_lambda_iter * m_ptc;
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
                const double conv_res = max_res; // [包02] 牛顿收敛判据仅依赖 raw_res，消除因行缩放过大导致的“假收敛”

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
                    A, b, eq_contribs, totalEq, params, linear_solver_cache);
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

                                    ADVar<3> krw, krg;
                                    CapRelPerm::kr_Mualem_vG<3>(Sw_ad, vg, rp, krw, krg);
                                    ADVar<3> pc = CapRelPerm::pc_vG<3>(Sw_ad, vg);

                                    nlohmann::json const_snap;
                                    const_snap["Sw"] = Sw_ad.val;
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

                            const std::string snap_path = "Test/Transient/Day6/" + caseName + "/fail_snapshot_step" +
                                std::to_string(step) + "_iter" + std::to_string(iter_used) + ".json";
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
                const double damp_scale = std::max(1.0e-12, dt_eff / dt_ref);
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
                            rel_update_inf);
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
                    if (dt <= params.dt_min && total_rollbacks > 20) throw std::runtime_error("[FAIL] dt reached lower bound with repeated rollback.");
                    if (total_rollbacks > 80) throw std::runtime_error("[FAIL] Max rollbacks exceeded.");
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

                auto export_vtk_snapshot = [&](const std::string& suffix) {
                    SyncStateToFieldManager(state, fm, mgr, sp_model, vg_cfg, rp_cfg);
                    const std::string fname = "Test/Transient/Day6/" + caseName + suffix;
                    if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                        PostProcess_2D(mgr, fm).ExportVTK(fname, t);
                    }
                    else {
                        PostProcess_3D(mgr, fm).ExportVTK(fname, t);
                    }
                    VerifyVtkExport(fname, N == 3);
                    if (suffix == "/final.vtk") final_vtk_exported = true;
                    vtk_export_count++;
                    std::cout << "    [VTK Export PASS] " << fname << "\n";
                    };

                const bool reached_target_end = (target_end_time_s > 0.0) && (t >= target_end_time_s - 1.0e-12);
                const double active_vtk_interval_s =
                    params.enable_two_stage_profile
                    ? (use_startup_stage ? params.startup_vtk_output_interval_s : params.long_vtk_output_interval_s)
                    : ((params.long_vtk_output_interval_s > 0.0) ? params.long_vtk_output_interval_s : params.startup_vtk_output_interval_s);

                if (active_vtk_interval_s > 0.0) {
                    if (!(next_vtk_output_time_s > 0.0)) {
                        next_vtk_output_time_s = t + active_vtk_interval_s;
                    }
                    if (t + 1.0e-12 >= next_vtk_output_time_s || reached_target_end) {
                        export_vtk_snapshot("/step_" + std::to_string(step) + ".vtk");
                        do {
                            next_vtk_output_time_s += active_vtk_interval_s;
                        } while (next_vtk_output_time_s <= t + 1.0e-12);
                    }
                }
                else if (step % 10 == 0 || step == params.max_steps || reached_target_end) {
                    export_vtk_snapshot((step == params.max_steps || reached_target_end)
                        ? "/final.vtk"
                        : ("/step_" + std::to_string(step) + ".vtk"));
                }
            }
        }
        if (!final_vtk_exported) {
            SyncStateToFieldManager(state, fm, mgr, sp_model, vg_cfg, rp_cfg);
            const std::string fname = "Test/Transient/Day6/" + caseName + "/final.vtk";
            if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                PostProcess_2D(mgr, fm).ExportVTK(fname, t);
            }
            else {
                PostProcess_3D(mgr, fm).ExportVTK(fname, t);
            }
            VerifyVtkExport(fname, N == 3);
            vtk_export_count++;
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

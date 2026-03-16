#pragma once

#include "../MeshManager.h"
#include "../3D_MeshManager.h"
#include "../BoundaryAssembler.h"
#include "../CapRelPerm_HD_AD.h"
#include "../FIM_TransientSupport.hpp"
#include "../UserDefineVarType.h"

#include <functional>
#include <type_traits>
#include <vector>

namespace FIM_Engine {

    enum class DiagLevel { Off, Summary, Hotspot, Forensic };
    enum class SolverRoute { FIM, IMPES };
    enum class LinearSolverType { SparseLU, BiCGSTAB, AMGCL };
    enum class SinglePhaseFluidModel { Water, CO2 };

    struct InitialConditions {
        double P_init = 2.0e5;
        double T_init = 300.0;
        double Sw_init = 0.2;
    };

    struct TransientStageProfile {
        double dt_max = 86400.0;
        int max_newton_iter = 8;
        double rel_res_tol = 1.0e-3;
        double rel_update_tol = 1.0e-6;

        bool enable_ptc = false;
        double ptc_lambda_init = 1.0;
        double ptc_lambda_decay = 0.5;
        double ptc_lambda_min = 0.0;

        bool enable_control_ramp = false;
        int control_ramp_steps = 5;
        double control_ramp_min = 0.2;
        bool control_ramp_apply_rate = true;
        bool control_ramp_apply_bhp = true;

        int dt_relres_iter_grow_hi = 8;
        int dt_relres_iter_neutral_hi = 14;
        int dt_relres_iter_soft_shrink_hi = 20;
        double dt_relres_grow_factor = 1.08;
        double dt_relres_neutral_factor = 1.00;
        double dt_relres_soft_shrink_factor = 0.98;
        double dt_relres_hard_shrink_factor = 0.92;
    };

    struct TransientSolverParams {
        int max_steps = 50;
        double dt_init = 1.0;
        double dt_min = 1e-4;
        double dt_max = 86400.0;
        double target_end_time_s = -1.0;

        bool enable_two_stage_profile = false;
        double startup_end_time_s = 0.0;
        TransientStageProfile startup_profile{};
        TransientStageProfile long_profile{};

        double startup_vtk_output_interval_s = -1.0;
        double long_vtk_output_interval_s = -1.0;
        double matrix_audit_coeff_tol = 1.0e-16;
        int matrix_audit_max_detail = 8;
        bool enable_matrix_audit = false;
        bool matrix_audit_strict = false;
        int matrix_audit_step = 1;
        int matrix_audit_iter = 1;
        double matrix_audit_eps = 1.0e-14;
        bool matrix_audit_require_nnc = true;
        bool matrix_audit_require_ff = true;

        int max_newton_iter = 8;
        double abs_res_tol = 1e-6;

        double stagnation_growth_tol = 0.995;
        double stagnation_abs_res_tol = 1.0e2;
        double stagnation_min_drop = 2.0e-3;
        double best_iter_growth_trigger = 1.5;
        int best_iter_guard_min_iter = 3;
        bool enable_best_iter_guard = false;
        bool enable_dt_floor_best = false;
        bool enable_stagnation_accept = false;
        bool enable_dt_floor_hold = false;

        double rel_res_tol = 1.0e-3;
        double rel_update_tol = 1.0e-6;

        bool enable_armijo_line_search = true;
        int armijo_max_backtracks = 8;
        double armijo_beta = 0.5;
        double armijo_c1 = 1.0e-4;
        bool enable_nonmonotone_line_search = true;
        int nonmonotone_window = 5;
        bool enable_ls_trace = true;
        bool enable_ls_base_check = true;
        double ls_base_mismatch_tol = 0.25;
        bool enable_controlled_accept_iter1 = true;
        int controlled_accept_min_armijo_reject = 6;
        double controlled_accept_relax = 2.05;
        double controlled_accept_max_rel_update = 0.20;
        int controlled_accept_max_per_step = 1;

        int ls_fail_rescue_threshold = 2;
        int ls_fail_rescue_max = 1;
        double ptc_rescue_boost = 5.0;
        double rollback_shrink_factor = 0.7;

        bool enable_ptc = false;
        double ptc_lambda_init = 1.0;
        double ptc_lambda_decay = 0.5;
        double ptc_lambda_min = 0.0;

        bool enable_control_ramp = false;
        int control_ramp_steps = 5;
        double control_ramp_min = 0.2;
        bool control_ramp_apply_rate = true;
        bool control_ramp_apply_bhp = true;
        bool enable_ramp_dt_protection = true;
        double ramp_dt_min_scale = 0.98;

        int dt_relres_iter_grow_hi = 8;
        int dt_relres_iter_neutral_hi = 14;
        int dt_relres_iter_soft_shrink_hi = 20;
        double dt_relres_grow_factor = 1.08;
        double dt_relres_neutral_factor = 1.00;
        double dt_relres_soft_shrink_factor = 0.98;
        double dt_relres_hard_shrink_factor = 0.92;

        bool enable_row_scaling = true;
        double row_scale_floor = 1.0;
        double row_scale_min = 1.0e-12;
        double row_scale_max = 1.0e+12;

        double max_dP = 1.0e4;
        double max_dT = 2.0;
        double max_dSw = 0.05;
        double min_alpha = 1.0e-8;
        bool enable_alpha_safe_two_phase = true;
        double sw_safe_eps = 1.0e-8;
        double sw_alpha_shrink = 0.98;
        bool clamp_state_to_eos_bounds = false;
        bool enforce_eos_domain = false;

        LinearSolverType lin_solver = LinearSolverType::AMGCL;
        double bicgstab_droptol = 1e-2;
        double well_source_sign = 1.0;
        Vector gravity_vector = Vector(0.0, 0.0, -9.81);

        double amgcl_tol = 1.0e-6;
        int amgcl_maxiter = 500;
        bool amgcl_use_fallback_sparselu = true;
        bool amgcl_log_on_failure = true;

        DiagLevel diag_level = DiagLevel::Summary;
        int diag_print_every_iter = 1;
        double diag_blowup_factor = 5.0;
        int diag_hot_repeat_iters = 3;
        double diag_hot_res_change_tol = 1e-2;
        int diag_clamp_trigger = 20;
        int diag_max_hot_conn = 5;
        int diag_max_clamp_dump = 10;
        double diag_flux_spike_factor = 10.0;
        double diag_eos_near_bound_ratio = 0.02;
        bool diag_incident_once_per_step = true;
    };

    template <typename T, typename = void>
    struct is_iterative_solver : std::false_type {};

    template <typename T>
    struct is_iterative_solver<T, std::void_t<decltype(std::declval<T>().iterations())>> : std::true_type {};

    template <typename MeshMgrType, typename FieldMgrType>
    struct TransientOptionalModules {
        std::function<void(MeshMgrType&, FieldMgrType&)> property_initializer;
        const BoundarySetting::BoundaryConditionManager* pressure_bc = nullptr;
        const BoundarySetting::BoundaryConditionManager* saturation_bc = nullptr;
        const BoundarySetting::BoundaryConditionManager* temperature_bc = nullptr;

        SinglePhaseFluidModel single_phase_fluid = SinglePhaseFluidModel::Water;
        CapRelPerm::VGParams vg_params = CapRelPerm::VGParams();
        CapRelPerm::RelPermParams rp_params = CapRelPerm::RelPermParams();

        std::function<void(const MeshMgrType&, const std::vector<Vector>&, int,
            std::vector<double>&, std::vector<double>&, std::vector<double>*)> state_initializer;
    };

} // namespace FIM_Engine

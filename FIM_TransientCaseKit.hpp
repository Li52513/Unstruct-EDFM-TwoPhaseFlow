/**
 * @file FIM_TransientCaseKit.hpp
 * @brief Reusable builders for transient FIM test cases (mesh-size init, BC presets, property modules)
 */
#pragma once

#include "FIM_TransientEngine.hpp"

#include "2D_PhysicalPropertiesManager.h"
#include "3D_PhysicalPropertiesManager.h"
#include "BoundaryConditionManager.h"
#include "MeshDefinitions.h"
#include "AABB.h"
#include "CapRelPerm_HD.h"

#include <string>
#include <utility>
#include <vector>
#include <cmath>
#include <algorithm>

namespace FIM_CaseKit {

using BoundaryConditionManager = BoundarySetting::BoundaryConditionManager;

struct PropertyPreset2D {
    SolidProperties_RockMatrix rock_bg{};
    SolidProperties_RockMatrix rock_region{};
    AABB rock_region_box{};
    bool enable_rock_region = true;

    SolidProperties_Frac frac{};

    FIM_Engine::SinglePhaseFluidModel single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::Water;
    CapRelPerm::VGParams vg{};
    CapRelPerm::RelPermParams rp{};
};

struct PropertyPreset3D {
    RockPropertyParams rock_bg{};
    RockPropertyParams rock_region{};
    BoundingBox3D rock_region_box{};
    bool enable_rock_region = true;

    FractureGlobalParams frac{};

    FIM_Engine::SinglePhaseFluidModel single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::Water;
    CapRelPerm::VGParams vg{};
    CapRelPerm::RelPermParams rp{};
};

inline PropertyPreset2D MakeDefaultPropertyPreset2D() {
    PropertyPreset2D p;

    p.rock_bg.phi_r = 0.18;
    p.rock_bg.kxx = 1.0e-13;
    p.rock_bg.kyy = 1.0e-13;
    p.rock_bg.kzz = 1.0e-13;
    p.rock_bg.compressibility = 1.0e-9;
    p.rock_bg.rho_r = 2600.0;
    p.rock_bg.cp_r = 1000.0;
    p.rock_bg.k_r = 2.0;

    p.rock_region = p.rock_bg;
    p.rock_region.kxx = 3.0e-13;
    p.rock_region.kyy = 3.0e-13;
    p.rock_region_box = AABB(Vector(3.0, 3.0, -1.0), Vector(7.0, 7.0, 1.0));

    p.frac.phi_f = 0.6;
    p.frac.permeability = 1.0e-11;
    p.frac.compressibility = 1.0e-9;
    p.frac.aperture = 1.0e-3;
    p.frac.rho_f = 2400.0;
    p.frac.cp_f = 900.0;
    p.frac.k_f = 2.0;

    p.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::Water;
    p.vg.alpha = 1.0e-5;
    p.vg.n = 2.0;
    p.vg.Swr = 0.0;
    p.vg.Sgr = 0.0;
    p.vg.Pc_max = 1.0e7;
    p.vg.Se_eps = 1.0e-3;
    p.rp.L = 0.5;

    return p;
}

inline PropertyPreset3D MakeDefaultPropertyPreset3D() {
    PropertyPreset3D p;

    p.rock_bg = RockPropertyParams(100.0, 100.0, 50.0, 0.16, 2600.0, 1000.0, 2.0, 1.0e-9);
    p.rock_region = RockPropertyParams(300.0, 300.0, 120.0, 0.18, 2600.0, 1000.0, 2.0, 1.0e-9);
    p.rock_region_box = BoundingBox3D{2.0, 8.0, 2.0, 8.0, 0.0, 10.0};

    p.frac = FractureGlobalParams(0.7, 2400.0, 900.0, 2.0, 0.1);

    p.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::Water;
    p.vg.alpha = 1.0e-5;
    p.vg.n = 2.0;
    p.vg.Swr = 0.0;
    p.vg.Sgr = 0.0;
    p.vg.Pc_max = 1.0e7;
    p.vg.Se_eps = 1.0e-3;
    p.rp.L = 0.5;

    return p;
}

inline int MatrixBlockCount(const MeshManager& mgr) {
    return mgr.getMatrixDOFCount();
}

inline int MatrixBlockCount(const MeshManager_3D& mgr) {
    return mgr.fracture_network().getSolverIndexOffset();
}

inline void InitFieldManager(const MeshManager& mgr, FieldManager_2D& fm) {
    size_t nNNC = 0;
    for (const auto& kv : mgr.getNNCTopologyMap()) {
        nNNC += kv.second.size();
    }

    size_t nFI = 0;
    for (const auto& frac : mgr.fracture_network().fractures) {
        if (frac.elements.size() > 1) nFI += frac.elements.size() - 1;
    }

    const size_t nFF = mgr.fracture_network().globalFFPts.size() * 6;
    fm.InitSizes(
        mgr.mesh().getGridCount(),
        mgr.fracture_network().getOrderedFractureElements().size(),
        nNNC,
        nFF,
        mgr.mesh().getFaces().size(),
        nFI);
}

inline void InitFieldManager(const MeshManager_3D& mgr, FieldManager_3D& fm) {
    const size_t nNNC = mgr.getInteractionPairs().size();
    const size_t nFI = mgr.fracture_network().getGlobalEdges().size();
    const size_t nFF = mgr.fracture_network().getOrderedFractureElements().size() * 4;
    fm.InitSizes(
        mgr.mesh().getGridCount(),
        mgr.fracture_network().getOrderedFractureElements().size(),
        nNNC,
        nFF,
        mgr.mesh().getFaces().size(),
        nFI);
}

// === 修改开始：更新 BuildSolverParams 的默认策略 (约第 159 行) ===
inline FIM_Engine::TransientSolverParams BuildSolverParams(bool twoPhase,
    int maxSteps = 50,
    double dtInit = 0.25) {
    FIM_Engine::TransientSolverParams params;
    params.well_source_sign = 1.0;
    params.max_steps = maxSteps;
    params.dt_init = dtInit;

    // Prefer AMGCL when available; otherwise fallback to SparseLU at compile-time.
#if FIM_HAS_AMGCL
    params.lin_solver = FIM_Engine::LinearSolverType::AMGCL;
#else
    params.lin_solver = FIM_Engine::LinearSolverType::SparseLU;
#endif

    // Day6 transient presets
    params.max_newton_iter = 16;
    params.max_dP = 2.0e5;
    params.max_dT = 5.0;
    params.max_dSw = twoPhase ? 0.08 : params.max_dSw;

    params.abs_res_tol = 1.0e-4;
    params.rel_res_tol = 2.0e-1;
    params.rel_update_tol = 1.0e-5;

    params.enable_nonmonotone_line_search = true;
    params.enable_control_ramp = true;
    params.control_ramp_steps = 10;
    params.control_ramp_min = 0.15;
    params.enable_ramp_dt_protection = true;
    params.ramp_dt_min_scale = 0.98;

    params.enable_ptc = true;
    params.ptc_lambda_init = 2.0;
    params.ptc_lambda_decay = 0.8;
    params.ptc_lambda_min = 0.05;

    params.dt_relres_iter_grow_hi = 8;
    params.dt_relres_iter_neutral_hi = 16;
    params.dt_relres_iter_soft_shrink_hi = 22;
    params.dt_relres_grow_factor = 1.05;
    params.dt_relres_neutral_factor = 1.00;
    params.dt_relres_soft_shrink_factor = 1.00;
    params.dt_relres_hard_shrink_factor = 0.95;

    // Turn off overly verbose diagnostics for acceptance runs
    params.enable_ls_trace = false;
    params.diag_level = FIM_Engine::DiagLevel::Off;
    params.stagnation_min_drop = 5.0e-3;
    params.stagnation_abs_res_tol = 5.0e6;

    return params;
}


enum class LongRunScenario {
    S2D_SP,
    S2D_TP,
    S3D_SP,
    S3D_TP
};

inline double DefaultStartupDays(LongRunScenario scenario) {
    switch (scenario) {
    case LongRunScenario::S2D_SP: return 7.0;
    case LongRunScenario::S2D_TP: return 14.0;
    case LongRunScenario::S3D_SP: return 21.0;
    case LongRunScenario::S3D_TP: return 30.0;
    default: return 14.0;
    }
}

inline FIM_Engine::TransientSolverParams BuildLongRunTemplate(
    LongRunScenario scenario,
    double horizon_years,
    double startup_days = -1.0) {

    const bool is3d = (scenario == LongRunScenario::S3D_SP || scenario == LongRunScenario::S3D_TP);
    const bool twoPhase = (scenario == LongRunScenario::S2D_TP || scenario == LongRunScenario::S3D_TP);
    const double years = std::max(1.0e-6, horizon_years);
    const double year_s = 365.0 * 86400.0;
    const double target_end_time_s = years * year_s;

    double dt_init = 5.0e-2;
    if (twoPhase && !is3d) dt_init = 1.0e-2;
    if (!twoPhase && is3d) dt_init = 1.0e-2;
    if (twoPhase && is3d) dt_init = 5.0e-3;

    const int max_steps = static_cast<int>(std::min(1000000.0, std::max(50000.0, years * 20000.0)));
    auto params = BuildSolverParams(twoPhase, max_steps, dt_init);

    params.enable_two_stage_profile = true;
    params.target_end_time_s = target_end_time_s;
    params.dt_min = 1.0e-6;

    const double startup_base_days = DefaultStartupDays(scenario);
    const double startup_days_eff = (startup_days > 0.0)
        ? startup_days
        : std::max(3.0, std::min(90.0, startup_base_days * std::sqrt(years)));
    params.startup_end_time_s = startup_days_eff * 86400.0;

    FIM_Engine::TransientStageProfile startup = params.startup_profile;
    FIM_Engine::TransientStageProfile longrun = params.long_profile;

    startup.dt_max = (!twoPhase && !is3d) ? 3600.0
        : (twoPhase && !is3d) ? 1800.0
        : (!twoPhase && is3d) ? 1200.0
        : 600.0;
    startup.max_newton_iter = (!twoPhase && !is3d) ? 16
        : (twoPhase && !is3d) ? 20
        : (!twoPhase && is3d) ? 20
        : 24;
    startup.rel_res_tol = 2.0e-1;
    startup.rel_update_tol = 1.0e-5;
    startup.enable_ptc = true;
    startup.ptc_lambda_init = twoPhase ? 2.0 : 1.5;
    startup.ptc_lambda_decay = 0.75;
    startup.ptc_lambda_min = 0.05;
    startup.enable_control_ramp = true;
    startup.control_ramp_steps = is3d ? 24 : 16;
    startup.control_ramp_min = twoPhase ? 0.08 : 0.10;
    startup.control_ramp_apply_rate = true;
    startup.control_ramp_apply_bhp = true;
    startup.dt_relres_iter_grow_hi = twoPhase ? 8 : 10;
    startup.dt_relres_iter_neutral_hi = twoPhase ? 12 : 14;
    startup.dt_relres_iter_soft_shrink_hi = twoPhase ? 18 : 20;
    startup.dt_relres_grow_factor = 1.08;
    startup.dt_relres_neutral_factor = 1.00;
    startup.dt_relres_soft_shrink_factor = 0.98;
    startup.dt_relres_hard_shrink_factor = 0.90;

    const double long_dt_scale = (years >= 10.0) ? 2.5 : 1.0;
    longrun.dt_max = ((!twoPhase && !is3d) ? 2.0 * 86400.0
        : (twoPhase && !is3d) ? 1.0 * 86400.0
        : (!twoPhase && is3d) ? 0.5 * 86400.0
        : 0.25 * 86400.0) * long_dt_scale;
    longrun.max_newton_iter = twoPhase ? (is3d ? 14 : 12) : (is3d ? 12 : 10);
    longrun.rel_res_tol = 2.0e-1;
    longrun.rel_update_tol = 1.0e-5;
    longrun.enable_ptc = true;
    longrun.ptc_lambda_init = 0.2;
    longrun.ptc_lambda_decay = 0.5;
    longrun.ptc_lambda_min = 0.0;
    longrun.enable_control_ramp = false;
    longrun.control_ramp_steps = startup.control_ramp_steps;
    longrun.control_ramp_min = startup.control_ramp_min;
    longrun.control_ramp_apply_rate = true;
    longrun.control_ramp_apply_bhp = true;
    longrun.dt_relres_iter_grow_hi = twoPhase ? 10 : 12;
    longrun.dt_relres_iter_neutral_hi = twoPhase ? 14 : 16;
    longrun.dt_relres_iter_soft_shrink_hi = twoPhase ? 20 : 24;
    longrun.dt_relres_grow_factor = 1.15;
    longrun.dt_relres_neutral_factor = 1.03;
    longrun.dt_relres_soft_shrink_factor = 1.00;
    longrun.dt_relres_hard_shrink_factor = 0.90;

    params.startup_profile = startup;
    params.long_profile = longrun;

    params.dt_max = longrun.dt_max;
    params.max_newton_iter = longrun.max_newton_iter;
    params.rel_res_tol = longrun.rel_res_tol;
    params.rel_update_tol = longrun.rel_update_tol;
    params.enable_ptc = longrun.enable_ptc;
    params.ptc_lambda_init = longrun.ptc_lambda_init;
    params.ptc_lambda_decay = longrun.ptc_lambda_decay;
    params.ptc_lambda_min = longrun.ptc_lambda_min;
    params.enable_control_ramp = longrun.enable_control_ramp;
    params.control_ramp_steps = longrun.control_ramp_steps;
    params.control_ramp_min = longrun.control_ramp_min;
    params.control_ramp_apply_rate = longrun.control_ramp_apply_rate;
    params.control_ramp_apply_bhp = longrun.control_ramp_apply_bhp;
    params.dt_relres_iter_grow_hi = longrun.dt_relres_iter_grow_hi;
    params.dt_relres_iter_neutral_hi = longrun.dt_relres_iter_neutral_hi;
    params.dt_relres_iter_soft_shrink_hi = longrun.dt_relres_iter_soft_shrink_hi;
    params.dt_relres_grow_factor = longrun.dt_relres_grow_factor;
    params.dt_relres_neutral_factor = longrun.dt_relres_neutral_factor;
    params.dt_relres_soft_shrink_factor = longrun.dt_relres_soft_shrink_factor;
    params.dt_relres_hard_shrink_factor = longrun.dt_relres_hard_shrink_factor;

    params.startup_vtk_output_interval_s = twoPhase ? 86400.0 : 43200.0;
    params.long_vtk_output_interval_s = twoPhase ? (years >= 10.0 ? 90.0 * 86400.0 : 15.0 * 86400.0)
                                                 : (years >= 10.0 ? 180.0 * 86400.0 : 30.0 * 86400.0);

    return params;
}
inline void ConfigurePressureBC2D(BoundaryConditionManager& bc, double pLeft, double pRight, double pFar) {
    bc.Clear();
    bc.SetDirichletBC(MeshTags::LEFT, pLeft);
    bc.SetDirichletBC(MeshTags::RIGHT, pRight);
    bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bc.SetRobinBC(MeshTags::TOP, 1.0e-8, pFar);
}

inline void ConfigureTemperatureBC2D(BoundaryConditionManager& bc, double tempFar) {
    bc.Clear();
    bc.SetDirichletBC(MeshTags::LEFT, tempFar);
    bc.SetDirichletBC(MeshTags::RIGHT, tempFar);
    bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bc.SetRobinBC(MeshTags::TOP, 1.0e-4, tempFar);
}

inline void ConfigureSaturationBC2D(BoundaryConditionManager& bc, double swLeft, double swRight) {
    bc.Clear();
    bc.SetDirichletBC(MeshTags::LEFT, swLeft);
    bc.SetDirichletBC(MeshTags::RIGHT, swRight);
    bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bc.SetRobinBC(MeshTags::TOP, 1.0e-8, swRight);
}

inline void ConfigurePressureBC3D(BoundaryConditionManager& bc, double pLeft, double pRight, double pFar) {
    bc.Clear();
    bc.SetDirichletBC(MeshTags::LEFT, pLeft);
    bc.SetDirichletBC(MeshTags::RIGHT, pRight);
    bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bc.SetRobinBC(MeshTags::TOP, 1.0e-8, pFar);
    bc.SetNeumannBC(MeshTags::TAG_FRONT, 0.0);
    bc.SetRobinBC(MeshTags::TAG_BACK, 1.0e-8, pFar);
}

inline void ConfigureTemperatureBC3D(BoundaryConditionManager& bc, double tempFar) {
    bc.Clear();
    bc.SetDirichletBC(MeshTags::LEFT, tempFar);
    bc.SetDirichletBC(MeshTags::RIGHT, tempFar);
    bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bc.SetRobinBC(MeshTags::TOP, 1.0e-4, tempFar);
    bc.SetNeumannBC(MeshTags::TAG_FRONT, 0.0);
    bc.SetRobinBC(MeshTags::TAG_BACK, 1.0e-4, tempFar);
}

inline void ConfigureSaturationBC3D(BoundaryConditionManager& bc, double swLeft, double swRight) {
    bc.Clear();
    bc.SetDirichletBC(MeshTags::LEFT, swLeft);
    bc.SetDirichletBC(MeshTags::RIGHT, swRight);
    bc.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bc.SetRobinBC(MeshTags::TOP, 1.0e-8, swRight);
    bc.SetNeumannBC(MeshTags::TAG_FRONT, 0.0);
    bc.SetRobinBC(MeshTags::TAG_BACK, 1.0e-8, swRight);
}

inline FIM_Engine::TransientOptionalModules<MeshManager, FieldManager_2D>
BuildModules2D(const PropertyPreset2D& preset,
               BoundaryConditionManager* pBC,
               BoundaryConditionManager* tBC,
               BoundaryConditionManager* sBC) {
    FIM_Engine::TransientOptionalModules<MeshManager, FieldManager_2D> modules;

    modules.property_initializer = [preset](MeshManager& mgr, FieldManager_2D& fm) {
        const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();

        PhysicalPropertiesManager_2D propMgr(mgr, fm, pCfg, tCfg);

        std::vector<std::pair<std::string, std::pair<AABB, SolidProperties_RockMatrix>>> rockRegions;
        if (preset.enable_rock_region) {
            rockRegions.emplace_back("rock_region_1", std::make_pair(preset.rock_region_box, preset.rock_region));
        }

        propMgr.InitRockProperties(preset.rock_bg, rockRegions);
        propMgr.InitFractureProperties(preset.frac);
    };

    modules.pressure_bc = pBC;
    modules.temperature_bc = tBC;
    modules.saturation_bc = sBC;
    modules.single_phase_fluid = preset.single_phase_fluid;
    modules.vg_params = preset.vg;
    modules.rp_params = preset.rp;
    return modules;
}

inline FIM_Engine::TransientOptionalModules<MeshManager, FieldManager_2D>
BuildModules2D(BoundaryConditionManager* pBC,
               BoundaryConditionManager* tBC,
               BoundaryConditionManager* sBC) {
    return BuildModules2D(MakeDefaultPropertyPreset2D(), pBC, tBC, sBC);
}

inline FIM_Engine::TransientOptionalModules<MeshManager_3D, FieldManager_3D>
BuildModules3D(const PropertyPreset3D& preset,
               BoundaryConditionManager* pBC,
               BoundaryConditionManager* tBC,
               BoundaryConditionManager* sBC) {
    FIM_Engine::TransientOptionalModules<MeshManager_3D, FieldManager_3D> modules;

    modules.property_initializer = [preset](MeshManager_3D& mgr, FieldManager_3D& fm) {
        const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();

        PhysicalPropertiesManager_3D propMgr(mgr, fm, pCfg, tCfg);

        std::vector<std::pair<std::string, std::pair<BoundingBox3D, RockPropertyParams>>> rockRegions;
        if (preset.enable_rock_region) {
            rockRegions.emplace_back("rock_region_1", std::make_pair(preset.rock_region_box, preset.rock_region));
        }

        propMgr.InitRockProperties(preset.rock_bg, rockRegions);
        propMgr.InitFractureProperties(preset.frac);
    };

    modules.pressure_bc = pBC;
    modules.temperature_bc = tBC;
    modules.saturation_bc = sBC;
    modules.single_phase_fluid = preset.single_phase_fluid;
    modules.vg_params = preset.vg;
    modules.rp_params = preset.rp;
    return modules;
}

inline FIM_Engine::TransientOptionalModules<MeshManager_3D, FieldManager_3D>
BuildModules3D(BoundaryConditionManager* pBC,
               BoundaryConditionManager* tBC,
               BoundaryConditionManager* sBC) {
    return BuildModules3D(MakeDefaultPropertyPreset3D(), pBC, tBC, sBC);
}

} // namespace FIM_CaseKit



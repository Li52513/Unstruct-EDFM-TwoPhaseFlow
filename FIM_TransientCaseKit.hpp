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

#include <string>
#include <utility>
#include <vector>

namespace FIM_CaseKit {

using BoundaryConditionManager = BoundarySetting::BoundaryConditionManager;

struct PropertyPreset2D {
    SolidProperties_RockMatrix rock_bg{};
    SolidProperties_RockMatrix rock_region{};
    AABB rock_region_box{};
    bool enable_rock_region = true;

    SolidProperties_Frac frac{};
};

struct PropertyPreset3D {
    RockPropertyParams rock_bg{};
    RockPropertyParams rock_region{};
    BoundingBox3D rock_region_box{};
    bool enable_rock_region = true;

    FractureGlobalParams frac{};
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

    return p;
}

inline PropertyPreset3D MakeDefaultPropertyPreset3D() {
    PropertyPreset3D p;

    p.rock_bg = RockPropertyParams(100.0, 100.0, 50.0, 0.16, 2600.0, 1000.0, 2.0, 1.0e-9);
    p.rock_region = RockPropertyParams(300.0, 300.0, 120.0, 0.18, 2600.0, 1000.0, 2.0, 1.0e-9);
    p.rock_region_box = BoundingBox3D{2.0, 8.0, 2.0, 8.0, 0.0, 10.0};

    p.frac = FractureGlobalParams(0.7, 2400.0, 900.0, 2.0, 0.1);

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

inline FIM_Engine::TransientSolverParams BuildSolverParams(bool twoPhase,
                                                           int maxSteps = 50,
                                                           double dtInit = 0.25) {
    FIM_Engine::TransientSolverParams params;
    params.well_source_sign = -1.0;
    params.max_steps = maxSteps;
    params.dt_init = dtInit;
    params.lin_solver = FIM_Engine::LinearSolverType::SparseLU;
    // Day6 transient presets: relax startup stiffness from BC + BHP while
    // keeping convergence guardrails in place.
    params.max_newton_iter = 12;
    params.max_dP = 2.0e5;
    params.max_dT = 5.0;
    params.max_dSw = twoPhase ? 0.08 : params.max_dSw;
    params.stagnation_min_drop = 5.0e-3;
    params.stagnation_abs_res_tol = 5.0e6;
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
    return modules;
}

inline FIM_Engine::TransientOptionalModules<MeshManager_3D, FieldManager_3D>
BuildModules3D(BoundaryConditionManager* pBC,
               BoundaryConditionManager* tBC,
               BoundaryConditionManager* sBC) {
    return BuildModules3D(MakeDefaultPropertyPreset3D(), pBC, tBC, sBC);
}

} // namespace FIM_CaseKit

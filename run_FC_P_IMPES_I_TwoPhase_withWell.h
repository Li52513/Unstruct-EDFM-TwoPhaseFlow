#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "PhysicalPropertiesManager_TwoPhase.h"
#include "FC_P_IMPES_Iteration_Loop.h"
#include "WellConfig_TwoPhase.h"

// New entry for FC-IMPES-I (keeps legacy IMPES untouched for A/B comparison).
inline int run_FC_P_IMPES_I_TwoPhase_WellCase()
{
    // ---------------- 0. Mesh + registries ----------------
    const double lengthX = 100;
    const double lengthY = 100;
    const double lengthZ = 0.0;

    const int sectionNumX = 50;
    const int sectionNumY = 50;
    const int sectionNumZ = 0;

    const bool usePrism = true;
    const bool useQuadBase = false;

    std::cout << "\n ==========Mesh: BUILDING (FC-IMPES-I) ==========\n";
    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OverRelaxed);
    std::cout << " \n ==========Mesh :BUILDING COMPLETED ==========\n";

    FieldRegistry reg;
    FaceFieldRegistry freg;

    // ---------------- 1. Initialize primary fields ----------------
    std::cout << "--- FC-IMPES-I: initializing primary fields ---\n";

    VGParams vg_params;
    vg_params.Swr = 0.1;
    vg_params.Sgr = 0.15;
    vg_params.alpha = 1e-4;
    vg_params.n = 2.0;

    RelPermParams rp_params;
    rp_params.L = 0.5;

    InitFields ic;
    ic.p_w0 = 1e7;
    ic.dp_wdx = 0.0;
    ic.dp_wdy = 0.0;
    ic.dp_wdz = 0.0;
    ic.s_w = 0.2;

    IMPES_Iteration::PressureEquation_String P_Eq;
    IMPES_Iteration::SaturationEquation_String S_Eq;

    Initializer::createPrimaryFields_TwoPhase_HT_IMPES(mgr.mesh(), reg, P_Eq.pressure_field, S_Eq.saturation, "T");
    Initializer::fillBaseDistributions_TwoPhase_HT_IMPES(mgr.mesh(), reg, ic);

    ensureTransientFields_scalar(mgr.mesh(), reg, P_Eq.pressure_field, P_Eq.pressure_old_field, P_Eq.pressure_prev_field);
    ensureTransientFields_scalar(mgr.mesh(), reg, S_Eq.saturation, S_Eq.saturation_old, S_Eq.saturation_prev);
    ensureTransientFields_scalar(mgr.mesh(), reg, "p_g", "p_g_old", "p_g_prev");

    // ---------------- 2. Rock properties ----------------
    PhysicalPropertiesManager ppm;
    ppm.classifyRockRegionsByGeometry(mgr, {}, Cell::RegionType::Medium);
    ppm.UpdateRockProperties(mgr, reg, P_Eq.pressure_field, "T");

    // ---------------- 3. Pressure boundary conditions ----------------
    PressureBC::Registry pbcReg;
    {
        std::vector<PressureBC::FaceID> bnd;
        bnd.reserve(mgr.mesh().getFaces().size());
        for (const auto& F : mgr.mesh().getFaces())
        {
            if (!F.isBoundary()) continue;
            bnd.push_back(F.id);
        }
        // default: no-flow Neumann (dp/dn = 0)
        pbcReg.addNeumannGrad(bnd, 0.0);
    }
    PressureBCAdapter PbcA(pbcReg);

    // ---------------- 4. Wells (same config as legacy case) ----------------
    std::vector<WellConfig_TwoPhase> wells_cfg;
    {
        WellConfig_TwoPhase inj;
        inj.name = "INJ";
        inj.role = WellDOF_TwoPhase::Role::Injector;
        inj.mode = WellDOF_TwoPhase::Mode::Rate;
        inj.target = 1; // kg/s total mass rate
        inj.Tin = 300.0;
        inj.s_w_bh = 1.0; // pure water injection
        inj.geom.pos = Vector{ 0.1 * lengthX, 0.1 * lengthY, 0.0 };
        inj.geom.rw = 0.1;
        inj.geom.skin = 0.0;
        inj.geom.H = 10;
        inj.geom.perfRadius = 0.0;
        inj.geom.maxHitCells = 5;
        wells_cfg.push_back(inj);
    }
    {
        WellConfig_TwoPhase prod;
        prod.name = "PROD";
        prod.role = WellDOF_TwoPhase::Role::Producer;
        prod.mode = WellDOF_TwoPhase::Mode::Pressure;
        prod.target = 8e6;
        prod.Tin = 300.0;
        prod.s_w_bh = 1.0;
        prod.geom.pos = Vector{ 0.9 * lengthX, 0.9 * lengthY, 0.0 };
        prod.geom.rw = 0.1;
        prod.geom.skin = 0.0;
        prod.geom.H = 10;
        prod.geom.perfRadius = 0.0;
        prod.geom.maxHitCells = 5;
        wells_cfg.push_back(prod);
    }

    build_masks_and_WI_for_all(mgr, reg, wells_cfg);

    const int Ncells = static_cast<int>(mgr.mesh().getCells().size());
    std::vector<WellDOF_TwoPhase> wells_dof;
    register_well_dofs_for_all_TwoPhase(Ncells, wells_cfg, wells_dof);

    // ---------------- 5. Saturation config ----------------
    IMPES_Iteration::SaturationTransportConfig satCfg;
    satCfg.VG_Parameter.vg_params = vg_params;
    satCfg.VG_Parameter.relperm_params = rp_params;
    satCfg.dS_max = 0.15;
    satCfg.CFL_safety = 0.9;
    satCfg.time_control_scheme = IMPES_Iteration::SatTimeControlScheme::SimpleCFL;
    satCfg.time_integration_scheme = IMPES_Iteration::SatTimeIntegrationScheme::ExplicitEuler;
	satCfg.primary_phase = IMPES_Iteration::SaturationPrimaryPhase::Water;

    // ---------------- 6. Pressure controls ----------------
    IMPES_Iteration::PressureSolveControls pCtrl;
    pCtrl.assembly.VG_Parameter.vg_params = vg_params;
    pCtrl.assembly.VG_Parameter.relperm_params = rp_params;
    pCtrl.assembly.enable_buoyancy = false;
    pCtrl.assembly.gradient_smoothing = 6;
    pCtrl.assembly.gravity = Vector{ 0.0, 0.0, 0.0 };

    pCtrl.max_outer = 1000;
    pCtrl.tol_abs = 1e6;
    pCtrl.tol_rel = 1e-4;
    pCtrl.under_relax = 1.0;
    pCtrl.verbose = true;

    // ---------------- 7. Flux naming container ----------------
    IMPES_Iteration::FluxSplitConfig fluxCfg;
    fluxCfg.boundary_inflow_fw = 1;
    fluxCfg.enforce_boundary_inflow_fw = true;
    fluxCfg.flux_sign_epsilon = 1e-10;
    fluxCfg.pressure_bc = &PbcA;

    // FC-IMPES-I flux/operator config
    FC_P_IMPES_I::PhaseOperatorConfig fccfg;
    fccfg.tag_total = "fc_case11_total";
    fccfg.tag_w = "fc_case11_w";
    fccfg.tag_g_p = "fc_case11_g_p";
    fccfg.tag_g_pc = "fc_case11_g_pc";
    fccfg.mf_total = "mf_total_fc";
    fccfg.mf_w = "mf_w_fc";
    fccfg.mf_g = "mf_g_fc";
    fccfg.mf_g_p = "mf_g_p_fc";
    fccfg.mf_g_pc = "mf_g_pc_fc";

    // ---------------- 8. Run ----------------
    const int    nSteps = 1000;
    double       dt_initial = 1e-5;

    FC_P_IMPES_I::TimeStepControl timeCtrl;
    timeCtrl.dt_min = 1e-5;
    timeCtrl.dt_max = 100;
    timeCtrl.grow_factor = 100;

    const int writeEveryP = 10;
    const int writeEverySw = 10;

    const std::string outPrefixP = "./Postprocess_Data/IMPES_Iteration_Test/FC_P_IMPES_ExplicitEuler/p_fc_impes_ps_withwell_CASE2/p_ps";
    const std::string outPrefixSw = "./Postprocess_Data/IMPES_Iteration_Test/FC_P_IMPES_ExplicitEuler/s_fc_impes_ps_withwell_CASE2/s_ps";
    const int snapshotEveryCsv = 10;
    const std::string snapshotPrefix = "./Postprocess_Data/csv_snapshots/FC_P_IMPES_ExplicitEuler/ps_state_withwell_fc_impes_i_CASE2";

    std::cout << "--- FC-IMPES-I: start transient run ---\n";
    const bool ok = FC_P_IMPES_I::runTransient_FC_P_IMPES_I(
        mgr, reg, freg,
        PbcA,
        wells_dof,
        nSteps,
        dt_initial,
        pCtrl,
        satCfg,
        fluxCfg,
        writeEveryP,
        writeEverySw,
        outPrefixP,
        outPrefixSw,
        snapshotEveryCsv,
        snapshotPrefix,
        fccfg,
        timeCtrl);

    if (!ok)
    {
        std::cerr << "[FC-IMPES-I] transient run failed.\n";
        return EXIT_FAILURE;
    }

    std::cout << "[FC-IMPES-I] finished successfully.\n";
    return EXIT_SUCCESS;
}

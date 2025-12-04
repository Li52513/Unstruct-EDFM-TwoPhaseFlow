#pragma once
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <iostream>

#include "0_PhysicalParametesCalculateandUpdata.h"
#include "IMPES_Iteration_Loop.h"
#include "WellConfig_TwoPhase.h"
#include "TwoPhaseWells_StrictRate.h"

namespace
{
    template <typename T>
    inline T clampValue(T value, T lo, T hi)
    {
        return std::max(lo, std::min(hi, value));
    }
}

namespace WellDetail
{
    template <typename T>
    inline T clampValue(T value, T lo, T hi)
    {
        return std::max(lo, std::min(hi, value));
    }
}

/**
 * @brief Waterflood-style IMPES run: pure-water injector displaces reservoir CO2.
 *
 * Scenario summary:
 *  - Domain size extended to 240 m (x) ¡Á 60 m (y) so the flood front has travel distance.
 *  - Grid: 120 ¡Á 30 quadrilateral sections, orthogonal correction enabled.
 *  - Boundaries: no-flow Neumann on all sides; flow is entirely well-driven.
 *  - Initial water saturation is only slightly above Swr, meaning the reservoir is mainly CO2.
 *  - Injector near the left boundary (30 m, 30 m) delivers 60 kg/s of almost pure water.
 *  - Producer on the right boundary (210 m, 30 m) is pressure-controlled at 8.8 MPa,
 *    so production is initially CO2 until the water front arrives.
 */
int run_IMPES_Iteration_TwoPhase_WithWell()
{
    // ---------- 0. Geometry and mesh ----------
    const double lengthX = 240.0;
    const double lengthY = 60.0;
    const double lengthZ = 0.0;

    const int sectionNumX = 120;
    const int sectionNumY = 30;
    const int sectionNumZ = 0;

    const bool usePrism = true;
    const bool useQuadBase = false;

    std::cout << "\n===== [CASE] IMPES two-phase with rate-controlled wells =====\n";
    std::cout << "--- IMPES: building mesh (240m x 60m, 120x30 sections) ---\n";

    MeshManager mgr(lengthX, lengthY, lengthZ,
        sectionNumX, sectionNumY, sectionNumZ,
        usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    FieldRegistry     reg;
    FaceFieldRegistry freg;

    // ---------- 1. Primary field initialization ----------
    VGParams vg_params;
    vg_params.Swr = 0.12;
    vg_params.Sgr = 0.1;
    vg_params.alpha = 2.5e-4;
    vg_params.n = 2.0;

    RelPermParams rp_params;
    rp_params.L = 0.4;

    InitFields ic;
    ic.p_w0 = 1.0e7;
    ic.dp_wdx = 0.0;
    ic.dp_wdy = 0.0;
    ic.dp_wdz = 0.0;
    // Start just above Swr so the formation mainly contains CO2.
    const double sw_background = WellDetail::clampValue(vg_params.Swr + 0.01, vg_params.Swr + 1e-5, 1.0 - vg_params.Sgr - 1e-5);
    ic.s_w = sw_background;

    Initializer::createPrimaryFields_TwoPhase_HT_IMPES(
        mgr.mesh(), reg, "p_w", "s_w", "T");
    Initializer::fillBaseDistributions_TwoPhase_HT_IMPES(
        mgr.mesh(), reg, ic);

    ensureTransientFields_scalar(mgr.mesh(), reg, "p_w", "p_w_old", "p_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "s_w", "s_w_old", "s_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "p_g", "p_g_old", "p_g_prev");

    // ---------- 2. Rock/fluid preparation ----------
    PhysicalPropertiesManager ppm;
    ppm.classifyRockRegionsByGeometry(mgr, {}, Cell::RegionType::Medium);
    ppm.UpdateMatrixRockAt(mgr, reg, "p_w", "T");

    // ---------- 3. Boundary conditions (no-flow) ----------
    const auto& bfaces = mgr.boundaryFaces();
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient noFlow{ 0.0, 1.0, 0.0 };
    PressureBC::setBoxBCs2D(pbc_pw, bfaces, noFlow, noFlow, noFlow, noFlow);
    PressureBCAdapter PbcA{ pbc_pw };

    // ---------- 4. Well configurations ----------
    std::vector<WellConfig_TwoPhase> wells_cfg;

    WellConfig_TwoPhase inj;
    inj.name = "INJ_A";
    inj.role = WellDOF_TwoPhase::Role::Injector;
    inj.mode = WellDOF_TwoPhase::Mode::RateWithPressureBias;
    inj.target = 60.0; // kg/s
    inj.weak_pressure_target = 1.15e7; // Pa
    inj.weak_pressure_weight = 5e-8;
    inj.Tin = 330.0;
    inj.s_w_bh = 0.999; // virtually pure water
    inj.mu_w_inj = 6.5e-4;
    inj.mu_g_inj = 3.5e-5;
    inj.rho_w_inj = 960.0;
    inj.rho_g_inj = 550.0;
    inj.cp_w_inj = 4200.0;
    inj.cp_g_inj = 1400.0;
    inj.geom.name = inj.name;
    inj.geom.pos = Vector{ 30.0, 30.0, 0.0 };
    inj.geom.rw = 0.1;
    inj.geom.skin = 0.0;
    inj.geom.H = 30.0;
    inj.geom.perfRadius = 7.5;
    inj.geom.maxHitCells = 6;
    inj.pm_2p.reFactor = 0.28;
    inj.pm_2p.fallbackKh = 2.0e-14;
    wells_cfg.push_back(inj);

    WellConfig_TwoPhase prod;
    prod.name = "PROD_B";
    prod.role = WellDOF_TwoPhase::Role::Producer;
    prod.mode = WellDOF_TwoPhase::Mode::Pressure;
    prod.target = 8.8e6; // Pa
    prod.Tin = 315.0;
    prod.s_w_bh = 0.05;
    prod.mu_w_inj = 5e-4;
    prod.mu_g_inj = 3e-5;
    prod.rho_w_inj = 950.0;
    prod.rho_g_inj = 600.0;
    prod.cp_w_inj = 4200.0;
    prod.cp_g_inj = 1400.0;
    prod.geom.name = prod.name;
    prod.geom.pos = Vector{ 210.0, 30.0, 0.0 };
    prod.geom.rw = 0.1;
    prod.geom.skin = 0.0;
    prod.geom.H = 30.0;
    prod.geom.perfRadius = 7.5;
    prod.geom.maxHitCells = 6;
    prod.pm_2p.reFactor = 0.28;
    prod.pm_2p.fallbackKh = 2.0e-14;
    wells_cfg.push_back(prod);

    build_masks_and_WI_for_all(mgr, reg, wells_cfg);

    const int Ncells = static_cast<int>(mgr.mesh().getCells().size());
    std::vector<WellDOF_TwoPhase> wells_dof;
    register_well_dofs_for_all_TwoPhase(Ncells, wells_cfg, wells_dof);

    // ---------- 5. Saturation transport configuration ----------
    IMPES_Iteration::SaturationTransportConfig satCfg;
    satCfg.saturation = "s_w";
    satCfg.saturation_old = "s_w_old";
    satCfg.saturation_prev = "s_w_prev";
    satCfg.water_source_field = "well_source_mass_w";
    satCfg.VG_Parameter.vg_params = vg_params;
    satCfg.VG_Parameter.relperm_params = rp_params;
    satCfg.time_control_scheme = IMPES_Iteration::SatTimeControlScheme::RedondoLike;
    satCfg.CFL_safety = 0.6;
    satCfg.dS_max = 0.05;

    // ---------- 6. Pressure solver configuration ----------
    IMPES_Iteration::PressureSolveControls pCtrl;
    pCtrl.max_outer = 2;
    pCtrl.tol_abs = 1e3;
    pCtrl.tol_rel = 1e-4;
    pCtrl.under_relax = 1;
    pCtrl.verbose = true;
    pCtrl.assembly.enable_buoyancy = false;
    pCtrl.assembly.gradient_smoothing = 1;
    pCtrl.assembly.gravity = Vector{ 0.0, 0.0, 0.0 };
    pCtrl.assembly.dirichlet_zero_flux_tol = 0.0;

    // ---------- 7. Mass-flux splitting configuration ----------
    IMPES_Iteration::FluxSplitConfig fluxCfg;
    fluxCfg.rho_water = TwoPhase::Water().rho_tag;
    fluxCfg.rho_gas = TwoPhase::CO2().rho_tag;
    fluxCfg.pressure_bc = &PbcA;

    // ---------- 8. Initialize Pc/p_g/old property layers ----------
    {
        TwoPhase::updateTwoPhasePropertiesAtTimeStep(mgr, reg, "s_w", vg_params, rp_params);

        auto p_w = reg.get<volScalarField>("p_w");
        auto p_g = reg.get<volScalarField>("p_g");
        auto p_w_old = reg.get<volScalarField>("p_w_old");
        auto p_w_prev = reg.get<volScalarField>("p_w_prev");
        auto p_g_old = reg.get<volScalarField>("p_g_old");
        auto p_g_prev = reg.get<volScalarField>("p_g_prev");
        auto Pc = reg.get<volScalarField>(TwoPhase::Auxiliaryparameters().Pc_tag);
        const auto& cells = mgr.mesh().getCells();
        const auto& id2idx = mgr.mesh().getCellId2Index();
        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const double pg = (*p_w)[i] + (*Pc)[i];
            (*p_g)[i] = pg;
            (*p_g_old)[i] = pg;
            (*p_g_prev)[i] = pg;
            (*p_w_old)[i] = (*p_w)[i];
            (*p_w_prev)[i] = (*p_w)[i];
        }

        TwoPhase::updateWaterBasicPropertiesAtStep(mgr, reg, "p_w", "T");
        TwoPhase::updateCO2BasicPropertiesAtStep(mgr, reg, "p_g", "T");
        TwoPhase::copyBasicPropertiesToOldLayer(reg);
    }

    // ---------- 9. Time-stepping controls ----------
    const int    nSteps = 4500;
    double       dt_initial = 2.5e-3;

    const int writeEveryP = 1;
    const int writeEverySw = 1;

    const std::string outPrefixP = "./Postprocess_Data/IMPES_Iteration_Test/Case7/p_impes_ps_revised/p_ps";
    const std::string outPrefixSw = "./Postprocess_Data/IMPES_Iteration_Test/Case7/s_impes_ps_revised/s_ps";
    const int snapshotEveryCsv = 1;
    const std::string snapshotPrefix = "./Postprocess_Data/csv_snapshots/Case7/waterflood_ps_state";
    std::vector<std::string> snapshotFields = { "p_w", "s_w" };

    std::cout << "--- IMPES: start transient run (two-phase with wells) ---\n";

    const bool ok = IMPES_Iteration::runTransient_IMPES_Iteration
    (
        mgr,
        reg,
        freg,
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
        snapshotFields
    );

    if (!ok)
    {
        std::cerr << "[CASE] IMPES transient run with wells failed.\n";
        return EXIT_FAILURE;
    }

    auto reportWellRates = [&]()
    {
        auto modeToString = [](WellDOF_TwoPhase::Mode mode)
        {
            switch (mode)
            {
            case WellDOF_TwoPhase::Mode::Pressure: return "Pressure";
            case WellDOF_TwoPhase::Mode::Rate: return "Rate";
            case WellDOF_TwoPhase::Mode::RateWithPressureBias: return "Rate+BHP";
            default: return "Unknown";
            }
        };
        std::cout << "\n[CASE] Well mass-flow summary (kg/s):\n";
        std::vector<double> Mw_cell, Mg_cell;
        for (const auto& well : wells_dof)
        {
            if (!FVM::TwoPhaseWellsStrict::compute_well_phase_mass_rates_strict(
                mgr,
                reg,
                well,
                satCfg.VG_Parameter.vg_params,
                satCfg.VG_Parameter.relperm_params,
                pCtrl.assembly.pressure_field,
                Mw_cell,
                Mg_cell))
            {
                std::cerr << "  * " << well.name << ": unable to evaluate rates.\n";
                continue;
            }
            const double totalMw = std::accumulate(Mw_cell.begin(), Mw_cell.end(), 0.0);
            const double totalMg = std::accumulate(Mg_cell.begin(), Mg_cell.end(), 0.0);
            const double totalMass = totalMw + totalMg;
            std::cout << "  * " << well.name
                << " [" << (well.role == WellDOF_TwoPhase::Role::Injector ? "Injector" : "Producer")
                << ", " << modeToString(well.mode)
                << "] Mw=" << totalMw << ", Mg=" << totalMg << ", total=" << totalMass << "\n";
        }
    };

    reportWellRates();

    std::cout << "[CASE] IMPES wells case finished successfully.\n";
    return EXIT_SUCCESS;
}

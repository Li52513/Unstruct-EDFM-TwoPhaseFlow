#pragma once

#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#if __cplusplus >= 201703L
#include <filesystem>
#endif

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "PhysicalPropertiesManager_TwoPhase.h"

#include "Solver_TimeLoopUtils.h"
#include "FC_P_IMPES_PressureSolver.h"
#include "FluxSplitterandSolver.h"               // reuse FluxSplitConfig as a naming container
#include "SaturationTransportEqAssemblerandSolver.h"

#include "PostProcess_.h"
#include "IMPES_PostProcessIO.h"
#include "PostProcess_OutputReport.h"
#include "TwoPhaseWells_StrictRate.h"
#include "IMPES_Iteration_WellHelpers.h"
#include "IMPES_ConservationDiagnostics.h"

namespace FC_P_IMPES_I
{
    struct TimeStepControl
    {
        double dt_min = 1e-5;
        double dt_max = 100;
        double grow_factor = 1.1;
        double shrink_factor = 0.7;
        double safety_factor =10000;
        int    max_retries = 8;
    };

    inline bool runTransient_FC_P_IMPES_I
    (
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& Pbc,
        std::vector<WellDOF_TwoPhase>& wells,
        int                           nSteps,
        double                        dt_initial,
        IMPES_Iteration::PressureSolveControls& pressureCtrl,
        IMPES_Iteration::SaturationTransportConfig& satCfg,
        const IMPES_Iteration::FluxSplitConfig& fluxCfg_in,
        int                           writeEveryP = 0,
        int                           writeEverySw = 0,
        const std::string& outPrefixP = "",
        const std::string& outPrefixSw = "",
        int                           snapshotEveryCsv = 0,
        const std::string& snapshotPrefix = "",
        PhaseOperatorConfig fccfg = PhaseOperatorConfig{}
    )
    {
        if (nSteps <= 0)
        {
            std::cerr << "[FC-IMPES-I][Iteration] invalid nSteps.\n";
            return false;
        }
        if (!(dt_initial > 0.0))
        {
            std::cerr << "[FC-IMPES-I][Iteration] invalid initial dt.\n";
            return false;
        }

#if __cplusplus >= 201703L
        try
        {
            if (!outPrefixP.empty())
            {
                auto dirP = std::filesystem::path(outPrefixP).parent_path();
                if (!dirP.empty()) std::filesystem::create_directories(dirP);
            }
            if (!outPrefixSw.empty())
            {
                auto dirSw = std::filesystem::path(outPrefixSw).parent_path();
                if (!dirSw.empty()) std::filesystem::create_directories(dirSw);
            }
            if (!snapshotPrefix.empty())
            {
                auto dirSnap = std::filesystem::path(snapshotPrefix).parent_path();
                if (!dirSnap.empty()) std::filesystem::create_directories(dirSnap);
            }
        }
        catch (...)
        {
            std::cerr << "[FC-IMPES-I][Iteration] cannot create directories for Tecplot / CSV output.\n";
            return false;
        }
#endif

        // Use FluxSplitConfig only as a field-name container for downstream modules (sat/diagnostics/postprocess).
        IMPES_Iteration::FluxSplitConfig fluxCfg = fluxCfg_in;
        fluxCfg.water_mass_flux = fccfg.mf_w;
        fluxCfg.gas_mass_flux = fccfg.mf_g;

        // Optional: operator sanity check once (Discretize(sum) vs sum(Discretize)).
        {
            FaceCoeffDiffSummary diff;
            if (compareFaceCoeffsCombinedVsSummed(mgr, reg, freg, Pbc, fccfg, &diff))
            {
                std::cout << "[FC-IMPES-I][OpsCheck] max|a_total-(a_w+a_g)|=" << diff.max_abs_a
                    << " (face " << diff.face_id_max_a << "), "
                    << "max|s_total-(s_w+s_g)|=" << diff.max_abs_s
                    << " (face " << diff.face_id_max_s << ")\n";
            }
        }

        TimeStepControl Tsc;
        Mesh& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        auto p_w = reg.get<volScalarField>(pressureCtrl.assembly.pressure_field);
        auto s_w = reg.get<volScalarField>(satCfg.saturation);
        auto p_w_old = reg.get<volScalarField>(pressureCtrl.assembly.pressure_old_field);
        auto p_w_prev = reg.get<volScalarField>(pressureCtrl.assembly.pressure_prev_field);
        auto s_w_old = reg.get<volScalarField>(satCfg.saturation_old);
        if (!p_w || !s_w || !p_w_old || !p_w_prev || !s_w_old)
        {
            std::cerr << "[FC-IMPES-I][Iteration] missing pressure/saturation or old/prev fields.\n";
            return false;
        }

        TwoPhase::updateTwoPhasePropertiesAtTimeStep(
            mgr, reg,
            satCfg.saturation,
            satCfg.VG_Parameter.vg_params,
            satCfg.VG_Parameter.relperm_params);

        auto p_g = reg.get<volScalarField>(pressureCtrl.assembly.pressure_g);
        auto p_g_old = reg.get<volScalarField>("p_g_old");
        auto p_g_prev = reg.get<volScalarField>("p_g_prev");
        auto Pc = reg.get<volScalarField>(pressureCtrl.assembly.Pc_field);
        if (!p_g || !p_g_old || !p_g_prev || !Pc)
        {
            std::cerr << "[FC-IMPES-I][Iteration] missing p_g/Pc or old/prev.\n";
            return false;
        }

        // Initial p_g = p_w + Pc
        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            (*p_g)[i] = (*p_w)[i] + (*Pc)[i];
        }

        // Initialize basic properties and align *_old layers (avoid spurious ddt mass jump at step 1)
        TwoPhase::updateWaterBasicPropertiesAtStep(mgr, reg, pressureCtrl.assembly.pressure_field, "T");
        TwoPhase::updateCO2BasicPropertiesAtStep(mgr, reg, pressureCtrl.assembly.pressure_g, "T");
        TwoPhase::copyBasicPropertiesToOldLayer(reg);

        const std::string Qw_well_name = "Qw_well";
        const std::string Qg_well_name = "Qg_well";
        const std::string phaseMassBalanceCsvFile =
            snapshotPrefix.empty() ? "" : (snapshotPrefix + "_phase_mass_balance.csv");
        std::vector<std::string> snapshotFields;

        double time = 0.0;
        double dt = dt_initial;

        for (int stepId = 1; stepId <= nSteps; /* advanced on accept */)
        {
            dt = std::max(Tsc.dt_min, std::min(Tsc.dt_max, dt));
            const double t_next = time + dt;

            std::cout << "[FC-IMPES-I][Iteration] Time step " << stepId
                << " (t^n = " << time << " s, dt = " << dt
                << " s, t^{n+1} = " << t_next << " s) =====\n";

            auto rollback_after_failure = [&](const char* stage) -> bool
                {
                    auto revert_scalar = [&](const std::string& x_name,
                                             const std::string& x_old_name,
                                             const std::string& x_prev_name) -> void
                        {
                            copyField(reg, x_old_name, x_name);
                            if (!x_prev_name.empty())
                                copyField(reg, x_old_name, x_prev_name);
                        };

                    revert_scalar(
                        pressureCtrl.assembly.pressure_field,
                        pressureCtrl.assembly.pressure_old_field,
                        pressureCtrl.assembly.pressure_prev_field);
                    revert_scalar(
                        satCfg.saturation,
                        satCfg.saturation_old,
                        satCfg.saturation_prev);
                    revert_scalar(
                        pressureCtrl.assembly.pressure_g,
                        "p_g_old",
                        "p_g_prev");

                    TwoPhase::updateTwoPhasePropertiesAtTimeStep(
                        mgr, reg,
                        satCfg.saturation,
                        satCfg.VG_Parameter.vg_params,
                        satCfg.VG_Parameter.relperm_params);

                    for (const auto& c : cells)
                    {
                        if (c.id < 0) continue;
                        const size_t i = id2idx.at(c.id);
                        (*p_g)[i] = (*p_w)[i] + (*Pc)[i];
                    }

                    TwoPhase::updateWaterBasicPropertiesAtStep(mgr, reg, pressureCtrl.assembly.pressure_field, "T");
                    TwoPhase::updateCO2BasicPropertiesAtStep(mgr, reg, pressureCtrl.assembly.pressure_g, "T");
                    TwoPhase::copyBasicPropertiesToOldLayer(reg);

                    dt *= Tsc.shrink_factor;
                    if (dt < Tsc.dt_min)
                    {
                        std::cerr << "[FC-IMPES-I][Iteration] dt fell below dt_min after " << stage << ".\n";
                        return false;
                    }
                    return true;
                };

            if (!startTimeStep_scalar(
                mesh, reg,
                pressureCtrl.assembly.pressure_field,
                pressureCtrl.assembly.pressure_old_field,
                pressureCtrl.assembly.pressure_prev_field))
            {
                std::cerr << "[FC-IMPES-I][Iteration] startTimeStep for pressure failed.\n";
                return false;
            }
            if (!startTimeStep_scalar(
                mesh, reg,
                satCfg.saturation,
                satCfg.saturation_old,
                satCfg.saturation_prev))
            {
                std::cerr << "[FC-IMPES-I][Iteration] startTimeStep for saturation failed.\n";
                return false;
            }

            bool accept_step = true;

            double linRes_last = 0.0;
            int    linIters_last = 0;
            int    used_outer = 0;
            bool   p_converged = false;

            for (int it = 0; it < pressureCtrl.max_outer; ++it)
            {
                for (const auto& c : cells)
                {
                    if (c.id < 0) continue;
                    const size_t i = id2idx.at(c.id);
                    (*p_g)[i] = (*p_w)[i] + (*Pc)[i];
                }

                TwoPhase::updateWaterBasicPropertiesAtStep(mgr, reg, pressureCtrl.assembly.pressure_field, "T");
                TwoPhase::updateCO2BasicPropertiesAtStep(mgr, reg, pressureCtrl.assembly.pressure_g, "T");

                double dp_inf = 0.0;
                double linRes = 0.0;
                int    linIts = 0;

                const bool okP = solver_FC_IMPES_I_PressureEq(
                    mgr, reg, freg,
                    Pbc,
                    wells,
                    dt,
                    pressureCtrl,
                    dp_inf,
                    linRes,
                    linIts,
                    fccfg);

                if (!okP)
                {
                    std::cerr << "[FC-IMPES-I][PressureStep] solver_FC_IMPES_I_PressureEq failed at outer it=" << it << ".\n";
                    accept_step = false;
                    break;
                }

                used_outer = it + 1;
                linRes_last = linRes;
                linIters_last = linIts;

                if (pressureCtrl.verbose)
                {
                    std::cout << "[FC-IMPES-I][Pressure] outer " << used_outer
                        << " dp_inf=" << dp_inf
                        << " linRes=" << linRes_last
                        << " linIters=" << linIters_last
                        << " (URF=" << pressureCtrl.under_relax << ")\n";
                }

                if (dp_inf <= pressureCtrl.tol_abs)
                {
                    std::cout << "[FC-IMPES-I][Pressure] converged at outer iter " << it
                        << " with dp_inf=" << dp_inf << "\n";
                    p_converged = true;
                    break;
                }
            }

            if (!p_converged || !accept_step)
            {
                std::cerr << "[FC-IMPES-I][Iteration] Pressure not converged at time step "
                    << stepId << ", rollback and shrink dt.\n";
                if (!rollback_after_failure("pressure")) return false;
                continue;
            }

            // Wells -> phase source fields (kg/s per cell)
            if (!IMPES_Iteration::buildWaterAndGasSourceFieldsFromWells(
                mgr, reg, wells,
                satCfg.VG_Parameter.vg_params,
                satCfg.VG_Parameter.relperm_params,
                pressureCtrl.assembly.pressure_field,
                Qw_well_name,
                Qg_well_name))
            {
                std::cerr << "[FC-IMPES-I][Well] failed to build well phase sources.\n";
                return false;
            }
            satCfg.water_source_field = Qw_well_name;

            // Saturation explicit step (uses fluxCfg.water_mass_flux = mf_w_fc)
            IMPES_Iteration::SaturationStepStats satStats;
            bool sat_ok = false;
            switch (satCfg.time_integration_scheme)
            {
            case IMPES_Iteration::SatTimeIntegrationScheme::ExplicitEuler:
                sat_ok = IMPES_Iteration::advanceSaturationExplicit_Euler(mgr, reg, freg, satCfg, fluxCfg, dt, satStats);
                break;
            case IMPES_Iteration::SatTimeIntegrationScheme::HeunRK2:
                sat_ok = IMPES_Iteration::advanceSaturationHeun_RK2(mgr, reg, freg, satCfg, fluxCfg, dt, satStats);
                break;
            default:
                std::cerr << "[FC-IMPES-I][Saturation] unknown time integration scheme.\n";
                sat_ok = false;
                break;
            }

            if (!sat_ok)
            {
                std::cerr << "[FC-IMPES-I][Iteration] saturation update failed at step " << stepId << ".\n";
                if (!rollback_after_failure("saturation")) return false;
                continue;
            }

            std::cout << "[FC-IMPES-I][Saturation] step " << stepId
                << " (" << (satCfg.time_integration_scheme == IMPES_Iteration::SatTimeIntegrationScheme::HeunRK2 ? "Heun/RK2" : "Euler")
                << "): max_CFL=" << satStats.max_CFL
                << ", max_dS=" << satStats.max_dS
                << ", dt_suggest=" << satStats.suggested_dt << "\n";

            // Mass-balance diagnostics (per-phase residual + global mismatch)
            IMPES_Iteration::PhaseMassBalanceDiagnostics massDiag;
            {
                const bool okDiag = IMPES_Iteration::computePhaseMassBalanceDiagnostics(
                    mgr, reg, freg,
                    satCfg.saturation,
                    satCfg.saturation_old,
                    fluxCfg,
                    satCfg.water_source_field,
                    Qg_well_name,
                    dt,
                    "mass_res_w",
                    "mass_res_g",
                    massDiag);

                if (!okDiag)
                {
                    std::cerr << "[FC-IMPES-I][MassDiag] computePhaseMassBalanceDiagnostics failed at step " << stepId << ".\n";
                    reg.getOrCreate<volScalarField>("mass_res_w", cells.size(), 0.0);
                    reg.getOrCreate<volScalarField>("mass_res_g", cells.size(), 0.0);
                }
                else
                {
                    std::cout << "[FC-IMPES-I][MassDiag] step " << stepId
                        << " (t=" << t_next << " s)\n"
                        << "  Water: global_res=" << massDiag.global_res_w_rate
                        << " kg/s, max|res|=" << massDiag.max_abs_res_w_rate
                        << " (cell " << massDiag.max_abs_res_w_cellId << ")\n"
                        << "  CO2:   global_res=" << massDiag.global_res_g_rate
                        << " kg/s, max|res|=" << massDiag.max_abs_res_g_rate
                        << " (cell " << massDiag.max_abs_res_g_cellId << ")\n";

                    IMPES_Iteration::appendPhaseMassBalanceCSVRow(
                        phaseMassBalanceCsvFile,
                        stepId,
                        t_next,
                        dt,
                        massDiag);
                }
            }

            // Snapshot CSV (state fields)
            if (snapshotEveryCsv > 0 &&
                !snapshotPrefix.empty() &&
                (stepId % snapshotEveryCsv) == 0)
            {
                if (snapshotFields.empty())
                {
                    snapshotFields =
                    {
                        pressureCtrl.assembly.pressure_field,
                        satCfg.saturation,
                        "mass_res_w",
                        "mass_res_g"
                    };
                }

                if (!IMPES::Output::writeFieldSnapshotCSV(
                    snapshotPrefix,
                    stepId,
                    t_next,
                    mgr,
                    reg,
                    snapshotFields))
                {
                    std::cerr << "[FC-IMPES-I][Iteration] snapshot CSV export failed at step "
                              << stepId << ".\n";
                    return false;
                }
            }

            // Tecplot output (optional)
            const bool dumpP = !outPrefixP.empty() && ((writeEveryP <= 0) || (stepId % writeEveryP == 0));
            const bool dumpSw = !outPrefixSw.empty() && ((writeEverySw <= 0) || (stepId % writeEverySw == 0));
            if (dumpP)
            {
                const std::vector<Vector> gradP =
                    computeCellGradients_LSQ_with_GG(
                        mgr.mesh(), reg,
                        pressureCtrl.assembly.pressure_field.c_str(), 0);

                std::ostringstream fnP;
                fnP << outPrefixP << "_step_" << std::setw(5) << std::setfill('0')
                    << stepId << ".plt";

                const std::string faceFieldName = pressureCtrl.assembly.pressure_field + "_face_tmp";
                const bool okPltP = outputTecplot_cellToFaceToNode_BC(
                    mgr, reg, freg,
                    nullptr, &Pbc,
                    pressureCtrl.assembly.pressure_field,
                    faceFieldName,
                    &gradP,
                    fnP.str());

                if (!okPltP)
                {
                    std::cerr << "[FC-IMPES-I][Iteration] Tecplot export (pressure) failed at step "
                              << stepId << ".\n";
                    return false;
                }
            }
            if (dumpSw)
            {
                const std::vector<Vector> gradSw =
                    computeCellGradients_LSQ_with_GG(
                        mgr.mesh(), reg,
                        satCfg.saturation.c_str(), 0);

                std::ostringstream fnSw;
                fnSw << outPrefixSw << "_step_" << std::setw(5) << std::setfill('0')
                    << stepId << ".plt";

                const std::string faceFieldName = satCfg.saturation + "_face_tmp";
                const bool okPltSw = outputTecplot_cellToFaceToNode_BC(
                    mgr, reg, freg,
                    nullptr, nullptr,
                    satCfg.saturation,
                    faceFieldName,
                    &gradSw,
                    fnSw.str());

                if (!okPltSw)
                {
                    std::cerr << "[FC-IMPES-I][Iteration] Tecplot export (saturation) failed at step "
                              << stepId << ".\n";
                    return false;
                }
            }

            // Accept step: update properties to the new time layer (for next step)
            TwoPhase::updateTwoPhasePropertiesAtTimeStep(
                mgr, reg,
                satCfg.saturation,
                satCfg.VG_Parameter.vg_params,
                satCfg.VG_Parameter.relperm_params);

            for (const auto& c : cells)
            {
                if (c.id < 0) continue;
                const size_t i = id2idx.at(c.id);
                (*p_g)[i] = (*p_w)[i] + (*Pc)[i];
            }

            TwoPhase::updateWaterBasicPropertiesAtStep(mgr, reg, pressureCtrl.assembly.pressure_field, "T");
            TwoPhase::updateCO2BasicPropertiesAtStep(mgr, reg, pressureCtrl.assembly.pressure_g, "T");
            TwoPhase::copyBasicPropertiesToOldLayer(reg);

            time = t_next;
            ++stepId;

            dt = std::min(Tsc.dt_max, dt * Tsc.grow_factor);
            if (satStats.suggested_dt > 0.0)
                dt = std::min(dt, satStats.suggested_dt * Tsc.safety_factor);
        }

        return true;
    }
}

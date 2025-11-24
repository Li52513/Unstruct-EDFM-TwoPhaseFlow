#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "FVM_WellDOF_TwoPhase.h"
#include "FVM_WellCoupling_TwoPhase.h"
#include "InitConfig.h"
#include "IMPES_CommonUtils.h"
#include "IMPES_FluxSplitterandSolver.h"
#include "TwoPhaseWells_StrictRate.h"

namespace IMPES
{
    /**
     * @struct SaturationTransportConfig
     * @brief Configuration for the explicit water-saturation update.
     */
    struct SaturationTransportConfig
    {
        std::string saturation = "s_w";            // Current water saturation
        std::string saturation_old = "s_w_old";    // Previous time level saturation
        std::string saturation_prev = "s_w_prev";  // Outer-iteration reference (if any)
        std::string water_flux = "mf_w";           // Face water mass flux (kg/s)
        std::string porosity = "phi_r";            // Porosity
        std::string rho_water = "rho_w";           // Water density (kg/m3)
        std::string rho_gas = "rho_g";             // Gas density (kg/m3) if needed for volume balance
        std::string pressure = "p_w";              // Pressure field for well source calculations
        double thickness = 1.0;
        double Sw_min = 0.0;
        double Sw_max = 1.0;

        /// CFL 控制
        double CFL_limit = 0.5;        ///< 建议的 CFL 上限
        bool   track_cfl = true;       ///< 是否启用 CFL 统计/控制

        /// 饱和度抖动控制：限制单步内 |ΔS_w| 的最大值
        double dS_limit = 0.2;         ///< 允许的单步最大 |ΔS_w|，典型 0.1~0.2
        bool   track_saturation_jitter = false; ///< 是否基于 |ΔS_w| 控制时间步长

        VGParams vg_params;                        // VG parameters for well source computation
        RelPermParams rp_params;
        bool include_well_sources = true;
    };

    /**
     * @brief Report containing CFL info and (optionally) mass balance diagnostics.
     */
    struct SaturationTransportReport
    {
        double max_CFL = 0.0;  ///< 本时间步的最大 CFL
        double suggested_dt = 0.0;  ///< 同时考虑 CFL + |ΔS_w| 的推荐时间步
        double mass_balance_err = 0.0;///< 预留的质量平衡误差
        double max_dS = 0.0;  ///< 本时间步的最大 |ΔS_w|
    };

    namespace detail
    {
        /**
         * @brief Accumulate net mass flux for a cell from face fluxes (positive = inflow).
         */
        inline double accumulateCellFlux(
            const Mesh& mesh,
            const FieldRegistry& reg,
            const FaceFieldRegistry& freg,
            const std::string& flux_name,
            int cellIndex,
            double thickness)
        {
            auto mf = freg.get<faceScalarField>(flux_name);
            if (!mf) return 0.0;

            const auto& cells = mesh.getCells();
            const auto& faces = mesh.getFaces();
            if (cellIndex < 0 || cellIndex >= static_cast<int>(cells.size()))
            {
                std::cerr << "[IMPES][Sat] invalid cell index " << cellIndex
                          << " while accumulating flux for field '" << flux_name << "'.\n";
                return 0.0;
            }

            const Cell& c = cells[cellIndex];
            double acc = 0.0;
            bool faceIdxWarning = false;
            bool fluxIdxWarning = false;
            for (int faceId : c.CellFaceIDs)
            {
                const int faceIdxSigned = faceId - 1; // Cell stores 1-based face ids
                if (faceIdxSigned < 0 || faceIdxSigned >= static_cast<int>(faces.size()))
                {
                    if (!faceIdxWarning)
                    {
                        std::cerr << "[IMPES][Sat] invalid face index (1-based id) " << faceId
                                  << " encountered while computing flux for cell " << c.id << ".\n";
                        faceIdxWarning = true;
                    }
                    continue;
                }
                const size_t faceIdx = static_cast<size_t>(faceIdxSigned);
                if (faceIdx >= mf->data.size())
                {
                    if (!fluxIdxWarning)
                    {
                    std::cerr << "[IMPES][Sat] face field '" << flux_name
                              << "' size (" << mf->data.size()
                              << ") is smaller than requested index " << faceIdx
                              << " (for face id " << faceId << ").\n";
                        fluxIdxWarning = true;
                    }
                    continue;
                }
                const Face& f = faces[faceIdx];
                const double mdot = (*mf)[faceIdx];
                if (f.ownerCell == c.id)
                {
                    acc -= mdot;
                }
                else if (f.neighborCell == c.id)
                {
                    acc += mdot;
                }
            }
            return acc * thickness;
        }

    } // namespace detail

    /**
     * @brief Explicit advance of water saturation using phase mass fluxes.
     * Computes CFL and suggested dt if requested; can include well source terms.
     */
    inline bool advanceSaturationExplicit
    (
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const std::vector<WellDOF_TwoPhase>& wells,
        double dt,
        const SaturationTransportConfig& cfg,
        SaturationTransportReport* report = nullptr
    )
    {
        if (!(dt > 0.0))
        {
            std::cerr << "[IMPES][Sat] invalid dt.\n";
            return false;
        }

        auto Sw = reg.get<volScalarField>(cfg.saturation.c_str());
        auto Sw_old = reg.get<volScalarField>(cfg.saturation_old.c_str());
        if (!Sw || !Sw_old)
        {
            std::cerr << "[IMPES][Sat] missing saturation fields.\n";
            return false;
        }

        auto phi = reg.get<volScalarField>(cfg.porosity.c_str());
        auto rho_w = reg.get<volScalarField>(cfg.rho_water.c_str());
        if (!phi || !rho_w)
        {
            std::cerr << "[IMPES][Sat] missing porosity or rho_w fields.\n";
            return false;
        }
        auto p_field = reg.get<volScalarField>(cfg.pressure.c_str());
        if (cfg.include_well_sources && !wells.empty() && !p_field)
        {
            std::cerr << "[IMPES][Sat] missing pressure field '" << cfg.pressure << "' for well coupling.\n";
            return false;
        }

        Mesh& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        const size_t nSw = Sw->data.size();
        const size_t nSwOld = Sw_old->data.size();
        const size_t nPhi = phi->data.size();
        const size_t nRhoW = rho_w->data.size();

        std::vector<double> wellSources(Sw->data.size(), 0.0);
        if (cfg.include_well_sources && !wells.empty())
        {
            for (const auto& well : wells)
            {
                FVM::TwoPhaseWellsStrict::build_saturation_well_sources_strict(
                    mgr, reg, wells,
                    cfg.vg_params, cfg.rp_params,
                    cfg.pressure,   // "p_w"
                    wellSources);
            }
        }

        double maxCFL = 0.0;
        double max_dS = 0.0;
        const double thickness = std::max(cfg.thickness, 1e-12);
        bool missingCellIdLogged = false;
        bool indexMismatchLogged = false;
        bool fieldSizeMismatchLogged = false;

        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            auto it = id2idx.find(c.id);
            if (it == id2idx.end())
            {
                if (!missingCellIdLogged)
                {
                    std::cerr << "[IMPES][Sat] unable to locate local index for cell id "
                              << c.id << ". Subsequent cells without indices will be skipped.\n";
                    missingCellIdLogged = true;
                }
                continue;
            }

            const size_t i = static_cast<size_t>(it->second);
            if (i >= cells.size() || i >= wellSources.size())
            {
                if (!indexMismatchLogged)
                {
                    std::cerr << "[IMPES][Sat] cell index mismatch while accessing well sources "
                              << "(cell id " << c.id << ", local index " << i << ").\n";
                    indexMismatchLogged = true;
                }
                continue;
            }
            if (i >= nSw || i >= nSwOld || i >= nPhi || i >= nRhoW)
            {
                if (!fieldSizeMismatchLogged)
                {
                    std::cerr << "[IMPES][Sat] field size mismatch (cell id " << c.id
                              << ", local index " << i << "). Sizes => Sw:" << nSw
                              << " Sw_old:" << nSwOld << " phi:" << nPhi
                              << " rho_w:" << nRhoW << ".\n";
                    fieldSizeMismatchLogged = true;
                }
                continue;
            }

            const double V = std::max(c.volume * thickness, 1e-30);
            const double phi_i = clampValue((*phi)[i], 0.0, 1.0);
            const double rho_w_i = std::max((*rho_w)[i], 0.0);
            double massFlux = detail::accumulateCellFlux(mesh, reg, freg, cfg.water_flux, static_cast<int>(i), thickness);
            massFlux += wellSources[i];

            const double accumulation = phi_i * rho_w_i * V;
            const double deltaS = (massFlux * dt) / std::max(accumulation, 1e-30);
            (*Sw)[i] = clampValue((*Sw_old)[i] + deltaS, cfg.Sw_min, cfg.Sw_max);

            if (cfg.track_cfl)
            {
                const double CFL = std::abs(massFlux) * dt / std::max(accumulation, 1e-30);
                maxCFL = std::max(maxCFL, CFL);
            }

            if (cfg.track_saturation_jitter)
            {
                max_dS = std::max(max_dS, std::abs(deltaS));
            }
        }

        if (report)
        {
            report->max_CFL = maxCFL;
            report->max_dS = max_dS;

            double dt_suggest = dt;

            // 基于 CFL 的建议
            if (cfg.track_cfl && maxCFL > 1e-12)
            {
                const double dt_cfl = dt * cfg.CFL_limit / maxCFL;
                dt_suggest = std::min(dt_suggest, dt_cfl);
            }

            // 基于饱和度抖动的建议
            if (cfg.track_saturation_jitter && max_dS > 1e-12)
            {
                const double factor = cfg.dS_limit / max_dS;
                const double dt_sat = dt * factor;
                dt_suggest = std::min(dt_suggest, dt_sat);
            }

            report->suggested_dt = dt_suggest;
        }
        return true;
    }
} // namespace IMPES

namespace IMPES_revised
{
    /**
     * @brief Configuration for revised explicit water-saturation update
     *        using total mass flux + phase-splitting.
     *
     * This wraps the existing IMPES::SaturationTransportConfig and
     * IMPES::FluxSplitConfig so that:
     *  - total mass flux (from the revised pressure solve) is split
     *    into phase fluxes (mf_w, mf_g);
     *  - the explicit saturation step then reuses the original
     *    advanceSaturationExplicit implementation in a consistent way.
     */
    struct SaturationTransportConfig
    {
        // --- saturation update settings (mirrors IMPES::SaturationTransportConfig) ---
        std::string saturation        = "s_w";
        std::string saturation_old    = "s_w_old";
        std::string saturation_prev   = "s_w_prev";
        std::string porosity          = "phi_r";
        std::string rho_water         = "rho_w";
        std::string rho_gas           = "rho_g";
        std::string pressure          = "p_w";
        double thickness              = 1.0;
        double Sw_min                 = 0.0;
        double Sw_max                 = 1.0;
        /// CFL 控制
        double CFL_limit = 0.5;        ///< 建议的 CFL 上限
        bool   track_cfl = true;       ///< 是否启用 CFL 统计/控制
        /// 饱和度抖动控制：限制单步内 |ΔS_w| 的最大值
        double dS_limit = 0.2;         ///< 允许的单步最大 |ΔS_w|，典型 0.1~0.2
        bool   track_saturation_jitter = false; ///< 是否基于 |ΔS_w| 控制时间步长
        VGParams vg_params;
        RelPermParams rp_params;
        bool include_well_sources     = true;

        // --- flux-splitting settings (wraps IMPES::FluxSplitConfig) ---
        bool   recompute_phase_flux   = true;          // if true, call splitTwoPhaseMassFlux first
        std::string total_mass_flux   = "mf_total";    // input: from IMPES_revised pressure step
        std::string water_mass_flux   = "mf_w";        // output / saturation input
        std::string gas_mass_flux     = "mf_g";        // optional: phase-2 mass flux
        std::string fractional_flow_face = "fw_face";  // optional: face fw
        std::string lambda_water      = "lambda_w";
        std::string lambda_gas        = "lambda_g";
        std::string capillary_correction_flux = "mf_capillary_corr"; // optional: added to mf_w, subtracted from mf_g
        std::string gravity_correction_flux = "mf_gravity_corr";     // optional: added to mf_w, subtracted from mf_g
        double min_lambda             = 1e-30;
        double flux_sign_epsilon      = 1e-15;
    };

    using SaturationTransportReport = IMPES::SaturationTransportReport;

    /**
     * @brief Revised explicit saturation advance:
     *        1) optionally split total mass flux into phase fluxes;
     *        2) delegate to IMPES::advanceSaturationExplicit.
     */
    inline bool advanceSaturationExplicit(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const std::vector<WellDOF_TwoPhase>& wells,
        double dt,
        const SaturationTransportConfig& cfg,
        SaturationTransportReport* report = nullptr)
    {
        if (!(dt > 0.0))
        {
            std::cerr << "[IMPES_revised][Sat] invalid dt.\n";
            return false;
        }

        // 1) Optional: split total mass flux into phase fluxes, to obtain mf_w.
        if (cfg.recompute_phase_flux)
        {
            IMPES::FluxSplitConfig fcfg;
            fcfg.total_mass_flux          = cfg.total_mass_flux;
            fcfg.water_mass_flux          = cfg.water_mass_flux;
            fcfg.gas_mass_flux            = cfg.gas_mass_flux;
            fcfg.fractional_flow_face     = cfg.fractional_flow_face;
            fcfg.lambda_water             = cfg.lambda_water;
            fcfg.lambda_gas               = cfg.lambda_gas;
            fcfg.rho_water                = cfg.rho_water;
            fcfg.rho_gas                  = cfg.rho_gas;
            fcfg.saturation               = cfg.saturation;
            fcfg.capillary_correction_flux = cfg.capillary_correction_flux;
            fcfg.gravity_correction_flux   = cfg.gravity_correction_flux;
            fcfg.min_lambda               = cfg.min_lambda;
            fcfg.flux_sign_epsilon        = cfg.flux_sign_epsilon;
            fcfg.pressure_bc              = nullptr; // can be supplied externally if needed

            if (!IMPES::splitTwoPhaseMassFlux(mgr, reg, freg, fcfg))
            {
                std::cerr << "[IMPES_revised][Sat] flux splitting failed.\n";
                return false;
            }
        }

        // 2) Build base config for the original saturation updater.
        IMPES::SaturationTransportConfig baseCfg;
        baseCfg.saturation       = cfg.saturation;
        baseCfg.saturation_old   = cfg.saturation_old;
        baseCfg.saturation_prev  = cfg.saturation_prev;
        baseCfg.water_flux       = cfg.water_mass_flux;
        baseCfg.porosity         = cfg.porosity;
        baseCfg.rho_water        = cfg.rho_water;
        baseCfg.rho_gas          = cfg.rho_gas;
        baseCfg.pressure         = cfg.pressure;
        baseCfg.thickness        = cfg.thickness;
        baseCfg.Sw_min           = cfg.Sw_min;
        baseCfg.Sw_max           = cfg.Sw_max;
        baseCfg.CFL_limit        = cfg.CFL_limit;
        baseCfg.track_cfl        = cfg.track_cfl;
		baseCfg.dS_limit         = cfg.dS_limit;
		baseCfg.track_saturation_jitter = cfg.track_saturation_jitter;
        baseCfg.vg_params        = cfg.vg_params;
        baseCfg.rp_params        = cfg.rp_params;
        baseCfg.include_well_sources = cfg.include_well_sources;

        // 3) Delegate to the existing explicit saturation step.
        return IMPES::advanceSaturationExplicit(mgr, reg, freg, wells, dt, baseCfg, report);
    }

} // namespace IMPES_revised

namespace IMPES_Iteration
{
    struct SaturationTransportConfig
    {
        std::string saturation = "s_w";
        std::string saturation_old = "s_w_old";
        std::string saturation_prev = "s_w_prev";
        VGParams vg_params;                        
        RelPermParams rp_params;
    };
}

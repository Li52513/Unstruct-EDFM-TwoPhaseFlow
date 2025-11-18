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
        double thickness = 1.0;                    // Layer thickness for 2D plane models
        double Sw_min = 0.0;
        double Sw_max = 1.0;
        double CFL_limit = 0.5;                    // Recommended CFL threshold
        bool track_cfl = true;
        VGParams vg_params;                        // VG parameters for well source computation
        RelPermParams rp_params;
        bool include_well_sources = true;
    };

    /**
     * @brief Report containing CFL info and (optionally) mass balance diagnostics.
     */
    struct SaturationTransportReport
    {
        double max_CFL = 0.0;
        double suggested_dt = 0.0;
        double mass_balance_err = 0.0;
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
                double p_bh = well.target; // TODO: replace with solved BHP for rate-controlled wells.
                FVM::TwoPhaseWellCoupling::add_well_source_to_saturation_rhs(
                    wellSources,
                    mgr,
                    reg,
                    well,
                    p_bh,
                    *p_field,
                    cfg.vg_params,
                    cfg.rp_params);
            }
        }

        double maxCFL = 0.0;
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
        }

        if (report)
        {
            report->max_CFL = maxCFL;
            if (cfg.track_cfl && maxCFL > 1e-12)
            {
                report->suggested_dt = dt * cfg.CFL_limit / maxCFL;
            }
            else
            {
                report->suggested_dt = dt;
            }
        }
        return true;
    }
} // namespace IMPES


```cpp`r`n#pragma once

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

namespace IMPES
{
    template<typename T>
    inline T clampValue(T v, T lo, T hi)
    {
        if (v < lo) return lo;
        if (v > hi) return hi;
        return v;
    }
    /**
     * @brief 楗卞拰搴︽樉寮?鍗婇殣寮忔洿鏂伴厤缃?     */
    struct SaturationTransportConfig
    {
        std::string saturation = "s_w";            
        std::string porosity = "phi_r";            
        std::string rho_gas = "rho_g";             
        double thickness = 1.0;                    
        double Sw_min = 0.0;
        double Sw_max = 1.0;
        double CFL_limit = 0.5;                    
        VGParams vg_params;                        
        RelPermParams rp_params;
        bool include_well_sources = true;
    };

    struct SaturationTransportReport
    {
        double max_CFL = 0.0;
        double suggested_dt = 0.0;
        double mass_balance_err = 0.0;
    };

    namespace detail
    {
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
            const auto& id2idx = mesh.getCellId2Index();

            const Cell& c = cells[cellIndex];
            double acc = 0.0;
            for (int fid : c.CellFaceIDs)
            {
                const Face& f = faces[fid];
                const double mdot = (*mf)[fid];
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

    } 

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

        std::vector<double> wellSources(cells.size(), 0.0);
        if (cfg.include_well_sources && !wells.empty())
        {
            for (const auto& well : wells)
            {
                double p_bh = well.target; 
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

        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
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
} 
`r`n```

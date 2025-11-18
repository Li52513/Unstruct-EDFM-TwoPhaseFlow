```C++
#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "BCAdapter.h"
#include "DiffusionCentral.h"
#include "ConvectionUpwind_Flux.h"
#include "Solver_TimeLoopUtils.h"
#include "TimeTerm_IMPES.h"
#include "Solver_AssemblerCOO.h"
#include "FVM_WellCoupling_TwoPhase.h"

namespace IMPES
{
    template<typename T>
    inline T clampValue(T v, T lo, T hi)
    {
        if (v < lo) return lo;
        if (v > hi) return hi;
        return v;
    }

    struct PressureAssemblyConfig
    {
        std::string operator_tag = "p_impes";
        std::string pressure_old_field = "p_w_old";
        std::string rho_total_field = "rho_t";
        std::string drho_water_field = "drho_w_dp";
        std::string drho_gas_field = "drho_g_dp";
        std::string temperature_eval_field = "T";
        std::string kxx_field = "kxx";
        std::string kyy_field = "kyy";
        std::string kzz_field = "kzz";
        std::string total_mass_flux_name = "mf_total";
        std::string total_vol_flux_name = "Qf_total";
        std::string total_velocity_name = "ufn_total";
        double rock_compressibility = 0.0;
        Vector gravity = { 0.0, -9.81, 0.0 };
        bool enable_buoyancy = true;
        bool auto_update_total_density = true;
        int gradient_smoothing = 0;
    }


    struct PressureAssemblyResult
    {
            SparseSystemCOO system;
            std::shared_ptr<faceScalarField> total_vol_flux;
            std::shared_ptr<faceScalarField> total_face_velocity;
    };

        namespace detail
        {
            inline bool ensureTotalDensityField(
                MeshManager& mgr,
                FieldRegistry& reg,
                const PressureAssemblyConfig& cfg)
            {
                if (!cfg.auto_update_total_density || cfg.rho_total_field.empty()) return true;

                auto sw = reg.get<volScalarField>(cfg.saturation_field);
                auto rho_w = reg.get<volScalarField>(cfg.rho_w_field);
                auto rho_g = reg.get<volScalarField>(cfg.rho_g_field);
                if (!sw || !rho_w || !rho_g) return reg.has(cfg.rho_total_field);

                const size_t n = mgr.mesh().getCells().size();
                if (!cfg.rho_total_old_field.empty())
                {
                    auto existing = reg.get<volScalarField>(cfg.rho_total_field);
                    if (existing)
                    {
                        auto rho_old = reg.getOrCreate<volScalarField>(cfg.rho_total_old_field, n, 0.0);
                        rho_old->data = existing->data;
                    }
                }

                auto rho_t = reg.getOrCreate<volScalarField>(cfg.rho_total_field, n, 0.0);
                std::shared_ptr<volScalarField> drho_total;
                auto drho_w = cfg.drho_water_field.empty() ? nullptr : reg.get<volScalarField>(cfg.drho_water_field);
                auto drho_g = cfg.drho_gas_field.empty() ? nullptr : reg.get<volScalarField>(cfg.drho_gas_field);
                if (!cfg.drho_dp_field.empty() && drho_w && drho_g)
                {
                    drho_total = reg.getOrCreate<volScalarField>(cfg.drho_dp_field, n, 0.0);
                    std::fill(drho_total->data.begin(), drho_total->data.end(), 0.0);
                }

                const auto& cells = mgr.mesh().getCells();
                const auto& id2idx = mgr.mesh().getCellId2Index();
                for (const auto& c : cells)
                {
                    if (c.id < 0) continue;
                    const size_t i = id2idx.at(c.id);
                    const double sw_i = clampValue((*sw)[i], 0.0, 1.0);
                    const double rho_w_i = std::max((*rho_w)[i], 0.0);
                    const double rho_g_i = std::max((*rho_g)[i], 0.0);
                    (*rho_t)[i] = sw_i * rho_w_i + (1.0 - sw_i) * rho_g_i;

                    if (drho_total && drho_w && drho_g)
                    {
                        const double dw = (*drho_w)[i];
                        const double dg = (*drho_g)[i];
                        (*drho_total)[i] = sw_i * dw + (1.0 - sw_i) * dg;
                    }
                }
                return true;
            }

            inline void ensureFieldExists(
                MeshManager& mgr,
                FieldRegistry& reg,
                const std::string& name,
                double init_value = 0.0)
            {
                if (name.empty() || reg.has(name)) return;
                const size_t n = mgr.mesh().getCells().size();
                reg.getOrCreate<volScalarField>(name, n, init_value);
            }

            inline bool buildTotalFluxFields(
                MeshManager& mgr,
                FieldRegistry& reg,
                FaceFieldRegistry& freg,
                const PressureBCAdapter& pbc,
                const PressureAssemblyConfig& cfg,
                const OperatorFieldNames& nm)
            {
                return FVM::Convection::buildFlux_Darcy_Mass(
                    mgr, reg, freg,
                    nm.a_f_diff, nm.s_f_diff,
                    cfg.pressure_field,
                    cfg.rho_total_field,
                    cfg.total_mass_flux_name,
                    cfg.total_vol_flux_name,
                    cfg.total_velocity_name,
                    &pbc,
                    cfg.enable_buoyancy);
            }

            inline bool updatePhaseDensityAndDerivative(
                MeshManager& mgr,
                FieldRegistry& reg,
                const PressureAssemblyConfig& cfg,
                const std::string& phase,
                const std::string& rho_field,
                const std::string& drho_field)
            {
                if (rho_field.empty() || drho_field.empty()) return true;
                return computeRhoAndDrhoDpAt(
                    mgr, reg,
                    cfg.pressure_field,
                    cfg.temperature_eval_field,
                    phase,
                    rho_field,
                    drho_field);
            }
        }

        inline bool assemblePressureSystem(
            MeshManager& mgr,
            FieldRegistry& reg,
            FaceFieldRegistry& freg,
            const PressureBCAdapter& pbc,
            const std::vector<WellDOF_TwoPhase>& wells,
            double dt,
            const PressureAssemblyConfig& cfg,
            PressureAssemblyResult& result)
        {
            if (!(dt > 0.0))
            {
                std::cerr << "[IMPES][Pressure] invalid time step.\n";
                return false;
            }

            if (!detail::updatePhaseDensityAndDerivative(mgr, reg, cfg, "water", cfg.rho_w_field, cfg.drho_water_field))
            {
                std::cerr << "[IMPES][Pressure] failed to update water density/derivative.\n";
                return false;
            }
            if (!detail::updatePhaseDensityAndDerivative(mgr, reg, cfg, "co2", cfg.rho_g_field, cfg.drho_gas_field))
            {
                std::cerr << "[IMPES][Pressure] failed to update CO2 density/derivative.\n";
                return false;
            }

            if (!detail::ensureTotalDensityField(mgr, reg, cfg))
            {
                std::cerr << "[IMPES][Pressure] failed to build total density field.\n";
                return false;
            }

            detail::ensureFieldExists(mgr, reg, cfg.rho_total_old_field);
            detail::ensureFieldExists(mgr, reg, cfg.drho_dp_field);

            if (!reg.has(cfg.lambda_total_field))
            {
                std::cerr << "[IMPES][Pressure] missing total mobility field '" << cfg.lambda_total_field << "'.\n";
                return false;
            }

            const OperatorFieldNames nm = makeNames(cfg.operator_tag);


            ​

                cfg.temperature_eval_field.empty() ? cfg.pressure_field : cfg.temperature_eval_field;

            if (!TimeTerm_IMPES_Pressure(
                mgr, reg, dt,
                cfg.phi_field,
                cfg.pressure_old_field,
                cfg.pressure_field,
                cfg.rho_total_old_field,
                cfg.rho_total_field,
                cfg.drho_dp_field,
                cfg.rock_compressibility,
                nm.a_time,
                nm.b_time))
            {
                std::cerr << "[IMPES][Pressure] IMPES time-term assembly failed.\n";
                return false;
            }

            mobility_tokens.emplace_back("kxx:" + cfg.kxx_field);
            mobility_tokens.emplace_back("kyy:" + cfg.kyy_field);
            mobility_tokens.emplace_back("kzz:" + cfg.kzz_field);
            mobility_tokens.emplace_back(cfg.lambda_total_field);

            if (!FVM::Diffusion::build_FaceCoeffs_Central(
                mgr, reg, freg,
                nm.a_f_diff, nm.s_f_diff,
                cfg.pressure_field,
                mobility_tokens,
                cfg.rho_total_field,
                FVM::Diffusion::RhoFaceMethod::Linear,
                cfg.gravity,
                pbc,
                cfg.enable_buoyancy,
                cfg.gradient_smoothing))
            {
                std::cerr << "[IMPES][Pressure] diffusion operator build failed.\n";
                return false;
            }


            SparseSystemCOO sys;
            if (!assemble_COO(mgr, reg, freg, "ddt+diffusion", nm, &sys))
            {
                std::cerr << "[IMPES][Pressure] assemble_COO failed.\n";
                return false;
            }


            int N = 0;
            result.cell_lid = buildUnknownMap(mgr.mesh(), N);
            if (static_cast<size_t>(N) != result.cell_lid.size())
            {
                std::cerr << "[IMPES][Pressure] inconsistent unknown map.\n";
                return false;
            }

            for (const auto& well : wells)
            {
                FVM::TwoPhaseWellCoupling::couple_well_to_pressure_equation(
                    sys, mgr, reg, well, result.cell_lid);
            }

            if (!cfg.total_mass_flux_name.empty())
            {
                if (!detail::buildTotalFluxFields(mgr, reg, freg, pbc, cfg, nm))
                {
                    std::cerr << "[IMPES][Pressure] failed to compute total mass flux.\n";
                    return false;
                }
                result.total_mass_flux = freg.get<faceScalarField>(cfg.total_mass_flux_name);
                result.total_vol_flux = freg.get<faceScalarField>(cfg.total_vol_flux_name);
                result.total_face_velocity = freg.get<faceScalarField>(cfg.total_velocity_name);
            }
            else
            {
                result.total_mass_flux.reset();
                result.total_vol_flux.reset();
                result.total_face_velocity.reset();
            }

            result.system = std::move(sys);
            return true;
        }

    
}


```


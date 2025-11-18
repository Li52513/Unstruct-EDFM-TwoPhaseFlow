#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <memory>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "FaceSignMask.hpp"
#include "BCAdapter.h"
#include "IMPES_CommonUtils.h"

namespace IMPES
{
    /**
     * @brief Configuration for two-phase mass-flux splitting.
     *
     * Uses total Darcy mass flux and cell mobilities to split into phase mass fluxes,
     * with optional fractional-flow output on faces. Additional capillary/gravity
     * correction fluxes can be injected via optional face fields.
     */
    struct FluxSplitConfig
    {
        std::string total_mass_flux = "mf_total";      // Total mass flux from pressure step
        std::string water_mass_flux = "mf_w";          // Output: water-phase mass flux
        std::string gas_mass_flux = "mf_g";            // Output: gas-phase (CO2) mass flux
        std::string fractional_flow_face = "fw_face";  // Optional: water fractional flow on faces (can be empty)
        std::string lambda_water = "lambda_w";         // Cell field: water mobility
        std::string lambda_gas = "lambda_g";           // Cell field: gas mobility
        std::string rho_water;                         // Optional: water density for mass-based fractional flow
        std::string rho_gas;                           // Optional: gas density for mass-based fractional flow
        std::string saturation = "s_w";                // Cell field: saturation (sanity/debug)
        std::string capillary_correction_flux;         // Optional: face field (kg/s) added to water, subtracted from gas
        std::string gravity_correction_flux;           // Optional: face field (kg/s) added to water, subtracted from gas
        double min_lambda = 1e-30;                     // Safeguard to avoid divide-by-zero
        double flux_sign_epsilon = 1e-15;              // Threshold to detect face flux direction
        const PressureBCAdapter* pressure_bc = nullptr; // Optional: BC adapter used to identify sealed Neumann faces
        double no_flow_a_epsilon = 1e-30;              // |a| <= eps -> treat as Neumann
        double no_flow_c_epsilon = 1e-30;              // |c| <= eps -> zero flux Neumann
    };

    /**
     * @brief Outputs of flux splitting, including sign-mask update info.
     */
    struct FluxSplitResult
    {
        std::shared_ptr<faceScalarField> mf_w;
        std::shared_ptr<faceScalarField> mf_g;
        std::shared_ptr<faceScalarField> fw_face;
        FaceSignUpdateInfo signInfo;
        std::vector<int> flippedFaces;
    };

    /**
     * @brief Split total mass flux into water/CO2 phase fluxes using upwind mobilities.
     * @param mgr    Mesh manager (topology/connectivity).
     * @param reg    Cell field registry (mobilities, saturation).
     * @param freg   Face field registry (total and phase fluxes).
     * @param cfg    Splitting configuration.
     * @param mask   Optional: face sign mask to reuse/update upwind orientation.
     * @param result Optional: output pointers to phase fluxes and sign-mask info.
     */
    inline bool splitTwoPhaseMassFlux(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const FluxSplitConfig& cfg,
        FaceSignMask* mask = nullptr,
        FluxSplitResult* result = nullptr)
    {
        auto mf_total = freg.get<faceScalarField>(cfg.total_mass_flux.c_str());
        if (!mf_total)
        {
            std::cerr << "[IMPES][Flux] missing total mass flux field '" << cfg.total_mass_flux << "'.\n";
            return false;
        }

        auto lambda_w = reg.get<volScalarField>(cfg.lambda_water.c_str());
        auto lambda_g = reg.get<volScalarField>(cfg.lambda_gas.c_str());
        if (!lambda_w || !lambda_g)
        {
            std::cerr << "[IMPES][Flux] missing lambda fields '" << cfg.lambda_water
                << "' or '" << cfg.lambda_gas << "'.\n";
            return false;
        }

        auto rho_w = cfg.rho_water.empty() ? nullptr : reg.get<volScalarField>(cfg.rho_water.c_str());
        auto rho_g = cfg.rho_gas.empty() ? nullptr : reg.get<volScalarField>(cfg.rho_gas.c_str());
        if ((!cfg.rho_water.empty() && !rho_w) || (!cfg.rho_gas.empty() && !rho_g))
        {
            std::cerr << "[IMPES][Flux] requested density fields '" << cfg.rho_water
                << "' or '" << cfg.rho_gas << "' not found.\n";
            return false;
        }

        auto cap_corr = cfg.capillary_correction_flux.empty()
            ? nullptr
            : freg.get<faceScalarField>(cfg.capillary_correction_flux.c_str());
        if (!cfg.capillary_correction_flux.empty() && !cap_corr)
        {
            std::cerr << "[IMPES][Flux] capillary correction flux field '"
                << cfg.capillary_correction_flux << "' not found.\n";
            return false;
        }

        auto grav_corr = cfg.gravity_correction_flux.empty()
            ? nullptr
            : freg.get<faceScalarField>(cfg.gravity_correction_flux.c_str());
        if (!cfg.gravity_correction_flux.empty() && !grav_corr)
        {
            std::cerr << "[IMPES][Flux] gravity correction flux field '"
                << cfg.gravity_correction_flux << "' not found.\n";
            return false;
        }

        auto mf_w = freg.getOrCreate<faceScalarField>(cfg.water_mass_flux.c_str(), mf_total->data.size(), 0.0);
        auto mf_g = freg.getOrCreate<faceScalarField>(cfg.gas_mass_flux.c_str(), mf_total->data.size(), 0.0);
        std::shared_ptr<faceScalarField> fw_face;
        if (!cfg.fractional_flow_face.empty())
        {
            fw_face = freg.getOrCreate<faceScalarField>(cfg.fractional_flow_face.c_str(), mf_total->data.size(), 0.0);
        }

        Mesh& mesh = mgr.mesh();
        const auto& faces = mesh.getFaces();
        const auto& id2idx = mesh.getCellId2Index();

        const auto clamp_no_flow_boundary = [&](const Face& F) -> bool
        {
            if (!cfg.pressure_bc) return false;
            if (!F.isBoundary()) return false;
            double a = 0.0, b = 0.0, c = 0.0;
            if (!cfg.pressure_bc->getABC(F.id, a, b, c)) return false;
            const double aTol = std::max(cfg.no_flow_a_epsilon, 0.0);
            const double cTol = std::max(cfg.no_flow_c_epsilon, 0.0);
            if (std::abs(a) > aTol) return false;
            return std::abs(c) <= cTol;
        };

        for (const auto& F : faces)
        {
            const int iF = F.id - 1;
            if (clamp_no_flow_boundary(F))
            {
                (*mf_w)[iF] = 0.0;
                (*mf_g)[iF] = 0.0;
                if (fw_face)
                {
                    (*fw_face)[iF] = 0.0;
                }
                continue;
            }

            const double flux = (*mf_total)[iF];
            double fw_value = 0.0;
            double mfW = 0.0;
            double mfG = 0.0;

            if (std::abs(flux) > cfg.flux_sign_epsilon)
            {
                const int P = F.ownerCell;
                const int N = F.neighborCell;
                int upCell = -1;
                if (flux >= 0.0 || N < 0)
                {
                    upCell = (P >= 0) ? static_cast<int>(id2idx.at(P)) : -1;
                }
                else
                {
                    upCell = (N >= 0) ? static_cast<int>(id2idx.at(N)) : (P >= 0 ? static_cast<int>(id2idx.at(P)) : -1);
                }

                if (upCell >= 0)
                {
                    const double lamW = std::max((*lambda_w)[upCell], 0.0);
                    const double lamG = std::max((*lambda_g)[upCell], 0.0);
                    double lamW_eff = lamW;
                    double lamG_eff = lamG;
                    if (rho_w && rho_g)
                    {
                        const double rhoW = std::max((*rho_w)[upCell], 0.0);
                        const double rhoG = std::max((*rho_g)[upCell], 0.0);
                        lamW_eff *= rhoW;
                        lamG_eff *= rhoG;
                    }
                    const double denom = std::max(lamW_eff + lamG_eff, cfg.min_lambda);
                    fw_value = clampValue(lamW_eff / denom, 0.0, 1.0);
                }
                else
                {
                    fw_value = 0.0;
                }

                mfW = fw_value * flux;
                mfG = flux - mfW;
            }

            double correction = 0.0;
            if (cap_corr)
            {
                correction += (*cap_corr)[iF];
            }
            if (grav_corr)
            {
                correction += (*grav_corr)[iF];
            }

            (*mf_w)[iF] = mfW + correction;
            (*mf_g)[iF] = mfG - correction;
            if (fw_face)
            {
                (*fw_face)[iF] = fw_value;
            }
        }

        FaceSignUpdateInfo info{};
        if (mask)
        {
            std::vector<int> flippedTmp;
            auto& flipped = result ? result->flippedFaces : flippedTmp;
            info = updateFaceSignMask_fromFlux(*mf_total, cfg.flux_sign_epsilon, *mask, flipped);
        }

        if (result)
        {
            result->mf_w = mf_w;
            result->mf_g = mf_g;
            result->fw_face = fw_face;
            result->signInfo = info;
            if (!mask)
            {
                result->flippedFaces.clear();
            }
        }

        return true;
    }
} // namespace IMPES

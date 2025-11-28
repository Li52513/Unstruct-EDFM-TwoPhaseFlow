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

#include "SolverContrlStrName.h"
#include "0_PhysicalParametesCalculateandUpdata.h"

namespace IMPES_Iteration
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
        FluxSplitConfig_String Fsc;
        std::string water_mass_flux =               Fsc.water_mass_flux;          // Output: water-phase mass flux
        std::string gas_mass_flux =                 Fsc.gas_mass_flux;            // Output: gas-phase (CO2) mass flux
        std::string fractional_flow_face =          Fsc.fractional_flow_face;  // Optional: water fractional flow on faces (can be empty)

        std::string rho_water =                     Fsc.rho_water;                         // Optional: water density for mass-based fractional flow
        std::string rho_gas =                       Fsc.rho_gas;                           // Optional: gas density for mass-based fractional flow

        std::string capillary_correction_flux =     Fsc.capillary_correction_flux;         // Optional: face field (kg/s) added to water, subtracted from gas
        std::string gravity_correction_flux =       Fsc.gravity_correction_flux;           // Optional: face field (kg/s) added to water, subtracted from gas

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
     * @brief Physically consistent split of total mass flux into water/CO2 phase mass fluxes.
     *
     * Total pure-pressure mass flux F_press is given by PressureAssemblyConfig().total_mass_flux_name,
     * which corresponds to
     *   F_press = -k (lambda_w*rho_w + lambda_g*rho_g) grad(p_w)
     *
     * Capillary and gravity correction fluxes are given as:
     *   F_cap,total  = -k lambda_g rho_g grad(p_c)               (gas-phase capillary flux)
     *   F_grav,total = k (lambda_w*rho_w^2 + lambda_g*rho_g^2) g (sum of phase gravity fluxes)
     *
     * This function reconstructs phase mass fluxes as:
     *   F_w = fw_mass * F_press + F_grav,w
     *   F_g = (1-fw_mass)*F_press + F_cap,total + F_grav,g
     *
     * where gravity split is based on lambda * rho^2 weighting:
     *   F_grav,w = F_grav,total * (lambda_w*rho_w^2)/(lambda_w*rho_w^2 + lambda_g*rho_g^2)
     *   F_grav,g = F_grav,total - F_grav,w
     *
     * so that F_w + F_g = F_press + F_cap,total + F_grav,total.
     */
    inline bool splitTwoPhaseMassFlux(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const FluxSplitConfig& cfg,
        FaceSignMask* mask = nullptr,
        FluxSplitResult* result = nullptr
    )
    {
        auto mf_total = freg.get<volScalarField>(PressureEquation_String().total_mass_flux_name);
        if (!mf_total)
        {
            std::cerr << "[IMPES][Flux] missing total mass flux field '" << PressureEquation_String().total_mass_flux_name << "'.\n";
            return false;
        }

        auto lambda_w = reg.get<volScalarField>(TwoPhase::Water().lambda_w_tag);
        auto lambda_g = reg.get<volScalarField>(TwoPhase::CO2().lambda_g_tag);
        if (!lambda_w || !lambda_g)
        {
            std::cerr << "[IMPES][Flux] missing lambda fields '" << TwoPhase::Water().lambda_w_tag
                << "' or '" << TwoPhase::CO2().lambda_g_tag << "'.\n";
            return false;
        }

        auto rho_w = cfg.rho_water.empty() ? nullptr : reg.get<volScalarField>(TwoPhase::Water().rho_tag);
        auto rho_g = cfg.rho_gas.empty() ? nullptr : reg.get<volScalarField>(TwoPhase::CO2().rho_tag);

        if ((!cfg.rho_water.empty() && !rho_w) || (!cfg.rho_gas.empty() && !rho_g))
        {
            std::cerr << "[IMPES][Flux] requested density fields '" << cfg.rho_water
                << "' or '" << cfg.rho_gas << "' not found.\n";
            return false;
        }

        auto cap_corr = cfg.capillary_correction_flux.empty() ? nullptr : freg.get<faceScalarField>(PressureEquation_String().capillary_correction_flux_name);
        if (!cfg.capillary_correction_flux.empty() && !cap_corr)
        {
            std::cerr << "[IMPES][Flux] capillary correction flux field '"
                << cfg.capillary_correction_flux << "' not found.\n";
            return false;
        }

        auto grav_corr = cfg.gravity_correction_flux.empty() ? nullptr : freg.get<faceScalarField>(PressureEquation_String().gravity_correction_flux_name);
        if (!cfg.gravity_correction_flux.empty() && !grav_corr)
        {
            std::cerr << "[IMPES][Flux] gravity correction flux field '"
                << cfg.gravity_correction_flux << "' not found.\n";
            return false;
        }

        auto mf_w = freg.getOrCreate<faceScalarField>(cfg.water_mass_flux, mf_total->data.size(), 0.0);
        auto mf_g = freg.getOrCreate<faceScalarField>(cfg.gas_mass_flux, mf_total->data.size(), 0.0);
        std::shared_ptr<faceScalarField> fw_face;
        if (!cfg.fractional_flow_face.empty())
        {
            fw_face = freg.getOrCreate<faceScalarField>(cfg.fractional_flow_face, mf_total->data.size(), 0.0);
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
            if (iF < 0 || iF >= static_cast<int>(mf_total->data.size()))
                continue;

            //如果是无流边界面则通量为0
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

            const double flux_press = (*mf_total)[iF]; // F_press on this face

            // ---- 6.1 Determine upwind cell for fractional flow & weighting ---- //
            const int P = F.ownerCell;
            const int N = F.neighborCell;
            int upCell = -1;
            if (std::abs(flux_press) > cfg.flux_sign_epsilon)
            {
                if (flux_press >= 0.0 || N < 0) //对于内部网格面，为流出该网格；对于边界面即（N<0) 取P单元的值，即一阶迎风 (这里所有边界面均取内部p网格
                {
                    upCell = (P >= 0) ? static_cast<int>(id2idx.at(P)) : -1;
                }
                else  //似乎进行不到这个，因为边界面都在上面取了
                {
                    upCell = (N >= 0) ? static_cast<int>(id2idx.at(N)): (P >= 0 ? static_cast<int>(id2idx.at(P)) : -1);
                }
            }
            else
            {
                // Nearly zero pressure flux: fall back to owner (or neighbor)
                if (P >= 0)      upCell = static_cast<int>(id2idx.at(P)); 
                else if (N >= 0) upCell = static_cast<int>(id2idx.at(N)); //几乎也不存在这个工况 因为所有面都有P单元
            }
            double fw_mass = 0.0;
            double lamW = 0.0, lamG = 0.0;
            double rhoW = 1.0, rhoG = 1.0;

            if (upCell >= 0)
            {
                lamW = std::max((*lambda_w)[upCell], 0.0);
                lamG = std::max((*lambda_g)[upCell], 0.0);
                rhoW = std::max((*rho_w)[upCell], 0.0);
                rhoG = std::max((*rho_g)[upCell], 0.0);

                const double lamW_eff = lamW * rhoW;
                const double lamG_eff = lamG * rhoG;
                const double denom = std::max(lamW_eff + lamG_eff, cfg.min_lambda);

                fw_mass = clampValue(lamW_eff / denom, 0.0, 1.0);
            }
            else
            {
                fw_mass = 0.0;
            }

            //---- 6.2 Pressure-driven part: split F_press ----//
            double Fw_press = fw_mass * flux_press;
            double Fg_press = flux_press - Fw_press;

            // ---- 6.3 Gravity part: split F_grav,total based on lambda*rho^2 ---- //
            double F_grav_total = 0.0;
            if (grav_corr) F_grav_total = (*grav_corr)[iF];

            double Fw_grav = 0.0;
            double Fg_grav = 0.0;
            if (std::abs(F_grav_total) > 0)
            {
                if (upCell >= 0)
                {
                    const double w_w = lamW * rhoW * rhoW;
                    const double w_g = lamG * rhoG * rhoG;
                    const double denom_wg = std::max(w_w + w_g, cfg.min_lambda);

                    Fw_grav = F_grav_total * (w_w / denom_wg);
                    Fg_grav = F_grav_total - Fw_grav;
                }
                else
                {
                    // No valid upCell: fallback split
                    cout << "Face:" << F.id << "has no upCell";
                    Fw_grav = F_grav_total * fw_mass;
                    Fg_grav = F_grav_total - Fw_grav;
                }
            }
            // ---- 6.4 Capillary part: entire F_cap,total belongs to gas phase ---- //
            double F_cap_total = 0.0;
            if (cap_corr) F_cap_total = (*cap_corr)[iF];
            const double Fw_cap = 0.0;          // water has no direct Pc term in standard formulation
            const double Fg_cap = F_cap_total;  // gas takes full capillary flux

            // ---- 6.5 Compose phase mass fluxes ---- //
            const double Fw = Fw_press + Fw_grav;
            const double Fg = Fg_press + Fg_grav + Fg_cap;

            (*mf_w)[iF] = Fw;
            (*mf_g)[iF] = Fg;
            if (fw_face) (*fw_face)[iF] = fw_mass;
        }
        // ----- 7) Optional: update FaceSignMask using total pressure flux ----- //
        FaceSignUpdateInfo info{};
        std::vector<int> flippedTmp;
        if (mask)
        {
            auto& flipped = result ? result->flippedFaces : flippedTmp;
            info = updateFaceSignMask_fromFlux(
                *mf_total, cfg.flux_sign_epsilon, *mask, flipped);
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

}// namespace IMPES_Iteration
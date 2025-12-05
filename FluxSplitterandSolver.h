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
#include "FaceMassRateCalculate.h"

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
        double flux_sign_epsilon = 1e-1;              // Threshold to detect face flux direction

        // Inflow boundary handling: override fw on boundary inflow faces (N<0, flux<0)
        bool   enforce_boundary_inflow_fw = true;      // enable/disable override
        double boundary_inflow_fw = 1.0;               // default inlet fractional flow (mass-based), 1 = pure water

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
     * \brief 将总质量通量 + 毛细/重力质量通量拆分为水/气两相质量通量。
     *
     * 输入：
     *   - fluxNames.total_mass_flux:   F_press     (总压力驱动质量通量)
     *   - fluxNames.capillary_mass_flux: F_cap,total  (毛细质量通量，全部归气相)
     *   - fluxNames.gravity_mass_flux:   F_grav,total (重力质量通量 = F_grav,w + F_grav,g)
     *
     *   - λ_w, λ_g: 由 FieldRegistry 中的 TwoPhase::Water().lambda_w_tag /
     *               TwoPhase::CO2().lambda_g_tag 提供
     *   - ρ_w, ρ_g: 由 cfg.rho_water / cfg.rho_gas 指定的 cell 场提供
     *
     * 拆分公式（对每个面）：
     *
     *   1) 质量型分相函数（在 upwind 单元上计算）：
     *      fw_mass = (λ_w ρ_w) / (λ_w ρ_w + λ_g ρ_g)
     *
     *   2) 压力驱动部分：
     *      Fw_press = fw_mass * F_press
     *      Fg_press = F_press - Fw_press
     *
     *   3) 重力部分：
     *      F_grav,total 给定；按 λ ρ² 权重拆分：
     *      w_w = λ_w ρ_w²,  w_g = λ_g ρ_g²
     *      Fw_grav = F_grav,total * w_w / (w_w + w_g)
     *      Fg_grav = F_grav,total - Fw_grav
     *
     *   4) 毛细部分：
     *      F_cap,total 全部归气相：
     *      Fw_cap = 0
     *      Fg_cap = F_cap,total
     *
     *   5) 最终相通量：
     *      F_w = Fw_press + Fw_grav
     *      F_g = Fg_press + Fg_grav + Fg_cap
     *
     * 同时：
     *   - 对 no-flow Neumann 边界 (a≈0,b≠0,c≈0) 直接强制 F_w=F_g=0；
     *   - `boundaryInflow` 通过 F_press<0 && N<0 判定，若启用
     *     cfg.enforce_boundary_inflow_fw，则用 boundary_inflow_fw 覆盖 fw_mass。
     */
     /**
   * \brief 将总质量通量 + 毛细/重力质量通量拆分为水/气两相质量通量。
   *
   * 输入：
   *   - fluxNames.total_mass_flux:   F_press     (总压力驱动质量通量)
   *   - fluxNames.capillary_mass_flux: F_cap,total  (毛细质量通量，全部归气相)
   *   - fluxNames.gravity_mass_flux:   F_grav,total (重力质量通量 = F_grav,w + F_grav,g)
   *
   *   - λ_w, λ_g: 由 FieldRegistry 中的 TwoPhase::Water().lambda_w_tag /
   *               TwoPhase::CO2().lambda_g_tag 提供
   *   - ρ_w, ρ_g: 由 cfg.rho_water / cfg.rho_gas 指定的 cell 场提供
   *
   * 拆分公式（对每个面）：
   *
   *   1) 质量型分相函数（在 upwind 单元上计算）：
   *      fw_mass = (λ_w ρ_w) / (λ_w ρ_w + λ_g ρ_g)
   *
   *   2) 压力驱动部分：
   *      Fw_press = fw_mass * F_press
   *      Fg_press = F_press - Fw_press
   *
   *   3) 重力部分：
   *      F_grav,total 给定；按 λ ρ² 权重拆分：
   *      w_w = λ_w ρ_w²,  w_g = λ_g ρ_g²
   *      Fw_grav = F_grav,total * w_w / (w_w + w_g)
   *      Fg_grav = F_grav,total - Fw_grav
   *
   *   4) 毛细部分：
   *      F_cap,total 全部归气相：
   *      Fw_cap = 0
   *      Fg_cap = F_cap,total
   *
   *   5) 最终相通量：
   *      F_w = Fw_press + Fw_grav
   *      F_g = Fg_press + Fg_grav + Fg_cap
   *
   * 同时：
   *   - 对 no-flow Neumann 边界 (a≈0,b≠0,c≈0) 直接强制 F_w=F_g=0；
   *   - `boundaryInflow` 通过 F_press<0 && N<0 判定，若启用
   *     cfg.enforce_boundary_inflow_fw，则用 boundary_inflow_fw 覆盖 fw_mass。
   */
    inline bool splitTwoPhaseMassFlux(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const FluxSplitConfig& cfg,
        const FaceMassRateConfig& fluxNames,
        FaceSignMask* mask = nullptr,
        FluxSplitResult* result = nullptr
    )
    {
        // ===== 1) 取总质量通量 & 可选修正通量 ===== //
        auto mf_total = freg.get<faceScalarField>(fluxNames.total_mass_flux.c_str());
        if (!mf_total)
        {
            std::cerr << "[IMPES][Flux] missing total mass flux field '"
                << fluxNames.total_mass_flux << "'.\n";
            return false;
        }

        std::shared_ptr<faceScalarField> cap_corr = nullptr;
        if (!fluxNames.capillary_mass_flux.empty())
        {
            cap_corr = freg.get<faceScalarField>(fluxNames.capillary_mass_flux.c_str());
            if (!cap_corr)
            {
                std::cerr << "[IMPES][Flux] capillary mass flux field '"
                    << fluxNames.capillary_mass_flux << "' not found.\n";
                return false;
            }
        }

        std::shared_ptr<faceScalarField> grav_corr = nullptr;
        if (!fluxNames.gravity_mass_flux.empty())
        {
            grav_corr = freg.get<faceScalarField>(fluxNames.gravity_mass_flux.c_str());
            if (!grav_corr)
            {
                std::cerr << "[IMPES][Flux] gravity mass flux field '"
                    << fluxNames.gravity_mass_flux << "' not found.\n";
                return false;
            }
        }

        // ===== 2) 取 λ_w, λ_g, ρ_w, ρ_g ===== //
        auto lambda_w = reg.get<volScalarField>(TwoPhase::Water().lambda_w_tag);
        auto lambda_g = reg.get<volScalarField>(TwoPhase::CO2().lambda_g_tag);
        if (!lambda_w || !lambda_g)
        {
            std::cerr << "[IMPES][Flux] missing lambda fields '"
                << TwoPhase::Water().lambda_w_tag << "' or '"
                << TwoPhase::CO2().lambda_g_tag << "'.\n";
            return false;
        }

        if (cfg.rho_water.empty() || cfg.rho_gas.empty())
        {
            std::cerr << "[IMPES][Flux] rho_water/rho_gas field names must be provided in FluxSplitConfig.\n";
            return false;
        }

        auto rho_w = reg.get<volScalarField>(cfg.rho_water.c_str());
        auto rho_g = reg.get<volScalarField>(cfg.rho_gas.c_str());
        if (!rho_w || !rho_g)
        {
            std::cerr << "[IMPES][Flux] density fields '" << cfg.rho_water
                << "' or '" << cfg.rho_gas << "' not found.\n";
            return false;
        }

        // ===== 3) 分配输出面场 ===== //
        auto mf_w = freg.getOrCreate<faceScalarField>(
            cfg.water_mass_flux.c_str(), mf_total->data.size(), 0.0);
        auto mf_g = freg.getOrCreate<faceScalarField>(
            cfg.gas_mass_flux.c_str(), mf_total->data.size(), 0.0);

        std::shared_ptr<faceScalarField> fw_face;
        if (!cfg.fractional_flow_face.empty())
        {
            fw_face = freg.getOrCreate<faceScalarField>(
                cfg.fractional_flow_face.c_str(), mf_total->data.size(), 0.0);
        }

        Mesh& mesh = mgr.mesh();
        const auto& faces = mesh.getFaces();
        const auto& id2idx = mesh.getCellId2Index();

        // ===== 4) 辅助：识别 no-flow Neumann 边界（与 buildFlux_Darcy_Mass 对齐） ===== //
        const auto clamp_no_flow_boundary = [&](const Face& F) -> bool
            {
                if (!cfg.pressure_bc) return false;
                if (!F.isBoundary())  return false;

                double a = 0.0, b = 0.0, c = 0.0;
                if (!cfg.pressure_bc->getABC(F.id, a, b, c))
                    return false;

                const double aTol = std::max(cfg.no_flow_a_epsilon, 0.0);
                const double cTol = std::max(cfg.no_flow_c_epsilon, 0.0);

                const bool isZeroGradNeumann =
                    (std::abs(a) <= aTol && std::abs(b) > aTol);
                const bool isZeroFlux =
                    (std::abs(c) <= cTol);

                return isZeroGradNeumann && isZeroFlux;
            };

        // ===== 5) 主循环：逐面拆分 ===== //
        for (const auto& F : faces)
        {
            const int iF = F.id - 1;
            if (iF < 0 || iF >= static_cast<int>(mf_total->data.size()))
                continue;

            // 5.0 无流边界：直接置零
            if (clamp_no_flow_boundary(F))
            {
                (*mf_w)[iF] = 0.0;
                (*mf_g)[iF] = 0.0;
                if (fw_face) (*fw_face)[iF] = 0.0;
                continue;
            }

            const double flux_press = (*mf_total)[iF]; // F_press on this face

            // 5.1 upwind 单元：由总压力通量符号决定
            const int P = F.ownerCell;
            const int N = F.neighborCell;

            const bool boundaryInflow =
                (flux_press < -cfg.flux_sign_epsilon) && (N < 0);

            int upCell = -1;
            if (std::abs(flux_press) > cfg.flux_sign_epsilon)
            {
                if (flux_press >= 0.0 || N < 0)
                {
                    // 正向或边界：up = owner
                    upCell = (P >= 0) ? static_cast<int>(id2idx.at(P)) : -1;
                }
                else
                {
                    // 负向内部：up = neighbor
                    if (N >= 0) upCell = static_cast<int>(id2idx.at(N));
                    else if (P >= 0) upCell = static_cast<int>(id2idx.at(P));
                }
            }
            else
            {
                // 几乎零通量：fallback 到 owner（或 neighbor）
                if (P >= 0)      upCell = static_cast<int>(id2idx.at(P));
                else if (N >= 0) upCell = static_cast<int>(id2idx.at(N));
            }

            // 5.2 在 upwind 单元上计算质量型分相函数 fw_mass
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

            // 边界 inflow：可选地用 boundary_inflow_fw 覆盖 fw_mass
            if (boundaryInflow && cfg.enforce_boundary_inflow_fw)
            {
                fw_mass = clampValue(cfg.boundary_inflow_fw, 0.0, 1.0);
            }

            // 5.3 压力驱动部分拆分
            double Fw_press = fw_mass * flux_press;
            double Fg_press = flux_press - Fw_press;

            // 5.4 重力部分（如果提供了 grav_corr）
            double F_grav_total = 0.0;
            if (grav_corr)
                F_grav_total = (*grav_corr)[iF];

            double Fw_grav = 0.0;
            double Fg_grav = 0.0;
            if (std::abs(F_grav_total) > 0.0)
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
                    // 没有有效 upCell：fallback，用 fw_mass 拆分
                    std::cout << "[IMPES][Flux] Face " << F.id
                        << " has no valid upCell, gravity split fallback.\n";
                    Fw_grav = F_grav_total * fw_mass;
                    Fg_grav = F_grav_total - Fw_grav;
                }
            }

            // 5.5 毛细部分（全部归气相）
            double F_cap_total = 0.0;
            if (cap_corr)
                F_cap_total = (*cap_corr)[iF];

            const double Fw_cap = 0.0;
            const double Fg_cap = F_cap_total;

            // 5.6 最终两相质量通量
            const double Fw = Fw_press + Fw_grav;
            const double Fg = Fg_press + Fg_grav + Fg_cap;

            (*mf_w)[iF] = Fw;
            (*mf_g)[iF] = Fg;
            if (fw_face)
                (*fw_face)[iF] = fw_mass;
        }

        // ===== 6) 可选：用总压力通量更新 FaceSignMask ===== //
        FaceSignUpdateInfo info{};
        std::vector<int> flippedTmp;
        if (mask)
        {
            auto& flipped = (result ? result->flippedFaces : flippedTmp);
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
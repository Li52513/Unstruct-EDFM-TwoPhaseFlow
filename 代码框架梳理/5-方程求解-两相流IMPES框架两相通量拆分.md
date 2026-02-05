```C++
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
#include "PhysicalPropertiesManager_TwoPhase.h"
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

        double min_lambda = 1e-30;                     // Safeguard to avoid divide-by-zero
        double flux_sign_epsilon = 1e-1;              // Threshold to detect face flux direction

        // Inflow boundary handling: override fw on boundary inflow faces (N<0, flux<0)
        bool   enforce_boundary_inflow_fw = true;      // enable/disable override
        double boundary_inflow_fw = 1.0;               // default inlet fractional flow (mass-based), 1 = pure water

        const PressureBCAdapter* pressure_bc = nullptr; // Optional: BC adapter used to identify sealed Neumann faces
        double no_flow_a_epsilon = 1e-30;              // |a| <= eps -> treat as Neumann
        double no_flow_c_epsilon = 1e-30;              // |c| <= eps -> zero flux Neumann

        // ===== B: Phase-potential upwinding switch (A/B) =====
        bool use_phase_potential_upwind = false;   ///< false=A(原有), true=B(相势能迎风)
        // 相势能所需场名（默认取你 PressureEquation_String 的命名）
		std::string pg_field = PressureEquation_String().pressure_g;  ///< CO2 相压力场名
		std::string pw_field = PressureEquation_String().pressure_field;  ///< 水相压力场名
		std::string pc_field = PressureEquation_String().Pc_field;    ///< 毛细压力场名

        Vector gravity = { 0.0, 0.0, 0.0 };    ///< 重力加速度矢量（用于相势能计算）
        double min_rho_for_potential = 1e-12;      ///< rho 下限，避免 0 导致势能判别失效
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
    namespace detail 
    {
        /// \brief 计算点积（Vector 无 dot 时的兜底）
        inline double dotVec(const Vector& a, const Vector& b)
        {
            return a.m_x * b.m_x + a.m_y * b.m_y + a.m_z * b.m_z;
        }

        /// \brief 通过通量符号选择唯一 upCell（A 的核心；B 在边界/失败时也 fallback 用它）
        inline int selectUpCellByFlux(
            const Face& F, double flux, double eps,
            const std::unordered_map<int, int>& id2idx)
        {
            const int P = F.ownerCell;
            const int N = F.neighborCell;

			int upCell = -1;
            if (std::abs(flux) > eps)
            {
                if (flux >= 0.0 || N < 0)
                {
                    upCell = (P >= 0) ? static_cast<int>(id2idx.at(P)) : -1;
                }
                else
                {
                    if (N >= 0) upCell = static_cast<int>(id2idx.at(N));
                    else if (P >= 0) upCell = static_cast<int>(id2idx.at(P));
                }
            }
            else
            {
                if (P >= 0)      upCell = static_cast<int>(id2idx.at(P));
                else if (N >= 0) upCell = static_cast<int>(id2idx.at(N));
            }
            return upCell;
        }

        /**
         * \brief 通过相势能差选择某一相的上游单元（B 的核心）
         *
         * \f[
         * \psi_\alpha = p_\alpha - \rho_{\alpha,avg}\,(\mathbf g\cdot \mathbf x)
         * \Rightarrow \Delta\psi = (p_N-p_P) - \rho_{avg}\,\mathbf g\cdot(\mathbf x_N-\mathbf x_P)
         * \f]
         * 若 \f$\Delta\psi<0\f$，则 \f$\psi_P>\psi_N\f$，流向 P->N，上游为 P。
         */

        inline int selectUpCellByPhasePotential(
            const Face& F,
            const Mesh& mesh,
            const std::unordered_map<int, int>& id2idx,
            double pP, double pN,
            double rhoP, double rhoN,
            const Vector& g,
            double min_rho
            )
        {
            const int Pid = F.ownerCell;
            const int Nid = F.neighborCell;
            if (Pid < 0 || Nid < 0) return -1;

            const int iP = static_cast<int>(id2idx.at(Pid));
            const int iN = static_cast<int>(id2idx.at(Nid));

            const auto& cells = mesh.getCells();
            const Vector xP = cells[static_cast<size_t>(iP)].center;
            const Vector xN = cells[static_cast<size_t>(iN)].center;

            const Vector dX = Vector(xN.m_x - xP.m_x, xN.m_y - xP.m_y, xN.m_z - xP.m_z);
            const double gDotdX = dotVec(g, dX);

            const double rP = std::max(rhoP, min_rho);
			const double rN = std::max(rhoN, min_rho);
            const double rho_avg = 0.5 * (rP + rN);

            const double dpsi = (pN - pP) - rho_avg * gDotdX; // psi_N - psi_P
            return (dpsi < 0.0) ? iP : iN;
        }

        /// \brief 简单闭合检查：max |(Fw+Fg)-(Fpress+Fgrav+Fcap)|
        inline void diagnoseFluxSplitClosure(
            const faceScalarField& mf_w,
            const faceScalarField& mf_g,
            const faceScalarField& mf_press,
            const faceScalarField* mf_grav,
            const faceScalarField* mf_cap,
            const std::string& label
        )
        {
            const size_t nF = mf_press.data.size();
            double maxAbs = 0.0;
            for (size_t i = 0; i < nF; ++i)
            {
                const double rhs = mf_press.data[i]
                    + (mf_grav ? mf_grav->data[i] : 0.0)
                    + (mf_cap ? mf_cap->data[i] : 0.0);
                const double lhs = mf_w.data[i] + mf_g.data[i];
                maxAbs = std::max(maxAbs, std::abs(lhs - rhs));
            }
            std::cout << "[FluxSplitClosure][" << label << "] max |(Fw+Fg)-(Fpress+Fgrav+Fcap)| = "
                << maxAbs << "\n";
        }
	} // namespace detail

    inline bool splitTwoPhaseMassFlux_PotentialCalcute(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const FluxSplitConfig& cfg,
        const FaceMassRateConfig& fluxNames,
        FaceSignMask* mask = nullptr,
        FluxSplitResult* result = nullptr
    )
    {
        // ===== 1) 取总压力通量 & 修正通量 ===== //
		auto mf_press = freg.get<faceScalarField>(fluxNames.total_mass_flux);
        if (!mf_press)
        {
            std::cerr << "[IMPES][Flux] missing total mass flux field '"
                << fluxNames.total_mass_flux << "'.\n";
            return false;
        }
        std::shared_ptr<faceScalarField> cap_corr = nullptr;
        if (!fluxNames.capillary_mass_flux.empty())
        {
            cap_corr = freg.get<faceScalarField>(fluxNames.capillary_mass_flux);
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
            grav_corr = freg.get<faceScalarField>(fluxNames.gravity_mass_flux);
            if (!grav_corr)
            {
                std::cerr << "[IMPES][Flux] gravity mass flux field '"
                    << fluxNames.gravity_mass_flux << "' not found.\n";
                return false;
            }
        }
        // ===== 2) 取 λ 与 ρ（cell 场） ===== //
        auto lambda_w = reg.get<volScalarField>(PhysicalProperties_string::Water().lambda_w_tag);
		auto lambda_g = reg.get<volScalarField>(PhysicalProperties_string::CO2().lambda_g_tag);
		auto rho_w = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_tag);
		auto rho_g = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_tag);
        if (!lambda_w || !lambda_g || !rho_w || !rho_g)
        {
            std::cerr << "[IMPES][Flux] missing lambda/rho fields.\n";
            return false;
        }

        // ===== 3) B 模式所需压力场（仅当 use_phase_potential_upwind=true 才强制） ===== //
        std::shared_ptr<volScalarField> pw = nullptr, pc = nullptr, pg = nullptr;
        if (cfg.use_phase_potential_upwind)
        {
            pw = reg.get<volScalarField>(cfg.pw_field);
            if (!pw)
            {
                std::cerr << "[IMPES][Flux-B] missing pw field '" << cfg.pw_field << "'.\n";
                return false;
            }

            // pg 优先
            if (!cfg.pg_field.empty())
                pg = reg.get<volScalarField>(cfg.pg_field);

        }
        // ===== 4) 分配输出面场 ===== //
		auto mf_w = freg.getOrCreate<faceScalarField>(cfg.water_mass_flux, mf_press->data.size(), 0.0);
        auto mf_g = freg.getOrCreate<faceScalarField>(cfg.gas_mass_flux, mf_press->data.size(), 0.0);
        std::shared_ptr<faceScalarField> fw_face = nullptr;
        if (!cfg.fractional_flow_face.empty())
        {
            fw_face = freg.getOrCreate<faceScalarField>(
                cfg.fractional_flow_face, mf_press->data.size(), 0.0);
        }
        Mesh& mesh = mgr.mesh();
        const auto& faces = mesh.getFaces();
        const auto& id2idx = mesh.getCellId2Index();
        // ===== 5) no-flow Neumann 边界识别（与你 buildFlux_Darcy_Mass 对齐） ===== //
        const auto clamp_no_flow_boundary = [&](const Face& F) -> bool
        {
                if (!cfg.pressure_bc) return false;
                if (!F.isBoundary())  return false;

                double a = 0.0, b = 0.0, c = 0.0;
                if (!cfg.pressure_bc->getABC(F.id, a, b, c)) return false;
                const double aTol = std::max(cfg.no_flow_a_epsilon, 0.0);
                const double cTol = std::max(cfg.no_flow_c_epsilon, 0.0);
                const bool isZeroGradNeumann = (std::abs(a) <= aTol && std::abs(b) > aTol);
                const bool isZeroFlux = (std::abs(c) <= cTol);
                return isZeroGradNeumann && isZeroFlux;
        };
        // ===== 6) 主循环：逐面拆分 ===== //
        for (const auto& F : faces)
        {
            const int iF = F.id - 1;
            if (iF < 0 || iF >= static_cast<int>(mf_press->data.size())) continue;

            // 6.0 无流边界：强制 0（cap/gravity 在边界本来就是 0）
            if (clamp_no_flow_boundary(F))
            {
                (*mf_w)[iF] = 0.0;
                (*mf_g)[iF] = 0.0;
                if (fw_face) (*fw_face)[iF] = 0.0;
                continue;
            }

            const double F_press = (*mf_press)[iF];

            // 传统方式A 的唯一上游：由 F_press 决定（B 失败/边界时 fallback）
			const int upPress = detail::selectUpCellByFlux(F, F_press, cfg.flux_sign_epsilon, id2idx);
            // 边界 inflow：仍以“总压力驱动通量”判定，与你现有入口组分控制一致
            const int N = F.neighborCell;
            const bool boundaryInflow = (F_press < -cfg.flux_sign_epsilon) && (N < 0);
            int upW = upPress;
            int upG = upPress;

            // 6.1 B：相势能迎风（仅内部面可判）
            if (cfg.use_phase_potential_upwind && !F.isBoundary() && F.ownerCell >= 0 && F.neighborCell >= 0)
            {
                const int Pid = F.ownerCell;
                const int Nid = F.neighborCell;
                const int iP = static_cast<int>(id2idx.at(Pid));
                const int iN = static_cast<int>(id2idx.at(Nid));
                // water potential uses pw
                const double pwP = (*pw)[iP];
                const double pwN = (*pw)[iN];
                const double rwP = (*rho_w)[iP];
                const double rwN = (*rho_w)[iN];
                const int upW_try = detail::selectUpCellByPhasePotential(
                    F, mesh, id2idx, pwP, pwN, rwP, rwN,
                    cfg.gravity, cfg.min_rho_for_potential);

                // gas potential: pg (preferred) or pw (+pc if enabled) or pw (if not enabled)
                double pgP = pwP, pgN = pwN;
                if (pg)
                {
                    pgP = (*pg)[iP];
                    pgN = (*pg)[iN];
                }
                else {
                    pgP = pwP ;
                    pgN = pwN ;
                }
				const double rgP = (*rho_g)[iP];
                const double rgN = (*rho_g)[iN];
                const int upG_try = detail::selectUpCellByPhasePotential(
                    F, mesh, id2idx, pgP, pgN, rgP, rgN,
                    cfg.gravity, cfg.min_rho_for_potential);

                if (upW_try >= 0) upW = upW_try;
                if (upG_try >= 0) upG = upG_try;
            }
            // 6.2 A/B：取用于计算 fw_mass 的 (λ,ρ)
           // A: upW=upG=upPress
           // B: upW/upG 各自独立
            double lamW = 0.0, lamG = 0.0;
            double rhoW = 1.0, rhoG = 1.0;
            if (upW >= 0)
            {
				lamW = std::max((*lambda_w)[upW], cfg.min_lambda);
                rhoW = std::max((*rho_w)[upW], 0.0);
            }
            if (upG >= 0)
            {
                lamG = std::max((*lambda_g)[upG], 0.0);
                rhoG = std::max((*rho_g)[upG], 0.0);
            }
            // 质量型分相函数（A/B 差异核心之一：λg,ρg 可来自不同上游）
            double fw_mass = 0.0;
            {
                const double lamW_eff = lamW * rhoW;
                const double lamG_eff = lamG * rhoG;
                const double denom = std::max(lamW_eff + lamG_eff, cfg.min_lambda);
                fw_mass = clampValue(lamW_eff / denom, 0.0, 1.0);
            }
            // 边界 inflow：覆盖 fw（保持你现有入口控制）
            if (boundaryInflow && cfg.enforce_boundary_inflow_fw)
            {
                fw_mass = clampValue(cfg.boundary_inflow_fw, 0.0, 1.0);  //可以同给定cfg.boundary_inflow_fw的值来决定从边界上流入的是什么
            }
            // 6.3 压力驱动部分拆分
            const double Fw_press = fw_mass * F_press;
            const double Fg_press = F_press - Fw_press;
            // 6.4 重力部分：按 λ ρ² 权重拆分（A/B 差异之二：权重用各自上游）
            double F_grav_total = 0.0;
            if (grav_corr) F_grav_total = (*grav_corr)[iF];
            double Fw_grav = 0.0;
            double Fg_grav = 0.0;
            if (std::abs(F_grav_total) > 0.0)
            {
                const double w_w = lamW * rhoW * rhoW;
                const double w_g = lamG * rhoG * rhoG;
                const double denom_wg = std::max(w_w + w_g, cfg.min_lambda);
                Fw_grav = F_grav_total * (w_w / denom_wg);
                Fg_grav = F_grav_total - Fw_grav;
            }
            // 6.5 毛细部分：全部归气相（与你当前 FaceMassRate/Pressure 设计一致）
            double F_cap_total = 0.0;
            if (cap_corr) F_cap_total = (*cap_corr)[iF];
            const double Fw_cap = 0.0;
            const double Fg_cap = F_cap_total;
            // 6.6 最终两相质量通量
            const double Fw = Fw_press + Fw_grav + Fw_cap;
            const double Fg = Fg_press + Fg_grav + Fg_cap;
            (*mf_w)[iF] = Fw;
            (*mf_g)[iF] = Fg;
            if (fw_face) (*fw_face)[iF] = fw_mass;
		} // for faces
          // ===== 7) 可选：用 F_press 更新 FaceSignMask（保持你现有逻辑） ===== //
        FaceSignUpdateInfo info{};
        std::vector<int> flippedTmp;
        if (mask)
        {
            auto& flipped = (result ? result->flippedFaces : flippedTmp);
            info = updateFaceSignMask_fromFlux(*mf_press, cfg.flux_sign_epsilon, *mask, flipped);
        }

        if (result)
        {
            result->mf_w = mf_w;
            result->mf_g = mf_g;
            result->fw_face = fw_face;
            result->signInfo = info;
            if (!mask) result->flippedFaces.clear();
        }

        // ===== 8) 闭合检查 +（可选）cell net flux 诊断 ===== //
        detail::diagnoseFluxSplitClosure(
            *mf_w, *mf_g, *mf_press,
            grav_corr.get(), cap_corr.get(),
            cfg.use_phase_potential_upwind ? "B_phasePotentialUpwind" : "A_fluxPressUpwind");
        return true;
    }
```


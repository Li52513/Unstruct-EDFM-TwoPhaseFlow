#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"

// ===================== Timeterm =====================
namespace FVM {
    namespace Timeterm {

        template<class ABC>
        inline std::vector<char>
            mark_strong_BC_cells(MeshManager& mgr, const ABC& bc_adapter)
        {
            auto& mesh = mgr.mesh();
            const auto& faces = mesh.getFaces();
            const auto& id2idx = mesh.getCellId2Index();

            std::vector<char> mask(mesh.getCells().size(), 0);

            for (const auto& F : faces) {
                if (!F.isBoundary()) continue;

                double a = 0, b = 0, c = 0;
                if (!bc_adapter.getABC(F.id, a, b, c)) continue;
                if (std::abs(a) <= 1e-30) continue;      // 不是 Dirichlet/Robin，就不“钉”

                const int P = F.ownerCell;
                if (P >= 0) mask[id2idx.at(P)] = 1;      // 只需置 1，不要叠加
            }
            return mask;
        }


        // 全隐时间项：能量存储 C_eff * T
        // a_C = (V/dt)*C_eff^*，b_C = (V/dt)*C_eff^* * T_old
        // strongBCmask（可选）：对 mask[i]!=0 的“强边界单元”（a!=0）不写入时间项
        inline bool TimeTerm_FullyImplicit_SinglePhase_Temperature(
            MeshManager& mgr,
            FieldRegistry& reg,
            double dt,
            const std::string& Ceff_name,     // C_eff
            const std::string& T_old_name,    // T^n
            const std::string& aC_name,       // 输出：对角
            const std::string& bC_name
        ) {
            if (dt <= 0.0) { std::cerr << "[TimeTerm_T] invalid dt.\n"; return false; }

            auto& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const auto& id2idx = mesh.getCellId2Index();

            auto Ceff = reg.get<volScalarField>(Ceff_name.c_str());
            auto Told = reg.get<volScalarField>(T_old_name.c_str());
            if (!Ceff || !Told) {
                std::cerr << "[TimeTerm_T] missing fields: " << Ceff_name << " or " << T_old_name << "\n";
                return false;
            }

            auto aC = reg.getOrCreate<volScalarField>(aC_name.c_str(), cells.size(), 0.0);
            auto bC = reg.getOrCreate<volScalarField>(bC_name.c_str(), cells.size(), 0.0);
            std::fill(aC->data.begin(), aC->data.end(), 0.0);
            std::fill(bC->data.begin(), bC->data.end(), 0.0);

            const double inv_dt = 1.0 / dt;
            const double epsC = 0.0;

            for (const auto& c : cells) {
                if (c.id < 0) continue;
                const size_t i = id2idx.at(c.id);

                const double V = std::max(0.0, c.volume);
                const double Ceff_star = std::max(epsC, (*Ceff)[i]);
                const double T_old = (*Told)[i];

                // 常规时间项
                double a = V * inv_dt * Ceff_star;
                double b = V * inv_dt * Ceff_star * T_old;

                (*aC)[i] = a;
                (*bC)[i] = b;
            }
            return true;
        }

        inline bool TimeTerm_FullyImplicit_SinglePhase_Temperature_revised(
            MeshManager& mgr,
            FieldRegistry& reg,
            double dt,
            const std::string& Ceff_star_name,   // C_eff^*  (current iterate)
            const std::string& Ceff_old_name,    // C_eff^n  (frozen at step start)
            const std::string& T_old_name,       // T^n      (frozen at step start)
            const std::string& aC_name,          // output diagonal
            const std::string& bC_name           // output RHS constant
        ) {
            if (dt <= 0.0 || !std::isfinite(dt)) {
                std::cerr << "[TimeTerm_T] invalid dt.\n";
                return false;
            }

            auto& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const auto& id2idx = mesh.getCellId2Index();

            auto Cstar = reg.get<volScalarField>(Ceff_star_name.c_str());
            auto Cold = reg.get<volScalarField>(Ceff_old_name.c_str());
            auto Told = reg.get<volScalarField>(T_old_name.c_str());
            if (!Cstar || !Cold || !Told) {
                std::cerr << "[TimeTerm_T] missing fields: "
                    << Ceff_star_name << " or " << Ceff_old_name
                    << " or " << T_old_name << "\n";
                return false;
            }

            auto aC = reg.getOrCreate<volScalarField>(aC_name.c_str(), cells.size(), 0.0);
            auto bC = reg.getOrCreate<volScalarField>(bC_name.c_str(), cells.size(), 0.0);
            std::fill(aC->data.begin(), aC->data.end(), 0.0);
            std::fill(bC->data.begin(), bC->data.end(), 0.0);

            const double inv_dt = 1.0 / dt;
            const double epsC = 0.0; // 可设为小正数避免极端病态

            for (const auto& c : cells) {
                if (c.id < 0) continue;
                const size_t i = id2idx.at(c.id);

                const double V = std::max(0.0, c.volume);
                const double C_s = std::max(epsC, (*Cstar)[i]); // C_eff^*
                const double C_n = std::max(epsC, (*Cold)[i]); // C_eff^n
                const double T_n = (*Told)[i];                  // T^n

                // Backward Euler + Picard:
                // (V/dt) * ( C^* T^{n+1} - C^n T^n )
                (*aC)[i] = V * inv_dt * C_s;
                (*bC)[i] = V * inv_dt * C_n * T_n;
            }
            return true;
        }


        inline bool TimeTerm_FullyImplicit_SinglePhase_Flow(
            MeshManager& mgr,
            FieldRegistry& reg,
            double dt,
            const std::string& c_phi_name,          // 常数孔隙度压缩性(1/Pa)，无则传 0（体场）
            // ——字段名（全部已存在的体场）——
            const std::string& phi_name,            // φ^n
            const std::string& p_old_name,          // p^n
            const std::string& rho_old_name,        // ρ^n
            const std::string& p_lin_name,          // p^⋆
            const std::string& rho_lin_name,        // ρ^⋆
            const std::string& drdp_lin_name,       // (∂ρ/∂p)^⋆
            // ——输出——
            const std::string& aC_name,             // 时间项对角
            const std::string& bC_name            // 时间项常数项
        )
        {
            auto& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const auto& id2idx = mesh.getCellId2Index();

            auto phi = reg.get<volScalarField>(phi_name);
            auto p_n = reg.get<volScalarField>(p_old_name);
            auto rho_n = reg.get<volScalarField>(rho_old_name);
            auto p_lin = reg.get<volScalarField>(p_lin_name);
            auto rho_lin = reg.get<volScalarField>(rho_lin_name);
            auto drdp_lin = reg.get<volScalarField>(drdp_lin_name);
            auto c_phi = reg.get<volScalarField>(c_phi_name);

            if (!phi || !p_n || !rho_n || !p_lin || !rho_lin || !drdp_lin || !c_phi) return false;

            auto aC = reg.getOrCreate<volScalarField>(aC_name, cells.size(), 0.0);
            auto bC = reg.getOrCreate<volScalarField>(bC_name, cells.size(), 0.0);
            std::fill(aC->data.begin(), aC->data.end(), 0.0);
            std::fill(bC->data.begin(), bC->data.end(), 0.0);

            for (const auto& c : cells) {
                if (c.id < 0) continue;
                const size_t i = id2idx.at(c.id);
                const double V = std::max(0.0, c.volume);
                const double ph = (*phi)[i];
                const double pn = (*p_n)[i];
                const double rhoN = std::max(0.0, (*rho_n)[i]);
                const double pStar = (*p_lin)[i];
                const double rhoStar = std::max(0.0, (*rho_lin)[i]);
                const double drdp = std::max(0.0, (*drdp_lin)[i]);
                const double cphi = (*c_phi)[i];

                const double a = (V / dt) * ph * (drdp + rhoStar * cphi);
                const double b = (V / dt) * (ph * rhoN
                    - ph * rhoStar
                    + ph * drdp * pStar
                    + ph * rhoStar * cphi * pStar);

                (*aC)[i] = a;
                (*bC)[i] = b;
            }
            return true;
        }

       
    } // namespace Timeterm
} // namespace FVM

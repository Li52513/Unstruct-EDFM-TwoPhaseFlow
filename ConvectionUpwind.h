#pragma once
#include <string>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "FieldAcessForDiscre.h"
#include "ConvectionUpwind_Flux.h"
#include "Solver_AssemblerCOO.h"


template<class ABC>
inline bool getDirichletFromABC(const ABC& bc, int faceId, double& phi_b) {
    double a = 0, b = 0, c = 0;
    if (!bc.getABC(faceId, a, b, c)) return false;
    if (std::abs(a) <= 1e-30)      return false; // Neumann/通量型条件不给定值
    phi_b = c / a; // Dirichlet 或 Robin 返回对应的边界值
    return true;
}

namespace FVM {
    namespace Convection
    {

        /**
         * @brief 一阶迎风离散系数装配（统一用法）
         *
         * 面通量 flux 可以是 "Qf" 或 "mf"，通过“迎风选取”来给出系数 W_up。
         * 面系数 β_f = flux_f * W_up；内点不写 RHS，边界如果是 Dirichlet 则写 RHS。
         *
         * @tparam PhiABC  提供 getABC(faceId,a,b,c) 的边界条件适配器，如 TemperatureBCAdapter / PressureBCAdapter / ...
         * @param var_key  被离散变量的标识，用于调试或查找
         * @param flux_face_name 面通量面场名，通常是 "Qf" 或 "mf"
         * @param upwind_scalar_fields 迎风携带的标量列表，如 {"cp_g"} 或 {"rho_g","cp_g"}
         * @param nm aPP_conv / aPN_conv / bP_conv 的字段名集合
         * @param phiBC 边界条件适配器（如 TemperatureBCAdapter）
         * @param scalarBC 可选：为某些标量提供边界值的回调
         * @param wprod_face_debug 可选：调试输出的字段名
         */

        template<class PhiABC>
        inline bool build_FaceCoeffs_Upwind
        (
            MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
            const std::string& var_key,                   // 例如 "T"
            const std::string& flux_face_name,            // "mf_g"（质量通量）或 "Qf_g"（体积通量）
            const std::vector<std::string>& upwind_scalar_fields, // 随流物性，如 {"cp_g"} 或 {"rho_g","cp_g"}
            const OperatorFieldNames& nm,                 // aPP_conv / aPN_conv / bP_conv
            const PhiABC& phiBC,
            const std::unordered_map<std::string, std::function<bool(int, double&)>>& scalarBC = {},
            const std::string& wprod_face_debug = ""
        )
        {
            Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
            auto& faces = const_cast<std::vector<Face>&>(mesh.getFaces());
            const auto& id2idx = mesh.getCellId2Index();

            // 面通量
            auto fluxF = freg.get<faceScalarField>(flux_face_name.c_str());
            if (!fluxF) { std::cerr << "[Upwind] Missing face flux: " << flux_face_name << "\n"; return false; }

            // 使用 mf（质量通量）时不要把 rho_* 也放进 upwind 列表，避免 double-count
            const bool usingMf = (flux_face_name.find("mf") != std::string::npos);
            for (const auto& nmChi : upwind_scalar_fields) {
                if (usingMf && nmChi.find("rho") != std::string::npos) {
                    std::cerr << "[Upwind] WARNING: '" << flux_face_name
                        << "' looks like MASS flux; remove '" << nmChi
                        << "' from upwind list to avoid ρ× twice.\n";
                }
            }

            // 随流物性体场
            std::vector<std::shared_ptr<volScalarField>> chi_ptrs;
            chi_ptrs.reserve(upwind_scalar_fields.size());
            for (const auto& nmChi : upwind_scalar_fields) {
                auto f = reg.get<volScalarField>(nmChi.c_str());
                if (!f) { std::cerr << "[Upwind] Missing vol field: " << nmChi << "\n"; return false; }
                chi_ptrs.emplace_back(f);
            }

            // 输出面系数
            auto aPP = freg.getOrCreate<faceScalarField>(nm.aPP_conv.c_str(), faces.size(), 0.0);
            auto aPN = freg.getOrCreate<faceScalarField>(nm.aPN_conv.c_str(), faces.size(), 0.0);
            auto bP = freg.getOrCreate<faceScalarField>(nm.bP_conv.c_str(), faces.size(), 0.0);
            std::fill(aPP->data.begin(), aPP->data.end(), 0.0);
            std::fill(aPN->data.begin(), aPN->data.end(), 0.0);
            std::fill(bP->data.begin(), bP->data.end(), 0.0);

            // 可选：导出上风物性乘积（调试）
            std::shared_ptr<faceScalarField> Wface;
            if (!wprod_face_debug.empty()) {
                Wface = freg.getOrCreate<faceScalarField>(wprod_face_debug.c_str(), faces.size(), 0.0);
                std::fill(Wface->data.begin(), Wface->data.end(), 0.0);
            }

            auto product_at_cell = [&](size_t i)->double {
                double w = 1.0;
                for (const auto& chi : chi_ptrs) w *= (*chi)[i];
                return w;
                };

            const double epsF = 1e-18; // 小通量裁剪

            for (const auto& F : faces)
            {
                const int iF = F.id - 1;
                const double flux = (*fluxF)[iF]; // 约定：flux>0 表 owner→neighbor（owner 出流）
                if (std::abs(flux) <= epsF) continue;

                const int P = F.ownerCell;
                const size_t iP = id2idx.at(P);

                if (!F.isBoundary())
                {
                    const int N = F.neighborCell;
                    const size_t iN = id2idx.at(N);
                    const size_t iUP = (flux >= 0.0) ? iP : iN;  // 出流→取 owner；入流→取 neighbor
                    const double Wup = product_at_cell(iUP);
                    if (Wface) (*Wface)[iF] = Wup;

                    const double beta = flux * Wup;
                    if (beta >= 0.0) {
                        (*aPP)[iF] = beta;   // P 出流
                        (*aPN)[iF] = 0.0;
                    }
                    else {
                        (*aPP)[iF] = 0.0;
                        (*aPN)[iF] = beta;   // 负值，装配器会给 N 对角镜像
                    }
                }
                else
                {
                    // —— 边界 —— //
                    if (flux > 0.0) {
                        // 出流：耗散项（把热量带出域） → 对角
                        const double WP = product_at_cell(iP);
                        if (Wface) (*Wface)[iF] = WP;
                        (*aPP)[iF] += flux * WP;
                    }
                    else {
                        // 入流：RHS（边界给定面值），不加对角
                        double Wb = 1.0;
                        for (size_t k = 0; k < upwind_scalar_fields.size(); ++k) {
                            const std::string& nmChi = upwind_scalar_fields[k];
                            auto it = scalarBC.find(nmChi);
                            double chi_b = 0.0;
                            if (it != scalarBC.end() && it->second(F.id, chi_b)) Wb *= chi_b;
                            else Wb *= (*chi_ptrs[k])[iP]; // 无物性BC → 用 owner 近似
                        }
                        if (Wface) (*Wface)[iF] = Wb;

                        double phi_b = 0.0;
                        if (getDirichletFromABC(phiBC, F.id, phi_b)) {(*bP)[iF] += (-flux) * Wb * phi_b; // 注意：-flux>0
                        }
                    }
                }
            }
            return true;
        }
    }
}

#pragma once
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "FieldAcessForDiscre.h"
#include "BCAdapter.h"
#include "Diff_TPFA_PermeabilityOperation.h"   // kEffAlong(ehat,kxx,kyy,kzz)
#include "Diff_TPFA_GradientsOperation.h"      // computeCellGradients_LSQ_with_GG

namespace FVM {
    namespace Diffusion {

        enum class RhoFaceMethod { Linear, Harmonic };

        // ========== 解析 tokens==========
        struct ParsedMobility {
            std::string iso; std::string kxx, kyy, kzz;
            std::vector<std::string> mul, div;
            std::string rho; // 若给出 → 开启“质量形”：把 ρ_f 乘到 a/s
        };

        inline ParsedMobility parseMobilityTokens(const std::vector<std::string>& tokens) {
            ParsedMobility pm;
            auto trim = [](std::string s)->std::string {
                while (!s.empty() && isspace((unsigned char)s.back())) s.pop_back();
                size_t i = 0; while (i < s.size() && isspace((unsigned char)s[i])) ++i;
                return s.substr(i);
                };
            auto starts = [&](const std::string& s, const char* p) { return s.rfind(p, 0) == 0; };

            for (auto t0 : tokens) {
                auto t = trim(t0);
                if (t.empty()) continue;

                if (starts(t, "iso:")) { pm.iso = t.substr(4); continue; }
                if (starts(t, "kxx:")) { pm.kxx = t.substr(4); continue; }
                if (starts(t, "kyy:")) { pm.kyy = t.substr(4); continue; }
                if (starts(t, "kzz:")) { pm.kzz = t.substr(4); continue; }
                if (starts(t, "rho:")) { pm.rho = t.substr(4); continue; }

                if (t[0] == '/') { pm.div.emplace_back(t.substr(1)); continue; }

                auto posPow = t.find("^-1");
                if (posPow != std::string::npos) { pm.div.emplace_back(t.substr(0, posPow)); continue; }

                pm.mul.emplace_back(t); // 默认乘上
            }
            return pm;
        }


        /**
         * @brief 字段名驱动的 TPFA 中心差分扩散系数生成（无 Provider）。
         *
         * mobility_tokens 示例：
         *  - 各向同性导热：{"iso:lambda_eff"}
         *  - 单相达西：{"kxx:kxx","kyy:kyy","kzz:kzz","/mu_g"}  （λ = k_ê / μ_g）
         *  - 两相（逐相调用累加）：{"kxx:kxx","kyy:kyy","kzz:kzz","kr_w","/mu_w"}
         *
         * @param a_vol   输出面场名（体积形 a_f）
         * @param s_vol   输出面场名（体积形 s_f，含 cross/buoy/BC）
         * @param x_name  用于梯度的标量体场名（如 "p_g" 或 "T"）
         * @param mobility_tokens 见上
         * @param rho_field  若非空且 enable_buoy=true，则参与浮力 s_buoy；为空则不算浮力
         * @param rho_method 面密度插值方法
         * @param g, bc, enable_buoy, gradSmoothIters 与你现有版本等价
         */
        template<class BCProvider>
        inline bool build_FaceCoeffs_Central(
            MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
            const std::string& a_name,                 // 输出：面系数 a_f
            const std::string& s_name,                 // 输出：面源项 s_f
            const std::string& x_name,                 // 自变量体场（压力 p 或 温度 T）
            const std::vector<std::string>& mobility_tokens, // 组合 λ 的 tokens（支持 rho: 开启“质量形”）
            const std::string& rho_field_for_buoy,     // 浮力线性化用的密度体场（若 tokens 未给 rho:，仍可用于浮力）
            RhoFaceMethod rho_method,                  // 面密度插值：Linear/Harmonic
            const Vector& g,                           // 重力加速度向量
            const BCProvider& bc,                      // getABC(face,a,b,c)
            bool enable_buoy = true,                   // 是否加入浮力项
            int  gradSmoothIters = 0                   // 单元梯度平滑迭代
        )
        {
            Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
            auto& faces = const_cast<std::vector<Face>&>(mesh.getFaces());
            const auto& id2idx = mesh.getCellId2Index();
            const auto& cells = mesh.getCells();

            auto aF = freg.getOrCreate<faceScalarField>(a_name.c_str(), faces.size(), 0.0);
            auto sF = freg.getOrCreate<faceScalarField>(s_name.c_str(), faces.size(), 0.0);

            // —— 解析 tokens —— //
            const ParsedMobility pm = parseMobilityTokens(mobility_tokens);
            const bool MASS_FORM = !pm.rho.empty();                     // 给了 rho: → 质量形
            const std::string rho_used = MASS_FORM ? pm.rho : rho_field_for_buoy;

            // —— 准备体场指针 —— //
            auto isoF = pm.iso.empty() ? nullptr : reg.get<volScalarField>(pm.iso);
            auto kxxF = pm.kxx.empty() ? nullptr : reg.get<volScalarField>(pm.kxx);
            auto kyyF = pm.kyy.empty() ? nullptr : reg.get<volScalarField>(pm.kyy);
            auto kzzF = pm.kzz.empty() ? nullptr : reg.get<volScalarField>(pm.kzz);
            auto rhoF = rho_used.empty() ? nullptr : reg.get<volScalarField>(rho_used);

            std::vector<const volScalarField*> mulFs; mulFs.reserve(pm.mul.size());
            for (const auto& nm : pm.mul) {
                auto f = reg.get<volScalarField>(nm);
                if (f) mulFs.push_back(f.get());
                else   std::cerr << "[Diffusion] warn: mul field '" << nm << "' not found, ignored.\n";
            }
            std::vector<const volScalarField*> divFs; divFs.reserve(pm.div.size());
            for (const auto& nm : pm.div) {
                auto f = reg.get<volScalarField>(nm);
                if (f) divFs.push_back(f.get());
                else   std::cerr << "[Diffusion] warn: div field '" << nm << "' not found, ignored.\n";
            }

            // 自变量梯度（p 或 T），供非正交修正
            const std::vector<Vector> grad = computeCellGradients_LSQ_with_GG(mesh, reg, x_name.c_str(), gradSmoothIters);

            // —— 小工具：按 ehat 投影得到 λ(cell,ehat) —— //
            const double eps_l = 1e-30, eps_d = 1e-14, eps_A = 1e-30;
            auto lamAlong = [&](int cellId, const Vector& ehat)->double {
                const size_t i = id2idx.at(cellId);
                double lam = 0.0;

                if (isoF) {
                    lam = (*isoF)[i];
                }
                else {
                    const double kx = kxxF ? (*kxxF)[i] : 0.0;
                    const double ky = kyyF ? (*kyyF)[i] : 0.0;
                    const double kz = kzzF ? (*kzzF)[i] : 0.0;
                    const double ex = ehat.m_x, ey = ehat.m_y, ez = ehat.m_z;
                    lam = kx * ex * ex + ky * ey * ey + kz * ez * ez;  // 对角各向异性；若需全张量可在此扩展
                }
                for (const auto* f : mulFs) lam *= (*f)[i];
                for (const auto* f : divFs) lam /= std::max((*f)[i], eps_l);
                return std::max(lam, eps_l);
                };

            // 面密度插值（Linear/Harmonic）
            auto rho_face = [&](double gamma, int P, int N)->double {
                if (!rhoF) return 1.0;
                const double rP = (*rhoF)[id2idx.at(P)];
                if (N < 0) return rP;  // 边界
                const double rN = (*rhoF)[id2idx.at(N)];
                if (rho_method == RhoFaceMethod::Harmonic) {
                    const double d = std::max(gamma / rP + (1.0 - gamma) / rN, 1e-30);
                    return 1.0 / d;
                }
                return (1.0 - gamma) * rP + gamma * rN; // 线性插值
                };

            // —— 主循环：每个面构造 a_f, s_f —— //
            for (const auto& F : faces) {
                const int iF = F.id - 1;

                const Vector A = F.vectorE + F.vectorT;
                const double Aabs = A.Mag();
                const double Eabs = F.vectorE.Mag();

                const int P = F.ownerCell;  if (P < 0) { (*aF)[iF] = 0; (*sF)[iF] = 0; continue; }
                const int N = F.neighborCell;

                // ehat 与 gamma（若几何未给 gamma，按投影求）
                double dON = std::max(F.ownerToNeighbor.Mag(), eps_d);
                Vector ehat = (dON > eps_d ? F.ownerToNeighbor / dON : F.normal);
                double gamma = F.f_linearInterpolationCoef;
                if (!(gamma > 0.0 && gamma < 1.0) && N >= 0) {
                    const Vector eON = F.ownerToNeighbor / dON;
                    const double s = (F.midpoint - cells[id2idx.at(P)].center) * eON;
                    gamma = std::max(0.0, std::min(1.0, s / dON));
                }
                else if (N < 0) {
                    gamma = 0.0; // 边界退化：仅 P 侧
                }

                // 体积形导通率
                const double lamP = lamAlong(P, ehat);
                double a_vol = 0.0, s_cross = 0.0, s_buoy = 0.0;

                if (N >= 0) {
                    const double lamN = lamAlong(N, ehat);
                    const double lam_f = 1.0 / std::max(gamma / lamP + (1.0 - gamma) / lamN, eps_l);

                    a_vol = lam_f * Eabs / dON;

                    // 非正交修正（体积形）
                    const Vector& gP = grad[id2idx.at(P)];
                    const Vector& gN = grad[id2idx.at(N)];
                    const Vector  gf = (1.0 - gamma) * gP + gamma * gN;
                    s_cross = lam_f * (gf * F.vectorT);

                    // 浮力（体积形，线性 ρ）
                    if (enable_buoy && rhoF) {
                        const double rho_lin = rho_face(gamma, P, N);
                        s_buoy = -lam_f * rho_lin * (g * A);
                    }
                }
                else {
                    // ===== 边界：P 侧闭合 + Robin 统一 =====
                    const double a_face0 = lamP * Eabs / dON;  // ★ 原始面对角（未缩减）
                    a_vol = a_face0;                           // 暂存，后面再做缩减

                    // 交叉项：仅用 P 侧梯度（未乘ρ；质量形统一在末尾处理）
                    const Vector& gP = grad[id2idx.at(P)];
                    s_cross = lamP * (gP * F.vectorT);

                    // 浮力项：体积形（未乘ρ；质量形统一在末尾处理）
                    if (enable_buoy && rhoF) {
                        const double rhoP = (*rhoF)[id2idx.at(P)];
                        s_buoy = -lamP * rhoP * (g * A);
                    }

                    // —— Robin（a,b,c 面积化）：对角缩减用 (1-α)，BC 源项用“未缩减的 a_face0 × β_B”
                    double a = 0.0, b = 0.0, c = 0.0;
                    double s_BC_unmass = 0.0;                  // ★ 先不乘ρ
                    if (bc.getABC(F.id, a, b, c)) {
                        const double den = a * std::max(Aabs, eps_A) + b * (Eabs / dON);
                        const double alpha_f = (den > eps_A) ? (b * (Eabs / dON) / den) : 0.0;
                        const double betaB_f = (std::abs(a) > eps_A && den > eps_A) ? ((c * Aabs) / den) : 0.0;

                        a_vol = a_face0 * (1.0 - alpha_f);     // ★ 只缩减对角
                        s_BC_unmass = a_face0 * betaB_f;       // ★ BC 注入用“未缩减的 a_face0”
                    }

                    // 把未乘ρ的 BC 注入暂存到 s_cross 里，后面统一乘ρ并做限幅拆分
                    s_cross += s_BC_unmass;                    // ★ 注意：稍后会把“BC部分”从限幅中剥离
                }

                // —— 质量形：把 ρ_f 乘进去 —— //
                // —— 将“BC 注入”与“交叉+浮力”拆分，限幅只作用后者 —— //
                // 这里把 s_cross 里可能混入的“BC 注入”拆出来：方法是重新计算 s_BC_unmass
                double s_BC_unmass = 0.0;
                if (N < 0) {
                    // 仅边界面才可能有 BC 注入；复用上面的 a,b,c 计算逻辑
                    double a = 0.0, b = 0.0, c = 0.0;
                    if (bc.getABC(F.id, a, b, c)) {
                        const double den = a * std::max(Aabs, eps_A) + b * (Eabs / dON);
                        const double betaB_f = (std::abs(a) > eps_A && den > eps_A) ? ((c * Aabs) / den) : 0.0;
                        const double a_face0 = lamP * Eabs / dON;
                        s_BC_unmass = a_face0 * betaB_f;   // 未缩减、未乘ρ
                    }
                }
                // s_nb_unmass = 非正交交叉 + 浮力（未乘ρ），BC 从中减掉，单独保留
                const double s_nb_unmass = (s_cross + s_buoy) - s_BC_unmass;

                // —— 质量形（若给了 rho: token）统一乘 ρ_f —— //
                const double rho_f = MASS_FORM ? rho_face(gamma, P, N) : 1.0;
                const double a_face = rho_f * a_vol;
                const double s_nb = rho_f * s_nb_unmass;     // 仅“交叉+浮力”
                const double s_BC = rho_f * s_BC_unmass;     // “BC 注入”不参与限幅

                // —— 限幅：只限“交叉+浮力”部分，幅度±0.5 a_face —— //
                const double s_nb_limited = std::max(-0.5 * a_face, std::min(0.5 * a_face, s_nb));
                const double s_total = s_nb_limited + s_BC;

                // —— 写回 —— //
                (*aF)[iF] = a_face;
                (*sF)[iF] = s_total;
            }

            return true;
        }


    }
}



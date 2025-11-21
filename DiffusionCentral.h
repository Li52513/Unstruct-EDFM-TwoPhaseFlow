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
#include "Diff_TPFA_PermeabilityOperation.h"
#include "Diff_TPFA_GradientsOperation.h"

namespace FVM 
{
    namespace Diffusion 
    {

        enum class RhoFaceMethod { Linear, Harmonic };

        // --- new: 工具函数 ----------------------------------------------------------
        inline double clamp01(double v) { return std::max(0.0, std::min(1.0, v)); }

        inline double orthogonalityLimiter(const Vector& faceNormal,
            const Vector& dir,
            double minLimiter = 0.05,
            double power = 2.0)
        {
            const double denom = std::max(faceNormal.Mag() * dir.Mag(), 1e-14);
            if (denom <= 1e-14) return minLimiter;
            const double cosine = clamp01(std::abs(faceNormal * dir) / denom);
            const double limiter = std::pow(cosine, power);
            return std::max(minLimiter, limiter);
        }
        // ---------------------------------------------------------------------------

        struct ParsedMobility 
        {
            std::string iso;                    // 各向同性迁移率参数
            std::string kxx, kyy, kzz;          // 各向异性迁移率参数（x,y,z方向）
            std::vector<std::string> mul, div;  // 乘法（除法）因子列表
            std::string rho;                    // 密度参数
            std::string rho_buoy;               // 浮力项的密度参数
        };

        inline ParsedMobility parseMobilityTokens(const std::vector<std::string>& tokens)
        {
            ParsedMobility pm;                  // 创建解析结果对象
            auto trim = [](std::string s)->std::string
                {
                    // 去除尾部空白字符
                    while (!s.empty() && isspace((unsigned char)s.back())) s.pop_back();
                    // 去除头部空白字符
                    size_t i = 0; while (i < s.size() && isspace((unsigned char)s[i])) ++i;
                    return s.substr(i);
                };
            // 定义辅助lambda函数：检查字符串是否以指定前缀开头
            auto starts = [&](const std::string& s, const char* p) { return s.rfind(p, 0) == 0; };

            for (auto t0 : tokens)
            {
                auto t = trim(t0);
                if (t.empty()) continue;

                if (starts(t, "iso:")) { pm.iso = t.substr(4); continue; }
                if (starts(t, "kxx:")) { pm.kxx = t.substr(4); continue; }
                if (starts(t, "kyy:")) { pm.kyy = t.substr(4); continue; }
                if (starts(t, "kzz:")) { pm.kzz = t.substr(4); continue; }
                if (starts(t, "rho_buoy:")) { pm.rho_buoy = t.substr(9); continue; }
                if (starts(t, "rho_coeff:")) { pm.rho = t.substr(10); continue; }
                if (starts(t, "rho:")) { pm.rho = t.substr(4); continue; }

                if (t[0] == '/') { pm.div.emplace_back(t.substr(1)); continue; }

                auto posPow = t.find("^-1");
                if (posPow != std::string::npos) { pm.div.emplace_back(t.substr(0, posPow)); continue; }

                pm.mul.emplace_back(t);
            }
            return pm;
        }

        template<class BCProvider>
        inline bool build_FaceCoeffs_Central(
            MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
            const std::string& a_name,
            const std::string& s_name,
            const std::string& x_name,
            const std::vector<std::string>& mobility_tokens,
            const std::string& rho_field_for_buoy,
            RhoFaceMethod rho_method,
            const Vector& g,
            const BCProvider& bc,
            bool enable_buoy =false,
            int  gradSmoothIters = 0)
        {
            Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
            auto& faces = const_cast<std::vector<Face>&>(mesh.getFaces());
            const auto& id2idx = mesh.getCellId2Index();
            const auto& cells = mesh.getCells();

            auto aF = freg.getOrCreate<faceScalarField>(a_name.c_str(), faces.size(), 0.0);
            auto sF = freg.getOrCreate<faceScalarField>(s_name.c_str(), faces.size(), 0.0);

            ParsedMobility pm = parseMobilityTokens(mobility_tokens);

            auto isoF = (!pm.iso.empty() ? reg.get<volScalarField>(pm.iso.c_str()) : nullptr);
            auto kxxF = (!pm.kxx.empty() ? reg.get<volScalarField>(pm.kxx.c_str()) : nullptr);
            auto kyyF = (!pm.kyy.empty() ? reg.get<volScalarField>(pm.kyy.c_str()) : nullptr);
            auto kzzF = (!pm.kzz.empty() ? reg.get<volScalarField>(pm.kzz.c_str()) : nullptr);

            std::vector<std::shared_ptr<volScalarField>> mulKeep, divKeep;   // 持有者
            std::vector<const volScalarField*>           mulFs, divFs;      // 裸指针

            for (const auto& nm : pm.mul)
            {
                auto ptr = reg.get<volScalarField>(nm.c_str());
                if (ptr) {
                    mulKeep.push_back(ptr);
                    mulFs.push_back(ptr.get());
                }
            }
            for (const auto& nm : pm.div)
            {
                auto ptr = reg.get<volScalarField>(nm.c_str());
                if (ptr) {
                    divKeep.push_back(ptr);
                    divFs.push_back(ptr.get());
                }
            }

            auto rhoF = (!pm.rho.empty() ? reg.get<volScalarField>(pm.rho.c_str()) : nullptr);
            if (!rhoF && !rho_field_for_buoy.empty())
                rhoF = reg.get<volScalarField>(rho_field_for_buoy.c_str());

            auto grad = computeCellGradients_LSQ_with_GG(mesh, reg, x_name.c_str(), gradSmoothIters);

            const double epsA = 1e-30;
            const double eps_d = 1e-12;
            const double eps_l = 1e-20;

            auto lamAlong = [&](int cellId, const Vector& ehat)->double 
                {
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
                        lam = kx * ex * ex + ky * ey * ey + kz * ez * ez;
                    }
                    for (const auto* f : mulFs) lam *= (*f)[i];
                    for (const auto* f : divFs) lam /= std::max((*f)[i], eps_l);
                    return std::max(lam, eps_l);
                };

            auto rho_face = [&](double gamma, int P, int N)->double 
                {
                    if (!rhoF) return 1.0;
                    const double rP = (*rhoF)[id2idx.at(P)];
                    if (N < 0) return rP;
                    const double rN = (*rhoF)[id2idx.at(N)];
                    if (rho_method == RhoFaceMethod::Harmonic) {
                        const double d = std::max(gamma / rP + (1.0 - gamma) / rN, 1e-30);
                        return 1.0 / d;
                    }
                    return (1.0 - gamma) * rP + gamma * rN;
                };

            for (const auto& F : faces) 
            {
                const int iF = F.id - 1;

                const Vector A = F.vectorE + F.vectorT;
                const double Aabs = A.Mag();
                const double Eabs = F.vectorE.Mag();

                const int P = F.ownerCell;
                if (P < 0) {
                    (*aF)[iF] = 0.0;
                    (*sF)[iF] = 0.0;
                    continue;
                }
                const int N = F.neighborCell;

                double dON = std::max(F.ownerToNeighbor.Mag(), eps_d);
                Vector ehat = (dON > eps_d) ? (F.ownerToNeighbor / dON) : F.normal;

                double gamma = F.f_linearInterpolationCoef;
                if (!(gamma > 0.0 && gamma < 1.0) && N >= 0) {
                    const Vector eON = F.ownerToNeighbor / dON;
                    const double s = (F.midpoint - cells.at(id2idx.at(P)).center) * eON;
                    gamma = std::max(0.0, std::min(1.0, s / dON));
                }
                else if (N < 0) 
                {
                    gamma = 0.0;
                }

                const double orthoLimiter = orthogonalityLimiter(F.normal, ehat);
                const double lamP = lamAlong(P, ehat);

                double a_vol = 0.0;
                double s_cross = 0.0;
                double s_buoy = 0.0;
                double s_BC_unmass = 0.0;
                double s_total_vol = 0.0;

                if (N >= 0) 
                {
                    // —— 内部面 —— //
                    const double lamN = lamAlong(N, ehat);
                    const double lam_f = 1.0 / std::max(gamma / lamP + (1.0 - gamma) / lamN, eps_l);
                    a_vol = lam_f * Eabs / dON;

                    const Vector& gP = grad.at(id2idx.at(P));
                    const Vector& gN = grad.at(id2idx.at(N));
                    const Vector  gf = (1.0 - gamma) * gP + gamma * gN;
                    s_cross = lam_f * (gf * F.vectorT);

                    if (enable_buoy && rhoF) 
                    {
                        const double rho_lin = rho_face(gamma, P, N);
                        s_buoy = -lam_f * rho_lin * (g * A);
                    }

                    s_total_vol = s_cross + s_buoy;
                }
                else 
                {
                    // —— 边界面 —— //
                    const double a_face0 = lamP * Eabs / dON;

                    if (enable_buoy && rhoF) { const double rhoP = (*rhoF)[id2idx.at(P)]; 
                    s_buoy = -lamP * rhoP * (g * A); } 
                    
                    double a = 0.0, b = 0.0, c = 0.0;
                    if (bc.getABC(F.id, a, b, c))
                    {
                        s_cross = 0; //因为当前不处理边界面为曲面的情况，所以A == E，不进行非正交修正

                        double alpha = b * (Eabs / dON);
                        alpha = alpha / (a * Aabs + alpha);

                        double beta = c * Aabs;  //因为当前不处理边界面为曲面的情况，所以A == E，不进行非正交修正 所以 gf * F.vectorT =0
                        beta = beta / (a * Aabs + b * (Eabs / dON));

                        a_vol = (1 - alpha) * a_face0;
                        s_BC_unmass = a_face0 * beta;

                    }
                    s_total_vol = s_BC_unmass + s_buoy;

                }
                const double rho_f = rho_face(gamma, P, N);
                (*aF)[iF] = rho_f * a_vol;
                (*sF)[iF] = s_total_vol * rho_f;
              
        }

            return true;
        }

        // 
        struct TwoPhaseDiffusionTemplate
        {
            // 面系数模板：允许扩散系数与浮力修正使用不同的密度场
            template<class BCProvider>
            static bool build_FaceCoeffs_Central(
                MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
                const std::string& a_name,
                const std::string& s_name,
                const std::string& x_name,
                const std::vector<std::string>& mobility_tokens,
                const std::string& rho_coeff_fallback,
                const std::string& rho_buoy_fallback,
                RhoFaceMethod rho_method,
                const Vector& g,
                const BCProvider& bc,
                bool enable_buoy = false,
                int gradSmoothIters = 0)
            {
                Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
                auto& faces = const_cast<std::vector<Face>&>(mesh.getFaces());
                const auto& id2idx = mesh.getCellId2Index();
                const auto& cells = mesh.getCells();

                auto aF = freg.getOrCreate<faceScalarField>(a_name.c_str(), faces.size(), 0.0);
                auto sF = freg.getOrCreate<faceScalarField>(s_name.c_str(), faces.size(), 0.0);

                ParsedMobility pm = parseMobilityTokens(mobility_tokens);

                auto isoF = (!pm.iso.empty() ? reg.get<volScalarField>(pm.iso.c_str()) : nullptr);
                auto kxxF = (!pm.kxx.empty() ? reg.get<volScalarField>(pm.kxx.c_str()) : nullptr);
                auto kyyF = (!pm.kyy.empty() ? reg.get<volScalarField>(pm.kyy.c_str()) : nullptr);
                auto kzzF = (!pm.kzz.empty() ? reg.get<volScalarField>(pm.kzz.c_str()) : nullptr);

                std::vector<std::shared_ptr<volScalarField>> mulKeep, divKeep;
                std::vector<const volScalarField*>           mulFs, divFs;

                for (const auto& nm : pm.mul)
                {
                    auto ptr = reg.get<volScalarField>(nm.c_str());
                    if (ptr) {
                        mulKeep.push_back(ptr);
                        mulFs.push_back(ptr.get());
                    }
                }
                for (const auto& nm : pm.div)
                {
                    auto ptr = reg.get<volScalarField>(nm.c_str());
                    if (ptr) {
                        divKeep.push_back(ptr);
                        divFs.push_back(ptr.get());
                    }
                }

                auto rhoCoeffF = (!pm.rho.empty() ? reg.get<volScalarField>(pm.rho.c_str()) : nullptr);
                auto rhoBuoyF = (!pm.rho_buoy.empty() ? reg.get<volScalarField>(pm.rho_buoy.c_str()) : nullptr);
                if (!rhoCoeffF && !rho_coeff_fallback.empty())
                    rhoCoeffF = reg.get<volScalarField>(rho_coeff_fallback.c_str());
                if (!rhoBuoyF && !rho_buoy_fallback.empty())
                    rhoBuoyF = reg.get<volScalarField>(rho_buoy_fallback.c_str());
                if (!rhoBuoyF)
                    rhoBuoyF = rhoCoeffF;

                auto grad = computeCellGradients_LSQ_with_GG(mesh, reg, x_name.c_str(), gradSmoothIters);

                const double epsA = 1e-30;
                const double eps_d = 1e-12;
                const double eps_l = 1e-20;

                auto lamAlong = [&](int cellId, const Vector& ehat)->double
                    {
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
                            lam = kx * ex * ex + ky * ey * ey + kz * ez * ez;
                        }
                        for (const auto* f : mulFs) lam *= (*f)[i];
                        for (const auto* f : divFs) lam /= std::max((*f)[i], eps_l);
                        return std::max(lam, eps_l);
                    };

                auto rho_face = [&](double gamma, int P, int N, const volScalarField* fld)->double
                    {
                        if (!fld) return 1.0;
                        const double rP = (*fld)[id2idx.at(P)];
                        if (N < 0) return rP;
                        const double rN = (*fld)[id2idx.at(N)];
                        if (rho_method == RhoFaceMethod::Harmonic) {
                            const double denom = std::max(gamma / std::max(rP, eps_l) + (1.0 - gamma) / std::max(rN, eps_l), 1e-30);
                            return 1.0 / denom;
                        }
                        return (1.0 - gamma) * rP + gamma * rN;
                    };

                const volScalarField* rhoCoeffRaw = rhoCoeffF ? rhoCoeffF.get() : nullptr;
                const volScalarField* rhoBuoyRaw = rhoBuoyF ? rhoBuoyF.get() : nullptr;

                for (const auto& F : faces)
                {
                    const int iF = F.id - 1;

                    const Vector A = F.vectorE + F.vectorT;
                    const double Aabs = A.Mag();
                    const double Eabs = F.vectorE.Mag();

                    const int P = F.ownerCell;
                    if (P < 0) {
                        (*aF)[iF] = 0.0;
                        (*sF)[iF] = 0.0;
                        continue;
                    }
                    const int N = F.neighborCell;

                    double dON = std::max(F.ownerToNeighbor.Mag(), eps_d);
                    Vector ehat = (dON > eps_d) ? (F.ownerToNeighbor / dON) : F.normal;

                    double gamma = F.f_linearInterpolationCoef;
                    if (!(gamma > 0.0 && gamma < 1.0) && N >= 0) {
                        const Vector eON = F.ownerToNeighbor / dON;
                        const double s = (F.midpoint - cells.at(id2idx.at(P)).center) * eON;
                        gamma = std::max(0.0, std::min(1.0, s / dON));
                    }
                    else if (N < 0)
                    {
                        gamma = 0.0;
                    }

                    const double orthoLimiter = orthogonalityLimiter(F.normal, ehat);
                    const double lamP = lamAlong(P, ehat);

                    double a_vol = 0.0;
                    double s_cross = 0.0;
                    double s_buoy = 0.0;
                    double s_BC_unmass = 0.0;

                    if (N >= 0)
                    {
                        const double lamN = lamAlong(N, ehat);
                        const double lam_f = 1.0 / std::max(gamma / lamP + (1.0 - gamma) / lamN, eps_l);
                        a_vol = lam_f * Eabs / dON;

                        const Vector& gP = grad.at(id2idx.at(P));
                        const Vector& gN = grad.at(id2idx.at(N));
                        const Vector  gf = (1.0 - gamma) * gP + gamma * gN;
                        s_cross = lam_f * (gf * F.vectorT);

                        if (enable_buoy)
                        {
                            s_buoy = -lam_f * (g * A);
                        }
                    }
                    else
                    {
                        const double a_face0 = lamP * Eabs / dON;

                        if (enable_buoy)
                        {
                            s_buoy = -lamP * (g * A);
                        }

                        double a = 0.0, b = 0.0, c = 0.0;
                        if (bc.getABC(F.id, a, b, c))
                        {
                            s_cross = 0.0;

                            double alpha = b * (Eabs / dON);
                            alpha = alpha / (a * Aabs + alpha);

                            double beta = c * Aabs;
                            beta = beta / (a * Aabs + b * (Eabs / dON));

                            a_vol = (1 - alpha) * a_face0;
                            s_BC_unmass = a_face0 * beta;
                        }
                    }

                    const double rhoCoeff_face = rho_face(gamma, P, N, rhoCoeffRaw);
                    const double rhoBuoy_face = rho_face(gamma, P, N, rhoBuoyRaw ? rhoBuoyRaw : rhoCoeffRaw);

                    (*aF)[iF] = rhoCoeff_face * a_vol;
                    const double source = rhoCoeff_face * (s_cross + s_BC_unmass)
                        + rhoCoeff_face * rhoBuoy_face * s_buoy;
                    (*sF)[iF] = source;
                }

                return true;
            }
        };

    } // namespace Diffusion
} // namespace FVM

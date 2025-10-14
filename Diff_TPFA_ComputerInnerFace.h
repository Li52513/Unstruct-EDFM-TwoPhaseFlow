#pragma once
#include <string>
#include <memory>
#include <tuple>
#include <utility>
#include <algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "FieldAcessForDiscre.h"
#include "Diff_TPFA_PermeabilityOperation.h"
#include "Diff_TPFA_GradientsOperation.h"
#include "Diff_TPFA_UpwindforGravityandDensity.h"
#include "Diff_TPFA_Mobility.h"

// 内部面：TPFA（达西/纯扩散皆可）
template<class MobilityProvider, class DensityPolicy>
inline void Diffusion_TPFA_InnerFace_SinglePhase(
    MeshManager& mgr,
    const FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const GravUpwind& gu,
    const MobilityProvider& mob,
    const DensityPolicy& rhoPol,
    const std::string& a_name = "a_f_Diff",
    const std::string& s_name = "s_f_Diff",
    const std::string& x_name = "p",
    bool enable_buoy = true,
    int  gradSmoothIters = 0,
    // 新增：是否把 ρ_up 乘进矩阵系数（默认 true 以保持你当前行为）
    bool rho_in_matrix = true
) {
    Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
    auto& faces = const_cast<std::vector<Face>&>(mesh.getFaces());
    const auto& id2idx = mesh.getCellId2Index();

    auto a_f_Diff = freg.getOrCreate<faceScalarField>(a_name, faces.size(), 0.0);
    auto s_f_Diff = freg.getOrCreate<faceScalarField>(s_name, faces.size(), 0.0);

    // 预计算 ∇x（LSQ + GG 兜底，可选平滑）
    const std::vector<Vector> grad =
        computeCellGradients_LSQ_with_GG(mesh, reg, x_name.c_str(), gradSmoothIters);

    const double eps_d = 1e-14, eps_mu = 1e-30, eps_l = 1e-30;

#pragma omp parallel for schedule(static)
    for (int f = 0; f < static_cast<int>(faces.size()); ++f) {
        const Face& F = faces[f];
        const int idxF = F.id - 1;

        if (F.isBoundary()) {
            (*a_f_Diff)[idxF] = 0.0;
            (*s_f_Diff)[idxF] = 0.0;
            continue;
        }

        const int Pid = F.ownerCell;
        const int Nid = F.neighborCell;
        const int iP = id2idx.at(Pid);
        const int iN = id2idx.at(Nid);

        const double dperp = std::max(F.ownerToNeighbor.Mag(), eps_d);
        const Vector  ehat = (dperp > 1e-14 ? F.ownerToNeighbor / dperp : F.normal);
        const double  Eabs = F.vectorE.Mag();

        // 线性插值系数 gamma（修正：一律用 id2idx 取中心）
        double gamma = F.f_linearInterpolationCoef;
        if (!(gamma > 0.0 && gamma < 1.0)) {
            const Vector& CO = mesh.getCells()[iP].center;
            const Vector& CN = mesh.getCells()[iN].center;
            const Vector  dON = CN - CO;
            const double  D = dON.Mag();
            if (D > eps_d) {
                const Vector eON = dON / D;
                const double s = (F.midpoint - CO) * eON; // owner→面心投影
                gamma = std::min(1.0, std::max(0.0, s / D));
            }
            else {
                gamma = 0.5;
            }
        }

        // mobility 沿 ehat 的等效导通率（调和到面）
        const double lamP = std::max(mob.mobilityAlong(mesh, reg, Pid, ehat), eps_l);
        const double lamN = std::max(mob.mobilityAlong(mesh, reg, Nid, ehat), eps_l);
        const double lam_f = 1.0 / std::max(gamma / lamP + (1.0 - gamma) / lamN, eps_l);

        // 压力/温度（对达西类取 p，纯扩散类取 T）
        const double p_p = cellScalar(reg, mesh, x_name.c_str(), Pid, 0.0);
        const double p_n = cellScalar(reg, mesh, x_name.c_str(), Nid, 0.0);

        const Vector& CP = mesh.getCells()[iP].center;
        const Vector& CN = mesh.getCells()[iN].center;

        // ρ：交错平均用于浮力；迎风用于可选地缩放矩阵（达西-质量流速化）
        const double rho_inter = rhoPol.rhoBar(mesh, reg, Pid, Nid, gamma);
        const double rho_up = rhoPol.rhoUp(mesh, reg, Pid, Nid, p_p, p_n, CP, CN, gu);

        // β_f = (optionally ρ_up) * λ_f
        const double beta_f = std::max((rho_in_matrix ? rho_up : 1.0) * lam_f, 0.0);

        // a_face = β_f * |E| / d_perp
        const double a_face = beta_f * Eabs / dperp;
        (*a_f_Diff)[idxF] = a_face;

        // 交叉项：+β_f (∇x)_f · T_f
        const Vector& gP = grad[iP];
        const Vector& gN = grad[iN];
        const Vector  grad_f = (1.0 - gamma) * gP + gamma * gN;
        const double  s_cross = beta_f * (grad_f * F.vectorT);

        // 浮力项（若启用）：−β_f ρ̄ (g·ê)|E|  − β_f ρ̄ g·T
        double s_buoy = 0.0;
        if (enable_buoy) {
            const double gdot = gu.g * ehat;
            const double s_buoy_e = -beta_f * rho_inter * gdot * Eabs;
            const double s_buoy_t = -beta_f * rho_inter * (gu.g * F.vectorT);
            s_buoy = s_buoy_e + s_buoy_t;
        }

        const double s_face_raw = s_cross + s_buoy;

        // 轻度限幅（防止极端非正定情况）
        const double s_face_lim = std::max(-0.5 * a_face, std::min(0.5 * a_face, s_face_raw));
        (*s_f_Diff)[idxF] = s_face_lim;
    }
}

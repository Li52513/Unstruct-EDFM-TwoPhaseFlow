#pragma once
#include <algorithm>
#include <cmath>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "Initializer.h"
#include "WaterPropertyTable.h"
#include "CO2PropertyTable.h"

struct TimeTermStats {
	size_t nCells = 0, ghosts = 0, oor = 0;
	double a_min = +1e300, a_max = -1e300;
	double b_min = +1e300, b_max = -1e300;
};

/**
 * @brief 单相压力方程：时间项系数/常数项累加（仅时间项）
 * 形式： (V/dt) * C* * p^{n+1}  进对角;   RHS += (V/dt) * C* * p^n
 * 其中 C* = (∂(φρ)/∂p)|_{lin}, 线性化点可选：默认 n 层；若给出迭代场，则取中点 (n + iter)/2  详见Solver_TimeLoopUtils.h
 */

// 统一 θ 版：压力时间项（只依赖在 eval 点的物性）
// 形式： (V/dt)*C* * p^{n+1} 进对角； RHS += (V/dt)*C* * p^n
// 其中 C* = ρ_eval * (∂φ/∂p)|const + φ * (∂ρ/∂p)|eval
inline bool TimeTerm_Theta_SinglePhase_Flow(
    MeshManager& mgr,
    FieldRegistry& reg,
    double dt,
    double c_phi,                       // 常数孔隙度压缩性(1/Pa)，无则传 0
    // ——字段名——
    const std::string& phi_name,        // φ^n（通常就用 φ，不随时变）
    const std::string& p_old_name,      // p^n
    const std::string& rho_old_name,    // ρ^n(p^n,T^n)
    const std::string& p_eval_name,     // p_eval = (1-θ)p^n + θ p^{k}（θ=1 为全隐）
    const std::string& T_eval_name,     // T_eval（若 EoS 需要温度）
    const std::string& rho_eval_name,   // ρ_eval(p_eval,T_eval)
    const std::string& drdp_eval_name,  // (∂ρ/∂p)_eval
    const std::string& aC_name,         // 输出：时间项对角线
    const std::string& bC_name          // 输出：时间项右端
) {
    auto& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();

    if (cells.empty() || !(dt > 0.0)) {
        std::cerr << "[TimeTerm][P] invalid mesh or dt.\n";
        return false;
    }

    // 取场
    auto phi = reg.get<volScalarField>(phi_name);
    auto p_old = reg.get<volScalarField>(p_old_name);
    auto rho_n = reg.get<volScalarField>(rho_old_name);

    auto p_eval = reg.get<volScalarField>(p_eval_name);
    auto T_eval = reg.get<volScalarField>(T_eval_name); // 可能不直接用，但保留一致性
    auto rho_eval = reg.get<volScalarField>(rho_eval_name);
    auto drdp_eval = reg.get<volScalarField>(drdp_eval_name);

    if (!phi || !p_old || !rho_n || !p_eval || !rho_eval || !drdp_eval) {
        std::cerr << "[TimeTerm][P] missing fields: "
            << (phi ? "" : "phi ") << (p_old ? "" : "p_old ")
            << (rho_n ? "" : "rho_old ") << (p_eval ? "" : "p_eval ")
            << (rho_eval ? "" : "rho_eval ") << (drdp_eval ? "" : "drdp_eval ")
            << "\n";
        return false;
    }

    auto aC = reg.getOrCreate<volScalarField>(aC_name, cells.size(), 0.0);
    auto bC = reg.getOrCreate<volScalarField>(bC_name, cells.size(), 0.0);
    std::fill(aC->data.begin(), aC->data.end(), 0.0);
    std::fill(bC->data.begin(), bC->data.end(), 0.0);

    // 逐单元装配
    for (const auto& c : cells) {
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);

        const double V = std::max(0.0, c.volume);
        const double ph_n = std::min(1.0, std::max(0.0, (*phi)[i])); // φ^n
        const double pn = (*p_old)[i];                              // p^n
        const double rhoN = std::max(0.0, (*rho_n)[i]);               // ρ^n

        const double pe = (*p_eval)[i];                             // p_eval
        const double rhoE = std::max(0.0, (*rho_eval)[i]);            // ρ_eval
        const double drdpE = std::max(0.0, (*drdp_eval)[i]);           // (∂ρ/∂p)_eval ≥ 0 护栏

        // a = (V/dt) * [ φ^n * (∂ρ/∂p)_eval + φ^n * ρ_eval * c_phi ]
        const double a = (V / dt) * (ph_n * drdpE + ph_n * rhoE * c_phi);

        // b = (V/dt) * [ φ^n ρ^n  - φ^n ρ_eval + φ^n (∂ρ/∂p)_eval * p_eval + φ^n ρ_eval c_phi * p^n ]
        const double b = (V / dt) * (ph_n * rhoN - ph_n * rhoE + ph_n * drdpE * pe + ph_n * rhoE * c_phi * pn);

        (*aC)[i] = a;
        (*bC)[i] = b;
    }
    return true;
}

// 标记贴有“强值”边界（Dirichlet/Robin，a!=0）的单元：owner cell 置 1
//template<class ABC>
//inline std::vector<char> mark_strong_BC_cells(MeshManager& mgr, const ABC& bc, double a_eps = 1e-30)
//{
//    auto& mesh = mgr.mesh();
//    const auto& faces = mesh.getFaces();
//    const auto& id2idx = mesh.getCellId2Index();
//
//    std::vector<char> isStrong(mesh.getCells().size(), 0);
//    for (const auto& F : faces) {
//        if (!F.isBoundary()) continue;
//        double a = 0, b = 0, c = 0;
//        if (bc.getABC(F.id, a, b, c) && std::abs(a) > a_eps) {        // a!=0 → 强制值（含 Robin）
//            const int P = F.ownerCell;
//            if (P >= 0) isStrong[id2idx.at(P)] = 1;                    // 标记贴面的所有者单元
//        }
//    }
//    return isStrong;
//}

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





// 全隐（θ=1）时间项：单相流（压力方程）
// 线性化参考由 *lin* 字段提供（通常为上一次外迭代 p^{k-1}）
// 若为常物性：rho_lin 常数、drdp_lin=0 即可。

//inline bool TimeTerm_FullyImplicit_SinglePhase_Flow(
//    MeshManager& mgr,
//    FieldRegistry& reg,
//    double dt,
//    const std::string& c_phi_name,                       // 常数孔隙度压缩性(1/Pa)，无则传 0
//    // ——字段名（全部是已存在的体场）——
//    const std::string& phi_name,        // φ^n（通常就用 φ，不随时变）
//    const std::string& p_old_name,      // p^n
//    const std::string& rho_old_name,    // ρ^n(p^n, T^n)
//    const std::string& p_lin_name,      // p^⋆：线性化参考（如 p_prev）
//    const std::string& rho_lin_name,    // ρ^⋆=ρ(p^⋆,T^⋆)
//    const std::string& drdp_lin_name,   // (∂ρ/∂p)^⋆（常物性可为 0）
//    // ——输出——
//    const std::string& aC_name,         // 输出：时间项对角线
//    const std::string& bC_name          // 输出：时间项右端
//)
//{
//    auto& mesh = mgr.mesh();
//    const auto& cells = mesh.getCells();
//    const auto& id2idx = mesh.getCellId2Index();
//
//    // 取输入体场
//	auto phi = reg.get<volScalarField>(phi_name); // φ^n 孔隙率
//	auto p_n = reg.get<volScalarField>(p_old_name); 	// p^n 压力旧时层
//	auto rho_n = reg.get<volScalarField>(rho_old_name);  // ρ^n 密度旧时层
//	auto p_lin = reg.get<volScalarField>(p_lin_name);    // p^⋆ 线性化参考点 使用上一迭代层的值
//	auto rho_lin = reg.get<volScalarField>(rho_lin_name); // 上一迭代出的密度场 ρ^⋆
//	auto drdp_lin = reg.get<volScalarField>(drdp_lin_name); // (∂ρ/∂p)^⋆
//	auto c_phi = reg.get< volScalarField>(c_phi_name); // 常数孔隙度压缩性
//
//    // 创建输出
//    auto aC = reg.getOrCreate<volScalarField>(aC_name, cells.size(), 0.0);
//    auto bC = reg.getOrCreate<volScalarField>(bC_name, cells.size(), 0.0);
//    std::fill(aC->data.begin(), aC->data.end(), 0.0);
//    std::fill(bC->data.begin(), bC->data.end(), 0.0);
//
//    // 装配
//       for (const auto& c : cells) {
//            if (c.id < 0) continue;
//            const size_t i = id2idx.at(c.id);
//
//            const double V = std::max(0.0, c.volume);
//            const double ph_n = (*phi)[i];     // φ^n
//            const double pn = (*p_n)[i];                           // p^n
//            const double rhoN = std::max(0.0, (*rho_n)[i]);        // ρ^n
//			const double c_phi_value = (*c_phi)[i];                     // 常数孔隙度压缩性
//
//            const double pStar = (*p_lin)[i];                      // p^⋆
//            const double rhoStar = std::max(0.0, (*rho_lin)[i]);   // ρ^⋆
//            const double drdpStar = std::max(0.0, (*drdp_lin)[i]); // (∂ρ/∂p)^⋆ (非负护栏)
//
//            // LHS 对角
//            const double a = (V / dt) * ph_n * (drdpStar + rhoStar * c_phi_value);
//
//            // RHS 常数项（按推导式）
//            const double b = (V / dt) * (ph_n * rhoN
//                - ph_n * rhoStar
//                + ph_n * drdpStar * pStar
//                + ph_n * rhoStar * c_phi_value * pn);
//
//            (*aC)[i] = a;
//            (*bC)[i] = b;
//       }
//    return true;
//
//}
// 全隐（θ=1）时间项：单相“压力”
// 可选 strongBCmask：若提供，mask[i]!=0 的单元不写入时间项（用于强边界）


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
    const std::string& bC_name,             // 时间项常数项
    // ——可选：强边界单元屏蔽——
    const std::vector<char>* strongBCmask = nullptr
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
        if (strongBCmask && (*strongBCmask)[i]) continue;   // ★ 强边界单元不写时间项

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
            + ph * rhoStar * cphi * pn);

        (*aC)[i] = a;
        (*bC)[i] = b;
    }
    return true;
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
    const std::string& bC_name,       // 输出：右端
    // 可选：对 Dirichlet 单元做强制钉扎
    const std::vector<char>* strong_mask_cells = nullptr,   // 1=强制
    const std::vector<double>* T_target_cells = nullptr,   // 目标 T_b（若空则用 T_old）
    double pin_weight = 0.0                                      // 钉扎强度（0 表示不用）
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

    const bool use_pin = (strong_mask_cells && pin_weight > 0.0);

    for (const auto& c : cells) {
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);

        const double V = std::max(0.0, c.volume);
        const double Ceff_star = std::max(epsC, (*Ceff)[i]);
        const double T_old = (*Told)[i];

        // 常规时间项
        double a = V * inv_dt * Ceff_star;
        double b = V * inv_dt * Ceff_star * T_old;

        // ——强制钉扎：把行改成 (a + W)*T_P = b + W*T_target —— //
        if (use_pin && (*strong_mask_cells)[i]) {
            const double Ttar = (T_target_cells ? (*T_target_cells)[i] : T_old);
            a += pin_weight;
            b += pin_weight * Ttar;
        }

        (*aC)[i] = a;
        (*bC)[i] = b;
    }
    return true;
}

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

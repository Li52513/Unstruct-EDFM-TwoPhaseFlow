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
inline bool TimeTerm_Euler_SinglePhase_Flow
(
    MeshManager& mgr, FieldRegistry& reg,
    double dt, double dphi_dp_const,
    const std::string& phi_name = "phi",
    const std::string& p_old_name = "p_w_old",
    const std::string& rho_lin_name = "rho_time",
    const std::string& drho_dp_name = "drho_dp_time",
    const std::string& aC_name = "aC_time_p",
    const std::string& bC_name = "bC_time_p"
)
{
    auto& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();
    if (dt <= 0.0 || cells.empty()) return false;

    auto phi = reg.get<volScalarField>(phi_name);
    auto p0 = reg.get<volScalarField>(p_old_name);
    auto rhoL = reg.get<volScalarField>(rho_lin_name);
    auto drdp = reg.get<volScalarField>(drho_dp_name);
    if (!phi || !p0 || !rhoL || !drdp)
    {
        std::cerr << "[TimeTerm][P] missing fields: phi/p_old/rho_time/drho_dp_time\n";
        return false;
    }

    auto aC = reg.getOrCreate<volScalarField>(aC_name, cells.size(), 0.0);
    auto bC = reg.getOrCreate<volScalarField>(bC_name, cells.size(), 0.0);
    std::fill(aC->data.begin(), aC->data.end(), 0.0);
    std::fill(bC->data.begin(), bC->data.end(), 0.0);

    for (const auto& c : cells)
    {
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);

        const double vol = std::max(c.volume, 0.0);
        const double ph = std::max(0.0, (*phi)[i]);
        const double pn = (*p0)[i];
        const double rL = (*rhoL)[i];
        const double dr = (*drdp)[i];

        const double Cstar = rL * dphi_dp_const + ph * dr;
        const double a = vol * Cstar / dt;
        const double b = a * pn;           // (V/dt) C* p^n

        (*aC)[i] = a;
        (*bC)[i] = b;
    }
    return true;
}


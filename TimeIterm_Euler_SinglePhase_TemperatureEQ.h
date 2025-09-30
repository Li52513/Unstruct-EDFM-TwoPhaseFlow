#pragma once
#include <algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"

inline bool TimeTerm_Euler_SinglePhase_Temperature
(
    MeshManager & mgr,
    FieldRegistry & reg,
    double dt,
    double Ceff_floor = 1e-12,
    // —— 字段名 —— //
	const std::string& phi_name = "phi", // 孔隙度
	const std::string& rho_r_name = "rho_r", // 基岩密度
	const std::string& cp_r_name = "cp_r",   // 基岩比热容
    const std::string & rho_f_name = "rho_w",   // 单相水：rho_w / 单相气：rho_g
    const std::string & cp_f_name = "cp_w",    // 单相水：cp_w  / 单相气：cp_g
    const std::string & T_name = "T",
    const std::string & T_old_name = "T_old",
    const std::string & aC_name = "aC_time_T",
    const std::string & bC_name = "bC_time_T"
)
{
    auto& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();
    if (dt <= 0.0 || cells.empty()) return false;

    auto phi = reg.get<volScalarField>(phi_name);
    auto rho_r = reg.get<volScalarField>(rho_r_name);
    auto cp_r = reg.get<volScalarField>(cp_r_name);
    auto rho_f = reg.get<volScalarField>(rho_f_name);
    auto cp_f = reg.get<volScalarField>(cp_f_name);
    auto T0 = reg.get<volScalarField>(T_old_name);
    if (!phi || !rho_r || !cp_r || !rho_f || !cp_f || !T0) {
        std::cerr << "[TimeTerm][T] missing field(s): phi/rho_r/cp_r/rho_f/cp_f/T_old\n";
        return false;
    }

    auto aC = reg.getOrCreate<volScalarField>(aC_name, cells.size(), 0.0);
    auto bC = reg.getOrCreate<volScalarField>(bC_name, cells.size(), 0.0);
    std::fill(aC->data.begin(), aC->data.end(), 0.0);
    std::fill(bC->data.begin(), bC->data.end(), 0.0);

    for (const auto& c : cells)
    {
        if (c.id < 0) continue; // ghost
        const size_t i = id2idx.at(c.id);

        const double vol = std::max(c.volume, 0.0);
        const double ph = std::max(0.0, (*phi)[i]);
        const double rr = std::max(0.0, (*rho_r)[i]);
        const double cpr = std::max(0.0, (*cp_r)[i]);
        const double rf = std::max(0.0, (*rho_f)[i]);
        const double cpf = std::max(0.0, (*cp_f)[i]);
        const double Tn = (*T0)[i];

        const double Ceff = std::max(Ceff_floor,
            (1.0 - ph) * rr * cpr + ph * rf * cpf);

        const double a = vol * Ceff / dt;
        const double b = vol * Ceff * Tn / dt;

        (*aC)[i] = a;
        (*bC)[i] = b;
    }
    return true;
}
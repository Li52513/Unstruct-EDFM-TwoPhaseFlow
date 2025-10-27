//==================== TimeIterm_Euler_SinglePhase_TemperatureEQ.h ====================
#pragma once
#include <algorithm>
#include <string>
#include <memory>
#include <iostream>
#include "MeshManager.h"
#include "FieldRegistry.h"

// 针对单相而言-混合体积热容：Ceff = (1-phi)*rho_r*cp_r + phi*rho_f*cp_f
inline double mixtureHeatCapacity(double phi, double rho_r, double cp_r, double rho_f, double cp_f)
{
    phi = std::min(1.0, std::max(0.0, phi));
    rho_r = std::max(0.0, rho_r);
    cp_r = std::max(0.0, cp_r);
    rho_f = std::max(0.0, rho_f);
    cp_f = std::max(0.0, cp_f);
    return (1.0 - phi) * rho_r * cp_r + phi * rho_f * cp_f;
}

/**
 * 统一风格温度时间项：
 *   a_i = (V_i/dt) * C_eval_i
 *   b_i = a_i * T^n_i
 * 其中 C_eval = (1-θ_T)·C(n) + θ_T·C(lin) //引入高阶格式的接口
 *
 * 约定：
 *   - 固体骨架(φ,ρ_r,cp_r)只取单一场（默认用当前注册场）
 *   - 流体物性分两套：n层(rho_f_n,cp_f_n) 与 lin点(rho_f_lin,cp_f_lin)
 */


// 统一 θ 版：温度时间项（物性已在 eval 场更新到各字段中）
// 形式： (V/dt)*Ceff_eval * T^{n+1} 进对角； RHS += (V/dt)*Ceff_eval * T^n
inline bool TimeTerm_Theta_SinglePhase_Temperature
(
    MeshManager& mgr, FieldRegistry& reg,
    double dt, double Ceff_floor,
    // fields
    const std::string& phi_name,
    const std::string& rho_r_name, const std::string& cp_r_name,
    const std::string& rho_f_name, const std::string& cp_f_name,
    const std::string& T_old_name,
    const std::string& aC_name, const std::string& bC_name,
    const std::string& Ceff_out = ""
) {
    auto& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();
    if (cells.empty() || !(dt > 0.0)) return false;

    auto phi = reg.get<volScalarField>(phi_name);
    auto rr = reg.get<volScalarField>(rho_r_name);
    auto cpr = reg.get<volScalarField>(cp_r_name);
    auto rf = reg.get<volScalarField>(rho_f_name);
    auto cpf = reg.get<volScalarField>(cp_f_name);
    auto Tn = reg.get<volScalarField>(T_old_name);
    if (!phi || !rr || !cpr || !rf || !cpf || !Tn) return false;

    auto aC = reg.getOrCreate<volScalarField>(aC_name, cells.size(), 0.0);
    auto bC = reg.getOrCreate<volScalarField>(bC_name, cells.size(), 0.0);
    std::fill(aC->data.begin(), aC->data.end(), 0.0);
    std::fill(bC->data.begin(), bC->data.end(), 0.0);

    std::shared_ptr<volScalarField> CeffField;
    if (!Ceff_out.empty())
        CeffField = reg.getOrCreate<volScalarField>(Ceff_out, cells.size(), 0.0);

    for (const auto& c : cells) {
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);

        const double V = std::max(0.0, c.volume);
        const double ph = (*phi)[i];
        const double rrI = std::max(0.0, (*rr)[i]);
        const double cprI = std::max(0.0, (*cpr)[i]);
        const double rfI = std::max(0.0, (*rf)[i]);
        const double cpfI = std::max(0.0, (*cpf)[i]);

        // 物性已在 eval 场更新到这些字段里，直接取即可
        double Ceff = mixtureHeatCapacity(ph, rrI, cprI, rfI, cpfI);
        Ceff = std::max(Ceff_floor, Ceff);

        const double a = V * Ceff / dt;
        const double b = a * (*Tn)[i];   // 右端用 T^n

        (*aC)[i] = a;
        (*bC)[i] = b;
        if (CeffField) (*CeffField)[i] = Ceff;
    }
    return true;
}






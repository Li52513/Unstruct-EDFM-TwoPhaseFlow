//==================== TimeIterm_Euler_SinglePhase_TemperatureEQ.h ====================
#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>
#include "MeshManager.h"
#include "FieldRegistry.h"

// 混合体积热容：Ceff = (1-phi)*rho_r*cp_r + phi*rho_f*cp_f
inline double mixtureHeatCapacity(double phi, double rr, double cpr, double rf, double cpf) {
    phi = std::min(1.0, std::max(0.0, phi));
    return (1.0 - phi) * rr * cpr + phi * rf * cpf;
}

/**
 * 全隐Euler的温度时间项离散：
 *   A_ii += V*Ceff/dt,  b_i += (V*Ceff/dt)*T^n
 * 可选增强：
 *   - use_time_lin=true：在时间线性化点 T_time_lin 评估/混合Ceff
 *   - theta：Ceff^{n+theta} = (1-theta)*Ceff(T^n) + theta*Ceff(T_lin)
 *   - Ceff_out：把使用的Ceff写成场，便于诊断
 */
inline bool TimeTerm_Euler_SinglePhase_Temperature(
    MeshManager& mgr,
    FieldRegistry& reg,
    double dt,
    double Ceff_floor = 1e-12,
    // 字段名（与原版保持一致的默认）
    const std::string& phi_name = "phi",
    const std::string& rho_r_name = "rho_r",
    const std::string& cp_r_name = "cp_r",
    const std::string& rho_f_name = "rho_w",   // CO2时传 "rho_g"
    const std::string& cp_f_name = "cp_w",    // CO2时传 "cp_g"
    const std::string& T_name = "T",
    const std::string& T_old_name = "T_old",
    const std::string& aC_name = "aC_time_T",
    const std::string& bC_name = "bC_time_T",
    // ——可选增强控制——
    bool   use_time_lin = false,
    const std::string& T_lin_name = "T_time_lin",
    double theta = 0.0,                 // 0=旧时层（与原实现等价）；0.5≈CN风格
    const std::string& Ceff_out = ""    // 非空则输出 Ceff 场
) {
    auto& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();
    if (cells.empty() || !(dt > 0.0)) {
        std::cerr << "[TimeTerm][T] invalid mesh or dt.\n";
        return false;
    }

    auto phi = reg.get<volScalarField>(phi_name);
    auto rho_r = reg.get<volScalarField>(rho_r_name);
    auto cp_r = reg.get<volScalarField>(cp_r_name);
    auto rho_f = reg.get<volScalarField>(rho_f_name);
    auto cp_f = reg.get<volScalarField>(cp_f_name);
    auto Tn = reg.get<volScalarField>(T_old_name);
    if (!phi || !rho_r || !cp_r || !rho_f || !cp_f || !Tn) {
        std::cerr << "[TimeTerm][T] missing phi/rho_r/cp_r/rho_f/cp_f/T_old\n";
        return false;
    }

    std::shared_ptr<volScalarField> Tlin;
    if (use_time_lin) {
        Tlin = reg.get<volScalarField>(T_lin_name);
        if (!Tlin) {
            std::cerr << "[TimeTerm][T] use_time_lin=true but '" << T_lin_name << "' not found.\n";
            return false;
        }
        theta = std::min(1.0, std::max(0.0, theta));
    }

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
        const double ph = std::min(1.0, std::max(0.0, (*phi)[i]));
        const double rr = std::max(0.0, (*rho_r)[i]);
        const double cpr = std::max(0.0, (*cp_r)[i]);
        const double rf = std::max(0.0, (*rho_f)[i]);
        const double cpf = std::max(0.0, (*cp_f)[i]);

        const double Cn = mixtureHeatCapacity(ph, rr, cpr, rf, cpf);
        const double Clin = use_time_lin ? mixtureHeatCapacity(ph, rr, cpr, rf, cpf) : Cn;
        double Ceff = (1.0 - theta) * Cn + theta * Clin;
        Ceff = std::max(Ceff_floor, Ceff);

        const double a = V * Ceff / dt;
        const double b = a * (*Tn)[i];   // 保持线性结构：RHS仍用 T^n

        (*aC)[i] = a;
        (*bC)[i] = b;
        if (CeffField) (*CeffField)[i] = Ceff;
    }
    return true;
}




//#pragma once
//#include <algorithm>
//#include "MeshManager.h"
//#include "FieldRegistry.h"
//
//inline bool TimeTerm_Euler_SinglePhase_Temperature
//(
//    MeshManager & mgr,
//    FieldRegistry & reg,
//    double dt,
//    double Ceff_floor = 1e-12,
//    // —— 字段名 —— //
//	const std::string& phi_name = "phi", // 孔隙度
//	const std::string& rho_r_name = "rho_r", // 基岩密度
//	const std::string& cp_r_name = "cp_r",   // 基岩比热容
//    const std::string & rho_f_name = "rho_w",   // 单相水：rho_w / 单相气：rho_g
//    const std::string & cp_f_name = "cp_w",    // 单相水：cp_w  / 单相气：cp_g
//    const std::string & T_name = "T",
//    const std::string & T_old_name = "T_old",
//    const std::string & aC_name = "aC_time_T",
//    const std::string & bC_name = "bC_time_T"
//)
//{
//    auto& mesh = mgr.mesh();
//    const auto& cells = mesh.getCells();
//    const auto& id2idx = mesh.getCellId2Index();
//    if (dt <= 0.0 || cells.empty()) return false;
//
//    auto phi = reg.get<volScalarField>(phi_name);
//    auto rho_r = reg.get<volScalarField>(rho_r_name);
//    auto cp_r = reg.get<volScalarField>(cp_r_name);
//    auto rho_f = reg.get<volScalarField>(rho_f_name);
//    auto cp_f = reg.get<volScalarField>(cp_f_name);
//    auto T0 = reg.get<volScalarField>(T_old_name);
//    if (!phi || !rho_r || !cp_r || !rho_f || !cp_f || !T0) {
//        std::cerr << "[TimeTerm][T] missing field(s): phi/rho_r/cp_r/rho_f/cp_f/T_old\n";
//        return false;
//    }
//
//    auto aC = reg.getOrCreate<volScalarField>(aC_name, cells.size(), 0.0);
//    auto bC = reg.getOrCreate<volScalarField>(bC_name, cells.size(), 0.0);
//    std::fill(aC->data.begin(), aC->data.end(), 0.0);
//    std::fill(bC->data.begin(), bC->data.end(), 0.0);
//
//    for (const auto& c : cells)
//    {
//        if (c.id < 0) continue; // ghost
//        const size_t i = id2idx.at(c.id);
//
//        const double vol = std::max(c.volume, 0.0);
//        const double ph = std::max(0.0, (*phi)[i]);
//        const double rr = std::max(0.0, (*rho_r)[i]);
//        const double cpr = std::max(0.0, (*cp_r)[i]);
//        const double rf = std::max(0.0, (*rho_f)[i]);
//        const double cpf = std::max(0.0, (*cp_f)[i]);
//        const double Tn = (*T0)[i];
//
//        const double Ceff = std::max(Ceff_floor,
//            (1.0 - ph) * rr * cpr + ph * rf * cpf);
//
//        const double a = vol * Ceff / dt;
//        const double b = vol * Ceff * Tn / dt;
//
//        (*aC)[i] = a;
//        (*bC)[i] = b;
//    }
//    return true;
//}
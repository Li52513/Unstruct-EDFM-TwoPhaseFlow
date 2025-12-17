#pragma once
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "PhysicalPropertiesManager_TwoPhase.h"
#include "SolverContrlStrName.h"

namespace IMPES_Iteration
{
    /**
    * @brief Assemble IMPES pressure accumulation (ddt) term.
    * Equation :
    A_{ii}^{\text{time}} = \frac{V_i}{\Delta t} C_{total}^\nu $
    $$ b_i^{\text{time}} = \frac{V_i}{\Delta t} \left[ \underbrace{(\phi \rho_t)^n}_{\text{上一时刻真实质量}} - \underbrace{(\phi \rho_t)^\nu}_{\text{当前猜测质量}} + \underbrace{C_{total}^\nu P^\nu}_{\text{导数项修正}} \right] $$
    * Computes diagonal/time coefficient aC and source bC for the implicit pressure
    */
    inline bool TimeTerm_IMPES_Pressure(
        MeshManager& mgr,           //网格信息
		FieldRegistry& reg,         //场注册表
		double dt,                  //时间步长
        const std::string& p_old_name,
        const std::string& p_eval_name,
        const std::string& aC_name,
        const std::string& bC_name)
    {
        if (dt <= 0.0)
        {
            std::cerr << "[IMPES][TimeTerm] invalid dt.\n";
            return false;
        }

        auto& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        auto p_old = reg.get<volScalarField>(p_old_name);
        auto p_eval = reg.get<volScalarField>(p_eval_name);
        auto phi = reg.get<volScalarField>(PhysicalProperties_string::Rock().phi_tag);
        auto s_w = reg.get<volScalarField>(SaturationEquation_String().saturation);

        auto rho_w_eval = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_tag);
        auto rho_w_old = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_old_tag);
        auto drho_w_dp = reg.get<volScalarField>(PhysicalProperties_string::Water().drho_w_dp_tag);

        auto rho_g_eval = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_tag);
        auto rho_g_old = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_old_tag);
        auto drho_g_dp = reg.get<volScalarField>(PhysicalProperties_string::CO2().drho_g_dp_tag);

        auto c_r = reg.get<volScalarField>(PhysicalProperties_string::Rock().c_r_tag);

        if (!p_old || !p_eval || !phi || !s_w || !rho_w_eval || !rho_w_old || !drho_w_dp || !rho_g_eval || !rho_g_old || !drho_g_dp || !c_r)
        {
            std::cerr << "[IMPES][TimeTerm] missing fields for pressure accumulation.\n";
            return false;
        }

        auto aC = reg.getOrCreate<volScalarField>(aC_name.c_str(), cells.size(), 0.0);
        auto bC = reg.getOrCreate<volScalarField>(bC_name.c_str(), cells.size(), 0.0);
        std::fill(aC->data.begin(), aC->data.end(), 0.0);
        std::fill(bC->data.begin(), bC->data.end(), 0.0);

        const double inv_dt = 1.0 / dt;

        // 饱和度夹取范围（可根据实际需要调整）
        const double Sw_min = 1e-8;
        const double Sw_max = 1.0 - 1e-8;

        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const double V = std::max(0.0, c.volume);
            const double phi_i = std::max(0.0, std::min(1.0, (*phi)[i]));
            double sw_i = (*s_w)[i];
            sw_i = std::max(Sw_min, std::min(Sw_max, sw_i));
            const double sg_i = 1.0 - sw_i;

            const double p_n = (*p_old)[i]; // previous time p_w^n
            const double p_star = (*p_eval)[i]; // current evaluated p_w

            const double rho_w_n = std::max(0.0, (*rho_w_old)[i]);
            const double rho_w_star = std::max(0.0, (*rho_w_eval)[i]);

            const double rho_g_n = std::max(0.0, (*rho_g_old)[i]);
            const double rho_g_star = std::max(0.0, (*rho_g_eval)[i]);

            const double drho_w_dp_i = std::max(0.0, (*drho_w_dp)[i]);
            const double drho_g_dp_i = std::max(0.0, (*drho_g_dp)[i]);

            // 计算总密度及其对压力的导数
            const double rho_t_n = sw_i * rho_w_n + sg_i * rho_g_n;
            const double rho_t_star = sw_i * rho_w_star + sg_i * rho_g_star;

            const double drho_tdp_star = sw_i * drho_w_dp_i + sg_i * drho_g_dp_i;

            // 取出岩石压缩系数
            const double cr = std::max(0.0, (*c_r)[i]);

            //总压缩系数
            const double C_tot = phi_i * (drho_tdp_star + rho_t_star * cr);
            // 组装时间项系数
            double a = V * inv_dt * C_tot;
            // b = V/dt * [ φ ρ_t^n - φ ρ_t^ν + φ (dρ_t/dP) P^ν + φ ρ_t^ν c_r P^n ]
            double b = V * inv_dt * (phi_i * rho_t_n - phi_i * rho_t_star + phi_i * drho_tdp_star * p_star + phi_i * rho_t_star * cr * p_star);

            (*aC)[i] = a;
            (*bC)[i] = b;
        }
        return true;

    }
}
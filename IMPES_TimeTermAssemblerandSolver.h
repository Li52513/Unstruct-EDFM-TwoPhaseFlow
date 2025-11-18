#pragma once

#include <string>
#include <iostream>
#include <algorithm>
#include <vector>

#include "MeshManager.h"
#include "FieldRegistry.h"

namespace IMPES
{
    /**
     * @brief Assemble IMPES pressure accumulation (ddt) term.
     *
     * Computes diagonal/time coefficient aC and source bC for the implicit pressure
     * step using porosity, density and d(rho)/dp. Supports an optional strong_mask
     * to skip Dirichlet cells.
     */
    inline bool TimeTerm_IMPES_Pressure(
        MeshManager& mgr,
        FieldRegistry& reg,
        double dt,
        const std::string& phi_name,
        const std::string& p_old_name,
        const std::string& p_eval_name,
        const std::string& rho_old_name,
        const std::string& rho_eval_name,
        const std::string& drho_dp_name,
        double rock_compressibility,
        const std::string& aC_name,
        const std::string& bC_name,
        const std::vector<char>* strong_mask = nullptr)
    {
        if (dt <= 0.0)
        {
            std::cerr << "[IMPES][TimeTerm] invalid dt.\n";
            return false;
        }

        auto& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        auto phi = reg.get<volScalarField>(phi_name);
        auto p_old = reg.get<volScalarField>(p_old_name);
        auto p_eval = reg.get<volScalarField>(p_eval_name);
        auto rho_old = reg.get<volScalarField>(rho_old_name);
        auto rho_eval = reg.get<volScalarField>(rho_eval_name);
        auto drho_dp = reg.get<volScalarField>(drho_dp_name);

        if (!phi || !p_old || !p_eval || !rho_old || !rho_eval || !drho_dp)
        {
            std::cerr << "[IMPES][TimeTerm] missing fields for pressure accumulation.\n";
            return false;
        }

        auto aC = reg.getOrCreate<volScalarField>(aC_name.c_str(), cells.size(), 0.0);
        auto bC = reg.getOrCreate<volScalarField>(bC_name.c_str(), cells.size(), 0.0);
        std::fill(aC->data.begin(), aC->data.end(), 0.0);
        std::fill(bC->data.begin(), bC->data.end(), 0.0);

        const double inv_dt = 1.0 / dt;

        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            if (strong_mask && (*strong_mask)[i]) continue;

            const double V = std::max(0.0, c.volume);
            const double phi_i = std::max(0.0, std::min(1.0, (*phi)[i]));
            const double rho_n = std::max(0.0, (*rho_old)[i]);
            const double rho_star = std::max(0.0, (*rho_eval)[i]);
            const double drdp = std::max(0.0, (*drho_dp)[i]);
            const double p_n = (*p_old)[i];
            const double p_star = (*p_eval)[i];

            const double a = V * inv_dt * phi_i * (drdp + rho_star * rock_compressibility);
            const double b = V * inv_dt * (phi_i * rho_n
                - phi_i * rho_star
                + phi_i * drdp * p_star
                + phi_i * rho_star * rock_compressibility * p_n);

            (*aC)[i] = a;
            (*bC)[i] = b;
        }

        return true;
    }

	/**
	 * @brief Assemble IMPES pressure accumulation (ddt) term for two-phase flow.
	 *
	 * Computes diagonal/time coefficient aC and source bC for the implicit pressure
	 * step using porosity, saturations, densities and compressibilities of both phases.
	 * Supports an optional strong_mask to skip Dirichlet cells.
     * 
     * Equation£º\phi\left(c_{r}+c_{w}s_{w}\rho_{w}+c_{g}s_{g}\rho_{g}\right)\frac{\partial p}{\partial t} 
	 * phi_name-> \phi
	 * sat_w_name -> s_{w}
	 * rho_w_name -> \rho_{w}
	 * rho_g_name -> \rho_{g}
	 * comp_w_name -> c_{w}
	 * comp_g_name -> c_{g}
	 * p_old_name -> p^{n}
	 * rock_compressibility -> c_{r}
	 * aC_name -> diagonal/time coefficient
	 * bC_name -> source term
	 */
    inline bool TimeTerm_IMPES_Pressure_new(
        MeshManager& mgr,
        FieldRegistry& reg,
        double dt,
        const std::string& phi_name,
        const std::string& sat_w_name,
        const std::string& rho_w_name,
        const std::string& rho_g_name,
        const std::string& comp_w_name,
        const std::string& comp_g_name,
        const std::string& p_old_name,
        double rock_compressibility,
        const std::string& aC_name,
        const std::string& bC_name,
        const std::vector<char>* strong_mask = nullptr
    )
    {
        if (dt <= 0.0)
        {
            std::cerr << "[IMPES][TimeTerm] invalid dt.\n";
            return false;
        }

        auto& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        auto phi = reg.get<volScalarField>(phi_name);
        auto sw = reg.get<volScalarField>(sat_w_name);
        auto rho_w = reg.get<volScalarField>(rho_w_name);
        auto rho_g = reg.get<volScalarField>(rho_g_name);
        auto comp_w = reg.get<volScalarField>(comp_w_name);
        auto comp_g = reg.get<volScalarField>(comp_g_name);
        auto p_old = reg.get<volScalarField>(p_old_name);

		if (!phi || !sw || !rho_w || !rho_g || !comp_w || !comp_g || !p_old)
		{
			std::cerr << "[IMPES][TimeTerm] missing fields for pressure accumulation.\n";
			return false;
		}


        auto aC = reg.getOrCreate<volScalarField>(aC_name.c_str(), cells.size(), 0.0);
        auto bC = reg.getOrCreate<volScalarField>(bC_name.c_str(), cells.size(), 0.0);
        std::fill(aC->data.begin(), aC->data.end(), 0.0);
        std::fill(bC->data.begin(), bC->data.end(), 0.0);

		const double inv_dt = 1.0 / dt;

        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            if (strong_mask && (*strong_mask)[i]) continue;

            const double V = std::max(0.0, c.volume);
            const double phi_i = std::max(0.0, std::min(1.0, (*phi)[i]));
            const double sw_i = std::max(0.0, std::min(1.0, (*sw)[i]));
            const double sg_i = std::max(0.0, std::min(1.0, 1.0 - sw_i));
            const double rho_w_i = std::max(0.0, (*rho_w)[i]);
            const double rho_g_i = std::max(0.0, (*rho_g)[i]);
            const double cw_i = std::max(0.0, (*comp_w)[i]);
            const double cg_i = std::max(0.0, (*comp_g)[i]);

            const double accum = phi_i * (rock_compressibility + sw_i * rho_w_i * cw_i + sg_i * rho_g_i * cg_i);
            const double a = V * accum * inv_dt;
            const double b = a * (*p_old)[i];

            (*aC)[i] = a;
            (*bC)[i] = b;
        }
        return true;
    }






} // namespace IMPES

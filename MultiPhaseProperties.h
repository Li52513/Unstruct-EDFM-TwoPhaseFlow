#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "MeshManager.h"      // For MeshManager and FractureNetwork
#include "FieldRegistry.h"      // For FieldRegistry and volScalarField
#include "InitConfig.h"         // Where VGParams and RelPermParams are defined
#include "CapRelPerm.h"         // Where vg_params_valid, pc_vG, etc. are defined


namespace multiPhase {
	
	/**
	 * @brief Updates all two-phase derived properties for the ROCK MATRIX cells.
	 *
	 * This function is the core of the property update step for the rock matrix. It takes the
	 * primary variable Sw, calculates all dependent capillary and relative permeability properties,
	 * and computes the effective mixed thermal properties for the energy equation.
	 * It is designed to be robust, with checks for valid inputs and parameters.
	 *
	 * @param mgr The MeshManager, providing access to the mesh topology.
	 * @param reg The FieldRegistry for the ROCK MATRIX, where all fields are stored.
	 * @param vg The Van Genuchten model parameters.
	 * @param rp The Mualem relative permeability model parameters.
	 * @return true if the update was successful, false otherwise.
	 */

	inline bool updateRockTwoPhaseProperties_IMPES(
		MeshManager& mgr,
		FieldRegistry& reg,
		const VGParams& vg,
		const RelPermParams& rp)
	{
		// ===== Step 1: Sanity Checks & Setup =====
		if (!vg_params_valid(vg)) {
			std::cerr << "ERROR [multiPhase::Rock]: Invalid Van Genuchten parameters provided." << std::endl;
			return false;
		}

		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();

		// ===== Step 2: Get All Required INPUT Fields =====
		// Primary variable:
		auto s_w = reg.get<volScalarField>("s_w");
		auto p_w = reg.get<volScalarField>("p_w");
		// Rock properties:
		auto phi_r = reg.get<volScalarField>("phi_r");
		auto rho_r = reg.get<volScalarField>("rho_r");
		auto cp_r = reg.get<volScalarField>("cp_r");
		auto lambda_r = reg.get<volScalarField>("lambda_r"); //导热系数
		// Water properties:
		auto rho_w = reg.get<volScalarField>("rho_w");
		auto mu_w = reg.get<volScalarField>("mu_w");
		auto cp_w = reg.get<volScalarField>("cp_w");
		auto k_w = reg.get<volScalarField>("k_w");
		auto c_w = reg.get<volScalarField>("c_w"); //水的压缩系数
		// CO2 properties:
		auto rho_g = reg.get<volScalarField>("rho_g");
		auto mu_g = reg.get<volScalarField>("mu_g");
		auto cp_g = reg.get<volScalarField>("cp_g");
		auto k_g = reg.get<volScalarField>("k_g");
		auto c_g = reg.get<volScalarField>("c_g"); //气体的压缩系数

		// Check if any pointer is null. This is a critical safety check.
		if (!s_w || !p_w || !phi_r || !rho_r || !cp_r || !lambda_r ||
			!rho_w || !mu_w || !cp_w || !k_w ||
			!rho_g || !mu_g || !cp_g || !k_g||!c_g||!c_w) {
			std::cerr << "ERROR [multiPhase::Rock]: Missing one or more required fields for property update." << std::endl;
			return false;
		}

		// ===== Step 3: Get or Create All OUTPUT Fields =====
		auto Pc = reg.getOrCreate<volScalarField>("Pc", n, 0.0);
		auto p_g = reg.getOrCreate<volScalarField>("p_g", n, 0.0); // Gas pressure
		auto k_rw = reg.getOrCreate<volScalarField>("k_rw", n, 0.0);
		auto k_rg = reg.getOrCreate<volScalarField>("k_rg", n, 0.0);
		auto dPc_dSw = reg.getOrCreate<volScalarField>("dPc_dSw", n, 0.0);
		// Mobilities (lambda = k_r / mu) are central to the pressure equation
		auto lambda_w = reg.getOrCreate<volScalarField>("lambda_w", n, 0.0);
		auto lambda_g = reg.getOrCreate<volScalarField>("lambda_g", n, 0.0);
		auto lambda_t = reg.getOrCreate<volScalarField>("lambda_t", n, 0.0);
		//Time term  Total density: ρ_t = S_w * ρ_w + S_g * ρ_g and Total Compressibility (c_t):c_t = c_r + (S_w * ρ_w * c_w + S_g * ρ_g * c_g) / (S_w * ρ_w + S_g * ρ_g)
		auto rho_t = reg.getOrCreate<volScalarField>("rho_t", n, 0.0);
		auto c_t = reg.getOrCreate<volScalarField>("c_t", n, 0.0);
		auto c_r = reg.get<volScalarField>("c_r"); // Rock compressibility

		// Effective thermal properties
		auto C_eff = reg.getOrCreate<volScalarField>("C_eff", n, 0.0);
		auto lambda_eff = reg.getOrCreate<volScalarField>("lambda_eff", n, 0.0);

		// ===== Step 4: Main Loop over all Rock Cells =====
#pragma omp parallel for schedule(static)
		for (int ic = 0; ic < static_cast<int>(cells.size()); ++ic)
		{
			const auto& cell = cells[ic];
			if (cell.id < 0) continue;
			const size_t i = mesh.getCellId2Index().at(cell.id);

			// --- 4.1. Get Sw and perform clipping (robustness) ---
			double sw_val = (*s_w)[i]; //主变量-已经在初场初始化的时候进行了初始化 
			sw_val = clamp(sw_val, vg.Swr, 1.0 - vg.Sgr); //计算有效饱和度
			(*s_w)[i] = sw_val; // Write the clipped value back for consistency

			// --- 4.2. Calculate Pc, RelPerm, and derivatives ---
			(*Pc)[i] = pc_vG(sw_val, vg);
			kr_Mualem_vG(sw_val, vg, rp, (*k_rw)[i], (*k_rg)[i]);
			(*dPc_dSw)[i] = dpc_dSw_vG(sw_val, vg);

			// --- 4.3. Update Gas Pressure (p_g = p_w + Pc) ---
			(*p_g)[i] = (*p_w)[i] + (*Pc)[i];

			// --- 4.4. Calculate Phase Mobilities (lambda = k_r / mu) ---
			(*lambda_w)[i] = (*k_rw)[i] / std::max((*mu_w)[i], kTiny);
			(*lambda_g)[i] = (*k_rg)[i] / std::max((*mu_g)[i], kTiny);
			(*lambda_t)[i] = (*lambda_w)[i] + (*lambda_g)[i];

			// --- 4.4. Calculate Total Density and Total Compressibility ---
			const double sg_val = 1.0 - sw_val;
			const double rho_w_val = (*rho_w)[i];
			const double rho_g_val = (*rho_g)[i];
			const double c_w_val = (*c_w)[i];
			const double c_g_val = (*c_g)[i];
			// Total density
			(*rho_t)[i] = sw_val * rho_w_val + sg_val * rho_g_val;
			// Total compressibility
			(*c_t)[i] = (*c_r)[i] + (sw_val * rho_w_val * c_w_val + sg_val * rho_g_val * c_g_val) / std::max((*rho_t)[i], kTiny);

			// --- 4.5. Calculate Effective Thermal Properties ---
			const double phi_val = (*phi_r)[i];
			const double sg_val = 1.0 - sw_val;

			// Effective volumetric heat capacity
			(*C_eff)[i] = (1.0 - phi_val) * (*rho_r)[i] * (*cp_r)[i] +
				phi_val * (sw_val * (*rho_w)[i] * (*cp_w)[i] +
					sg_val * (*rho_g)[i] * (*cp_g)[i]);

			// Effective thermal conductivity
			(*lambda_eff)[i] = (1.0 - phi_val) * (*lambda_r)[i] +
				phi_val * (sw_val * (*k_w)[i] +
					sg_val * (*k_g)[i]);
		}

		return true; // Indicate success
	}


	/**
	 * @brief Updates all two-phase derived properties for the FRACTURE elements.
	 *
	 * This function mirrors the logic of `updateRockTwoPhaseProperties` but operates on
	 * the fracture elements and their corresponding FieldRegistry (`reg_fr`).
	 *
	 * @param mgr The MeshManager, providing access to the fracture network.
	 * @param reg_fr The FieldRegistry for the FRACTURE network.
	 * @param vg The Van Genuchten model parameters (can be the same or different from rock).
	 * @param rp The Mualem relative permeability model parameters.
	 * @return true if the update was successful, false otherwise.
	 */
	inline bool updateFractureTwoPhaseProperties_IMPES(
		MeshManager& mgr,
		FieldRegistry& reg_fr,
		const VGParams& vg,
		const RelPermParams& rp)
	{
		// ===== Step 1: Sanity Checks & Setup =====
		if (!vg_params_valid(vg)) {
			std::cerr << "ERROR [multiPhase::Fracture]: Invalid Van Genuchten parameters provided." << std::endl;
			return false;
		}

		const auto& frNet = mgr.fracture_network();
		size_t n_elements = 0;
		for (const auto& frac : frNet.fractures) {
			n_elements += frac.elements.size();
		}
		if (n_elements == 0) return true; // Nothing to do

		// ===== Step 2: Get All Required INPUT Fields from reg_fr =====
		// Primary variables:
		auto s_w_fr = reg_fr.get<volScalarField>("fr_s_w");
		auto p_w_fr = reg_fr.get<volScalarField>("fr_p_w");
		// Fracture solid properties:
		auto phi_fr = reg_fr.get<volScalarField>("fr_phi");
		auto rho_fr = reg_fr.get<volScalarField>("fr_rho_r"); // Note: field names from your setup
		auto cp_fr = reg_fr.get<volScalarField>("fr_cp_r");
		auto lambda_fr = reg_fr.get<volScalarField>("fr_lambda_r");
		// Fluid properties in fracture:
		auto rho_w_fr = reg_fr.get<volScalarField>("fr_rho_w");
		auto mu_w_fr = reg_fr.get<volScalarField>("fr_mu_w");
		auto cp_w_fr = reg_fr.get<volScalarField>("fr_cp_w");
		auto k_w_fr = reg_fr.get<volScalarField>("fr_k_w");
		auto rho_g_fr = reg_fr.get<volScalarField>("fr_rho_g");
		auto mu_g_fr = reg_fr.get<volScalarField>("fr_mu_g");
		auto cp_g_fr = reg_fr.get<volScalarField>("fr_cp_g");
		auto k_g_fr = reg_fr.get<volScalarField>("fr_k_g");
		if (!s_w_fr || !p_w_fr || !phi_fr || !rho_fr || !cp_fr || !lambda_fr ||
			!rho_w_fr || !mu_w_fr || !cp_w_fr || !k_w_fr ||
			!rho_g_fr || !mu_g_fr || !cp_g_fr || !k_g_fr) {
			std::cerr << "ERROR [multiPhase::Fracture]: Missing one or more required fields for property update." << std::endl;
			return false;
		}

		// ===== Step 3: Get or Create All OUTPUT Fields in reg_fr =====
		auto Pc_fr = reg_fr.getOrCreate<volScalarField>("fr_Pc", n_elements, 0.0);
		auto p_g_fr = reg_fr.getOrCreate<volScalarField>("fr_p_g", n_elements, 0.0);
		auto k_rw_fr = reg_fr.getOrCreate<volScalarField>("fr_k_rw", n_elements, 0.0);
		auto k_rg_fr = reg_fr.getOrCreate<volScalarField>("fr_k_rg", n_elements, 0.0);
		auto lambda_w_fr = reg_fr.getOrCreate<volScalarField>("fr_lambda_w", n_elements, 0.0);
		auto lambda_g_fr = reg_fr.getOrCreate<volScalarField>("fr_lambda_g", n_elements, 0.0);
		auto lambda_t_fr = reg_fr.getOrCreate<volScalarField>("fr_lambda_t", n_elements, 0.0);
		auto C_eff_fr = reg_fr.getOrCreate<volScalarField>("fr_C_eff", n_elements, 0.0);
		auto lambda_eff_fr = reg_fr.getOrCreate<volScalarField>("fr_lambda_eff", n_elements, 0.0);

		size_t global_idx = 0;
		for (const auto& frac : frNet.fractures)
		{
			for (const auto& elem : frac.elements) 
			{
				// The index `global_idx` corresponds to the position in the volScalarField data vectors.

				// --- 4.1. Clipping ---
				double sw_val = (*s_w_fr)[global_idx];
				sw_val = clamp(sw_val, vg.Swr, 1.0 - vg.Sgr);
				(*s_w_fr)[global_idx] = sw_val;

				// --- 4.2. Pc & RelPerm ---
				(*Pc_fr)[global_idx] = pc_vG(sw_val, vg);
				kr_Mualem_vG(sw_val, vg, rp, (*k_rw_fr)[global_idx], (*k_rg_fr)[global_idx]);

				// --- 4.3. Gas Pressure ---
				(*p_g_fr)[global_idx] = (*p_w_fr)[global_idx] + (*Pc_fr)[global_idx];

				// --- 4.4. Mobilities ---
				(*lambda_w_fr)[global_idx] = (*k_rw_fr)[global_idx] / std::max((*mu_w_fr)[global_idx], kTiny);
				(*lambda_g_fr)[global_idx] = (*k_rg_fr)[global_idx] / std::max((*mu_g_fr)[global_idx], kTiny);
				(*lambda_t_fr)[global_idx] = (*lambda_w_fr)[global_idx] + (*lambda_g_fr)[global_idx];

				// --- 4.5. Effective Thermal Properties ---
				const double phi_val = (*phi_fr)[global_idx];
				const double sg_val = 1.0 - sw_val;

				(*C_eff_fr)[global_idx] = (1.0 - phi_val) * (*rho_fr)[global_idx] * (*cp_fr)[global_idx] +
					phi_val * (sw_val * (*rho_w_fr)[global_idx] * (*cp_w_fr)[global_idx] +
						sg_val * (*rho_g_fr)[global_idx] * (*cp_g_fr)[global_idx]);

				(*lambda_eff_fr)[global_idx] = (1.0 - phi_val) * (*lambda_fr)[global_idx] +
					phi_val * (sw_val * (*k_w_fr)[global_idx] +
						sg_val * (*k_g_fr)[global_idx]);

				global_idx++;
			}
		}

		return true; // Indicate success
	}



}
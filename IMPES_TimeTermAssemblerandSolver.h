#pragma once

#include <string>
#include <iostream>
#include <algorithm>
#include <vector>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "0_PhysicalParametesCalculateandUpdata.h"

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
} // namespace IMPES



namespace IMPES_revised
{
    /*
     * @brief Assemble IMPES pressure accumulation (ddt) term based on a total compressibility formulation.
     *
     * This function models the time derivative term of the pressure equation using the
     * total compressibility concept, which simplifies the linearization of the mass accumulation term.
     *
     * The underlying physical equation being discretized is the conservation of total mass:
     *      ∂(φ * M_t)/∂t + ... = 0 , where M_t = S_w*ρ_w + S_g*ρ_g
     *
     * This term is linearized with respect to pressure 'p' as follows:
     *      ∂(φ * M_t)/∂t  ≈  [∂(φ * M_t)/∂p] * (∂p/∂t)
     *
     * The result of this linearization is commonly written in a compact engineering form:
     *
     *      [∂(φ * M_t)/∂p] * (∂p/∂t)  ≈  φ * ρ_t * c_t * (∂p/∂t)
     *
     * where each component is defined as:
     *
     * 1. Total Density (ρ_t):
     *    The saturation-weighted average density of the fluid mixture.
     *
     *      ρ_t = S_w * ρ_w + S_g * ρ_g
     *
     * 2. Total Compressibility (c_t):
     *    The effective compressibility of the entire system (rock + fluids). It is the sum of
     *    the rock compressibility and the mass-fraction-weighted average of the fluid compressibilities.
     *
     *      c_t = c_r + (S_w * ρ_w * c_w + S_g * ρ_g * c_g) / (S_w * ρ_w + S_g * ρ_g)
     *
     * where:
     *    - φ:   Porosity
     *    - S_w: Water saturation
     *    - S_g: Gas saturation (S_g = 1 - S_w)
     *    - ρ_w: Water density
     *    - ρ_g: Gas density
     *    - c_r: Rock compressibility (can be constant or a field)
     *    - c_w: Water compressibility (often assumed constant)
     *    - c_g: Gas compressibility (a strong function of pressure, e.g., c_g ≈ 1/p)
     */
  //  inline bool TimeTerm_IMPES_Pressure(
  //      MeshManager& mgr,
  //      FieldRegistry& reg,
  //      double dt,
  //      const std::string& phi_name,
  //      const std::string& p_old_name,
  //      const std::string& p_name, // Note: This is unused in this specific formulation but kept for signature consistency.
  //      const std::string& rho_t_name,
  //      const std::string& s_w_name,
  //      const std::string& rho_w_name,
  //      const std::string& rho_g_name,
		//const std::string& c_t_name,
  //      const std::string& aC_name,
  //      const std::string& bC_name)
  //  {
  //      // --- 1. Basic Validation ---
  //      if (dt <= 0.0)
  //      {
  //          std::cerr << "[IMPES_revised][TimeTerm] Error: Invalid dt provided.\n";
  //          return false;
  //      }

  //      // --- 2. Get Mesh and Fields ---
  //      auto& mesh = mgr.mesh();
  //      const auto& cells = mesh.getCells();
  //      const auto& id2idx = mesh.getCellId2Index();

  //      // Retrieve all required fields from the registry
  //      auto phi = reg.get<volScalarField>(phi_name);
  //      auto p_old = reg.get<volScalarField>(p_old_name);
  //      auto rho_t = reg.get<volScalarField>(rho_t_name);
  //      auto s_w = reg.get<volScalarField>(s_w_name);
  //      auto rho_w = reg.get<volScalarField>(rho_w_name);
  //      auto rho_g = reg.get<volScalarField>(rho_g_name);
  //      auto c_t = reg.get<volScalarField>(c_t_name);

  //      // Robust check to ensure all fields exist
  //      if (!phi || !p_old || !rho_t || !s_w || !rho_w || !rho_g || !c_t )
  //      {
  //          std::cerr << "[IMPES_revised][TimeTerm] Error: One or more required fields are missing from the registry.\n";
  //          return false;
  //      }

  //      // --- 3. Initialize Output Fields ---
  //      auto aC = reg.getOrCreate<volScalarField>(aC_name, cells.size(), 0.0);
  //      auto bC = reg.getOrCreate<volScalarField>(bC_name, cells.size(), 0.0);
  //      std::fill(aC->data.begin(), aC->data.end(), 0.0);
  //      std::fill(bC->data.begin(), bC->data.end(), 0.0);

  //      const double inv_dt = 1.0 / dt;

  //      // --- 4. Loop Over Cells and Assemble Coefficients ---
  //      for (const auto& c : cells)
  //      {
  //          if (c.id < 0) continue; // Skip inactive cells
  //          const size_t i = id2idx.at(c.id);

  //          // --- Get local cell properties with physical clamping for robustness ---
  //          const double V = std::max(0.0, c.volume);
  //          const double phi_i = std::max(0.0, std::min(1.0, (*phi)[i]));

  //          // Per the note, s_w is already clipped for residuals, but we re-clip to [0,1] for safety.
  //          const double sw_i = std::max(0.0, std::min(1.0, (*s_w)[i]));
  //          const double sg_i = 1.0 - sw_i;

  //          const double rho_w_i = std::max(0.0, (*rho_w)[i]);
  //          const double rho_g_i = std::max(0.0, (*rho_g)[i]);

  //          /// calculated by Multiphase::updateRockTwoPhaseProperties_IMPES
  //          const double ct_i = std::max(0.0, (*c_t)[i]);
  //          const double rho_t_i = std::max(1e-9, (*rho_t)[i]);

  //          // --- Discretize the Time Term: V * [φ * ρ_t * c_t] * (p_new - p_old) / dt ---
  //          // The coefficient multiplying p_new goes into the matrix diagonal (aC).
  //          // The rest goes into the source vector (bC).

  //          const double accum_coeff = V * phi_i * rho_t_i * ct_i * inv_dt;

  //          const double a = accum_coeff;
  //          const double b = accum_coeff * (*p_old)[i];

  //          (*aC)[i] = a;
  //          (*bC)[i] = b;
  //      }

  //      return true;
  //  }



        /**
     * @brief Assemble IMPES pressure accumulation term for two-phase flow
     *        using total density ρ_t and total compressibility c_t.
     *
     * Discrete form mimics the single-phase TimeTerm_IMPES_Pressure in
     * 《时间项参考.md》, but with:
     *   ρ -> ρ_t = Sw ρ_w + Sg ρ_g
     *   dρ/dp -> dρ_t/dp = ρ_t (c_t - c_r)
     *
     * For each cell i:
     *   aC[i] = V_i / dt * φ_i * ( dρ_t/dp |_⋆ + ρ_t,⋆ * c_r,i )
     *   bC[i] = V_i / dt * ( φ_i ρ_t^n
     *                      - φ_i ρ_t,⋆
     *                      + φ_i dρ_t/dp |_⋆ p_w,⋆
     *                      + φ_i ρ_t,⋆ c_r,i p_w^n )
     *
     * where:
     *   - ⋆ denotes evaluation at current outer iteration (rho_eval, p_eval)
     *   - n denotes previous time level (rho_old, p_old).
     */
inline bool TimeTerm_IMPES_Pressure(
    MeshManager& mgr,
    FieldRegistry& reg,
    double dt,
    const std::string& phi_name,          ///< porosity φ_r
    const std::string& p_old_name,        ///< previous time p_w^n
    const std::string& p_eval_name,       ///< current evaluated p_w
    const std::string& rho_old_name,      ///< previous time ρ_t^n
    const std::string& rho_eval_name,     ///< current evaluated ρ_t
    const std::string& c_t_name,          ///< total compressibility c_t
    const std::string& c_r_name,          ///< rock compressibility c_r
    const std::string& aC_name,           ///< output: diagonal coefficients
    const std::string& bC_name,           ///< output: RHS contributions
    const std::vector<char>* strong_mask = nullptr)
    {
    if (dt <= 0.0) {
        std::cerr << "[IMPES_revised][TimeTerm] invalid dt.\n";
        return false;
    }

    Mesh& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();

    auto phi = reg.get<volScalarField>(phi_name);
    auto p_old = reg.get<volScalarField>(p_old_name);
    auto p_eval = reg.get<volScalarField>(p_eval_name);
    auto rho_old = reg.get<volScalarField>(rho_old_name);
    auto rho_eval = reg.get<volScalarField>(rho_eval_name);
    auto c_t = reg.get<volScalarField>(c_t_name);
    auto c_r = reg.get<volScalarField>(c_r_name);

    if (!phi || !p_old || !p_eval || !rho_old || !rho_eval || !c_t || !c_r) {
        std::cerr << "[IMPES_revised][TimeTerm] missing fields for pressure accumulation.\n";
        return false;
    }

    auto aC = reg.getOrCreate<volScalarField>(aC_name.c_str(), cells.size(), 0.0);
    auto bC = reg.getOrCreate<volScalarField>(bC_name.c_str(), cells.size(), 0.0);
    std::fill(aC->data.begin(), aC->data.end(), 0.0);
    std::fill(bC->data.begin(), bC->data.end(), 0.0);

    const double inv_dt = 1.0 / dt;

    for (const auto& c : cells) {
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);
        if (strong_mask && (*strong_mask)[i]) continue; // strong Dirichlet cell

        const double V = std::max(0.0, c.volume);
        const double phi_i = std::max(0.0, std::min(1.0, (*phi)[i]));
        const double rho_n = std::max(0.0, (*rho_old)[i]);
        const double rho_st = std::max(0.0, (*rho_eval)[i]);
        const double ct_i = std::max(0.0, (*c_t)[i]);
        const double cr_i = std::max(0.0, (*c_r)[i]);
        const double p_n = (*p_old)[i];
        const double p_st = (*p_eval)[i];

        // total density derivative: dρ_t/dp = ρ_t (c_t - c_r)
        const double ct_minus_cr = std::max(ct_i - cr_i, 0.0);
        const double drho_dp_tot = rho_st * ct_minus_cr;

        const double a = V * inv_dt * phi_i * (drho_dp_tot + rho_st * cr_i);
        const double b = V * inv_dt * (
            phi_i * rho_n
            - phi_i * rho_st
            + phi_i * drho_dp_tot * p_st
            + phi_i * rho_st * cr_i * p_n
            );

        (*aC)[i] = a;
        (*bC)[i] = b;
    }

    return true;
    }
} // namespace IMPES_revised


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
        MeshManager& mgr,
        FieldRegistry& reg,
        double dt,
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
        auto phi = reg.get<volScalarField>(TwoPhase::Rock().phi_tag);
        auto s_w = reg.get<volScalarField>("s_w");

		auto rho_w_eval = reg.get<volScalarField>(TwoPhase::Water().rho_tag);
		auto rho_w_old = reg.get<volScalarField>(TwoPhase::Water().rho_old_tag);
        auto drho_w_dp = reg.get<volScalarField>(TwoPhase::Water().drho_w_dp_tag);

        auto rho_g_eval = reg.get<volScalarField>(TwoPhase::CO2().rho_tag);
		auto rho_g_old = reg.get<volScalarField>(TwoPhase::CO2().rho_old_tag);
		auto drho_g_dp = reg.get<volScalarField>(TwoPhase::CO2().drho_g_dp_tag);

        auto c_r = reg.get<volScalarField>(TwoPhase::Rock().c_r_tag);

		if (!p_old || !p_eval || !phi || !s_w || !rho_w_eval || !rho_w_old || !drho_w_dp || !rho_g_eval || !rho_g_old || !drho_g_dp||!c_r)
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
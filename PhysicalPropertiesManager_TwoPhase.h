#pragma once
#include <algorithm>
#include <cmath>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "InitConfig.h"					// Where VGParams and RelPermParams are defined
#include "CapRelPerm_Simple.h"         // Simplified Pc/kr model (linear Pc + Corey kr)
#include "Solver_TimeLoopUtils.h"		// CopyField here
#include "SolverContrlStrName.h"
#include "PhysicalPropertiesManager.h"
#include "CO2_Properties.h"
#include "Water_Properties.h"
#include "RockSolidProperties.h"

namespace TwoPhase
{
	inline double clamp_local(double v, double lo, double hi) { return std::max(lo, std::min(hi, v)); }
	/*
	* To ensure the CO2 properties fields exist in the FieldRegistry
	*/
	inline void ensureCO2inRockFields(FieldRegistry& reg, std::size_t n)
	{
		CO2_Prop::ensure_CO2Prop_Fields(reg, n);
	}
	/**
	* To ensure the Water properties fields exist in the FieldRegistry
	*/
	inline void ensureWaterinRockFields(FieldRegistry& reg, std::size_t n)
	{
		Water_Prop::ensure_WaterProp_Fields(reg, n);
	}
	/**
	* To ensure the Rock properties fields exist in the FieldRegistry
	*/
	inline void ensureRockFields(FieldRegistry& reg, std::size_t n)
	{
		rock::ensure_RockProp_Fields(reg, n);
	}
	/**
	* To ensure the Auxiliary parameters fields exist in the FieldRegistry
	*/
	inline void ensureAuxiliaryParametersFields(FieldRegistry& reg, std::size_t n)
	{
		// Ensure Auxiliary parameters fields exist 
		// If they do not exist, create them with default constant values
		
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::TwoPhase_case().Pc_tag, n, 0.0);// Capillary pressure
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::TwoPhase_case().dPc_dSw_tag, n, 0.0);// Derivative of capillary pressure to water saturation
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::TwoPhase_case().lambda_mass_tag, n, 0.0);// Mobility mass
	}
	/**
	*  To calculate the water basic properties at a certain step to get a output step  计算当前时刻水的基本物性参数   // old层的物性不需要计算  在进入迭代层前将计算后的基本物性参数复制存储到*_old中
	*/

	inline void computeWaterBasicPropertiesAtStep(MeshManager& mgr, FieldRegistry& reg, const std::string& p_w_inputfield, const std::string& T_inputfield)
	{

		Water_Prop::compute_water_properties_const(mgr, reg, p_w_inputfield, T_inputfield);
	}
	/**
	*  To calculate the CO2 basic properties at a certain step to get a output step  计算当前时刻CO2的基本物性参数   // old层的物性不需要计算  在进入迭代层前将计算后的基本物性参数复制存储到*_old中
	*/

	inline void computeCO2BasicPropertiesAtStep(MeshManager& mgr, FieldRegistry& reg, const std::string& p_g_inputfield, const std::string& T_inputfield)
	{

		CO2_Prop::compute_CO2_properties_const(mgr, reg, p_g_inputfield, T_inputfield);
	}
	/**
	* 
	* To calculate the Two-Phase properties of Water at a certain step,just calculate by sw
	* 
	*/
	inline void computerWaterandCO2TwoPhasePropertiesAtTimeStep(MeshManager& mgr, FieldRegistry& reg,const std::string& Sw_field, const VGParams& vg, const RelPermParams& rp)
	{
		// ===== Step 1: Sanity Checks & Setup =====
		if (!SimpleCapRelPerm::params_valid(vg)) 
		{
			std::cerr << "ERROR [multiPhase::Rock]: Invalid Van Genuchten parameters provided." << std::endl;
			return ;
		}

		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		ensureAuxiliaryParametersFields(reg, n); 

		// ===== Step 2: Get All Required INPUT Fields =====
		// Primary variable:
		auto s_w = reg.get<volScalarField>(Sw_field);
		if (!s_w)
		{
			std::cerr << "ERROR [multiPhase::Rock]: Missing one or more required fields for property update." << std::endl;
			return;
		}

		// ===== Step 3: Get or Create All OUTPUT Fields =====
		auto k_rw = reg.get<volScalarField>(PhysicalProperties_string::Water().k_rw_tag);
		auto k_rg = reg.get<volScalarField>(PhysicalProperties_string::CO2().k_rg_tag);
		auto Pc = reg.get<volScalarField>(PhysicalProperties_string::TwoPhase_case().Pc_tag);
		auto dPc_dSw = reg.get<volScalarField>(PhysicalProperties_string::TwoPhase_case().dPc_dSw_tag);

		for (size_t ic = 0; ic < n; ++ic)
		{
			const auto& c = cells[ic];
			if (c.id < 0) continue;
			const size_t i = mesh.getCellId2Index().at(c.id);
			// Relative permeability of water
			const double sw_val = (*s_w)[i];
			(*Pc)[i] = pc_vG(sw_val, vg);
			kr_Mualem_vG(sw_val, vg, rp, (*k_rw)[i], (*k_rg)[i]);
			(*dPc_dSw)[i] = dpc_dSw_vG(sw_val, vg);
		}
	}
	//=============================================================封装后的函数===================================================================//
	/**
	*  在进行时间步推进前，基于当前时刻的Sw计算得到两相属性参数，即 Pc, krw, krg, dPc_dSw 
	*/
	
	inline void updateTwoPhasePropertiesAtTimeStep(MeshManager& mgr, FieldRegistry& reg, const std::string& Sw_field, const VGParams& vg, const RelPermParams& rp)
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		ensureAuxiliaryParametersFields(reg, n); // Ensure Auxiliary parameters fields exist and get Pc and dPc_dSw
		ensureCO2inRockFields(reg, n); // Ensure CO2 fields exist and get
		ensureWaterinRockFields(reg, n); // Ensure Water fields exist and get
		computerWaterandCO2TwoPhasePropertiesAtTimeStep(mgr, reg, Sw_field, vg, rp);
	}

	/**
	*	基于已知的两相渗流参数，并根据当前的主变量压力和温度，计算当前时刻水的基本物性参数，在迭代IMPES中，除了在迭代层内，组装压力方程的系数矩阵前调用该函数进行计算，计算rho_w, mu_w, drho_wdp。需要注意的是，lambda_w也是在随着迭代发生变化的，因为lambda_w的计算涉及到黏度mu_w。
	* 还需要注意的是，这个函数只计算当前时刻的基本物性参数，不涉及old层的存储，old层的存储在进入迭代层前进行。 是迭代IMPES的典型特征
	*/

	inline void updateWaterBasicPropertiesAtStep(MeshManager& mgr, FieldRegistry& reg, const std::string& p_w_field, const std::string& T_field)
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		ensureWaterinRockFields(reg, n); // Ensure Water fields exist and get
		computeWaterBasicPropertiesAtStep(mgr, reg,  p_w_field, T_field);
	}

	inline void updateCO2BasicPropertiesAtStep(MeshManager& mgr, FieldRegistry& reg, const std::string& p_g_field, const std::string& T_field)
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		ensureCO2inRockFields(reg, n); // Ensure CO2 fields exist and get
		computeCO2BasicPropertiesAtStep(mgr, reg, p_g_field, T_field);
	}

	inline bool copyBasicPropertiesToOldLayer(FieldRegistry& reg)
	{
		// Water properties
		
		auto rho_w = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_tag);
		auto mu_w = reg.get<volScalarField>(PhysicalProperties_string::Water().mu_tag);
		auto cp_w = reg.get<volScalarField>(PhysicalProperties_string::Water().cp_tag);
		auto k_w = reg.get<volScalarField>(PhysicalProperties_string::Water().k_tag);
		auto c_w = reg.get<volScalarField>(PhysicalProperties_string::Water().c_w_tag);
		auto drho_w_dp = reg.get<volScalarField>(PhysicalProperties_string::Water().drho_w_dp_tag);

		auto rho_w_old = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_old_tag);
		auto mu_w_old = reg.get<volScalarField>(PhysicalProperties_string::Water().mu_old_tag);
		auto cp_w_old = reg.get<volScalarField>(PhysicalProperties_string::Water().cp_old_tag);
		auto k_w_old = reg.get<volScalarField>(PhysicalProperties_string::Water().k_old_tag);
		auto c_w_old = reg.get<volScalarField>(PhysicalProperties_string::Water().c_w_old_tag);
		auto drho_w_dp_old = reg.get<volScalarField>(PhysicalProperties_string::Water().drho_w_dp_old_tag);

		
		if (!rho_w || !mu_w || !drho_w_dp || !rho_w_old || !mu_w_old || !drho_w_dp_old || !cp_w || !k_w || !c_w || !cp_w_old || !k_w_old || !c_w_old)
		{
			std::cerr << "[TwoPhase][Water] Missing water property fields for copying to old layer.\n";
			return false;
		}
		copyField(reg, PhysicalProperties_string::Water().rho_tag, PhysicalProperties_string::Water().rho_old_tag);
		copyField(reg, PhysicalProperties_string::Water().mu_tag, PhysicalProperties_string::Water().mu_old_tag);
		copyField(reg, PhysicalProperties_string::Water().cp_tag, PhysicalProperties_string::Water().cp_old_tag);
		copyField(reg, PhysicalProperties_string::Water().k_tag, PhysicalProperties_string::Water().k_old_tag);
		copyField(reg, PhysicalProperties_string::Water().c_w_tag, PhysicalProperties_string::Water().c_w_old_tag);
		copyField(reg, PhysicalProperties_string::Water().drho_w_dp_tag, PhysicalProperties_string::Water().drho_w_dp_old_tag);

		// CO2 properties
		auto rho_g = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_tag);
		auto mu_g = reg.get<volScalarField>(PhysicalProperties_string::CO2().mu_tag);
		auto cp_g = reg.get<volScalarField>(PhysicalProperties_string::CO2().cp_tag);
		auto k_g = reg.get<volScalarField>(PhysicalProperties_string::CO2().k_tag);
		auto c_g = reg.get<volScalarField>(PhysicalProperties_string::CO2().c_g_tag);
		auto drho_g_dp = reg.get<volScalarField>(PhysicalProperties_string::CO2().drho_g_dp_tag);

		auto rho_g_old = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_old_tag);
		auto mu_g_old = reg.get<volScalarField>(PhysicalProperties_string::CO2().mu_old_tag);
		auto cp_g_old = reg.get<volScalarField>(PhysicalProperties_string::CO2().cp_old_tag);
		auto k_g_old = reg.get<volScalarField>(PhysicalProperties_string::CO2().k_old_tag);
		auto c_g_old = reg.get<volScalarField>(PhysicalProperties_string::CO2().c_g_old_tag);
		auto drho_g_dp_old = reg.get<volScalarField>(PhysicalProperties_string::CO2().drho_g_dp_old_tag);

		if (!rho_g || !mu_g || !drho_g_dp || !rho_g_old || !mu_g_old || !drho_g_dp_old || !cp_g || !k_g || !c_g || !cp_g_old || !k_g_old || !c_g_old)
		{
			std::cerr << "[TwoPhase][CO2] Missing CO2 property fields for copying to old layer.\n";
			return false;
		}
		copyField(reg, PhysicalProperties_string::CO2().rho_tag, PhysicalProperties_string::CO2().rho_old_tag);
		copyField(reg, PhysicalProperties_string::CO2().mu_tag, PhysicalProperties_string::CO2().mu_old_tag);
		copyField(reg, PhysicalProperties_string::CO2().cp_tag, PhysicalProperties_string::CO2().cp_old_tag);
		copyField(reg, PhysicalProperties_string::CO2().k_tag, PhysicalProperties_string::CO2().k_old_tag);
		copyField(reg, PhysicalProperties_string::CO2().c_g_tag, PhysicalProperties_string::CO2().c_g_old_tag);
		copyField(reg, PhysicalProperties_string::CO2().drho_g_dp_tag, PhysicalProperties_string::CO2().drho_g_dp_old_tag);

		// Rock properties (porosity only)
		auto phi = reg.get<volScalarField>(PhysicalProperties_string::Rock().phi_tag);
		auto phi_old = reg.get<volScalarField>(PhysicalProperties_string::Rock().phi_old_tag);
		if (!phi || !phi_old)
		{
			std::cerr << "[TwoPhase][Rock] Missing porosity fields for copying to old layer.\n";
			return false;
		}
		copyField(reg, PhysicalProperties_string::Rock().phi_tag, PhysicalProperties_string::Rock().phi_old_tag);

		return true;
	}

	/// \brief Update rock porosity `phi` from `phi_old` using c_r and pressure change.
	///
	/// Uses: phi^{n+1} = phi^{n} * exp(c_r * (p^{n+1} - p^{n})).
	inline bool updateRockPorosityFromCompressibilityAtStep(
		MeshManager& mgr,
		FieldRegistry& reg,
		const std::string& p_old_field,
		const std::string& p_field)
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const auto& id2idx = mesh.getCellId2Index();

		auto p_old = reg.get<volScalarField>(p_old_field);
		auto p = reg.get<volScalarField>(p_field);
		auto phi = reg.get<volScalarField>(PhysicalProperties_string::Rock().phi_tag);
		auto phi_old = reg.get<volScalarField>(PhysicalProperties_string::Rock().phi_old_tag);
		auto c_r = reg.get<volScalarField>(PhysicalProperties_string::Rock().c_r_tag);

		if (!p_old || !p || !phi || !phi_old || !c_r)
		{
			std::cerr << "[TwoPhase][Rock] Missing fields for porosity update (p/phi/c_r).\n";
			return false;
		}

		const double phi_min = 1e-12;
		const double phi_max = 1.0 - 1e-12;

		for (const auto& c : cells)
		{
			if (c.id < 0) continue;
			const size_t i = id2idx.at(c.id);

			const double phi_n = std::max((*phi_old)[i], phi_min);
			const double cr = (*c_r)[i];
			const double dp = (*p)[i] - (*p_old)[i];

			const double phi_np1 = clamp_local(phi_n * std::exp(cr * dp), phi_min, phi_max);
			(*phi)[i] = phi_np1;
		}

		return true;
	}

	inline void calculateteTotalMobilityField(MeshManager& mgr, FieldRegistry& reg)
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		auto lambda_mass = reg.getOrCreate<volScalarField>(PhysicalProperties_string::TwoPhase_case().lambda_mass_tag, n, 0.0);
		auto lambda_w = reg.get<volScalarField>(PhysicalProperties_string::Water().lambda_w_tag);
		auto lambda_g = reg.get<volScalarField>(PhysicalProperties_string::CO2().lambda_g_tag);
		auto rho_w = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_tag);
		auto rho_g = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_tag);
		if (!lambda_w || !lambda_g || !rho_w || !rho_g)
		{
			std::cerr << "[TwoPhase] Missing fields for calculating total mobility.\n";
			return;
		}
		for (size_t ic = 0; ic < n; ++ic)
		{
			const auto& c = cells[ic];
			const size_t i = mesh.getCellId2Index().at(c.id);
			(*lambda_mass)[i] = (*lambda_w)[i] * (*rho_w)[i] + (*lambda_g)[i] * (*rho_g)[i];
		}
	}


}

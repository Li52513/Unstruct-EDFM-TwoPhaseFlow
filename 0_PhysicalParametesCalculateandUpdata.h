#pragma once
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "InitConfig.h"         // Where VGParams and RelPermParams are defined
#include "CapRelPerm_Simple.h"         // Simplified Pc/kr model (linear Pc + Corey kr)
#include "Solver_TimeLoopUtils.h" // CopyField here


namespace TwoPhase
{
	inline double clamp_local(double v, double lo, double hi) { return std::max(lo, std::min(hi, v)); }

	/**
	*  To define the properties tag of water phase
	*/
	struct Water
	{
		// Basic properties
		std::string rho_tag = "rho_w";
		std::string mu_tag = "mu_w";
		std::string drho_w_dp_tag = "drho_wdp";		//used for Time term of Pressure equation in IMPES scheme to replace constant compressibility
		
		std::string rho_old_tag = "rho_w_old";
		std::string mu_old_tag = "mu_w_old";
		std::string drho_w_dp_old_tag = "drho_wdp_old";		
	
		// Two-phase properties
		std::string k_rw_tag = "k_rw";				//relative permeability of water
		std::string lambda_w_tag = "lambda_w";		//mobility of water
	};
	/**
	*  To define the properties tag of CO2 phase
	*/
	struct CO2
	{
		// Basic properties
		std::string rho_tag = "rho_g";
		std::string mu_tag = "mu_g";
		std::string drho_g_dp_tag = "drho_gdp";		//used for Time term of Pressure equation in IMPES scheme to replace constant compressibility

		std::string rho_old_tag = "rho_g_old";
		std::string mu_old_tag = "mu_g_old";
		std::string drho_g_dp_old_tag = "drho_gdp_old";

		// Two-phase properties
		std::string k_rg_tag = "k_rg";				//relative permeability of CO2
		std::string lambda_g_tag = "lambda_g";		//mobility of CO2
	};
	/**
	*  To define the properties tag of Rock
	*/
	struct Rock
	{
		// Bsic properties
		std::string phi_tag = "phi_r";				//porosity
		std::string c_r_tag = "c_r";				//rock compressibility
	};
	/**
	*  To define the properties tag of Two-Phase system
	*/
	struct Auxiliaryparameters
	{
		std::string Pc_tag = "Pc";					//capillary pressure
		std::string dPc_dSw_tag = "dPc_dSw";		//derivative of capillary pressure to water saturation
		std::string lambda_mass_tag = "lambda_mass";// equation : lambda_mass = lambda_w * rho_w + lambda_g * rho_g
	};

	//=====================================================================================================================================================//
	/*
	* To ensure the CO2 properties fields exist in the FieldRegistry
	*/
	inline void ensureCO2inRockFields(FieldRegistry& reg, std::size_t n)
	{
		// Ensure CO2 properties fields exist 
		// If they do not exist, create them with default constant values
		//Basic properties
		
		
		reg.getOrCreate<volScalarField>(CO2().rho_tag, n, 800.0);// Density of CO2 
		reg.getOrCreate<volScalarField>(CO2().mu_tag, n, 1.48e-5);// Viscosity of CO2
		reg.getOrCreate<volScalarField>(CO2().drho_g_dp_tag, n, 0.0);// Derivative of density of CO2 to pressure

		reg.getOrCreate<volScalarField>(CO2().rho_old_tag, n, 800.0);// Density of CO2 old
		reg.getOrCreate<volScalarField>(CO2().mu_old_tag, n, 1.48e-5);// Viscosity of CO2 old
		reg.getOrCreate<volScalarField>(CO2().drho_g_dp_old_tag, n, 0.0);// Derivative of density of CO2 to pressure old
		
		// Two-phase properties
		reg.getOrCreate<volScalarField>(CO2().k_rg_tag, n, 0.0);// Relative permeability of CO2
		reg.getOrCreate<volScalarField>(CO2().lambda_g_tag, n, 0.0);// Mobility of CO2
	}
	/**
	* To ensure the Water properties fields exist in the FieldRegistry
	*/
	inline void ensureWaterinRockFields(FieldRegistry& reg, std::size_t n)
	{
		// Ensure Water properties fields exist 
		// If they do not exist, create them with default constant values
		//Basic properties
		reg.getOrCreate<volScalarField>(Water().rho_tag, n, 1000.0);// Density of Water 
		reg.getOrCreate<volScalarField>(Water().mu_tag, n, 1.0e-3);// Viscosity of Water
		reg.getOrCreate<volScalarField>(Water().drho_w_dp_tag, n, 0.0);// Derivative of density of Water to pressure

		reg.getOrCreate<volScalarField>(Water().rho_old_tag, n, 1000.0);// Density of Water old
		reg.getOrCreate<volScalarField>(Water().mu_old_tag, n, 1.0e-3);// Viscosity of Water old
		reg.getOrCreate<volScalarField>(Water().drho_w_dp_old_tag, n, 0.0);// Derivative of density of Water to pressure old
		
		// Two-phase properties
		reg.getOrCreate<volScalarField>(Water().k_rw_tag, n, 0.0);// Relative permeability of Water
		reg.getOrCreate<volScalarField>(Water().lambda_w_tag, n, 0.0);// Mobility of Water
	}
	/**
	* To ensure the Rock properties fields exist in the FieldRegistry
	*/
	inline void ensureRockFields(FieldRegistry& reg, std::size_t n)
	{
		// Ensure Rock properties fields exist 
		// If they do not exist, create them with default constant values
		auto Rock_props = Rock();
		reg.getOrCreate<volScalarField>(Rock_props.phi_tag, n, 0.2);// Porosity
		reg.getOrCreate<volScalarField>(Rock_props.c_r_tag, n, 1.0e-9);// Rock compressibility
	}
	/**
	* To ensure the Auxiliary parameters fields exist in the FieldRegistry
	*/
	inline void ensureAuxiliaryParametersFields(FieldRegistry& reg, std::size_t n)
	{
		// Ensure Auxiliary parameters fields exist 
		// If they do not exist, create them with default constant values
		auto Auxi_params = Auxiliaryparameters();
		reg.getOrCreate<volScalarField>(Auxi_params.Pc_tag, n, 0.0);// Capillary pressure
		reg.getOrCreate<volScalarField>(Auxi_params.dPc_dSw_tag, n, 0.0);// Derivative of capillary pressure to water saturation
		reg.getOrCreate<volScalarField>(Auxi_params.lambda_mass_tag, n, 0.0);// Mobility mass
	}
	/**
	*  To calculate the water basic properties at a certain step to get a output step  计算当前时刻水的基本物性参数   // old层的物性不需要计算  在进入迭代层前将计算后的基本物性参数复制存储到*_old中
	*/

	inline void computeWaterBasicPropertiesAtStep(MeshManager& mgr, FieldRegistry& reg, const std::string& p_w_inputfield, const std::string& T_inputfield)
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		auto TF = reg.get<volScalarField>(T_inputfield);
		auto p_wF = reg.get<volScalarField>(p_w_inputfield);
		if (!TF || !p_wF)
		{
			std::cerr << "[TwoPhase][Water] Missing temperature or pressure field for water property calculation.\n";
			return;
		}
		auto rho_w = reg.get<volScalarField>(Water().rho_tag);
		auto mu_w = reg.get<volScalarField>(Water().mu_tag);
		auto drho_w_dp = reg.get<volScalarField>(Water().drho_w_dp_tag);
		auto lambda_w = reg.get<volScalarField>(Water().lambda_w_tag);
		auto k_rw = reg.get<volScalarField>(Water().k_rw_tag);
		if (!rho_w || !mu_w || !drho_w_dp)
		{
			std::cerr << "[TwoPhase][Water] Missing water property fields for calculation at time step.\n";
			return;
		}

		for (size_t ic = 0; ic < n; ++ic)
		{
			const auto& c = cells[ic];
			if (c.id < 0) continue;
			const size_t i = mesh.getCellId2Index().at(c.id);
			double T = (*TF)[i];
			double p_w = (*p_wF)[i];

			//当前基本物性参数参数为常数
			(*rho_w)[i] = 1000.0; //kg/m3
			(*mu_w)[i] = 1.0e-3; //Pa.s
			(*drho_w_dp)[i] = 0.0; //kg/m3.Pa
			
			//由于流度的计算涉及黏度，所以这里也需要更新流度,注意这里的相对渗透率采用的是在进入迭代层前计算好的值
			(*lambda_w)[i] = (*k_rw)[i] / std::max((*mu_w)[i], 1e-20);

		}
	}
	/**
	*  To calculate the CO2 basic properties at a certain step to get a output step  计算当前时刻CO2的基本物性参数   // old层的物性不需要计算  在进入迭代层前将计算后的基本物性参数复制存储到*_old中
	*/

	inline void computeCO2BasicPropertiesAtStep(MeshManager& mgr, FieldRegistry& reg, const std::string& p_g_inputfield, const std::string& T_inputfield)
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		auto TF = reg.get<volScalarField>(T_inputfield);
		auto p_gF = reg.get<volScalarField>(p_g_inputfield);
		if (!TF || !p_gF)
		{
			std::cerr << "[TwoPhase][CO2] Missing temperature or pressure field for water property calculation.\n";
			return;
		}
		auto rho_g = reg.get<volScalarField>(CO2().rho_tag);
		auto mu_g = reg.get<volScalarField>(CO2().mu_tag);
		auto drho_g_dp = reg.get<volScalarField>(CO2().drho_g_dp_tag);
		auto lambda_g = reg.get<volScalarField>(CO2().lambda_g_tag);
		auto k_rg = reg.get<volScalarField>(CO2().k_rg_tag);
		if (!rho_g || !mu_g || !drho_g_dp)
		{
			std::cerr << "[TwoPhase][CO2] Missing water property fields for calculation at time step.\n";
			return;
		}

		for (size_t ic = 0; ic < n; ++ic)
		{
			const auto& c = cells[ic];
			const size_t i = mesh.getCellId2Index().at(c.id);
			if (c.id < 0) continue;
			double T = (*TF)[i];
			double p_g = (*p_gF)[i];

			//当前基本物性参数参数为常数
			(*rho_g)[i] = 800.0; //kg/m3
			(*mu_g)[i] =	5e-3; //Pa.s
			(*drho_g_dp)[i] = 0.0; //kg/m3.Pa

			//由于流度的计算涉及黏度，所以这里也需要更新流度,注意这里的相对渗透率采用的是在进入迭代层前计算好的值
			(*lambda_g)[i] = (*k_rg)[i] / std::max((*mu_g)[i], 1e-20);

		}
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

		// ===== Step 2: Get All Required INPUT Fields =====
		// Primary variable:
		auto s_w = reg.get<volScalarField>(Sw_field);
		if (!s_w)
		{
			std::cerr << "ERROR [multiPhase::Rock]: Missing one or more required fields for property update." << std::endl;
			return;
		}

		// ===== Step 3: Get or Create All OUTPUT Fields =====
		auto k_rw = reg.get<volScalarField>(Water().k_rw_tag);
		auto k_rg = reg.get<volScalarField>(CO2().k_rg_tag);
		auto Pc = reg.get<volScalarField>(Auxiliaryparameters().Pc_tag);
		auto dPc_dSw = reg.get<volScalarField>(Auxiliaryparameters().dPc_dSw_tag);

		for (size_t ic = 0; ic < n; ++ic)
		{
			const auto& c = cells[ic];
			if (c.id < 0) continue;
			const size_t i = mesh.getCellId2Index().at(c.id);
			// Relative permeability of water
			double sw_val = (*s_w)[i];
			double se = SimpleCapRelPerm::Se_from_Sw_simple(sw_val, vg);
			//主变量-已经在初场初始化的时候进行了初始化 															// Capillary pressure (linear Pc = pc_entry*(1-Se))
			(*Pc)[i] = SimpleCapRelPerm::pc_linear(se, vg);																	// Capillary pressure (linear Pc = pc_entry*(1-Se))
			SimpleCapRelPerm::kr_corey(se, vg, SimpleCapRelPerm::RelPermConfig{}, (*k_rw)[i], (*k_rg)[i]);					// Relative permeabilities and derivative of capillary pressure (Corey + linear Pc)
			(*dPc_dSw)[i] = SimpleCapRelPerm::dpc_dSw_linear(se, vg);
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
		auto Water_props = Water();
		auto rho_w = reg.get<volScalarField>(Water_props.rho_tag);
		auto mu_w = reg.get<volScalarField>(Water_props.mu_tag);
		auto drho_w_dp = reg.get<volScalarField>(Water_props.drho_w_dp_tag);


		auto rho_w_old = reg.get<volScalarField>(Water_props.rho_old_tag);
		auto mu_w_old = reg.get<volScalarField>(Water_props.mu_old_tag);
		auto drho_w_dp_old = reg.get<volScalarField>(Water_props.drho_w_dp_old_tag);
		
		if (!rho_w || !mu_w || !drho_w_dp || !rho_w_old || !mu_w_old || !drho_w_dp_old)
		{
			std::cerr << "[TwoPhase][Water] Missing water property fields for copying to old layer.\n";
			return false;
		}
		copyField(reg, Water_props.rho_tag, Water_props.rho_old_tag);
		copyField(reg, Water_props.mu_tag, Water_props.mu_old_tag);
		copyField(reg, Water_props.drho_w_dp_tag, Water_props.drho_w_dp_old_tag);

		// CO2 properties
		auto CO2_props = CO2();
		auto rho_g = reg.get<volScalarField>(CO2_props.rho_tag);
		auto mu_g = reg.get<volScalarField>(CO2_props.mu_tag);
		auto drho_g_dp = reg.get<volScalarField>(CO2_props.drho_g_dp_tag);

		auto rho_g_old = reg.get<volScalarField>(CO2_props.rho_old_tag);
		auto mu_g_old = reg.get<volScalarField>(CO2_props.mu_old_tag);
		auto drho_g_dp_old = reg.get<volScalarField>(CO2_props.drho_g_dp_old_tag);

		if (!rho_g || !mu_g || !drho_g_dp || !rho_g_old || !mu_g_old || !drho_g_dp_old)
		{
			std::cerr << "[TwoPhase][CO2] Missing CO2 property fields for copying to old layer.\n";
			return false;
		}
		copyField(reg, CO2_props.rho_tag, CO2_props.rho_old_tag);
		copyField(reg, CO2_props.mu_tag, CO2_props.mu_old_tag);
		copyField(reg, CO2_props.drho_g_dp_tag, CO2_props.drho_g_dp_old_tag);
		return true;
	}

	inline void calculateteTotalMobilityField(MeshManager& mgr, FieldRegistry& reg)
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		auto lambda_mass = reg.getOrCreate<volScalarField>(Auxiliaryparameters().lambda_mass_tag, n, 0.0);
		auto lambda_w = reg.get<volScalarField>(Water().lambda_w_tag);
		auto lambda_g = reg.get<volScalarField>(CO2().lambda_g_tag);
		auto rho_w = reg.get<volScalarField>(Water().rho_tag);
		auto rho_g = reg.get<volScalarField>(CO2().rho_tag);
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

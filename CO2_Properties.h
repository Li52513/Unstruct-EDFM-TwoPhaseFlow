#pragma once
#include "PropertiesSummary.h" // 里面定义 CO2所需要的物性参数
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "SolverContrlStrName.h"
#include "CO2PropertyTable.h"

namespace CO2_Prop
{
	static constexpr double BASE_rho_g = 800;
	static constexpr double BASE_cp_g = 1200;
	static constexpr double BASE_k_g = 0.05;
	static constexpr double BASE_mu_g = 1.48e-5;
	static constexpr double BASE_Drho_Dp_g = 0;
	static constexpr double BASE_c_g = 0;

	inline bool ensure_CO2Prop_Fields(FieldRegistry& reg, std::size_t n)
	{
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().rho_tag, n, BASE_rho_g);				//rho，kg/m³
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().mu_tag, n, BASE_mu_g);				//mu，Pa·s
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().cp_tag, n, BASE_cp_g);				//Cp，J/(kg·K)
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().k_tag, n, BASE_k_g);					//k，W/(m·K)
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().drho_g_dp_tag, n, BASE_Drho_Dp_g);	//drho_dp，kg/(m³·Pa)
		reg.getOrCreate <volScalarField>(PhysicalProperties_string::CO2().c_g_tag, n, BASE_c_g);				//二氧化碳的可压缩系数，1/Pa

		//old time layer
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().rho_old_tag, n, BASE_rho_g);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().mu_old_tag, n, BASE_mu_g);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().cp_old_tag, n, BASE_cp_g);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().k_old_tag, n, BASE_k_g);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().drho_g_dp_old_tag, n, BASE_Drho_Dp_g);
		reg.getOrCreate <volScalarField>(PhysicalProperties_string::CO2().c_g_old_tag, n, BASE_c_g);
		return true;

	}
	inline bool compute_CO2_properties_const(MeshManager& mgr, FieldRegistry& reg, const std::string& p_g_field,const std::string& T_field)  //这里传入基岩的压力和温度 但是对于IMPES来说，不会更新这里的参数
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		auto TF = reg.get<volScalarField>(T_field);
		auto p_gF= reg.get<volScalarField>(p_g_field);

		auto rho_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_tag);
		auto cp_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().cp_tag);
		auto k_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().k_tag);
		auto mu_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().mu_tag);
		auto Drho_Dp_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().drho_g_dp_tag);
		auto c_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().c_g_tag);
		for (size_t ic = 0; ic < cells.size(); ++ic)
		{
			const auto& c = cells[ic];

			const size_t i = mesh.getCellId2Index().at(c.id);
			double T = (*TF)[i];
			double p_g = (*p_gF)[i];

			//当前传入Base参数，开展常物性测试，但是提供温度和压力接口
			(*rho_gF)[i] = BASE_rho_g;
			(*cp_gF)[i]  = BASE_cp_g;
			(*k_gF)[i]   = BASE_k_g;
			(*mu_gF)[i]  = BASE_mu_g;
			(*Drho_Dp_gF)[i] = BASE_Drho_Dp_g;
			(*c_gF)[i] = BASE_c_g; //假设可压缩系数为1/p
		}
		return true;
	}

	inline bool computer_CO2_propertier_Span_Wagner(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field)
	{
		// TO DO 基于Span-Wagner状态方程计算CO2的物性参数
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		ensure_CO2Prop_Fields (reg, n);
		auto TF = reg.get<volScalarField>(T_field);
		auto pF = reg.get<volScalarField>(p_field);
		if (!pF || !TF) {
			std::cerr << "[PPM][Fluid] missing fields '" << p_field << "' or '" << T_field << "'.\n";
			return false;
		}
		auto rho_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_tag);
		auto cp_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().cp_tag);
		auto k_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().k_tag);
		auto mu_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().mu_tag);

		auto gt = CO2PropertyTable::instance();
		for (size_t ic = 0; ic < cells.size(); ++ic)
		{
			const auto& c = cells[ic];
			if (c.id < 0) continue;
			const size_t i = mesh.getCellId2Index().at(c.id);
			double p = (*pF)[i], T = (*TF)[i];
			double rho = BASE_rho_g, mu = BASE_mu_g, cp = BASE_cp_g, k = BASE_k_g;
			try { const auto G = gt.getProperties(p, T); rho = G.rho; mu = G.mu; cp = G.cp; k = G.k; }
			catch (...) {}
			(*rho_gF)[i] = rho; (*mu_gF)[i] = mu; (*cp_gF)[i] = cp; (*k_gF)[i] = k;
		}
		return true;
	}

}



#pragma once
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

	///为常物性参数工况外部修改提供接口
	struct CO2Properties
	{
		double rho = 900;
		double mu = 1.48e-5;
		double cp = 1100;
		double k = 0.03;
		double dRho_dP = 0.0;
		double c = 0.0;
	};

	inline bool ensure_CO2Prop_Fields(FieldRegistry& reg, std::size_t n)
	{
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().rho_tag, n, BASE_rho_g);				//rho，kg/m³
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().mu_tag, n, BASE_mu_g);				//mu，Pa·s
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().cp_tag, n, BASE_cp_g);				//Cp，J/(kg·K)
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().k_tag, n, BASE_k_g);					//k，W/(m·K)
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().drho_g_dp_tag, n, BASE_Drho_Dp_g);	//drho_dp，kg/(m³·Pa)
		reg.getOrCreate <volScalarField>(PhysicalProperties_string::CO2().c_g_tag, n, BASE_c_g);				//二氧化碳的可压缩系数，1/Pa
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().lambda_g_tag, n, 0);					// Mobility of Water
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().k_rg_tag, n, 0);						// Relative permeability of Water

		//old time layer
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().rho_old_tag, n, BASE_rho_g);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().mu_old_tag, n, BASE_mu_g);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().cp_old_tag, n, BASE_cp_g);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().k_old_tag, n, BASE_k_g);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().drho_g_dp_old_tag, n, BASE_Drho_Dp_g);
		reg.getOrCreate <volScalarField>(PhysicalProperties_string::CO2().c_g_old_tag, n, BASE_c_g);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().lambda_g_old_tag, n, 0);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().k_rg_old_tag, n, 0);
		return true;

	}
	inline bool compute_CO2_properties_const(MeshManager& mgr, FieldRegistry& reg, const std::string& p_g_field,const std::string& T_field)  //这里传入基岩的压力和温度 但是对于IMPES来说，不会更新这里的参数
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		auto TF = reg.get<volScalarField>(T_field);
		auto p_gF= reg.get<volScalarField>(p_g_field);

		ensure_CO2Prop_Fields(reg, n);

		auto rho_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_tag);
		auto cp_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().cp_tag);
		auto k_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().k_tag);
		auto mu_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().mu_tag);
		auto Drho_Dp_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().drho_g_dp_tag);
		auto c_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().c_g_tag);
		auto k_rgF = reg.get<volScalarField>(PhysicalProperties_string::CO2().k_rg_tag);
		auto lambda_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().lambda_g_tag);

		for (size_t ic = 0; ic < cells.size(); ++ic)
		{
			const auto& c = cells[ic];

			const size_t i = mesh.getCellId2Index().at(c.id);
			double T = (*TF)[i];
			double p_g = (*p_gF)[i];

			
			(*rho_gF)[i] = CO2Properties().rho;
			(*cp_gF)[i]  = CO2Properties().cp;
			(*k_gF)[i]   = CO2Properties().k;
			(*mu_gF)[i]  = CO2Properties().mu;
			(*Drho_Dp_gF)[i] = CO2Properties().dRho_dP;
			(*c_gF)[i] = CO2Properties().c; 
			//流度计算
			(*lambda_gF)[i] = (*k_rgF)[i] / std::max((*mu_gF)[i], 1e-20);
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
		auto k_rgF = reg.get<volScalarField>(PhysicalProperties_string::CO2().k_rg_tag);
		auto lambda_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().lambda_g_tag);

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
			//流度计算
			(*lambda_gF)[i] = (*k_rgF)[i] / std::max((*mu_gF)[i], 1e-20);
		}
		return true;
	}

	// 计算CO2的压缩系数：drho/dp 在给定的 p/T 点处
	inline bool compute_CO2_DrhoDp
	(
		MeshManager& mgr, FieldRegistry& reg,
		const std::string& p_eval_name,
		const std::string& T_eval_name,
		double dp_rel = 1e-4, double dp_abs_min = 1.0
	)
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const auto& id2idx = mesh.getCellId2Index();
		if (cells.empty()) return false;

		auto pE = reg.get<volScalarField>(p_eval_name);
		auto TE = reg.get<volScalarField>(T_eval_name);
		if (!pE || !TE)
		{
			std::cerr << "[computeRhoAndDrhoDp] missing eval fields '"
				<< p_eval_name << "' or '" << T_eval_name << "'\n";
			return false;
		}
		auto dF = reg.getOrCreate<volScalarField>(PhysicalProperties_string::CO2().drho_g_dp_old_tag, cells.size(), 0.0);
		auto gt = CO2PropertyTable::instance();
		for (const auto& c : cells)
		{
			const size_t i = id2idx.at(c.id);
			double p = (*pE)[i], T = (*TE)[i];
			double rho_lin = 1000.0;
			try { rho_lin = gt.getProperties(p, T).rho; }
			catch (...) { /* 出界兜底：保持 rho_lin 默认 */ }
			// 对称差分：drho/dp ≈ [ρ(p+dp,T) - ρ(p-dp,T)] / (2dp)
			const double dpa = std::max(dp_abs_min, std::abs(p) * dp_rel);
			double rp = rho_lin, rm = rho_lin;
			try
			{
				double pp = p + dpa, pm = p - dpa;
				rp = gt.getProperties(pp, T).rho;
				rm = gt.getProperties(pm, T).rho;
			}
			catch (...) { /* 出界兜底：rp/rm 用 rho_lin */ }
			(*dF)[i] = (rp - rm) / std::max(2.0 * dpa, 1e-12);
		}
		return true;
	}

	// 计算当CO2为基岩中唯一流体时，的有效热物性参数
	inline bool computer_effectiveThermalProperties(MeshManager& mgr, FieldRegistry& reg)
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();

		ensure_CO2Prop_Fields(reg, n);

		// 岩石参数（必须存在）
		auto phiF = reg.get<volScalarField>(PhysicalProperties_string::Rock().phi_tag);
		auto rrF = reg.get<volScalarField>(PhysicalProperties_string::Rock().rho_tag);
		auto cprF = reg.get<volScalarField>(PhysicalProperties_string::Rock().cp_tag);
		auto lamrF = reg.get<volScalarField>(PhysicalProperties_string::Rock().lambda_tag);


		//流体物性参数
		auto rho_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_tag);
		auto cp_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().cp_tag);
		auto k_gF = reg.get<volScalarField>(PhysicalProperties_string::CO2().k_tag);

		// 有效热参数场
		auto Ceff = reg.getOrCreate<volScalarField>(PhysicalProperties_string::SinglePhase_case().C_eff_tag, n, 0.0);
		auto lame = reg.getOrCreate<volScalarField>(PhysicalProperties_string::SinglePhase_case().lambda_eff_tag, n, 0.0);
		for (size_t ic = 0; ic < cells.size(); ++ic)
		{
			const auto& c = cells[ic];
			if (c.id < 0) continue;
			const size_t i = mesh.getCellId2Index().at(c.id);

			(*Ceff)[i] = (1.0 - (*phiF)[i]) * (*rrF)[i] * (*cprF)[i] + (*phiF)[i] * (*rho_gF)[i] * (*cp_gF)[i];
			(*lame)[i] = (1.0 - (*phiF)[i]) * (*lamrF)[i] + (*phiF)[i] * (*k_gF)[i];
		}
		return true;
	}	

}



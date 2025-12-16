#pragma once
#include "PropertiesSummary.h" // 里面定义 水所需要的物性参数
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "SolverContrlStrName.h"
#include "WaterPropertyTable.h"

namespace Water_Prop
{
	static constexpr double BASE_rho_w = 980;
	static constexpr double BASE_cp_w = 4300;
	static constexpr double BASE_k_w = 0.05;
	static constexpr double BASE_mu_w = 3e-54;
	static constexpr double BASE_Drho_Dp_w = 0;
	static constexpr double BASE_c_w = 0; //压缩系数

	inline bool ensure_WaterProp_Fields(FieldRegistry& reg, std::size_t n)
	{
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().rho_tag, n, BASE_rho_w);				//rho，kg/m³
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().mu_tag, n, BASE_mu_w);				//mu，Pa·s
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().cp_tag, n, BASE_cp_w);				//Cp，J/(kg·K)
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().k_tag, n, BASE_k_w);					//k，W/(m·K)
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().drho_w_dp_tag, n, BASE_Drho_Dp_w);	//drho_dp，kg/(m³·Pa)
		reg.getOrCreate <volScalarField>(PhysicalProperties_string::Water().c_w_tag, n, BASE_c_w);				//水的可压缩系数，1/Pa

		//old time layer
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().rho_old_tag, n, BASE_rho_w);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().mu_old_tag, n, BASE_mu_w);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().cp_old_tag, n, BASE_cp_w);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().k_old_tag, n, BASE_k_w);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().drho_w_dp_old_tag, n, BASE_Drho_Dp_w);
		reg.getOrCreate <volScalarField>(PhysicalProperties_string::Water().c_w_old_tag, n, BASE_c_w);

		return true;
	}

	inline bool compute_water_properties_const(MeshManager& mgr, FieldRegistry& reg, const std::string& p_w_field, const std::string& T_field)  //这里传入基岩的压力和温度 但是对于IMPES来说，不会更新这里的参数
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		auto TF = reg.get<volScalarField>(T_field);
		auto p_wF = reg.get<volScalarField>(p_w_field);

		auto rho_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_tag);
		auto cp_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().cp_tag);
		auto k_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().k_tag);
		auto mu_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().mu_tag);
		auto Drho_Dp_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().drho_w_dp_tag);
		auto c_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().c_w_tag);

		for (size_t ic = 0; ic < cells.size(); ++ic)
		{
			const auto& c = cells[ic];

			const size_t i = mesh.getCellId2Index().at(c.id);
			double T = (*TF)[i];
			double p_w = (*p_wF)[i];

			//当前传入Base参数，开展常物性测试，但是提供温度和压力接口
			(*rho_wF)[i] = BASE_rho_w;
			(*cp_wF)[i] = BASE_cp_w;
			(*k_wF)[i] = BASE_k_w;
			(*mu_wF)[i] = BASE_mu_w;
			(*Drho_Dp_wF)[i] = BASE_Drho_Dp_w;
			(*c_wF)[i] = BASE_c_w;
		}
		return true;
	}

	//基于LAPWS模型，利用插值表计算水的物性参数
	inline bool computer_water_properties_LAPWS(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field)
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();

		ensure_WaterProp_Fields(reg, n);

		auto pF = reg.get<volScalarField>(p_field);
		auto TF = reg.get<volScalarField>(T_field);
		if (!pF || !TF) {
			std::cerr << "[PPM][Fluid] missing fields '" << p_field << "' or '" << T_field << "'.\n";
			return false;
		}
		auto rho_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_tag);
		auto cp_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().cp_tag);
		auto k_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().k_tag);
		auto mu_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().mu_tag);

		auto wt = WaterPropertyTable::instance();

		for (size_t ic = 0; ic < cells.size(); ++ic)
		{
			const auto& c = cells[ic];
			const size_t i = mesh.getCellId2Index().at(c.id);
			double p = (*pF)[i], T = (*TF)[i]; 
			double rho = 1000, mu = 1e-3, cp = 4200, k = 0.6;
			try { const auto W = wt.getProperties(p, T); rho = W.rho; mu = W.mu; cp = W.cp; k = W.k; }
			catch (...) {  }
			(*rho_wF)[i] = rho; (*mu_wF)[i] = mu; (*cp_wF)[i] = cp; (*k_wF)[i] = k;
		}
		return true;
	}

	// 在 (p_eval, T_eval) 处评估 ρ 与 ∂ρ/∂p
	inline bool compute_water_DrhoDpAt
	(
		MeshManager& mgr, FieldRegistry& reg,
		const std::string& p_eval_name,
		const std::string& T_eval_name,
		const std::string& drhodp_out ,
		double dp_rel = 1e-4, double dp_abs_min = 10.0
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
			std::cerr << "[computeRhoAndDrhoDpAt] missing eval fields '"
				<< p_eval_name << "' or '" << T_eval_name << "'\n";
			return false;
		}
		auto dF = reg.getOrCreate<volScalarField>(drhodp_out, cells.size(), 0.0);
		auto& wt = WaterPropertyTable::instance();
		for (const auto& c : cells)
		{
			const size_t i = id2idx.at(c.id);
			double p = (*pE)[i], T = (*TE)[i];
			double rho_lin = 1000.0;
			try { rho_lin = wt.getProperties(p, T).rho; }
			catch (...) { /* 出界兜底：保持 rho_lin 默认 */ }
			// 对称差分：drho/dp ≈ [ρ(p+dp,T) - ρ(p-dp,T)] / (2dp)
			const double dpa = std::max(dp_abs_min, std::abs(p) * dp_rel);
			double rp = rho_lin, rm = rho_lin;
			try 
			{
				double pp = p + dpa, pm = p - dpa;
				rp = wt.getProperties(pp, T).rho;
				rm = wt.getProperties(pm, T).rho;
			}
			catch (...) { /* 出界兜底：rp/rm 用 rho_lin */ }
			(*dF)[i] = (rp - rm) / std::max(2.0 * dpa, 1e-12);
		}
		return true;
	}

} //namespace Water_Prop


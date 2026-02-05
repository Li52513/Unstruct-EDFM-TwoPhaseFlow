#pragma once
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
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().lambda_w_tag, n, 0);					// Mobility of Water
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().k_rw_tag, n, 0);						// Relative permeability of Water

		//old time layer
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().rho_old_tag, n, BASE_rho_w);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().mu_old_tag, n, BASE_mu_w);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().cp_old_tag, n, BASE_cp_w);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().k_old_tag, n, BASE_k_w);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().drho_w_dp_old_tag, n, BASE_Drho_Dp_w);
		reg.getOrCreate <volScalarField>(PhysicalProperties_string::Water().c_w_old_tag, n, BASE_c_w);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().lambda_w_old_tag, n, 0);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().k_rw_old_tag, n, 0);

		return true;
	}
 	/// 计算基岩内部水的物性参数
	////常数物性参数工况
	inline bool compute_water_properties_const(MeshManager& mgr, FieldRegistry& reg, const std::string& p_w_field, const std::string& T_field)  //这里传入基岩的压力和温度 但是对于IMPES来说，不会更新这里的参数
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		auto TF = reg.get<volScalarField>(T_field);
		auto p_wF = reg.get<volScalarField>(p_w_field);

		ensure_WaterProp_Fields(reg, n);

		auto rho_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_tag);
		auto cp_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().cp_tag);
		auto k_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().k_tag);
		auto mu_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().mu_tag);
		auto Drho_Dp_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().drho_w_dp_tag);
		auto c_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().c_w_tag);
		auto lambda_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().lambda_w_tag);
		auto k_rwF = reg.get<volScalarField>(PhysicalProperties_string::Water().k_rw_tag);


		for (size_t ic = 0; ic < cells.size(); ++ic)
		{
			const auto& c = cells[ic];

			const size_t i = mesh.getCellId2Index().at(c.id);
			double T = (*TF)[i];
			double p_w = (*p_wF)[i];

			//当前传入Base参数，开展常物性测试，但是提供温度和压力接口
			(*rho_wF)[i] = WaterProperties().rho;
			(*cp_wF)[i] = WaterProperties().cp;
			(*k_wF)[i] = WaterProperties().k;
			(*mu_wF)[i] = WaterProperties().mu;
			(*Drho_Dp_wF)[i] = WaterProperties().dRho_dP;
			(*c_wF)[i] = WaterProperties().c;
			//由于流度的计算涉及黏度，所以这里也需要更新流度,注意这里的相对渗透率采用的是在进入迭代层前计算好的值
			(*lambda_wF)[i] = (*k_rwF)[i] / std::max((*mu_wF)[i], 1e-20);

		}
		return true;
	}
	////基于LAPWS模型，利用插值表计算水的物性参数
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
		auto lambda_wF = reg.get<volScalarField>(PhysicalProperties_string::Water().lambda_w_tag);
		auto k_rwF = reg.get<volScalarField>(PhysicalProperties_string::Water().k_rw_tag);

		auto& wt = WaterPropertyTable::instance();

		size_t oorCount = 0;
		for (size_t ic = 0; ic < cells.size(); ++ic)
		{
			const auto& c = cells[ic];
			if (c.id < 0) continue;

			const size_t i = mesh.getCellId2Index().at(c.id);
			double p = (*pF)[i], T = (*TF)[i];

			try {
				const auto W = wt.getProperties(p, T);
				(*rho_wF)[i] = W.rho;
				(*mu_wF)[i] = W.mu;
				(*cp_wF)[i] = W.cp;
				(*k_wF)[i] = W.k;
			}
			catch (...) {
				++oorCount;
				(*rho_wF)[i] = 1000;
				(*mu_wF)[i] = 1e-3;
				(*cp_wF)[i] = 4200;
				(*k_wF)[i] = 0.6;
			}

			(*lambda_wF)[i] = (*k_rwF)[i] / std::max((*mu_wF)[i], 1e-20);
		}
		if (oorCount) {
			std::cerr << "[WaterProps] OOR count = " << oorCount
				<< " (using fallback constants)\n";
		}
		return true;
	}
	/// 计算裂缝中的水的物性参数
	/////常数物性参数工况
	inline bool compute_Water_properties_const_inFracture(MeshManager& mgr, FieldRegistry& reg_fr, const std::string& pf_field, const std::string& Tf_field)
	{
		const FractureNetwork& frNet = mgr.fracture_network();
		//调用裂缝段索引
		const auto idx = buildFracElemIndex(frNet);
		const size_t ne = idx.total;
		if (!ne) {
			std::cout << "[PPM] No fracture elements. Skip InitializeFractureElementsProperties.\n";
			return false;
		}
		ensure_WaterProp_Fields(reg_fr, ne);
		// 取出计算裂缝物性参数需要的主变量场指针
		auto pfw = reg_fr.get<volScalarField>(pf_field);
		auto Tf = reg_fr.get<volScalarField>(Tf_field);
		// 取出裂缝物性参数场指针
		auto rho_wF = reg_fr.get<volScalarField>(PhysicalProperties_string::Water().rho_tag);				//rho，kg/m³
		auto cp_wF = reg_fr.get<volScalarField>(PhysicalProperties_string::Water().cp_tag);				//Cp，J/(kg·K)
		auto k_wF = reg_fr.get<volScalarField>(PhysicalProperties_string::Water().k_tag);					//k，W/(m·K)
		auto mu_wF = reg_fr.get<volScalarField>(PhysicalProperties_string::Water().mu_tag);				//mu，Pa·s
		auto Drho_Dp_wF = reg_fr.get<volScalarField>(PhysicalProperties_string::Water().drho_w_dp_tag); 	//drho_dp，kg/(m³·Pa)
		auto c_wF = reg_fr.get<volScalarField>(PhysicalProperties_string::Water().c_w_tag);				//二氧化碳的可压缩系数，1/Pa
		auto k_rwF = reg_fr.get<volScalarField>(PhysicalProperties_string::Water().k_rw_tag);				// Relative permeability of Water
		auto lambda_wF = reg_fr.get<volScalarField>(PhysicalProperties_string::Water().lambda_w_tag);		// Mobility of Water

		// 遍历所有裂缝段
		for (size_t f = 0; f < frNet.fractures.size(); ++f)
		{
			const auto& F = frNet.fractures[f];
			const size_t base = idx.offset[f];
			for (size_t e = 0; e < F.elements.size(); ++e)
			{
				const size_t g = base + e; //裂缝段全局索引

				double P = (*pfw)[g]; // 裂缝段水相压力，Pa
				double T = (*Tf)[g];  // 裂缝段温度，K
				(*rho_wF)[g] = WaterProperties().rho;
				(*cp_wF)[g] = WaterProperties().cp;
				(*k_wF)[g] = WaterProperties().k;
				(*mu_wF)[g] = WaterProperties().mu;
				(*Drho_Dp_wF)[g] = WaterProperties().dRho_dP;
				(*c_wF)[g] = WaterProperties().c;
				//流度计算，在计算之前需要计算相对渗透率，这里假设相对渗透率已经计算完成并存储在k_rgF中
				(*lambda_wF)[g] = (*k_rwF)[g] / std::max((*mu_wF)[g], 1e-20);
			}
		}
		return true;
	}
	/////基于LAPWS模型，利用插值表计算水的物性参数
	inline bool computer_water_properties_LAPWS_inFracture(MeshManager& mgr, FieldRegistry& reg_fr, const std::string& pf_field, const std::string& Tf_field)
	{
		const FractureNetwork& frNet = mgr.fracture_network();
		//调用裂缝段索引
		const auto idx = buildFracElemIndex(frNet);
		const size_t ne = idx.total;
		if (!ne) {
			std::cout << "[PPM] No fracture elements. Skip InitializeFractureElementsProperties.\n";
			return false;
		}
		ensure_WaterProp_Fields(reg_fr, ne);
		// 取出计算裂缝物性参数需要的主变量场指针
		auto pfw = reg_fr.get<volScalarField>(pf_field);
		auto Tf = reg_fr.get<volScalarField>(Tf_field);
		// 取出裂缝物性参数场指针
		auto rho_wF = reg_fr.get<volScalarField>(PhysicalProperties_string::Water().rho_tag);				//rho，kg/m³
		auto cp_wF = reg_fr.get<volScalarField>(PhysicalProperties_string::Water().cp_tag);				//Cp，J/(kg·K)
		auto k_wF = reg_fr.get<volScalarField>(PhysicalProperties_string::Water().k_tag);					//k，W/(m·K)
		auto mu_wF = reg_fr.get<volScalarField>(PhysicalProperties_string::Water().mu_tag);				//mu，Pa·s
		auto k_rwF = reg_fr.get<volScalarField>(PhysicalProperties_string::Water().k_rw_tag);				// Relative permeability of Water
		auto lambda_wF = reg_fr.get<volScalarField>(PhysicalProperties_string::Water().lambda_w_tag);		// Mobility of Water

		auto& gt = WaterPropertyTable::instance();
		// 遍历所有裂缝段
		size_t oorCount = 0;      // out-of-range or query failure
		for (size_t f = 0; f < frNet.fractures.size(); ++f)
		{
			const auto& F = frNet.fractures[f];
			const size_t base = idx.offset[f];
			for (size_t e = 0; e < F.elements.size(); ++e)
			{
				const size_t g = base + e;
				double P = (*pfw)[g]; // 裂缝段水相压力，Pa
				double T = (*Tf)[g];  // 裂缝段温度，K
				double rho = BASE_rho_w;
				double mu = BASE_mu_w;
				double cp = BASE_cp_w;
				double k = BASE_k_w;
				try {
					const auto G = gt.getProperties(P, T);
					rho = G.rho;
					mu = G.mu;
					cp = G.cp;
					k = G.k;
				}
				catch (...) {
					++oorCount; // keep fallback
				}
				(*rho_wF)[g] = rho;
				(*mu_wF)[g] = mu;
				(*cp_wF)[g] = cp;
				(*k_wF)[g] = k;
				// mobility/lambda with protection
				const double mu_eff = std::max(mu, 1e-20);
				(*lambda_wF)[g] = (*k_rwF)[g] / mu_eff;
			}
			if (oorCount > 0) {
				std::cerr << "[PPM][Fluid][water] property query failed/OOR in "
					<< oorCount << " elements in Fracture ID" << F.id << "(using fallback constants)\n";
			}
			
		}
		return true;
	}

	// 基岩内在 (p_eval, T_eval) 处评估 ρ 与 ∂ρ/∂p
	inline bool compute_water_DrhoDp
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
		auto dF = reg.getOrCreate<volScalarField>(PhysicalProperties_string::Water().drho_w_dp_tag, cells.size(), 0.0);
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
	//裂缝内在 (p_eval, T_eval) 处评估 ρ 与 ∂ρ/∂p
	inline bool compute_water_DrhoDp_inFracture
	(
		MeshManager& mgr, FieldRegistry& reg_fr,
		const std::string& p_eval_name,
		const std::string& T_eval_name,
		double dp_rel = 1e-4, double dp_abs_min = 1.0
	)
	{
		const FractureNetwork& frNet = mgr.fracture_network();
		//调用裂缝段索引
		const auto idx = buildFracElemIndex(frNet);
		const size_t ne = idx.total;
		if (!ne) {
			std::cout << "[PPM] No fracture elements. Skip compute_water_DrhoDp_inFracture.\n";
			return true;
		}
		auto pE = reg_fr.get<volScalarField>(p_eval_name);
		auto TE = reg_fr.get<volScalarField>(T_eval_name);
		if (!pE || !TE)
		{
			std::cerr << "[computeRhoAndDrhoDp_inFracture] missing eval fields '"
				<< p_eval_name << "' or '" << T_eval_name << "'\n";
			return false;
		}
		auto dF = reg_fr.getOrCreate<volScalarField>(PhysicalProperties_string::Water().drho_w_dp_tag, ne, 0.0);
		auto& wt = WaterPropertyTable::instance();
		// 遍历所有裂缝段
		for (size_t f = 0; f < frNet.fractures.size(); ++f)
		{
			const auto& F = frNet.fractures[f];
			const size_t base = idx.offset[f];
			for (size_t e = 0; e < F.elements.size(); ++e)
			{
				const size_t g = base + e; //裂缝段全局索引
				double p = (*pE)[g];
				double T = (*TE)[g];
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
				(*dF)[g] = (rp - rm) / std::max(2.0 * dpa, 1e-12);
			}
		}
		return true;
	}

	// 计算当水为基岩中唯一流体时，的有效热物性参数 是否为常数取决于compute_water_properties_const还是computer_water_properties_LAPWS
	inline bool computer_effectiveThermalProperties(MeshManager& mgr, FieldRegistry& reg)
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();

		ensure_WaterProp_Fields(reg, n);

		// 岩石参数（必须存在）
		auto phiF = reg.get<volScalarField>(PhysicalProperties_string::Rock().phi_tag);
		auto rrF = reg.get<volScalarField>(PhysicalProperties_string::Rock().rho_tag);
		auto cprF = reg.get<volScalarField>(PhysicalProperties_string::Rock().cp_tag);
		auto lamrF = reg.get<volScalarField>(PhysicalProperties_string::Rock().lambda_tag);


		//流体物性参数
		auto rho_gF = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_tag);
		auto cp_gF = reg.get<volScalarField>(PhysicalProperties_string::Water().cp_tag);
		auto k_gF = reg.get<volScalarField>(PhysicalProperties_string::Water().k_tag);

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

} //namespace Water_Prop


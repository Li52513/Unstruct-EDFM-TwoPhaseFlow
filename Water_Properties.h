#pragma once
#include "PropertiesSummary.h" // 里面定义 水所需要的物性参数
#include "MeshManager.h"
#include "FieldRegistry.h"

namespace Water_in_rock
{
	static constexpr double BASE_rho_w = 980;
	static constexpr double BASE_cp_w = 4300;
	static constexpr double BASE_k_w = 0.05;
	static constexpr double BASE_mu_w = 3e-54;
	static constexpr double BASE_Drho_Dp_w = 0;

	inline bool computeWATERinROCKProperties(MeshManager& mgr, FieldRegistry& reg, const std::string& p_w_field, const std::string& T_field)  //这里传入基岩的压力和温度 但是对于IMPES来说，不会更新这里的参数
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		auto TF = reg.get<volScalarField>(T_field);
		auto p_wF = reg.get<volScalarField>(p_w_field);

		auto rho_wF = reg.get<volScalarField>("rho_w");
		auto cp_wF = reg.get<volScalarField>("cp_w");
		auto k_wF = reg.get<volScalarField>("k_w");
		auto mu_wF = reg.get<volScalarField>("mu_w");
		auto Drho_Dp_wF = reg.get<volScalarField>("Drho_Dp_w");

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
		}
		return true;
	}
}

namespace Water_in_fracture
{
	static constexpr double BASE_rho_w = 980;
	static constexpr double BASE_cp_w = 4300;
	static constexpr double BASE_k_w = 0.05;
	static constexpr double BASE_mu_w = 3e-4;
	static constexpr double BASE_Drho_Dp_w = 0;

	inline bool computerWaterinFractureProperties(MeshManager& mgr, FieldRegistry& reg, FieldRegistry& reg_fr, const std::string& p_field_fr, const std::string& T_field_fr)
	{
		auto& mesh = mgr.mesh();
		// 统计裂缝段总数
		size_t Nseg = 0;
		for (auto& F : mgr.fracture_network().fractures) Nseg += F.elements.size();

		auto pF = reg_fr.get<volScalarField>(p_field_fr);
		auto TF = reg_fr.get<volScalarField>(T_field_fr);
		if (!pF || !TF) { std::cerr << "[PPM][FrFluid] missing fields.\n"; }

		auto fr_rho_w = reg_fr.get<volScalarField>("fr_rho_w");
		auto fr_mu_w = reg_fr.get<volScalarField>("fr_mu_w");
		auto fr_cp_w = reg_fr.get<volScalarField>("fr_cp_w");
		auto fr_k_w = reg_fr.get<volScalarField>("fr_k_w");
		auto fr_Drho_Dp_w = reg_fr.get<volScalarField>("fr_Drho_Dp_w");

		size_t gid = 0;
		for (auto& F : mgr.fracture_network().fractures)
		{
			for (auto& E : F.elements)
			{
				(*fr_rho_w)[gid] = BASE_rho_w;
				(*fr_mu_w)[gid] = BASE_mu_w;
				(*fr_cp_w)[gid] = BASE_cp_w;
				(*fr_k_w)[gid] = BASE_k_w;
				(*fr_Drho_Dp_w)[gid] = BASE_Drho_Dp_w;
				++gid;
			}
		}

		return true;

	}
}


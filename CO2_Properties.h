#pragma once
#include "PropertiesSummary.h" // 里面定义 CO2所需要的物性参数
#include "MeshManager.h"
#include "FieldRegistry.h"

namespace CO2_in_rock
{
	static constexpr double BASE_rho_g = 800;
	static constexpr double BASE_cp_g = 1200;
	static constexpr double BASE_k_g = 0.05;
	static constexpr double BASE_mu_g = 1.48e-5;
	static constexpr double BASE_Drho_Dp_g = 0;
	static constexpr double BASE_c_g = 0;

	inline bool computeCO2inROCKProperties(MeshManager& mgr, FieldRegistry& reg, const std::string& p_g_field,const std::string& T_field)  //这里传入基岩的压力和温度 但是对于IMPES来说，不会更新这里的参数
	{
		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		auto TF = reg.get<volScalarField>(T_field);
		auto p_gF= reg.get<volScalarField>(p_g_field);

		auto rho_gF = reg.get<volScalarField>("rho_g");
		auto cp_gF = reg.get<volScalarField>("cp_g");
		auto k_gF = reg.get<volScalarField>("k_g");
		auto mu_gF = reg.get<volScalarField>("mu_g");
		auto Drho_Dp_gF = reg.get<volScalarField>("Drho_Dp_g");
		auto c_gF = reg.get<volScalarField>("c_g");


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
}

namespace CO2_in_fracture
{
	static constexpr double BASE_rho_g = 800;
	static constexpr double BASE_cp_g = 1200;
	static constexpr double BASE_k_g = 0.05;
	static constexpr double BASE_mu_g = 1.48e-5;
	static constexpr double BASE_Drho_Dp_g = 0;

	inline bool computerCO2inFractureProperties(MeshManager& mgr, FieldRegistry& reg, FieldRegistry& reg_fr, const std::string& p_field_fr, const std::string& T_field_fr)
	{
		auto& mesh = mgr.mesh();
		// 统计裂缝段总数
		size_t Nseg = 0;
		for (auto& F : mgr.fracture_network().fractures) Nseg += F.elements.size();

		auto pF = reg_fr.get<volScalarField>(p_field_fr);
		auto TF = reg_fr.get<volScalarField>(T_field_fr);
		if (!pF || !TF) { std::cerr << "[PPM][FrFluid] missing fields.\n"; }

		auto fr_rho_g = reg_fr.get<volScalarField>("fr_rho_g");
		auto fr_mu_g = reg_fr.get<volScalarField>("fr_mu_g");
		auto fr_cp_g = reg_fr.get<volScalarField>("fr_cp_g");
		auto fr_k_g = reg_fr.get<volScalarField>("fr_k_g");
		auto fr_Drho_Dp_g = reg_fr.get<volScalarField>("fr_Drho_Dp_g");
		auto fr_c_w_g = reg_fr.get<volScalarField>("fr_c_g"); // CO2 可压缩系数

		size_t gid = 0;
		for (auto& F : mgr.fracture_network().fractures) 
		{
			for (auto& E : F.elements) 
			{
				(*fr_rho_g)[gid] = BASE_rho_g;
				(*fr_mu_g)[gid] = BASE_mu_g;
				(*fr_cp_g)[gid] = BASE_cp_g;
				(*fr_k_g)[gid] = BASE_k_g;
				(*fr_Drho_Dp_g)[gid] = BASE_Drho_Dp_g;
				double p_g = (*pF)[gid];
				(*fr_c_w_g)[gid] = 1 / p_g; //假设可压缩系数为1/p
				++gid;
			}
		}
		return true;

	}
}

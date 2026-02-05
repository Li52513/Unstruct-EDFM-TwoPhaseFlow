#pragma once
#include <iostream>  //
#include <limits>  //
#include "FieldRegistry.h"  //
#include "InitConfig.h" //
#include "CapRelPerm.h" //
#include "Mesh.h" //
#include "WaterPropertyTable.h" //
#include "CO2PropertyTable.h" //
#include "FracIndex.h"   // 

/********注释********/
//创建并初始化主变量场
//计算 p_c, p_g, k_rw, k_rg
//即时查询表合成 C_eff、λ_eff（不保留 ρ/μ/cp/λ 场）


//初始化诊断
struct InitDiagnostics
{
	double pmin = +1e99;
	double pmax = -1e99;
	double Tmin = +1e99;
	double Tmax = -1e99;
	double Swmin = +1e99;
	double Swmax = -1e99;
	size_t satClamped = 0;
	size_t propOOR = 0;
};


/********注释********/
//可实现功能： 
//1. 创建基岩内部代求变量场（p_w, S_w, T, p_c, p_g, kr_w, kr_g） - createPrimaryFields
//2. 根据 InitFields 线性分布初始化 p_w, S_w, T  - fillBaseDistributions
//3. 对 S_w 进行限幅 - enforceSaturationBounds
//4.调用VG模型和相对渗透率模型计算闭合关系 - computeClosure
//5.计算有效热物性参数 - computerEffectiveThermals
//6.初始化裂缝内部代求主变量场（pf_w, Sf_w, Tf） - initFracturePrimaries
//7.printDiag

struct Initializer
{

	static void createPrimaryFields(Mesh& mesh, FieldRegistry& reg, const std::string& PhysicalFieldName)
	{
		const size_t n = mesh.getCells().size();
		if (!reg.has(PhysicalFieldName)) reg.create<volScalarField>(PhysicalFieldName.c_str(), n, 0);
	}

template<class InitFields>
	static void fillBaseDistributions1(const Mesh& mesh, FieldRegistry& reg, const InitFields& init, const std::string& PhysicalFieldName)
	{
		auto PhysicalField = reg.get<volScalarField>(PhysicalFieldName.c_str());
		for (const auto& c : mesh.getCells())
		{
			const size_t i = mesh.getCellId2Index().at(c.id);
			(*PhysicalField)[i] = init.x0 + init.x_dx * c.center.m_x + init.x_dy * c.center.m_y + init.x_dz * c.center.m_z;
		}
	}

//裂缝段主变量初始化，以所在基岩单元的主变量为初始值
template<class FractureNetworkT>
    static void initFracturePrimaries(const Mesh& mesh, FractureNetworkT& frNet, FieldRegistry& reg_r, FieldRegistry& reg_fr)
    {
		//取出基岩主变量场指针
		auto p  = reg_r.get<volScalarField>("p_w");
        auto Sw = reg_r.get<volScalarField>("S_w");
        auto T  = reg_r.get<volScalarField>("T");
		const auto& id2idx = mesh.getCellId2Index();

		//调用裂缝段索引
		const auto idx = buildFracElemIndex(frNet);
		if (!idx.total) 
		{
			std::cout << "[Initializer] No fracture elements. Skip initFracturePrimaries.\n";
			return;
		}

		//确保裂缝主变量场存在

		auto pfw = reg_fr.getOrCreate<volScalarField>("pf_w", idx.total, 1e6); //裂缝段水相压力，Pa
		auto Sfw = reg_fr.getOrCreate<volScalarField>("Sf_w", idx.total, 0.9); //裂缝段水相饱和度，1
		auto Tf = reg_fr.getOrCreate<volScalarField>("Tf", idx.total, 303.15); //裂缝段温度，K

		for (size_t f = 0; f < frNet.fractures.size(); ++f)
		{
			const auto& F = frNet.fractures[f];
			const size_t base = idx.offset[f]; //本条裂缝的全局起点
			for (size_t e = 0; e < F.elements.size(); ++e)
			{
				const auto& elem = F.elements[e];
				const size_t g = base + e; //裂缝段全局索引

				if (elem.cellID >= 0)
				{
					auto it = id2idx.find(elem.cellID);
					if (it != id2idx.end())
					{
						const size_t i = it->second; //基岩单元内部索引 
						(*pfw)[g] = (*p)[i];
						(*Sfw)[g] = (*Sw)[i];
						(*Tf)[g] = (*T)[i];
						continue;
					}
				}
			}
		}
    }

	static void printDiag(const InitDiagnostics& d)
	{
		std::cout << "[Init] p_w in [" << d.pmin << "," << d.pmax << "], T in [" << d.Tmin << "," << d.Tmax
			<< "], S_w in [" << d.Swmin << "," << d.Swmax << "]\n";
		std::cout << "[Init] S_w clamped = " << d.satClamped << ", props OOR = " << d.propOOR << "\n";
	}

};


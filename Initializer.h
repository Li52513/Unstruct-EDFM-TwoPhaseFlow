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
	static void createPrimaryFields(Mesh& mesh, FieldRegistry& reg)
	{
		const size_t n = mesh.getCells().size();
		if (!reg.has("p_w")) reg.create <volScalarField>("p_w", n, 1e6);  //水相压力，Pa
		if (!reg.has("S_w")) reg.create <volScalarField>("S_w", n, 0.9);   //水相饱和度，1
		if (!reg.has("T")) reg.create <volScalarField>("T", n, 303.15);		//温度，K
		if (!reg.has("p_c")) reg.create<volScalarField>("p_c", n, 0.0);   //毛细压力，Pa
		if (!reg.has("p_g")) reg.create<volScalarField>("p_g", n, 1e5);   //气相压力，Pa		
		if (!reg.has("kr_w"))reg.create<volScalarField>("kr_w", n, 0.0);	//水相相对渗透率，1
		if (!reg.has("kr_g"))reg.create<volScalarField>("kr_g", n, 0.0);	//气相相对渗透率，1

	}

	static void createPrimaryFields_test_singlePhase_CO2_T_diffusion(Mesh& mesh, FieldRegistry& reg)
	{
		const size_t n = mesh.getCells().size();
		if (!reg.has("T")) reg.create <volScalarField>("T", n, 303.15);		//温度，K
		if (!reg.has("p_g")) reg.create<volScalarField>("p_g", n, 1e5);   //气相压力，Pa		
	}

	static void createPrimaryFields_singlePhase_CO2_T(Mesh& mesh, FieldRegistry& reg)
	{
		const size_t n = mesh.getCells().size();
		if (!reg.has("T")) reg.create <volScalarField>("T", n, 303.15);		//温度，K
	}

	static void createPrimaryFields_singlePhase_CO2_P(Mesh& mesh, FieldRegistry& reg)
	{
		const size_t n = mesh.getCells().size();
		if (!reg.has("p_g")) reg.create<volScalarField>("p_g", n, 1e5);   //气相压力，Pa		
	}



	static void fillBaseDistributions_test_singlePhase_CO2_T_diffusion(Mesh& mesh, FieldRegistry& reg, const InitFields& init)
	{
		auto p_g = reg.get<volScalarField>("p_g");  //气相压力，Pa
		auto T = reg.get<volScalarField>("T");		//温度，K
		for (const auto& c : mesh.getCells())
		{
			const size_t i = mesh.getCellId2Index().at(c.id); //获取网格单元的内部下标
			(*p_g)[i] = init.p0 + init.dpdx * c.center.m_x + init.dpdy * c.center.m_y + init.dpdz * c.center.m_z;
			(*T)[i] = init.T0 + init.dTdx * c.center.m_x + init.dTdy * c.center.m_y + init.dTdz * c.center.m_z;
		}
	}

	static void fillBaseDistributions_singlePhase_CO2_P(Mesh& mesh, FieldRegistry& reg, const InitFields& init)
	{
		auto p_g = reg.get<volScalarField>("p_g");  //气相压力，Pa
		for (const auto& c : mesh.getCells())
		{
			const size_t i = mesh.getCellId2Index().at(c.id); //获取网格单元的内部下标
			(*p_g)[i] = init.p0 + init.dpdx * c.center.m_x + init.dpdy * c.center.m_y + init.dpdz * c.center.m_z;
		}
	}

	static void fillBaseDistributions_singlePhase_CO2_T(Mesh& mesh, FieldRegistry& reg, const InitFields& init)
	{
		auto T = reg.get<volScalarField>("T");		//温度，K
		for (const auto& c : mesh.getCells())
		{
			const size_t i = mesh.getCellId2Index().at(c.id); //获取网格单元的内部下标
			(*T)[i] = init.T0 + init.dTdx * c.center.m_x + init.dTdy * c.center.m_y + init.dTdz * c.center.m_z;
		}
	}



	static void fillBaseDistributions(Mesh& mesh, FieldRegistry& reg, const InitFields& init)
	{
		auto p_w = reg.get<volScalarField>("p_w"); //水相压力，Pa
		auto S_w = reg.get<volScalarField>("S_w");  //水相饱和度，1
		auto T = reg.get<volScalarField>("T");		//温度，K
		for (const auto& c : mesh.getCells())
		{
			const size_t i = mesh.getCellId2Index().at(c.id); //获取网格单元的内部下标
			(*p_w)[i] = init.p0 + init.dpdx * c.center.m_x+ init.dpdy * c.center.m_y + init.dpdz * c.center.m_z;
			(*T)[i] = init.T0 + init.dTdx * c.center.m_x + init.dTdy * c.center.m_y + init.dTdz * c.center.m_z;
			(*S_w)[i] = init.sw0;

		}
	}

	static void enforceSaturationBounds(FieldRegistry& reg, const VGParams& vg, InitDiagnostics& d)  //对水相饱和度进行限幅
	{
		auto Sw = reg.get<volScalarField>("S_w");
		d.Swmin = +1e99; d.Swmax = -1e99;
		for (double& s : Sw->data) 
		{
			const double s0 = s;
			s = clamp(s, vg.Swr + 1e-6, 1.0 - vg.Sgr - 1e-6);
			if (s != s0) ++d.satClamped;
			d.Swmin = std::min(d.Swmin, s);
			d.Swmax = std::max(d.Swmax, s);
		}
	}

	static void computeClosure(Mesh& mesh, FieldRegistry& reg, const VGParams& vg, const RelPermParams& rp, InitDiagnostics& d) //计算
	{	

		///获取场变量指针，准备进行读写操作
		auto p_w = reg.get<volScalarField>("p_w"); //水相压力，Pa
		auto T = reg.get<volScalarField>("T");	//温度，K
		auto Sw = reg.get<volScalarField>("S_w"); //水相饱和度，1
		auto pc = reg.get<volScalarField>("p_c"); //毛细压力，Pa
		auto pg = reg.get<volScalarField>("p_g"); //气相压力，Pa
		auto krw = reg.get<volScalarField>("kr_w"); //水相相对渗透率，1
		auto krg = reg.get<volScalarField>("kr_g"); //气相相对渗透率，1

		///初始化诊断变量
		d.pmin = +1e99; d.pmax = -1e99; d.Tmin = +1e99; d.Tmax = -1e99;
		d.Swmin = +1e99; d.Swmax = -1e99;

		for (const auto& c : mesh.getCells())
		{
			//记录主变量的极值
			const size_t i = mesh.getCellId2Index().at(c.id); //获取网格单元的内部下标
			d.pmin = min(d.pmin, (*p_w)[i]);
			d.pmax = max(d.pmax, (*p_w)[i]);
			d.Tmin = min(d.Tmin, (*T)[i]);
			d.Tmax = max(d.Tmax, (*T)[i]);
			d.Swmin = std::min(d.Swmin, (*Sw)[i]);
			d.Swmax = std::max(d.Swmax, (*Sw)[i]);

			//计算毛细压力
			const double pc_i = pc_vG((*Sw)[i], vg);
			if (pc_i < 0.0) { ++d.propOOR; }
			(*pc)[i] = pc_i;
			(*pg)[i] = (*p_w)[i] + pc_i;

			//计算相对渗透率
			double kw, kg;
			kr_Mualem_vG((*Sw)[i], vg, rp, kw, kg);
			(*krw)[i] = kw; (*krg)[i] = kg;
		}
	}

	static inline void clampPT(double& p, double& T)
	{
		const double pMin = 5e6; // Pa 插值表压力下限.5MPa
		const double pMax = 7e7; // Pa 插值表压力上限,70MPa
		const double TMin = 293.15; // K 插值表温度下限,20℃
		const double TMax = 693.15; // K 插值表温度上限,200℃

		if (p < pMin)
		{
			p = pMin;
			cout << "[Warning] p clamped to " << pMin << " Pa\n";
		}
		else if (p > pMax)
		{
			p = pMax;
			cout << "[Warning] p clamped to " << pMax << " Pa\n";
		}
		if (T < TMin)
		{
			T = TMin;
			cout << "[Warning] T clamped to " << TMin << " K\n";
		}
		else if (T > TMax)
		{
			T = TMax;
			cout << "[Warning] T clamped to " << TMax << " K\n";
		}
	}


	static void computerEffectiveThermals(Mesh& mesh, FieldRegistry& reg, const RockDefaults& rock, InitDiagnostics& d)
	{
		//获取主变量场指针
		auto Sw = reg.get<volScalarField>("S_w");
		auto p_w = reg.get<volScalarField>("p_w");
		auto p_g = reg.get<volScalarField>("p_g");
		auto T = reg.get<volScalarField>("T");
		//获取岩石物性参数场指针
		auto rho_rF = reg.get<volScalarField>("rho_r");
		auto cp_rF = reg.get<volScalarField>("cp_r");
		auto lambda_rF = reg.get<volScalarField>("lambda_r");
		auto phiF = reg.get<volScalarField>("phi");
		//获取或创建有效热物性参数场指针
		auto C_eff = reg.getOrCreate<volScalarField>("C_eff", mesh.getCells().size(), 0.0);
		auto lambda_eff = reg.getOrCreate<volScalarField>("lambda_eff", mesh.getCells().size(), 0.0);

		//创建并初始化流体物性参数表
		auto& waterTable = WaterPropertyTable::instance();
		auto& co2Table = CO2PropertyTable::instance();
		
		//遍历网格单元，计算有效热物性参数
		for (const auto& c : mesh.getCells())
		{

			const size_t i = mesh.getCellId2Index().at(c.id); //获取网格单元的内部下标
			const double phi = phiF ? (*phiF)[i] : rock.phi;		// 孔隙度
			const double rho_r = rho_rF ? (*rho_rF)[i] : rock.rho_r;   // 基岩密度，kg/m³
			const double cp_r = cp_rF ? (*cp_rF)[i] : rock.cp_r;     // 基岩比热容，J/(kg·K)
			const double lambda_r = lambda_rF ? (*lambda_rF)[i] : rock.lambda_r; // 基岩导热系数，W/(m·K)
			double pw = (*p_w)[i], Tg = (*T)[i];
			double pg = (*p_g)[i];
			clampPT(pw, Tg); clampPT(pg, Tg);
			const double Swc = (*Sw)[i];
			const double Sgc = 1.0 - Swc;
	
			double rho_w = 1000, mu_w = 1e-3, cp_w = 4200, lam_w = 0.6;
			double rho_g = 600, mu_g = 1e-4, cp_g = 850, lam_g = 0.08;

			try {
				const auto W = waterTable.getProperties(pw, Tg);
				rho_w = W.rho; mu_w = W.mu; cp_w = W.cp; lam_w = W.k;
			}
			catch (...) { ++d.propOOR; /* 兜底值已就绪 */ }

			try {
				const auto G = co2Table.getProperties(pg, Tg);
				rho_g = G.rho; mu_g = G.mu; cp_g = G.cp; lam_g = G.k;
			}
			catch (...) { ++d.propOOR; }

			(void)mu_w; (void)mu_g; // 流动方程装配时会用
            (*C_eff)[i] = (1.0 - phi) * rho_r * cp_r + phi * (Swc * rho_w * cp_w + Sgc * rho_g * cp_g);
		    (*lambda_eff)[i] = (1.0 - phi) * lambda_r + phi * (Swc * lam_w + Sgc * lam_g);

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

		size_t ghostHit = 0; // 记录落在 ghost cell 的裂缝段数

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
				//能走到这里，说明裂缝段所在单元是 ghost cell
				++ghostHit;
			}
		}
		if (ghostHit > 0) 
		{
			std::cout << "[Initializer] initFracturePrimaries: ghost/invalid host hits = "
				<< ghostHit << "\n";
		}
    }

	static void printDiag(const InitDiagnostics& d)
	{
		std::cout << "[Init] p_w in [" << d.pmin << "," << d.pmax << "], T in [" << d.Tmin << "," << d.Tmax
			<< "], S_w in [" << d.Swmin << "," << d.Swmax << "]\n";
		std::cout << "[Init] S_w clamped = " << d.satClamped << ", props OOR = " << d.propOOR << "\n";
	}

};


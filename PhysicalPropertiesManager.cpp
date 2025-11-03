#include "PhysicalPropertiesManager.h"
#include "FractureNetwork.h"
#include "FractureTypes.h"
#include "MeshManager.h"
#include "WaterPropertyTable.h"
#include "CO2PropertyTable.h"
#include "RockSolidProperties.h"
#include "FractureSolidProperties.h"
#include "PhysicalProperties_CO2.h"
#include <cassert> // for assert
#include <algorithm>  // for swap
#include "Initializer.h"
#include <iomanip>


//***************************基岩与裂缝区域分类*************************************//

//用来判断点是否在凸多边形内的函数
static bool pointInConvexPolygon(const Vector& p, const vector<Vector>& verts)
{
	assert(verts.size() >= 3); // 至少需要三个顶点 NOTE: 这里假设 verts 是凸多边形的顶点；其次assert只能在debug模式下生效
	auto cross2d = [](const Vector& u, const Vector& v) 
		{
			return u.m_x * v.m_y - u.m_y * v.m_x;  // 计算二维向量的叉积
		};

	// 取第一个点作基准
	const Vector& a0 = verts[0];
	// 对每一条边 (a_i -> a_{i+1}) 检查同侧
	double sign0 = 0;
	for (size_t i = 0, n = verts.size(); i < n; ++i) {
		const Vector& a = verts[i];
		const Vector& b = verts[(i + 1) % n];
		Vector va = a - p;
		Vector vb = b - p;
		double c = cross2d(va, vb);
		if (i == 0) sign0 = c;
		else {
			// 如果出现异号，点就在多边形外
			if (c * sign0 < 0) return false;
		}
	}
	return true;
}

//***************************小工具：确保物性参数场存在*************************************//
static inline void ensureMatrixFluidFields(FieldRegistry& reg, std::size_t n)
{
	reg.getOrCreate<volScalarField>("rho_w", n, 1000.0);	//水的密度，kg/m³
	reg.getOrCreate<volScalarField>("mu_w", n, 1e-3);		//水的粘度，Pa·s	
	reg.getOrCreate<volScalarField>("cp_w", n, 4182.0);		//水的比热容，J/(kg·K)
	reg.getOrCreate<volScalarField>("k_w", n, 0.6);			//水的导热系数，W/(m·K)
	reg.getOrCreate<volScalarField>("rho_g", n, 1.98);		//二氧化碳的密度，kg/m³
	reg.getOrCreate<volScalarField>("mu_g", n, 1.48e-5);	//二氧化碳的粘度，Pa·s
	reg.getOrCreate<volScalarField>("cp_g", n, 846.0);		//二氧化碳的比热容，J/(kg·K)
	reg.getOrCreate<volScalarField>("k_g", n, 0.0146);		//二氧化碳的导热系数，W/(m·K)
	reg.getOrCreate<volScalarField>("Drho_Dp_w", n, 0.0);	//水的密度对压力的导数，kg/(m³·Pa)
	reg.getOrCreate<volScalarField>("Drho_Dp_g", n, 0.0);	//二氧化碳的密度对压力的导数，kg/(m³·Pa)
}

static inline void ensureFractureFluidFields(FieldRegistry& reg_fr, std::size_t ne)
{
	reg_fr.getOrCreate<volScalarField>("fr_rho_w", ne, 1000.0);
	reg_fr.getOrCreate<volScalarField>("fr_mu_w", ne, 1e-3);
	reg_fr.getOrCreate<volScalarField>("fr_cp_w", ne, 4182.0);
	reg_fr.getOrCreate<volScalarField>("fr_k_w", ne, 0.6);
	reg_fr.getOrCreate<volScalarField>("fr_rho_g", ne, 1.98);
	reg_fr.getOrCreate<volScalarField>("fr_mu_g", ne, 1.48e-5);
	reg_fr.getOrCreate<volScalarField>("fr_cp_g", ne, 846.0);
	reg_fr.getOrCreate<volScalarField>("fr_k_g", ne, 0.0146);
}

static inline void ensureRockFields(FieldRegistry& reg, size_t n)
{
	reg.getOrCreate<volScalarField>("phi", n, 0.15);
	reg.getOrCreate<volScalarField>("kxx", n, 1e-14);
	reg.getOrCreate<volScalarField>("kyy", n, 1e-14);
	reg.getOrCreate<volScalarField>("kzz", n, 1e-14);
	reg.getOrCreate<volScalarField>("rho_r", n, 2650.0);
	reg.getOrCreate<volScalarField>("cp_r", n, 1000.0);
	reg.getOrCreate<volScalarField>("lambda_r", n, 2.5); 
	reg.getOrCreate<volScalarField>("c_phi", n, 1e-12); //孔隙度可压缩性
}

static inline void ensureFracRockFields(FieldRegistry& reg_fr, size_t ne) //确保裂缝物性参数场存在，若不存在则创建并赋默认值
{
	reg_fr.getOrCreate<volScalarField>("fr_phi", ne, 1); //裂隙孔隙度
	reg_fr.getOrCreate<volScalarField>("fr_k_t", ne, 1e-12);    // 切向等效渗透率
	reg_fr.getOrCreate<volScalarField>("fr_k_n", ne, 1e-16);    // 法向等效渗透率
	reg_fr.getOrCreate<volScalarField>("fr_rho_r", ne, 2650.0); //裂缝密度，kg/m³
	reg_fr.getOrCreate<volScalarField>("fr_cp_r", ne, 1000.0); //裂缝比热容，J/(kg·K)
	reg_fr.getOrCreate<volScalarField>("fr_lambda_r", ne, 2.5); //裂缝导热系数，W/(m·K)
	reg_fr.getOrCreate<volScalarField>("fr_aperture", ne, 1e-3); //裂缝开度
}

static inline void ensureRockPrimaryFields(FieldRegistry& reg, size_t n, const InitFields& init)
{
	reg.getOrCreate<volScalarField>("p_w", n, init.p0); //基岩初始水相压力，Pa
	reg.getOrCreate<volScalarField>("T", n, init.T0);   //基岩初始温度，K
	reg.getOrCreate<volScalarField>("S_w", n, init.sw0); //基岩初始水相饱和度，1
}

inline void ensureFracPrimaryFields(FieldRegistry& freg, size_t ne)
{
	freg.getOrCreate<volScalarField>("pf_w", ne, 1.0e6);
	freg.getOrCreate<volScalarField>("Sf_w", ne, 0.90);
	freg.getOrCreate<volScalarField>("Tf", ne, 303.15);
}

//常物性测试 test_constProperties_singlePhase_CO2_T_diffusion
void PhysicalPropertiesManager::RockProperties_test_constProperties_singlePhase_CO2_T_diffusion(MeshManager& mgr, FieldRegistry& reg)
{
	auto& mesh = mgr.mesh();
	const auto& cells = mesh.getCells();
	const size_t n = cells.size();
	ensureRockFields(reg, n);

	auto phi_r = reg.get<volScalarField>("phi"); //孔隙度
	auto rho_r = reg.get<volScalarField>("rho_r"); //基岩密度，kg/m³
	auto cp_r = reg.get<volScalarField>("cp_r"); //基岩比热容，J/(kg·K)
	auto lam_r = reg.get<volScalarField>("lambda_r"); //基岩导热系数，W/(m·K)

	for (size_t ic = 0; ic < cells.size(); ++ic)
	{
		const auto& cell = cells[ic];
		if (cell.id < 0) continue;
		const size_t i = mesh.getCellId2Index().at(cell.id);
		(*phi_r)[i] = 0.15;
		(*rho_r)[i] = 2650.0;
		(*cp_r)[i] = 1000.0;
		(*lam_r)[i] = 3;
	}
}



void PhysicalPropertiesManager::CO2Properties_test_constProperties_singlePhase_CO2_T_diffusion(MeshManager& mgr, FieldRegistry& reg)
{
	auto& mesh = mgr.mesh();
	const auto& cells = mesh.getCells();
	const size_t n = cells.size();
	ensureMatrixFluidFields(reg, n);

	auto rho_gF = reg.get<volScalarField>("rho_g");
	auto cp_gF = reg.get<volScalarField>("cp_g");
	auto k_gF = reg.get<volScalarField>("k_g");
	for (size_t ic = 0; ic < cells.size(); ++ic)
	{
		const auto& c = cells[ic];
		if (c.id < 0) continue;
		const size_t i = mesh.getCellId2Index().at(c.id);
		(*rho_gF)[i] = 800.0;
		(*cp_gF)[i] = 1200.0;
		(*k_gF)[i] = 0.05;
	}
}



void PhysicalPropertiesManager::ComputeEffectiveThermalProperties_test_constProperties_singlePhase_CO2_T_diffusion(MeshManager& mgr, FieldRegistry& reg)
{
	auto& mesh = mgr.mesh();
	const auto& cells = mesh.getCells();
	const size_t n = cells.size();
	//ensureMatrixFluidFields(reg, n);
	//ensureRockFields(reg, n);

	// 岩石参数（必须存在）
	auto phiF = reg.get<volScalarField>("phi");
	auto rrF = reg.get<volScalarField>("rho_r");
	auto cprF = reg.get<volScalarField>("cp_r");
	auto lamrF = reg.get<volScalarField>("lambda_r");

	//流体物性参数
	auto rho_gF = reg.get<volScalarField>("rho_g");
	auto cp_gF = reg.get<volScalarField>("cp_g");
	auto k_gF = reg.get<volScalarField>("k_g");
	// 有效热参数场
	auto Ceff = reg.getOrCreate<volScalarField>("C_eff", n, 0.0);
	auto lame = reg.getOrCreate<volScalarField>("lambda_eff", n, 0.0);

	for (size_t ic = 0; ic < cells.size(); ++ic)
	{
		const auto& c = cells[ic];
		if (c.id < 0) continue;
		const size_t i = mesh.getCellId2Index().at(c.id);

		(*Ceff)[i] = (1.0 - (*phiF)[i]) * (*rrF)[i] * (*cprF)[i] + (*phiF)[i] * (*rho_gF)[i] * (*cp_gF)[i];
		(*lame)[i] = (1.0 - (*phiF)[i]) * (*lamrF)[i] + (*phiF)[i] * (*k_gF)[i];
	}
}

void PhysicalPropertiesManager::ComputeEffectiveThermalProperties_constProperties_singlePhase_CO2_T_H(MeshManager& mgr, FieldRegistry& reg)
{
	auto& mesh = mgr.mesh();
	const auto& cells = mesh.getCells();
	const size_t n = cells.size();
	//ensureMatrixFluidFields(reg, n);
	//ensureRockFields(reg, n);

	// 岩石参数（必须存在）
	auto phiF = reg.get<volScalarField>("phi");
	auto rrF = reg.get<volScalarField>("rho_r");
	auto cprF = reg.get<volScalarField>("cp_r");
	auto lamrF = reg.get<volScalarField>("lambda_r");

	//流体物性参数
	auto rho_gF = reg.get<volScalarField>("rho_g");
	auto cp_gF = reg.get<volScalarField>("cp_g");
	auto k_gF = reg.get<volScalarField>("k_g");
	// 有效热参数场
	auto Ceff = reg.getOrCreate<volScalarField>("C_eff", n, 0.0);
	auto lame = reg.getOrCreate<volScalarField>("lambda_eff", n, 0.0);

	for (size_t ic = 0; ic < cells.size(); ++ic)
	{
		const auto& c = cells[ic];
		if (c.id < 0) continue;
		const size_t i = mesh.getCellId2Index().at(c.id);

		(*Ceff)[i] = (1.0 - (*phiF)[i]) * (*rrF)[i] * (*cprF)[i] + (*phiF)[i] * (*rho_gF)[i] * (*cp_gF)[i];
		(*lame)[i] = (1.0 - (*phiF)[i]) * (*lamrF)[i] + (*phiF)[i] * (*k_gF)[i];
	}
}


//固体常物性，流体变物性测试案例，取COMSOL计算模型
//流体相参数
void PhysicalPropertiesManager::CO2Properties_test_varProperties_singlePhase_CO2_T_diffusion(MeshManager& mgr, FieldRegistry& reg,  const std::string& T_field)
{
	auto& mesh = mgr.mesh();
	const auto& cells = mesh.getCells();
	const size_t n = cells.size();
	ensureMatrixFluidFields(reg, n);
	auto TF = reg.get<volScalarField>(T_field);

	auto rho_gF = reg.get<volScalarField>("rho_g");
	auto cp_gF = reg.get<volScalarField>("cp_g");
	auto k_gF = reg.get<volScalarField>("k_g");


	//取消并行计算
	for (size_t ic = 0; ic < cells.size(); ++ic)
	{
		const auto& c = cells[ic];

		const size_t i = mesh.getCellId2Index().at(c.id);
		double T = (*TF)[i]; 

			const double rho = CO2::rho_CO2_kg_m3(T);// ρ(T) = 539.7 / T
			const double mu = CO2::mu_CO2_Pa_s(T);   // 220–1000 K 多项式
			const double cp = CO2::cp_mass_J_kgK(T);   // 293–3000 K 分段多项式（质量比热）
			const double k = CO2::k_W_mK(T);    // 220–3273 K 分段多项式

			(*rho_gF)[i] = rho;
			(*cp_gF)[i] = cp;
			(*k_gF)[i] = k;

	}
}

void PhysicalPropertiesManager::ComputeEffectiveThermalProperties_test_varProperties_singlePhase_CO2_T_diffusion(MeshManager& mgr, FieldRegistry& reg, const std::string& Tf_field)
{
	auto& mesh = mgr.mesh();
	const auto& cells = mesh.getCells();
	const size_t n = cells.size();
	//ensureMatrixFluidFields(reg, n);
	//ensureRockFields(reg, n);

	// 岩石参数（必须存在）
	auto phiF = reg.get<volScalarField>("phi");
	auto rrF = reg.get<volScalarField>("rho_r");
	auto cprF = reg.get<volScalarField>("cp_r");
	auto lamrF = reg.get<volScalarField>("lambda_r");

	//流体物性参数
	auto rho_gF = reg.get<volScalarField>("rho_g");
	auto cp_gF = reg.get<volScalarField>("cp_g");
	auto k_gF = reg.get<volScalarField>("k_g");
	// 有效热参数场
	auto Ceff = reg.getOrCreate<volScalarField>("C_eff", n, 0.0);
	auto lame = reg.getOrCreate<volScalarField>("lambda_eff", n, 0.0);

	for (size_t ic = 0; ic < cells.size(); ++ic)
	{
		const auto& c = cells[ic];
		if (c.id < 0) continue;
		const size_t i = mesh.getCellId2Index().at(c.id);

		(*Ceff)[i] = (1.0 - (*phiF)[i]) * (*rrF)[i] * (*cprF)[i] + (*phiF)[i] * (*rho_gF)[i] * (*cp_gF)[i];
		(*lame)[i] = (1.0 - (*phiF)[i]) * (*lamrF)[i] + (*phiF)[i] * (*k_gF)[i];
	}
}



//对流项 常物性测试
void PhysicalPropertiesManager::RockProperties_test_constProperties_singlePhase_CO2(MeshManager& mgr, FieldRegistry& reg)
{
	auto& mesh = mgr.mesh();
	const auto& cells = mesh.getCells();
	const size_t n = cells.size();
	ensureRockFields(reg, n);

	auto phi_r = reg.get<volScalarField>("phi"); //孔隙度
	auto rho_r = reg.get<volScalarField>("rho_r"); //基岩密度，kg/m³
	auto cp_r = reg.get<volScalarField>("cp_r"); //基岩比热容，J/(kg·K)
	auto lam_r = reg.get<volScalarField>("lambda_r"); //基岩导热系数，W/(m·K)
	auto c_phi = reg.get<volScalarField>("c_phi"); //孔隙度可压缩性，1/Pa
	auto k_xx = reg.get<volScalarField>("kxx"); //基岩渗透率，m²
	auto k_yy = reg.get<volScalarField>("kyy");	
	auto k_zz = reg.get<volScalarField>("kzz");

	for (size_t ic = 0; ic < cells.size(); ++ic)
	{
		const auto& cell = cells[ic];
		if (cell.id < 0) continue;
		const size_t i = mesh.getCellId2Index().at(cell.id);
		(*phi_r)[i] = 0.15;
		(*rho_r)[i] = 2650.0;
		(*cp_r)[i] = 1000.0;
		(*lam_r)[i] = 3;
		(*c_phi)[i] = 0;
		(*k_xx)[i] = 1e-14;
		(*k_yy)[i] = 1e-14;
		(*k_zz)[i] = 1e-14;
	}
}

void PhysicalPropertiesManager::CO2Properties_test_constProperties_singlePhase_CO2 (MeshManager& mgr, FieldRegistry& reg)
{
	auto& mesh = mgr.mesh();
	const auto& cells = mesh.getCells();
	const size_t n = cells.size();
	ensureMatrixFluidFields(reg, n);

	auto rho_gF = reg.get<volScalarField>("rho_g");
	auto cp_gF = reg.get<volScalarField>("cp_g");
	auto k_gF = reg.get<volScalarField>("k_g");
	auto mu_gF = reg.get<volScalarField>("mu_g");
	auto Drho_Dp_gF = reg.get<volScalarField>("Drho_Dp_g");

	for (size_t ic = 0; ic < cells.size(); ++ic)
	{
		const auto& c = cells[ic];
		if (c.id < 0) continue;
		const size_t i = mesh.getCellId2Index().at(c.id);
		(*rho_gF)[i] = 800.0;
		(*cp_gF)[i] = 1200.0;
		(*k_gF)[i] = 0.05;
		(*mu_gF)[i] = 1.48e-5;
		(*Drho_Dp_gF)[i] = 0;
	}
}











void PhysicalPropertiesManager::classifyRockRegionsByGeometry(MeshManager& mgr, const vector<RegionGeometry>& regionGeoms, Cell::RegionType defaultRegion = Cell::RegionType::Medium)
{
	auto& mesh = mgr.mesh();
	for (auto& cell : mesh.getCells())
	{
		if (cell.id < 0) continue; 
		Cell::RegionType region = defaultRegion; // 默认区域
		for (auto const& rg : regionGeoms)
		{
			if (pointInConvexPolygon(cell.center, rg.vertices))
			{
				region = rg.type; // 如果点在凸多边形内，设置为对应的区域类型
				break; // 找到匹配的区域后跳出循环
			}
		}
		cell.region = region; // 更新单元的区域类型
	}
}

void PhysicalPropertiesManager::classifyFractureElementsByGeometry(MeshManager& mgr, int fracID, const Vector& regionStart, const Vector& regionEnd, FractureElementType insideType, FractureElementType outsideType)
{

	// 取出裂缝网络
	auto& fn = mgr.fracture_network();
	if (fracID < 0 || fracID >= (int)fn.fractures.size()) return;
	auto& frac = fn.fractures[fracID];

	// 1) 计算整条裂缝的方向向量和长度平方 L2
	Vector f0 = frac.start;
	Vector f1 = frac.end;
	Vector dir{ f1.m_x - f0.m_x,
				 f1.m_y - f0.m_y,
				 f1.m_z - f0.m_z };  //方向向量
	double L2 = dir.m_x * dir.m_x + dir.m_y * dir.m_y + dir.m_z * dir.m_z; //长度平方
	if (L2 <= 0) return;

	// 2) 计算 regionStart, regionEnd 在这条裂缝上的归一化参数 t0, t1
	auto paramOf = [&](const Vector& P)
		{
			Vector v{ P.m_x - f0.m_x,
					  P.m_y - f0.m_y,
					  P.m_z - f0.m_z };
			double dot = v.m_x * dir.m_x + v.m_y * dir.m_y + v.m_z * dir.m_z;
			return dot / L2;  // 0→1 范围
		};
	double t0 = paramOf(regionStart);
	double t1 = paramOf(regionEnd);
	if (t1 < t0) std::swap(t0, t1);

	// 3) 对每段，根据它的 param0/param1 中点决定标签
	for (auto& elem : frac.elements) {
		double mid = 0.5 * (elem.param0 + elem.param1);
		if (mid >= t0 && mid <= t1)
			elem.type = insideType;
		else
			elem.type = outsideType;
	}

}

void PhysicalPropertiesManager::classifyRockRegions(MeshManager& mgr, double poroLow, double poroHigh, double permLow, double permHigh)    //输入低渗孔隙率、渗透率阈值&高渗孔隙率、渗透率，现在这个函数功能还不能自适应
{
	auto& mesh = mgr.mesh();
	for (auto& cell : mesh.getCells())
	{
		if (cell.id < 0) continue; // 跳过 Ghost Cell
		const auto& solidproperties = cell.SolidMaterialProps;
		if (solidproperties.porosity <= poroLow && solidproperties.permeability <= permLow)
		{
			cell.region = Cell::RegionType::Low;
		}
		else if (solidproperties.porosity >= poroHigh && solidproperties.permeability >= permHigh)
		{
			cell.region = Cell::RegionType::High;
		}

		else
			cell.region = Cell::RegionType::Medium;

	
	}
}

void PhysicalPropertiesManager::classifyFractureElements(MeshManager& mgr,double permThreshold)
{
	auto& fn = const_cast<FractureNetwork&>(mgr.fracture_network());
	for (auto& frac : fn.fractures) 
	{
		for (auto& elem : frac.elements)
		{
			if (elem.solidProps.permeability< permThreshold)
			{
				elem.type = FractureElementType::Blocking;
			}
			else
			{
                elem.type = FractureElementType::Conductive;
			}
		}

	}
}






//***************************基岩物性参数计算与赋值**************************//
//固相参数
void PhysicalPropertiesManager::UpdateMatrixRockAt(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field)
{
	auto& mesh = mgr.mesh();
	const auto& cells = mesh.getCells();
	const size_t n = cells.size();

	ensureRockFields(reg, n);

	auto pF = reg.get<volScalarField>(p_field);
	auto TF = reg.get<volScalarField>(T_field);
	if (!pF || !TF)
	{
		std::cerr << "[PPM][Rock] missing fields '" << p_field << "' or '" << T_field << "'.\n";
		return;
	}
	auto phi_r = reg.get<volScalarField>("phi"); //孔隙度
	auto kxx_r = reg.get<volScalarField>("kxx"); //渗透率 xx方向
	auto kyy_r = reg.get<volScalarField>("kyy"); //渗透率 yy方向
	auto kzz_r = reg.get<volScalarField>("kzz"); //渗透率 zz方向
	auto rho_r = reg.get<volScalarField>("rho_r"); //基岩密度，kg/m³
	auto cp_r = reg.get<volScalarField>("cp_r"); //基岩比热容，J/(kg·K)
	auto lam_r = reg.get<volScalarField>("lambda_r"); //基岩导热系数，W/(m·K)
	
	//取消并行计算
	for (size_t ic = 0; ic < cells.size(); ++ic)
	{
		const auto& cell = cells[ic];
		if (cell.id < 0) continue;
		const size_t i = mesh.getCellId2Index().at(cell.id);
		double P = (*pF)[i], T = (*TF)[i];
		//Initializer::clampPT(P, T);
		const auto sp = rock::computeSolidProperties(cell.region, P, T);
		(*phi_r)[i] = sp.porosity;
		(*kxx_r)[i] = sp.permeability;
		(*kyy_r)[i] = sp.permeability;
		(*kzz_r)[i] = sp.permeability;
		(*rho_r)[i] = sp.rho_s;
		(*cp_r)[i] = sp.cp_s;
		(*lam_r)[i] = sp.k_s;

	}


//#pragma omp parallel for schedule(static)
//	for (int ic = 0; ic < static_cast<int>(cells.size()); ++ic) {
//		const auto& cell = cells[ic];
//		if (cell.id < 0) continue;
//		const size_t i = mesh.getCellId2Index().at(cell.id);
//		double P = (*pF)[i], T = (*TF)[i];
//		Initializer::clampPT(P, T);
//		const auto sp = rock::computeSolidProperties(cell.region, P, T);
//		(*phi_r)[i] = sp.porosity;
//		(*kxx_r)[i] = sp.permeability;
//		(*kyy_r)[i] = sp.permeability;
//		(*kzz_r)[i] = sp.permeability;
//		(*rho_r)[i] = sp.rho_s;
//		(*cp_r)[i] = sp.cp_s;
//		(*lam_r)[i] = sp.k_s;
//	}
}







//流体相参数
void PhysicalPropertiesManager::UpdateMatrixFluidAt(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field, const std::string& phase)
{
	auto& mesh = mgr.mesh();
	const auto& cells = mesh.getCells();
	const size_t n = cells.size();

	ensureMatrixFluidFields(reg, n);

	auto pF = reg.get<volScalarField>(p_field);
	auto TF = reg.get<volScalarField>(T_field);
	if (!pF || !TF) {
		std::cerr << "[PPM][Fluid] missing fields '" << p_field << "' or '" << T_field << "'.\n";
		return;
	}
	auto rho_wF = reg.get<volScalarField>("rho_w");
	auto mu_wF = reg.get<volScalarField>("mu_w");
	auto cp_wF = reg.get<volScalarField>("cp_w");
	auto k_wF = reg.get<volScalarField>("k_w");

	auto rho_gF = reg.get<volScalarField>("rho_g");
	auto mu_gF = reg.get<volScalarField>("mu_g");
	auto cp_gF = reg.get<volScalarField>("cp_g");
	auto k_gF = reg.get<volScalarField>("k_g");

	auto wt = WaterPropertyTable::instance();
	//auto gt = CO2PropertyTable::instance();

	std::size_t oor = 0, ghost = 0;
	const bool doW = (phase == "water" || phase == "both");
	const bool doG = (phase == "CO2" || phase == "both");

	//取消并行计算
	for (size_t ic = 0; ic < cells.size(); ++ic)
	{
		const auto& c = cells[ic];
		if (c.id < 0) { ++ghost; continue; }
		const size_t i = mesh.getCellId2Index().at(c.id);
		double p = (*pF)[i], T = (*TF)[i]; //Initializer::clampPT(p, T);
		if (doW) {
			double rho = 1000, mu = 1e-3, cp = 4200, k = 0.6;
			try { const auto W = wt.getProperties(p, T); rho = W.rho; mu = W.mu; cp = W.cp; k = W.k; }
			catch (...) { ++oor; }
			(*rho_wF)[i] = rho; (*mu_wF)[i] = mu; (*cp_wF)[i] = cp; (*k_wF)[i] = k;
		}
		if (doG)
		{
			/*double rho = 1.98, mu = 1.48e-5, cp = 846, k = 0.0146;
			try { const auto G = gt.getProperties(p, T); rho = G.rho; mu = G.mu; cp = G.cp; k = G.k; }
			catch (...) { ++oor; }
			(*rho_gF)[i] = rho; (*mu_gF)[i] = mu; (*cp_gF)[i] = cp; (*k_gF)[i] = k;*/
			const double rho = CO2::rho_CO2_kg_m3(T);// ρ(T) = 539.7 / T
			const double mu = CO2::mu_CO2_Pa_s(T);   // 220–1000 K 多项式
			const double cp = CO2::cp_mass_J_kgK(T);   // 293–3000 K 分段多项式（质量比热）
			const double k = CO2::k_W_mK(T);    // 220–3273 K 分段多项式

			(*rho_gF)[i] = rho;
			(*mu_gF)[i] = mu;
			(*cp_gF)[i] = cp;
			(*k_gF)[i] = k;

		}
	}

//#pragma omp parallel for schedule(static)
//	for (int ic = 0; ic < static_cast<int>(cells.size()); ++ic) {
//		const auto& c = cells[ic];
//		if (c.id < 0) { ++ghost; continue; }
//		const size_t i = mesh.getCellId2Index().at(c.id);
//		double p = (*pF)[i], T = (*TF)[i]; Initializer::clampPT(p, T);
//		if (doW) {
//			double rho = 1000, mu = 1e-3, cp = 4200, k = 0.6;
//			try { const auto W = wt.getProperties(p, T); rho = W.rho; mu = W.mu; cp = W.cp; k = W.k; }
//			catch (...) { ++oor; }
//			(*rho_wF)[i] = rho; (*mu_wF)[i] = mu; (*cp_wF)[i] = cp; (*k_wF)[i] = k;
//		}
//		if (doG) 
//		{
//			/*double rho = 1.98, mu = 1.48e-5, cp = 846, k = 0.0146;
//			try { const auto G = gt.getProperties(p, T); rho = G.rho; mu = G.mu; cp = G.cp; k = G.k; }
//			catch (...) { ++oor; }
//			(*rho_gF)[i] = rho; (*mu_gF)[i] = mu; (*cp_gF)[i] = cp; (*k_gF)[i] = k;*/
//			const double rho = CO2::rho_CO2_kg_m3(T);// ρ(T) = 539.7 / T
//			const double mu = CO2::mu_CO2_Pa_s(T);   // 220–1000 K 多项式
//			const double cp = CO2::cp_mass_J_kgK(T);   // 293–3000 K 分段多项式（质量比热）
//			const double k = CO2::k_W_mK(T);    // 220–3273 K 分段多项式
//
//			(*rho_gF)[i] = rho;
//			(*mu_gF)[i] = mu;
//			(*cp_gF)[i] = cp;
//			(*k_gF)[i] = k;
//
//		}
//	}

	if (ghost) std::cout << "[PPM][Fluid] skipped ghost=" << ghost << "\n";
	if (oor)   std::cout << "[PPM][Fluid] table OOR hits=" << oor << "\n";
}

// ====== 矩阵：单相和多相 C_eff / lambda_eff ======
void PhysicalPropertiesManager::ComputeMatrixEffectiveThermalsAt( MeshManager& mgr, FieldRegistry& reg,const std::string& p_field, const std::string& T_field,const std::string& phase, double Ceff_floor)
{
	auto& mesh = mgr.mesh();
	const auto& cells = mesh.getCells();
	const auto& id2 = mesh.getCellId2Index();
	const size_t n = cells.size();

	// 岩石参数（必须存在）
	auto phiF = reg.get<volScalarField>("phi");
	auto rrF = reg.get<volScalarField>("rho_r");
	auto cprF = reg.get<volScalarField>("cp_r");
	auto lamrF = reg.get<volScalarField>("lambda_r");
	if (!phiF || !rrF || !cprF || !lamrF) {
		std::cerr << "[Thermal] missing rock fields: phi/rho_r/cp_r/lambda_r.\n";
		return;
	}

	// 流体参数（优先用场；若缺失就兜底查表）
	auto rwF = reg.get<volScalarField>("rho_w");
	auto cwF = reg.get<volScalarField>("cp_w");
	auto kwF = reg.get<volScalarField>("k_w");

	auto rgF = reg.get<volScalarField>("rho_g");
	auto cgF = reg.get<volScalarField>("cp_g");
	auto kgF = reg.get<volScalarField>("k_g");

	// 主变量（兜底查表用）
	auto pF = reg.get<volScalarField>(p_field);
	auto TF = reg.get<volScalarField>(T_field);
	if (!pF || !TF) {
		std::cerr << "[Thermal] missing p/T fields: '" << p_field << "' or '" << T_field << "'.\n";
		return;
	}

	// 饱和度（两相才用）
	auto SwF = reg.get<volScalarField>("S_w");

	// 结果场
	auto Ceff = reg.getOrCreate<volScalarField>("C_eff", n, 0.0);
	auto lame = reg.getOrCreate<volScalarField>("lambda_eff", n, 0.0);


	// 物性表（仅在场缺失时兜底）
	auto wt = WaterPropertyTable::instance();
	auto gt = CO2PropertyTable::instance();

	// 规范化 phase（不区分大小写）
	auto tolower_str = [](std::string s) { for (auto& c : s) c = char(::tolower(c)); return s; };
	const std::string ph = tolower_str(phase);
	const bool doW = (ph == "water" || ph == "both");
	const bool doG = (ph == "co2" || ph == "both");

	//取消并行计算
	for (size_t ic = 0; ic < cells.size(); ++ic)
	{
		const auto& c = cells[ic];
		if (c.id < 0) continue;
		const size_t i = id2.at(c.id);

		// 岩石
		const double phi = std::min(1.0, std::max(0.0, (*phiF)[i]));	//孔隙度
		const double rr = std::max(0.0, (*rrF)[i]);						//基岩密度	
		const double cpr = std::max(0.0, (*cprF)[i]);					//基岩比热容
		const double lamr = std::max(0.0, (*lamrF)[i]);					//基岩导热系数

		// p/T（仅用于查表兜底）
		double p = (*pF)[i], T = (*TF)[i];
		//Initializer::clampPT(p, T);

		// 水相
		double rw = 0.0, cw = 0.0, kw = 0.0;
		if (doW) {
			if (rwF && cwF && kwF) { rw = std::max(0.0, (*rwF)[i]); cw = std::max(0.0, (*cwF)[i]); kw = std::max(0.0, (*kwF)[i]); }
			else { // 场缺失→查表
				try { const auto W = wt.getProperties(p, T); rw = W.rho; cw = W.cp; kw = W.k; }
				catch (...) { rw = 1000; cw = 4200; kw = 0.6; }
			}
		}

		// 气相（CO2）
		double rg = 0.0, cg = 0.0, kg = 0.0;
		if (doG) {
			if (rgF && cgF && kgF) { rg = std::max(0.0, (*rgF)[i]); cg = std::max(0.0, (*cgF)[i]); kg = std::max(0.0, (*kgF)[i]); }
			else {
				rg = CO2::rho_CO2_kg_m3(T);
				cg = CO2::cp_mass_J_kgK(T);
				kg = CO2::k_W_mK(T);
			}
		}

		// 饱和度：两相用；若没有 S_w 则默认退化为单相水（Sw=1）
		double Sw = 1.0;
		if (ph == "both") {
			if (SwF) Sw = std::min(1.0, std::max(0.0, (*SwF)[i]));
			else     Sw = 1.0; // 无 S_w → 退化单相水
		}
		else if (ph == "water") {
			Sw = 1.0;
		}
		else if (ph == "co2") {
			Sw = 0.0;
		}
		const double Sg = 1.0 - Sw;

		// 有效体积热容 C_eff
		double Cfluid = 0.0;
		if (doW && doG)      Cfluid = phi * (Sw * rw * cw + Sg * rg * cg);
		else if (doW)        Cfluid = phi * (rw * cw);
		else if (doG)        Cfluid = phi * (rg * cg);
		const double C = (1.0 - phi) * rr * cpr + Cfluid;
		(*Ceff)[i] = std::max(Ceff_floor, C);

		// 有效导热系数 λ_eff（体积分数线性混合；需要更复杂模型可再替换）
		double kfluid = 0.0;
		if (doW && doG)      kfluid = Sw * kw + Sg * kg;
		else if (doW)        kfluid = kw;
		else if (doG)        kfluid = kg;
		(*lame)[i] = (1.0 - phi) * lamr + phi * std::max(0.0, kfluid);

	}

//#pragma omp parallel for schedule(static)
//	for (int ic = 0; ic < static_cast<int>(cells.size()); ++ic) 
//	{
//		const auto& c = cells[ic];
//		if (c.id < 0) continue;
//		const size_t i = id2.at(c.id);
//
//		// 岩石
//		const double phi = std::min(1.0, std::max(0.0, (*phiF)[i]));	//孔隙度
//		const double rr = std::max(0.0, (*rrF)[i]);						//基岩密度	
//		const double cpr = std::max(0.0, (*cprF)[i]);					//基岩比热容
//		const double lamr = std::max(0.0, (*lamrF)[i]);					//基岩导热系数
//
//		// p/T（仅用于查表兜底）
//		double p = (*pF)[i], T = (*TF)[i];
//		Initializer::clampPT(p, T);
//
//		// 水相
//		double rw = 0.0, cw = 0.0, kw = 0.0;
//		if (doW) {
//			if (rwF && cwF && kwF) { rw = std::max(0.0, (*rwF)[i]); cw = std::max(0.0, (*cwF)[i]); kw = std::max(0.0, (*kwF)[i]); }
//			else { // 场缺失→查表
//				try { const auto W = wt.getProperties(p, T); rw = W.rho; cw = W.cp; kw = W.k; }
//				catch (...) { rw = 1000; cw = 4200; kw = 0.6; }
//			}
//		}
//
//		// 气相（CO2）
//		double rg = 0.0, cg = 0.0, kg = 0.0;
//		if (doG) {
//			if (rgF && cgF && kgF) { rg = std::max(0.0, (*rgF)[i]); cg = std::max(0.0, (*cgF)[i]); kg = std::max(0.0, (*kgF)[i]); }
//			else {
//				rg = CO2::rho_CO2_kg_m3(T);
//				cg = CO2::cp_mass_J_kgK(T);
//				kg = CO2::k_W_mK(T);
//			}
//		}
//
//		// 饱和度：两相用；若没有 S_w 则默认退化为单相水（Sw=1）
//		double Sw = 1.0;
//		if (ph == "both") {
//			if (SwF) Sw = std::min(1.0, std::max(0.0, (*SwF)[i]));
//			else     Sw = 1.0; // 无 S_w → 退化单相水
//		}
//		else if (ph == "water") {
//			Sw = 1.0;
//		}
//		else if (ph == "co2") {
//			Sw = 0.0;
//		}
//		const double Sg = 1.0 - Sw;
//
//		// 有效体积热容 C_eff
//		double Cfluid = 0.0;
//		if (doW && doG)      Cfluid = phi * (Sw * rw * cw + Sg * rg * cg);
//		else if (doW)        Cfluid = phi * (rw * cw);
//		else if (doG)        Cfluid = phi * (rg * cg);
//		const double C = (1.0 - phi) * rr * cpr + Cfluid;
//		(*Ceff)[i] = std::max(Ceff_floor, C);
//
//		// 有效导热系数 λ_eff（体积分数线性混合；需要更复杂模型可再替换）
//		double kfluid = 0.0;
//		if (doW && doG)      kfluid = Sw * kw + Sg * kg;
//		else if (doW)        kfluid = kw;
//		else if (doG)        kfluid = kg;
//		(*lame)[i] = (1.0 - phi) * lamr + phi * std::max(0.0, kfluid);
//	}
}

// ====== 裂缝：固相（这里先复用你原有接口逻辑，如果后续需要按 p/T 变化再扩展）======

void PhysicalPropertiesManager::UpdateFractureRockAt (MeshManager& mgr, FieldRegistry& reg_fr, FieldRegistry& reg, const std::string& pf_field, const std::string& Tf_field)
{
	const FractureNetwork& frNet = mgr.fracture_network(); // 取出裂缝网络
	//调用裂缝段索引
	const auto idx = buildFracElemIndex(frNet);
	const size_t ne = idx.total;

	if (!ne) {
		std::cout << "[PPM] No fracture elements. Skip InitializeFractureElementsProperties.\n";
		return;
	}

	//确保裂缝主变量场&裂缝固相场存在
	ensureFracPrimaryFields(reg_fr, ne);
	ensureFracRockFields(reg_fr, ne);

	// 取出计算裂缝物性参数需要的主变量场指针
	auto pfw = reg_fr.get<volScalarField>(pf_field);
	auto Tf = reg_fr.get<volScalarField>(Tf_field);
	// 取出裂缝物性参数场指针
	auto fr_phi = reg_fr.get<volScalarField>("fr_phi");
	auto fr_k_t = reg_fr.get<volScalarField>("fr_k_t");
	auto fr_k_n = reg_fr.get<volScalarField>("fr_k_n");
	auto fr_rho_r = reg_fr.get<volScalarField>("fr_rho_r");
	auto fr_cp_r = reg_fr.get<volScalarField>("fr_cp_r");
	auto fr_lam_r = reg_fr.get<volScalarField>("fr_lambda_r");
	auto fr_b = reg_fr.get<volScalarField>("fr_aperture");

	// 遍历所有裂缝段
	for (size_t f = 0; f < frNet.fractures.size(); ++f)
	{
		const auto& F = frNet.fractures[f];
		const size_t base = idx.offset[f]; //本条裂缝的全局起点
		for (size_t e = 0; e < F.elements.size(); ++e)
		{
			const size_t g = base + e; //裂缝段全局索引
			const auto& elem = F.elements[e]; //为了取出 elem.type

			double P = (*pfw)[g]; // 裂缝段水相压力，Pa
			double T = (*Tf)[g];  // 裂缝段温度，K

			const auto sp = fracture::computeSolidProperties(elem.type, P, T);
			(*fr_phi)[g] = sp.porosity;
			(*fr_k_t)[g] = sp.permeability; //切向等效渗透率
			(*fr_k_n)[g] = sp.permeability * 1e-4; //法向等效渗透率，假设比切向小4个数量级
			(*fr_rho_r)[g] = sp.rho_s;
			(*fr_cp_r)[g] = sp.cp_s;
			(*fr_lam_r)[g] = sp.k_s;
			(*fr_b)[g] = sp.aperture;
		}
	}
}

// ====== 裂缝：流体 ======
void PhysicalPropertiesManager::UpdateFractureFluidAt( MeshManager& mgr, FieldRegistry& reg, FieldRegistry& reg_fr,const std::string& p_field_fr, const std::string& T_field_fr,const std::string& phase)
{
	auto& mesh = mgr.mesh();
	// 统计裂缝段总数
	size_t Nseg = 0;
	for (auto& F : mgr.fracture_network().fractures) Nseg += F.elements.size();
	ensureFractureFluidFields(reg_fr, Nseg);

	auto pF = reg_fr.get<volScalarField>(p_field_fr);
	auto TF = reg_fr.get<volScalarField>(T_field_fr);
	if (!pF || !TF) { std::cerr << "[PPM][FrFluid] missing fields.\n"; return; }

	auto fr_rho_w = reg_fr.get<volScalarField>("fr_rho_w");
	auto fr_mu_w = reg_fr.get<volScalarField>("fr_mu_w");
	auto fr_cp_w = reg_fr.get<volScalarField>("fr_cp_w");
	auto fr_k_w = reg_fr.get<volScalarField>("fr_k_w");

	auto fr_rho_g = reg_fr.get<volScalarField>("fr_rho_g");
	auto fr_mu_g = reg_fr.get<volScalarField>("fr_mu_g");
	auto fr_cp_g = reg_fr.get<volScalarField>("fr_cp_g");
	auto fr_k_g = reg_fr.get<volScalarField>("fr_k_g");

	auto wt = WaterPropertyTable::instance();
	//auto gt = CO2PropertyTable::instance();

	const bool doW = (phase == "water" || phase == "both");
	const bool doG = (phase == "CO2" || phase == "both");

	size_t gid = 0;
	for (auto& F : mgr.fracture_network().fractures) {
		for (auto& E : F.elements) {
			double p = (*pF)[gid], T = (*TF)[gid]; //Initializer::clampPT(p, T);
			if (doW) {
				double rho = 1000, mu = 1e-3, cp = 4200, k = 0.6;
				try { const auto W = wt.getProperties(p, T); rho = W.rho; mu = W.mu; cp = W.cp; k = W.k; }
				catch (...) {}
				(*fr_rho_w)[gid] = rho; (*fr_mu_w)[gid] = mu; (*fr_cp_w)[gid] = cp; (*fr_k_w)[gid] = k;
			}
			if (doG)
			{
				const double rho = CO2::rho_CO2_kg_m3(T);// ρ(T) = 539.7 / T
				const double mu = CO2::mu_CO2_Pa_s(T);   // 220–1000 K 多项式
				const double cp = CO2::cp_mass_J_kgK(T);   // 293–3000 K 分段多项式（质量比热）
				const double k = CO2::k_W_mK(T);    // 220–3273 K 分段多项式
				(*fr_rho_g)[gid] = rho;
				(*fr_mu_g)[gid] = mu;
				(*fr_cp_g)[gid] = cp;
				(*fr_k_g)[gid] = k;
			}
			++gid;
		}
	}
}

// ====== 裂缝：单相有效热（如需对裂缝解温度时用；若暂不解可忽略）======
void PhysicalPropertiesManager::ComputeFractureEffectiveThermalsAt( MeshManager& mgr, FieldRegistry& reg_fr, FieldRegistry& reg,const std::string& p_field_fr, const std::string& T_field_fr,const std::string& phase, double Ceff_floor)
{
	// 如果裂缝也解温度，可在 reg_fr 里维护 fr_C_eff / fr_lambda_eff；
	// 此处略——和矩阵版本一致，把岩石参数换成裂缝等效属性即可。
	(void)mgr; (void)reg_fr; (void)reg; (void)p_field_fr; (void)T_field_fr;
	(void)phase; (void)Ceff_floor;
}


//************向后兼容就接口，保持你原有的初始化逻辑************//
//void PhysicalPropertiesManager::InitializeRockMatrixProperties(MeshManager& mgr, FieldRegistry& reg)
//{
//	UpdateMatrixRockAt(mgr, reg, "p_w", "T");
//}
//
//void PhysicalPropertiesManager::UpdateMatrixProperties(MeshManager& mgr, FieldRegistry& reg) 
//{
//	UpdateMatrixRockAt(mgr, reg, "p_w", "T"); // 或者 "p_w_prev","T_prev" 视你的时序
//}
//
//void PhysicalPropertiesManager::InitializeMatrixFluidProperties( MeshManager& mgr, FieldRegistry& reg, const VGParams&) 
//{
//	UpdateMatrixFluidAt(mgr, reg, "p_w", "T", "both"); // 原先你会同时填 w/g
//}
//
//void PhysicalPropertiesManager::UpdateMatrixFluidProperties(MeshManager& mgr, FieldRegistry& reg, const VGParams&) 
//{
//	UpdateMatrixFluidAt(mgr, reg, "p_w", "T", "both");
//}
//
//void PhysicalPropertiesManager::InitializeFractureElementsProperties(MeshManager& mgr, FieldRegistry& reg_fr)
//{
//
//	const FractureNetwork& frNet = mgr.fracture_network(); // 取出裂缝网络
//	//调用裂缝段索引
//	const auto idx = buildFracElemIndex(frNet);
//	const size_t ne = idx.total;
//
//	if (!ne) {
//		std::cout << "[PPM] No fracture elements. Skip InitializeFractureElementsProperties.\n";
//		return;
//	}
//
//	//确保裂缝主变量场&裂缝固相场存在
//	ensureFracPrimaryFields(reg_fr, ne);
//	ensureFracRockFields(reg_fr, ne);
//
//	// 取出计算裂缝物性参数需要的主变量场指针
//	auto pfw = reg_fr.get<volScalarField>("pf_w");
//	auto Tf = reg_fr.get<volScalarField>("Tf");
//	// 取出裂缝物性参数场指针
//	auto fr_phi = reg_fr.get<volScalarField>("fr_phi");
//	auto fr_k_t = reg_fr.get<volScalarField>("fr_k_t");
//	auto fr_k_n = reg_fr.get<volScalarField>("fr_k_n");
//	auto fr_rho_r = reg_fr.get<volScalarField>("fr_rho_r");
//	auto fr_cp_r = reg_fr.get<volScalarField>("fr_cp_r");
//	auto fr_lam_r = reg_fr.get<volScalarField>("fr_lambda_r");
//	auto fr_b = reg_fr.get<volScalarField>("fr_aperture");
//
//	// 遍历所有裂缝段
//	for (size_t f = 0; f < frNet.fractures.size(); ++f)
//	{
//		const auto& F = frNet.fractures[f];
//		const size_t base = idx.offset[f]; //本条裂缝的全局起点
//		for (size_t e = 0; e < F.elements.size(); ++e)
//		{
//			const size_t g = base + e; //裂缝段全局索引
//			const auto& elem = F.elements[e]; //为了取出 elem.type
//
//			double P = (*pfw)[g]; // 裂缝段水相压力，Pa
//			double T = (*Tf)[g];  // 裂缝段温度，K
//
//			const auto sp = fracture::computeSolidProperties(elem.type, P, T);
//			(*fr_phi)[g] = sp.porosity;
//			(*fr_k_t)[g] = sp.permeability; //切向等效渗透率
//			(*fr_k_n)[g] = sp.permeability * 1e-4; //法向等效渗透率，假设比切向小4个数量级
//			(*fr_rho_r)[g] = sp.rho_s;
//			(*fr_cp_r)[g] = sp.cp_s;
//			(*fr_lam_r)[g] = sp.k_s;
//			(*fr_b)[g] = sp.aperture;
//		}
//	}
//}
//
//void PhysicalPropertiesManager::UpdateFractureSolidProperties(MeshManager& mgr, FieldRegistry& reg_fr) {
//	InitializeFractureElementsProperties(mgr, reg_fr);
//}
//
//void PhysicalPropertiesManager::InitializeFractureFluidProperties(MeshManager& mgr, FieldRegistry& reg, FieldRegistry& reg_fr, const VGParams&) {
//	UpdateFractureFluidAt(mgr, reg, reg_fr, "pf_w", "Tf", "both");
//}
//
//void PhysicalPropertiesManager::UpdateFractureFluidProperties(MeshManager& mgr, FieldRegistry& reg, FieldRegistry& reg_fr, const VGParams&) {
//	UpdateFractureFluidAt(mgr, reg, reg_fr, "pf_w", "Tf", "both");
//}




//----------------------------------输出&调试----------------------------------------//

// 小工具：安全读取场值（不存在则给默认）
static inline double getOr(const FieldRegistry& reg,
	const std::string& name,
	std::size_t i,
	double def = std::numeric_limits<double>::quiet_NaN())
{
	auto f = reg.get<volScalarField>(name);
	return f ? (*f)[i] : def;
}

static inline double getFrOr(const FieldRegistry& reg_fr,
	const std::string& name,
	std::size_t g,
	double def = std::numeric_limits<double>::quiet_NaN())
{
	auto f = reg_fr.get<volScalarField>(name);
	return f ? (*f)[g] : def;
}



void PhysicalPropertiesManager::debugPrintProperties(MeshManager& mgr,
	const FieldRegistry& reg,
	const FieldRegistry& reg_fr,
	std::size_t maxPrint) const
{
	auto& mesh = mgr.mesh();
	const auto& cells = mesh.getCells();
	const auto& id2idx = mesh.getCellId2Index();

	std::cout << std::fixed << std::setprecision(6);

	// ===================== 基岩单元 =====================
	std::cout << "\n=== Matrix Cells (up to " << maxPrint << ") ===\n";
	std::size_t printed = 0;
	for (const auto& c : cells) {
		if (c.id < 0) continue;     // ghost
		const auto it = id2idx.find(c.id);
		if (it == id2idx.end()) continue;
		const std::size_t i = it->second;

		if (printed++ >= maxPrint) { std::cout << "... (truncated)\n"; break; }

		// 主变量
		double pw = getOr(reg, "p_w", i);
		double Sw = getOr(reg, "S_w", i);
		double T = getOr(reg, "T", i);
		double pc = getOr(reg, "p_c", i);
		double pg = getOr(reg, "p_g", i);
		double krw = getOr(reg, "kr_w", i);
		double krg = getOr(reg, "kr_g", i);

		// 固相（各向异性渗透率用 kxx/kyy/kzz）
		double phi = getOr(reg, "phi", i);
		double kxx = getOr(reg, "kxx", i);
		double kyy = getOr(reg, "kyy", i);
		double kzz = getOr(reg, "kzz", i);
		double rho_r = getOr(reg, "rho_r", i);
		double cp_r = getOr(reg, "cp_r", i);
		double lam_r = getOr(reg, "lambda_r", i);

		// 流体（水/CO2）
		double rho_w = getOr(reg, "rho_w", i);
		double mu_w = getOr(reg, "mu_w", i);
		double cp_w = getOr(reg, "cp_w", i);
		double k_w = getOr(reg, "k_w", i);

		double rho_g = getOr(reg, "rho_g", i);
		double mu_g = getOr(reg, "mu_g", i);
		double cp_g = getOr(reg, "cp_g", i);
		double k_g = getOr(reg, "k_g", i);

		// 有效热
		double Ceff = getOr(reg, "C_eff", i);
		double lame = getOr(reg, "lambda_eff", i);

		std::cout
			<< "Cell " << c.id
			<< " | region=" << static_cast<int>(c.region)
			<< " | center=(" << c.center.m_x << "," << c.center.m_y << "," << c.center.m_z << ")\n"
			<< "  Primaries: p_w=" << pw << " Pa, S_w=" << Sw
			<< ", T=" << T << " K, p_c=" << pc << " Pa, p_g=" << pg
			<< ", kr_w=" << krw << ", kr_g=" << krg << "\n"
			<< "  Rock: phi=" << phi
			<< ", kxx=" << kxx << ", kyy=" << kyy << ", kzz=" << kzz
			<< ", rho_r=" << rho_r << " kg/m^3"
			<< ", cp_r=" << cp_r << " J/(kg·K)"
			<< ", lambda_r=" << lam_r << " W/(m·K)\n"
			<< "  Fluids-W: rho_w=" << rho_w << " kg/m^3"
			<< ", mu_w=" << mu_w << " Pa·s"
			<< ", cp_w=" << cp_w << " J/(kg·K)"
			<< ", k_w=" << k_w << " W/(m·K)\n"
			<< "  Fluids-CO2: rho_g=" << rho_g << " kg/m^3"
			<< ", mu_g=" << mu_g << " Pa·s"
			<< ", cp_g=" << cp_g << " J/(kg·K)"
			<< ", k_g=" << k_g << " W/(m·K)\n"
			<< "  Effective: C_eff=" << Ceff << " J/(m^3·K)"
			<< ", lambda_eff=" << lame << " W/(m·K)\n";
	}

	// ===================== 裂缝段 =====================
	const auto& frNet = mgr.fracture_network();
	const auto idx = buildFracElemIndex(frNet);
	std::cout << "\n=== Fracture Elements (total " << idx.total
		<< ", showing up to " << maxPrint << ") ===\n";

	std::size_t g = 0, printedFr = 0;
	for (std::size_t f = 0; f < frNet.fractures.size(); ++f) {
		const auto& F = frNet.fractures[f];
		for (std::size_t e = 0; e < F.elements.size(); ++e, ++g) {
			if (printedFr++ >= maxPrint) { std::cout << "... (truncated)\n"; goto DONE_FR; }

			const auto& elem = F.elements[e];

			// 裂缝主变量
			double pfw = getFrOr(reg_fr, "pf_w", g);
			double Sfw = getFrOr(reg_fr, "Sf_w", g);
			double Tf = getFrOr(reg_fr, "Tf", g);

			// 裂缝固相
			double fr_phi = getFrOr(reg_fr, "fr_phi", g);
			double fr_kt = getFrOr(reg_fr, "fr_k_t", g);
			double fr_kn = getFrOr(reg_fr, "fr_k_n", g);
			double fr_rho = getFrOr(reg_fr, "fr_rho_r", g);
			double fr_cp = getFrOr(reg_fr, "fr_cp_r", g);
			double fr_lam = getFrOr(reg_fr, "fr_lambda_r", g);
			double fr_b = getFrOr(reg_fr, "fr_aperture", g);

			// 裂缝流体
			double fr_rho_w = getFrOr(reg_fr, "fr_rho_w", g);
			double fr_mu_w = getFrOr(reg_fr, "fr_mu_w", g);
			double fr_cp_w = getFrOr(reg_fr, "fr_cp_w", g);
			double fr_k_w = getFrOr(reg_fr, "fr_k_w", g);

			double fr_rho_g = getFrOr(reg_fr, "fr_rho_g", g);
			double fr_mu_g = getFrOr(reg_fr, "fr_mu_g", g);
			double fr_cp_g = getFrOr(reg_fr, "fr_cp_g", g);
			double fr_k_g = getFrOr(reg_fr, "fr_k_g", g);

			std::cout
				<< "Frac " << frNet.fractures[f].id
				<< " | Elem local=" << e << ", global=" << g
				<< " | hostCell=" << elem.cellID
				<< " | type=" << static_cast<int>(elem.type) << "\n"
				<< "  Primaries: pf_w=" << pfw << " Pa, Sf_w=" << Sfw
				<< ", Tf=" << Tf << " K\n"
				<< "  Rock(fr): phi=" << fr_phi
				<< ", k_t=" << fr_kt << ", k_n=" << fr_kn
				<< ", rho_r=" << fr_rho << " kg/m^3"
				<< ", cp_r=" << fr_cp << " J/(kg·K)"
				<< ", lambda_r=" << fr_lam << " W/(m·K)"
				<< ", aperture=" << fr_b << " m\n"
				<< "  Fluids(fr)-W: rho_w=" << fr_rho_w << " kg/m^3"
				<< ", mu_w=" << fr_mu_w << " Pa·s"
				<< ", cp_w=" << fr_cp_w << " J/(kg·K)"
				<< ", k_w=" << fr_k_w << " W/(m·K)\n"
				<< "  Fluids(fr)-CO2: rho_g=" << fr_rho_g << " kg/m^3"
				<< ", mu_g=" << fr_mu_g << " Pa·s"
				<< ", cp_g=" << fr_cp_g << " J/(kg·K)"
				<< ", k_g=" << fr_k_g << " W/(m·K)\n";
		}
	}
DONE_FR:
	std::cout << std::endl;
}


// ======================= 裂缝：初始化/更新 流体物性场 =======================//


void PhysicalPropertiesManager:: MatrixFluidPropertiesTest(const double& T, const double& P)
{
	auto wt = WaterPropertyTable::instance();
	auto gt = CO2PropertyTable::instance();

	try {
		const auto W = wt.getProperties(P, T);
		std::cout << "Water @ " << P << " Pa, " << T << " K: "
			<< " rho=" << W.rho << " kg/m^3"
			<< ", mu=" << W.mu << " Pa·s"
			<< ", cp=" << W.cp << " J/(kg·K)"
			<< ", k=" << W.k << " W/(m·K)\n";
	}
	catch (...) {
		std::cout << "[Error] WaterPropertyTable: OOR for P=" << P << " Pa, T=" << T << " K\n";
	}

	try {
		const auto G = gt.getProperties(P, T);
		std::cout << "CO2 @ " << P << " Pa, " << T << " K: "
			<< " rho=" << G.rho << " kg/m^3"
			<< ", mu=" << G.mu << " Pa·s"
			<< ", cp=" << G.cp << " J/(kg·K)"
			<< ", k=" << G.k << " W/(m·K)\n";


	}
	catch (...) {
		std::cout << "[Error] CO2PropertyTable: OOR for P=" << P << " Pa, T=" << T << " K\n";
	}
}



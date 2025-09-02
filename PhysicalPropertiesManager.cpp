#include "PhysicalPropertiesManager.h"
#include "FractureNetwork.h"
#include "FractureTypes.h"
#include "MeshManager.h"
#include "WaterPropertyTable.h"
#include "CO2PropertyTable.h"
#include "RockSolidProperties.h"
#include "FractureSolidProperties.h"
#include <cassert> // for assert
#include <algorithm>  // for swap
#include "Initializer.h"


//-——————————————-——分类———————————————————//

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

//—————————————————初始化/—————————————————//  
//Step1:构造物性参数场，若场已存在则不变

static inline void ensureRockFields(FieldRegistry& reg, size_t n)  //确保基岩物性参数场存在，若不存在则创建并赋默认值
{
	reg.getOrCreate<volScalarField>("phi", n, 0.15); //默认中等孔隙度
	reg.getOrCreate<volScalarField>("kxx", n, 1e-14); //默认中等渗透率 xx方向
	reg.getOrCreate<volScalarField>("kyy", n, 1e-14); //默认中等渗透率 yy方向
	reg.getOrCreate<volScalarField>("kzz", n, 1e-14); //默认中等渗透率 zz方向
	reg.getOrCreate<volScalarField>("rho_r", n, 2650.0); //基岩密度，kg/m³
	reg.getOrCreate<volScalarField>("cp_r", n, 1000.0); //基岩比热容，J/(kg·K)
	reg.getOrCreate<volScalarField>("lambda_r", n, 2.5); //基岩导热系数，W/(m·K)
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

//Step2:根据单元的区域类型和初始的 P 和 T 计算物性参数，并赋值给场
void PhysicalPropertiesManager::InitializeRockMatrixProperties(MeshManager& mgr, FieldRegistry& reg)
{
	auto& mesh = mgr.mesh();
	const size_t n = mesh.getCells().size();
	ensureRockFields(reg, n); //确保基岩物性参数场存在，若不存在则创建并赋默认值

	auto phi_r = reg.get<volScalarField>("phi");
	auto kxx_r = reg.get<volScalarField>("kxx");
	auto kyy_r = reg.get<volScalarField>("kyy");
	auto kzz_r = reg.get<volScalarField>("kzz");
	auto rho_r = reg.get<volScalarField>("rho_r");
	auto cp_r = reg.get<volScalarField>("cp_r");
	auto lambda_r = reg.get<volScalarField>("lambda_r");

	//从场中读取水相压力和温度
	auto p_w = reg.get<volScalarField>("p_w");
	auto T = reg.get<volScalarField>("T");

	for (const auto& cell : mesh.getCells())
	{
		if (cell.id < 0) continue; // 跳过 Ghost Cell
		const size_t i = mesh.getCellId2Index().at(cell.id); //获取网格单元的内部下标
		const double P = (*p_w)[i];   // 用湿润相压力作为固相 p 输入（如果需要全压，可替换为 p_g 或加常量差）
		const double Tc = (*T)[i];
		//将场中的物性参数赋值给单元
		const auto sp = rock::computeSolidProperties(cell.region, P, Tc);
		(*phi_r)[i] = sp.porosity;
		(*kxx_r)[i] = sp.permeability;
		(*kyy_r)[i] = sp.permeability;
		(*kzz_r)[i] = sp.permeability;
		(*rho_r)[i] = sp.rho_s;
		(*cp_r)[i] = sp.cp_s;
		(*lambda_r)[i] = sp.k_s;
	}
}

void PhysicalPropertiesManager::InitializeFractureElementsProperties(MeshManager& mgr, FieldRegistry& reg_fr)
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
	auto pfw = reg_fr.get<volScalarField>("pf_w");
	auto Tf = reg_fr.get<volScalarField>("Tf");
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

	//auto& fn = const_cast<FractureNetwork&>(mgr.fracture_network());
	//auto& waterProps = WaterPropertyTable::instance();
	//auto& CO2Props = CO2PropertyTable::instance();

	//for (auto& frac : fn.fractures)
	//{
	//	for (auto& elem : frac.elements)
	//	{
	//		double P = elem.p_fr;
	//		double T = elem.T_fr;
	//		elem.solidProps = fracture::computeSolidProperties(elem.type, P, T);
	//		// … fluid/gas …
	//		elem.fluidProps = waterProps.getProperties(P, T);
	//		elem.gasProps = CO2Props.getProperties(P, T);
	//	}
	//}

	//
	//auto& waterTable = WaterPropertyTable::instance();
	//auto& co2Table = CO2PropertyTable::instance(); //目前采用的SpanWagner方法 ，如果需要采用其他方法，将D:\dataBase_MatrialProperties中其他方法的,txt文件传进去

	//for (auto& cell : mesh.getCells())
	//{
	//	if (cell.id < 0) continue;  // 跳过 Ghost Cell

	//	// 取单元自身的初始压力和温度（调用初始化模块后应已赋值）
	//	double P = cell.pressure;
	//	double T = cell.temperature;

	//	// 1) 固相：按 cell.region、P、T 计算
	//	cell.SolidMaterialProps = rock::computeSolidProperties(cell.region, P, T);

	//	// 2) 水相 & CO₂ 相：表格中插值
	//	cell.WaterMaterialProps = waterTable.getProperties(P, T);  
	//	cell.CO2MaterialProps = co2Table.getProperties(P, T);	
//	}
//}



//----------------------------------更新----------------------------------------//
void PhysicalPropertiesManager::UpdateMatrixProperties(MeshManager& mgr, FieldRegistry& reg_r)
{
	InitializeRockMatrixProperties(mgr, reg_r);
	
	//auto& mesh = mgr.mesh();
	//// 1）根据当前压力P 和当前温度T 计算基岩和其内部的水相和CO2相物性参数
	//auto waterProps = WaterPropertyTable::instance();
	//auto CO2Props = CO2PropertyTable::instance();
		//for (auto& cell : mesh.getCells())
	//{
	//	if (cell.id < 0) continue; // 跳过 Ghost Cell
		//	double P = cell.pressure;
	//	double T = cell.temperature;
		//	cell.SolidMaterialProps = rock::computeSolidProperties(cell.region, P, T);
	//	// … fluid/gas …
	//	cell.WaterMaterialProps = waterProps.getProperties(P, T);
	//	cell.CO2MaterialProps = CO2Props.getProperties(P, T);
	//}
}

void PhysicalPropertiesManager::UpdateFractureSolidProperties(MeshManager& mgr, FieldRegistry& reg_fr)
{
	InitializeFractureElementsProperties(mgr, reg_fr);
	//auto waterProps = WaterPropertyTable::instance();
	//auto CO2Props = CO2PropertyTable::instance();
	//auto& fn = mgr.fracture_network();
	//for (auto& frac : fn.fractures) {
	//	for (auto& elem : frac.elements) 
	//	{
	//		double P = elem.p_fr;
	//		double T = elem.T_fr;
	//		elem.solidProps = fracture::computeSolidProperties(elem.type, P, T);
	//		// … fluid/gas …
	//		elem.fluidProps = waterProps.getProperties(P, T);
	//		elem.gasProps = CO2Props.getProperties(P, T);
	//	}
	//}
}

//----------------------------------输出&调试----------------------------------------//
void PhysicalPropertiesManager::debugPrintProperties( MeshManager& mgr) const
{
	//// 打印基岩单元
	//std::cout << "\n=== Rock Matrix Cells Properties ===\n";
	// auto& mesh = mgr.mesh();
	//for ( auto& cell : mesh.getCells()) {
	//	if (cell.id < 0) continue;  // ghost
	//	std::cout
	//		<< "Cell " << cell.id
	//		<< " | region=" << static_cast<int>(cell.region)   //0- Low, 1- Medium, 2- High
	//		<< " | P=" << cell.pressure
	//		<< " | T=" << cell.temperature
	//		<< " | poro=" << cell.SolidMaterialProps.porosity
	//		<< " | perm=" << cell.SolidMaterialProps.permeability
	//		<< " | mu_water=" << cell.WaterMaterialProps.mu
	//		<< " | rho_CO2=" << cell.CO2MaterialProps.rho
	//		<< "\n";
	//}

	//// 打印每条裂缝的每个裂缝段
	//std::cout << "\n=== Fracture Elements Properties ===\n";
	//const auto& fn = mgr.fracture_network();
	//for (const auto& frac : fn.fractures) {
	//	std::cout << "Fracture " << frac.id << ":\n";
	//	for (const auto& elem : frac.elements) {
	//		std::cout
	//			<< "  Elem " << elem.id
	//			<< " | cell=" << elem.cellID
	//			<< " | type=" << static_cast<int>(elem.type)    //  //0- Blocking, 1- Conductive
	//			<< " | P=" << elem.p_fr
	//			<< " | T=" << elem.T_fr
	//			<< " | poro=" << elem.solidProps.porosity
	//			<< " | perm=" << elem.solidProps.permeability
	//			<< " | mu_water=" << elem.fluidProps.mu
	//			<< " | rho_CO2=" << elem.gasProps.rho
	//			<< "\n";
	//	}
	//}
	//std::cout << std::endl;
}

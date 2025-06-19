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
	for (auto& cell : mesh.cells)
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
	for (auto& cell : mesh.cells)
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

//—————————————————初始化/—————————————————//   需要接入初始值，现在还存的是pressure
void PhysicalPropertiesManager::InitializeRockMatrixProperties(MeshManager& mgr)
{
	auto& mesh = mgr.mesh();
	auto& waterTable = WaterPropertyTable::instance();
	auto& co2Table = CO2PropertyTable::instance(); //目前采用的SpanWagner方法 ，如果需要采用其他方法，将D:\dataBase_MatrialProperties中其他方法的,txt文件传进去

	for (auto& cell : mesh.cells)
	{
		if (cell.id < 0) continue;  // 跳过 Ghost Cell

		// 取单元自身的初始压力和温度（调用初始化模块后应已赋值）
		double P = cell.pressure;
		double T = cell.temperature;

		// 1) 固相：按 cell.region、P、T 计算
		cell.SolidMaterialProps = rock::computeSolidProperties(cell.region, P, T);

		// 2) 水相 & CO₂ 相：表格中插值
		cell.WaterMaterialProps = waterTable.getProperties(P, T);  
		cell.CO2MaterialProps = co2Table.getProperties(P, T);	
	}
}

void PhysicalPropertiesManager::InitializeFractureElementsProperties(MeshManager& mgr)
{

	auto& fn = const_cast<FractureNetwork&>(mgr.fracture_network());
	auto& waterProps = WaterPropertyTable::instance();
	auto& CO2Props = CO2PropertyTable::instance();

	for (auto& frac : fn.fractures) 
	{
		for (auto& elem : frac.elements)
		{
			double P = elem.p_fr;
			double T = elem.T_fr;
			elem.solidProps = fracture::computeSolidProperties(elem.type, P, T);
			// … fluid/gas …
			elem.fluidProps = waterProps.getProperties(P, T);
			elem.gasProps = CO2Props.getProperties(P, T);
		}
	}
}

//----------------------------------更新----------------------------------------//
void PhysicalPropertiesManager::UpdateMatrixProperties(MeshManager& mgr) 
{
	
	
	// 1）根据当前压力P 和当前温度T 计算基岩和其内部的水相和CO2相物性参数
	auto waterProps = WaterPropertyTable::instance();
	auto CO2Props = CO2PropertyTable::instance();

	for (auto& cell : mgr.mesh().cells) 
	{
		if (cell.id < 0) continue; // 跳过 Ghost Cell

		double P = cell.pressure;
		double T = cell.temperature;

		cell.SolidMaterialProps = rock::computeSolidProperties(cell.region, P, T);
		// … fluid/gas …
		cell.WaterMaterialProps = waterProps.getProperties(P, T);
		cell.CO2MaterialProps = CO2Props.getProperties(P, T);
	}
}

void PhysicalPropertiesManager::UpdateFractureSolidProperties(MeshManager& mgr) 
{
	auto waterProps = WaterPropertyTable::instance();
	auto CO2Props = CO2PropertyTable::instance();
	auto& fn = mgr.fracture_network();
	for (auto& frac : fn.fractures) {
		for (auto& elem : frac.elements) 
		{
			double P = elem.p_fr;
			double T = elem.T_fr;
			elem.solidProps = fracture::computeSolidProperties(elem.type, P, T);
			// … fluid/gas …
			elem.fluidProps = waterProps.getProperties(P, T);
			elem.gasProps = CO2Props.getProperties(P, T);
		}
	}
}

//----------------------------------输出&调试----------------------------------------//
void PhysicalPropertiesManager::debugPrintProperties( MeshManager& mgr) const
{
	// 打印基岩单元
	std::cout << "\n=== Rock Matrix Cells Properties ===\n";
	const auto& mesh = mgr.mesh();
	for (const auto& cell : mesh.cells) {
		if (cell.id < 0) continue;  // ghost
		std::cout
			<< "Cell " << cell.id
			<< " | region=" << static_cast<int>(cell.region)   //0- Low, 1- Medium, 2- High
			<< " | P=" << cell.pressure
			<< " | T=" << cell.temperature
			<< " | poro=" << cell.SolidMaterialProps.porosity
			<< " | perm=" << cell.SolidMaterialProps.permeability
			<< " | mu_water=" << cell.WaterMaterialProps.mu
			<< " | rho_CO2=" << cell.CO2MaterialProps.rho
			<< "\n";
	}

	// 打印每条裂缝的每个裂缝段
	std::cout << "\n=== Fracture Elements Properties ===\n";
	const auto& fn = mgr.fracture_network();
	for (const auto& frac : fn.fractures) {
		std::cout << "Fracture " << frac.id << ":\n";
		for (const auto& elem : frac.elements) {
			std::cout
				<< "  Elem " << elem.id
				<< " | cell=" << elem.cellID
				<< " | type=" << static_cast<int>(elem.type)    //  //0- Blocking, 1- Conductive
				<< " | P=" << elem.p_fr
				<< " | T=" << elem.T_fr
				<< " | poro=" << elem.solidProps.porosity
				<< " | perm=" << elem.solidProps.permeability
				<< " | mu_water=" << elem.fluidProps.mu
				<< " | rho_CO2=" << elem.gasProps.rho
				<< "\n";
		}
	}
	std::cout << std::endl;
}

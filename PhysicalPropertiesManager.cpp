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
#include <iomanip>


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

template<class FieldT>
std::shared_ptr<FieldT>ensureSize(FieldRegistry& reg_fr, const std::string& name, std::size_t n, double defVal)
{
	auto f = reg_fr.get<FieldT>(name);
	if (!f)return reg_fr.getOrCreate<FieldT>(name, n, defVal);
	if (f->data.size()!= n) f->data.resize(n, defVal);
	return f;
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



//----------------------------------更新----------------------------------------//
void PhysicalPropertiesManager::UpdateMatrixProperties(MeshManager& mgr, FieldRegistry& reg_r)
{
	InitializeRockMatrixProperties(mgr, reg_r);
}

void PhysicalPropertiesManager::UpdateFractureSolidProperties(MeshManager& mgr, FieldRegistry& reg_fr)
{
	InitializeFractureElementsProperties(mgr, reg_fr);
}

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




//--------------------------流体物性参数赋值----------------------------------------//
inline void ensureMatrixFluidFields(FieldRegistry& reg, std::size_t n)
{
	reg.getOrCreate<volScalarField>("rho_w", n, 1000.0); //水密度，kg/m³
	reg.getOrCreate<volScalarField>("mu_w", n, 1e-3);    //水粘度，Pa·s
	reg.getOrCreate<volScalarField>("cp_w", n, 4182.0); //水比热容，J/(kg·K)
	reg.getOrCreate<volScalarField>("k_w", n, 0.6);   //水导热系数，W/(m·K)

	reg.getOrCreate<volScalarField>("rho_g", n, 1.98);   //CO2密度，kg/m³
	reg.getOrCreate<volScalarField>("mu_g", n, 1.48e-5); //CO2粘度，Pa·s
	reg.getOrCreate<volScalarField>("cp_g", n, 846.0);   //CO2比热容，J/(kg·K)
	reg.getOrCreate<volScalarField>("k_g", n, 0.0146); //CO2导热系数，W/(m·K)
}

inline void ensureFractureFluidFields(FieldRegistry& reg_fr, std::size_t ne)
{
	reg_fr.getOrCreate<volScalarField>("fr_rho_w", ne, 1000.0); //水密度，kg/m³
	reg_fr.getOrCreate<volScalarField>("fr_mu_w", ne, 1e-3);    //水粘度，Pa·s
	reg_fr.getOrCreate<volScalarField>("fr_cp_w", ne, 4182.0); //水比热容，J/(kg·K)
	reg_fr.getOrCreate<volScalarField>("fr_k_w", ne, 0.6);   //水导热系数，W/(m·K)
	reg_fr.getOrCreate<volScalarField>("fr_rho_g", ne, 1.98);   //CO2密度，kg/m³
	reg_fr.getOrCreate<volScalarField>("fr_mu_g", ne, 1.48e-5); //CO2粘度，Pa·s
	reg_fr.getOrCreate<volScalarField>("fr_cp_g", ne, 846.0);   //CO2比热容，J/(kg·K)
	reg_fr.getOrCreate<volScalarField>("fr_k_g", ne, 0.0146); //CO2导热系数，W/(m·K)
}

void PhysicalPropertiesManager::InitializeMatrixFluidProperties(MeshManager& mgr, FieldRegistry& reg, const VGParams& vg)
{
	auto& mesh = mgr.mesh();
	const auto& cells = mesh.getCells();
	const size_t n = cells.size();

	ensureMatrixFluidFields(reg, n); //确保基岩流体物性参数场存在，若不存在则创建并赋默认值

	//获取主变量场指针
	auto p_w = reg.get<volScalarField>("p_w");
	auto p_g = reg.get<volScalarField>("p_g");
	auto T = reg.get<volScalarField>("T");

	//获取物性场指针
	auto rho_wF = reg.get<volScalarField>("rho_w");
	auto mu_wF = reg.get<volScalarField>("mu_w");
	auto cp_wF = reg.get<volScalarField>("cp_w");
	auto k_wF = reg.get<volScalarField>("k_w");

	auto rho_gF = reg.get<volScalarField>("rho_g");
	auto mu_gF = reg.get<volScalarField>("mu_g");
	auto cp_gF = reg.get<volScalarField>("cp_g");
	auto k_gF = reg.get<volScalarField>("k_g");

	auto wt = WaterPropertyTable::instance();
	auto gt = CO2PropertyTable::instance();

	std::size_t oor = 0, ghost = 0;

	for (const auto& c : cells)
	{
		if (c.id < 0) { ++ghost; continue; } // 跳过 Ghost Cell

		const std::size_t i = mesh.getCellId2Index().at(c.id); //获取网格单元的内部下标
		double pw = (*p_w)[i], pg = (*p_g)[i], Tf = (*T)[i];
		Initializer::clampPT(pw, Tf); Initializer::clampPT(pg, Tf);

		double rho_w = 1000, mu_w = 1e-3, cp_w = 4200, k_w = 0.6;
		double rho_g = 600, mu_g = 1e-4, cp_g = 850, k_g = 0.08;

		try {
			const auto W = wt.getProperties(pw, Tf);
			rho_w = W.rho; mu_w = W.mu; cp_w = W.cp; k_w = W.k;
		}
		catch (...) { ++oor; }

		try {
			const auto G = gt.getProperties(pg, Tf);
			rho_g = G.rho; mu_g = G.mu; cp_g = G.cp; k_g = G.k;
		}
		catch (...) { ++oor; }

		(*rho_wF)[i] = rho_w; (*mu_wF)[i] = mu_w; (*cp_wF)[i] = cp_w; (*k_wF)[i] = k_w;
		(*rho_gF)[i] = rho_g; (*mu_gF)[i] = mu_g; (*cp_gF)[i] = cp_g; (*k_gF)[i] = k_g;

		if (ghost) std::cout << "[PPM] InitializeMatrixFluidProperties: skipped ghost cells = " << ghost << "\n";
		if (oor)   std::cout << "[PPM] InitializeMatrixFluidProperties: property-table OOR hits = " << oor << "\n";
	}
 
}

void PhysicalPropertiesManager::UpdateMatrixFluidProperties(MeshManager& mgr, FieldRegistry& reg, const VGParams& vg)
{
	InitializeMatrixFluidProperties(mgr, reg, vg);
}


// ======================= 裂缝：初始化/更新 流体物性场 =======================//

void  PhysicalPropertiesManager:: InitializeFractureFluidProperties(MeshManager& mgr, FieldRegistry& reg, FieldRegistry& reg_fr, const VGParams& vg)
{
	const auto& fr_Net = mgr.fracture_network(); // 取出裂缝网络
	const auto idx = buildFracElemIndex(fr_Net); // 调用裂缝段索引
	const size_t ne = idx.total; // 裂缝段总数
	if (!ne) {
		std::cout << "[PPM] No fracture elements. Skip InitializeFractureFluidProperties.\n";
		return;
	}

	// 获取裂缝主变量场（先前由 Initializer::initFracturePrimaries 写好）
	auto pf_w = ensureSize<volScalarField>(reg_fr, "pf_w", ne, 1.0e6);
	auto Sf_w = ensureSize<volScalarField>(reg_fr, "Sf_w", ne, 0.90);
	auto Tf = ensureSize<volScalarField>(reg_fr, "Tf", ne, 303.15);

	// 获取裂缝流体物性场
	ensureFractureFluidFields(reg_fr, ne);
	auto fr_rho_w = reg_fr.get<volScalarField>("fr_rho_w");
	auto fr_mu_w = reg_fr.get<volScalarField>("fr_mu_w");
	auto fr_cp_w = reg_fr.get<volScalarField>("fr_cp_w");
	auto fr_k_w = reg_fr.get<volScalarField>("fr_k_w");

	auto fr_rho_g = reg_fr.get<volScalarField>("fr_rho_g");
	auto fr_mu_g = reg_fr.get<volScalarField>("fr_mu_g");
	auto fr_cp_g = reg_fr.get<volScalarField>("fr_cp_g");
	auto fr_k_g = reg_fr.get<volScalarField>("fr_k_g");

	auto wt = WaterPropertyTable::instance();
	auto gt = CO2PropertyTable::instance();

	std::size_t oor = 0;

	for (size_t f = 0; f < fr_Net.fractures.size(); ++f)
	{
		const auto& F = fr_Net.fractures[f];
		const size_t base = idx.offset[f]; //本条裂缝的全局起点
		for (std::size_t e = 0; e < F.elements.size(); ++e)
		{
			const size_t g = base + e; //裂缝段全局索引

			double pw = (*pf_w)[g];
			double T = (*Tf)[g];

			// fracture gas pressure：用 vG 由 Sf_w → pc_f，再 pg_f = pw + pc_f
			const double pc_f = pc_vG((*Sf_w)[g], vg);
			double pg = pw + pc_f;

			Initializer::clampPT(pw, T); Initializer::clampPT(pg, T);

			double rho_w = 1000, mu_w = 1e-3, cp_w = 4200, k_w = 0.6;
			double rho_g = 600, mu_g = 1e-4, cp_g = 850, k_g = 0.08;

			try {
				const auto W = wt.getProperties(pw, T);
				rho_w = W.rho; mu_w = W.mu; cp_w = W.cp; k_w = W.k;
			}
			catch (...) { ++oor; }

			try {
				const auto G = gt.getProperties(pg, T);
				rho_g = G.rho; mu_g = G.mu; cp_g = G.cp; k_g = G.k;
			}
			catch (...) { ++oor; }
			(*fr_rho_w)[g] = rho_w; (*fr_mu_w)[g] = mu_w; (*fr_cp_w)[g] = cp_w; (*fr_k_w)[g] = k_w;
			(*fr_rho_g)[g] = rho_g; (*fr_mu_g)[g] = mu_g; (*fr_cp_g)[g] = cp_g; (*fr_k_g)[g] = k_g;
		}
		if (oor) std::cout << "[PPM] InitializeFractureFluidProperties: property-table OOR hits = " << oor << "\n";
	}
	
}

void  PhysicalPropertiesManager::UpdateFractureFluidProperties(MeshManager& mgr, FieldRegistry& reg, FieldRegistry& reg_fr, const VGParams& vg)
{
	InitializeFractureFluidProperties(mgr, reg, reg_fr, vg);
}


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
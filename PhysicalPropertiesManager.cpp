#include "PhysicalPropertiesManager.h"
#include "RockSolidProperties.h"
#include "Water_Properties.h"
#include "CO2_Properties.h"

#include "FractureSolidProperties.h "
#include "RockSolidProperties.h"

#include "WaterPropertyTable.h"
#include "CO2PropertyTable.h"
#include "CO2_Properties.h"
#include "Water_Properties.h"

#include "MeshManager.h"
#include "FractureNetwork.h"
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

//基岩区域分类
void PhysicalPropertiesManager::classifyRockRegionsByGeometry(MeshManager& mgr, const vector<RegionGeometry>& regionGeoms = {}, Cell::RegionType defaultRegion = Cell::RegionType::Medium)
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


//**********************基岩物性参数更新，调用RockSolidProperties.h**************************//
//
// brief :调用rock命名空间中的函数实现基岩物性参数注册与计算更新
//
bool PhysicalPropertiesManager::UpdateRockProperties(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field)
{
	return rock::computer_rock_properties(mgr, reg, p_field, T_field);
}

//********************基岩内部水相物性参数更新**************************************//
bool PhysicalPropertiesManager::UpdateWaterProperties(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field)
{
	return Water_Prop::compute_water_properties_const(mgr, reg, p_field, T_field);
}

//********************基岩内部二氧化碳物性参数更新**************************************//
bool PhysicalPropertiesManager::UpdateCO2Properties(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field)
{
	return CO2_Prop::compute_CO2_properties_const(mgr, reg, p_field, T_field);
}

//********************基岩内部单相有效热物性参数更新**************************************//

// brief: 单一水相
bool PhysicalPropertiesManager::ComputeEffectiveThermalProperties_Water(MeshManager& mgr, FieldRegistry& reg)
{
	return Water_Prop::computer_effectiveThermalProperties(mgr, reg);
}

// brief: 单一二氧化碳相
bool PhysicalPropertiesManager::ComputeEffectiveThermalProperties_CO2(MeshManager& mgr, FieldRegistry& reg)
{
	return CO2_Prop::computer_effectiveThermalProperties(mgr, reg);
}

//*********************基岩内部单相密度对压力导数更新**************************************//
// brief: 单一水相
bool PhysicalPropertiesManager::ComputeDrho_dp_Water(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field)
{
	return Water_Prop::compute_water_DrhoDp(mgr, reg, p_field, T_field);
}

// brief: 单一二氧化碳相
bool PhysicalPropertiesManager::ComputeDrho_dp_CO2(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field)
{
	return CO2_Prop::compute_CO2_DrhoDp(mgr, reg, p_field, T_field);
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


inline void ensureFracPrimaryFields(FieldRegistry& freg, size_t ne)
{
	freg.getOrCreate<volScalarField>("pf_w", ne, 1.0e6);
	freg.getOrCreate<volScalarField>("Sf_w", ne, 0.90);
	freg.getOrCreate<volScalarField>("Tf", ne, 303.15);
}

 


//**********************裂缝物性参数注册与赋值******************************//
//注册
static inline void ensureFracRockFields(FieldRegistry& reg_fr, size_t ne) //确保裂缝物性参数场存在，若不存在则创建并赋默认值
{
	reg_fr.getOrCreate<volScalarField>("fr_phi_r", ne, 1); //裂隙孔隙度
	reg_fr.getOrCreate<volScalarField>("fr_k_t", ne, 1e-12);    // 切向等效渗透率
	reg_fr.getOrCreate<volScalarField>("fr_k_n", ne, 1e-16);    // 法向等效渗透率
	reg_fr.getOrCreate<volScalarField>("fr_rho_r", ne, 2650.0); //裂缝密度，kg/m³
	reg_fr.getOrCreate<volScalarField>("fr_cp_r", ne, 1000.0); //裂缝比热容，J/(kg·K)
	reg_fr.getOrCreate<volScalarField>("fr_lambda_r", ne, 2.5); //裂缝导热系数，W/(m·K)
	reg_fr.getOrCreate<volScalarField>("fr_aperture", ne, 1e-3); //裂缝开度
}

//裂缝区域划分
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

//裂缝固相物性参数注册与计算更新
void PhysicalPropertiesManager::UpdateFractureRockAt(MeshManager& mgr, FieldRegistry& reg, FieldRegistry& reg_fr, const std::string& pf_field, const std::string& Tf_field)
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
	ensureFracRockFields(reg_fr, ne);

	// 取出计算裂缝物性参数需要的主变量场指针
	auto pfw = reg_fr.get<volScalarField>(pf_field);
	auto Tf = reg_fr.get<volScalarField>(Tf_field);
	// 取出裂缝物性参数场指针
	auto fr_phi = reg_fr.get<volScalarField>("fr_phi_r");
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
			(*fr_phi)[g] = sp.phi_f;
			(*fr_k_t)[g] = sp.permeability; //切向等效渗透率
			(*fr_k_n)[g] = sp.permeability * 1e-4; //法向等效渗透率，假设比切向小4个数量级
			(*fr_rho_r)[g] = sp.rho_f;
			(*fr_cp_r)[g] = sp.cp_f;
			(*fr_lam_r)[g] = sp.k_f;
			(*fr_b)[g] = sp.aperture;
		}
	}
}
//****************************************************************************//



//*********************裂缝中CO2的物性参数名称注册与赋值****************************//
//注册
static inline void ensureCO2inFractureFields(FieldRegistry& reg_fr, std::size_t ne)
{
	reg_fr.getOrCreate<volScalarField>("fr_rho_g", ne, 1.98);
	reg_fr.getOrCreate<volScalarField>("fr_mu_g", ne, 1.48e-5);
	reg_fr.getOrCreate<volScalarField>("fr_cp_g", ne, 846.0);
	reg_fr.getOrCreate<volScalarField>("fr_k_g", ne, 0.0146);
	reg_fr.getOrCreate<volScalarField>("fr_Drho_Dp_g", ne, 0.0);	//二氧化碳的密度对压力的导数，kg/(m³·Pa)
	reg_fr.getOrCreate <volScalarField>("fr_c_g", ne, 0.0);      //二氧化碳的可压缩系数，1/Pa
}


//*********************裂缝中水的物性参数名称注册与赋值****************************//
static inline void ensureWaterinFractureFields(FieldRegistry& reg_fr, std::size_t ne)
{
	reg_fr.getOrCreate<volScalarField>("fr_rho_w", ne, 1.98);
	reg_fr.getOrCreate<volScalarField>("fr_mu_w", ne, 1.48e-5);
	reg_fr.getOrCreate<volScalarField>("fr_cp_w", ne, 846.0);
	reg_fr.getOrCreate<volScalarField>("fr_k_w", ne, 0.0146);
	reg_fr.getOrCreate<volScalarField>("fr_Drho_Dp_w", ne, 0.0);	//二氧化碳的密度对压力的导数，kg/(m³·Pa)
	reg_fr.getOrCreate <volScalarField>("fr_c_w", ne, 0.0);		//水的可压缩系数，1/Pa
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

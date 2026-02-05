#pragma once
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

//基岩区域分类
void PhysicalPropertiesManager::classifyRockRegionsByGeometry(MeshManager& mgr, const vector<rock::RegionGeometry>& regionGeoms = {}, Cell::RegionType defaultRegion = Cell::RegionType::Medium)
{
	rock::classifyRockRegionsByGeometry(mgr, regionGeoms, defaultRegion);
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

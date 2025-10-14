#pragma once
#include <string> // std::string
#include <memory> // std::shared_ptr
#include <tuple> // std::tie
#include <utility> // std::forward
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "BCAdapter.h"
#include "FieldAcessForDiscre.h"
#include "Diff_TPFA_PermeabilityOperation.h"
#include "Diff_TPFA_GradientsOperation.h"
#include "Diff_TPFA_UpwindforGravityandDensity.h"
#include "Diff_TPFA_Mobility.h"
#include "Diff_TPFA_ComputerBoundaryFace.h"		
#include "Diff_TPFA_ComputerInnerFace.h"
#include "Diff_TPFA_DebugPrint.h"
#include "TemperatureBCAdapter.h"


/*****************************水单相达西流扩散项离散封装*********************************/

////Step1: 水单相达西流内部面扩散项离散
inline void DiffusionIterm_singlePhase_DarcyFlow_water_TPFA_Innerface
(
	MeshManager& mgr,
	const FieldRegistry& reg,
	FaceFieldRegistry& freg,
	const GravUpwind& gu,
	double k_default, // 各向同性渗透率
	const std::string& a_name = "a_f_Diff_p_w",
	const std::string& s_name = "s_f_Diff_p_w",
	const std::string& x_name = "p_w",
	bool enable_buoy = true,
	int gradSmoothIters = 0)
{
	DarcyWaterMobility_singlePhase mob_w(k_default);
	UpwindDensityByPotential_water rhoPol_w(true);
	Diffusion_TPFA_InnerFace_SinglePhase(mgr, reg, freg, gu, mob_w, rhoPol_w,
		a_name, s_name, x_name, /*buoy*/true, /*gradSmooth*/0,
		/*rho_in_matrix*/ false);
	debugPrintInnerFaces_TPFA(mgr, reg, freg, gu, mob_w, rhoPol_w,
		a_name.c_str(), s_name.c_str(), x_name.c_str(),
		/*enable_buoy*/true, /*gradSmoothIters*/0,
		/*max_to_print*/100);
	
}

////Step2: 水单相达西流边界面扩散项离散
inline void DiffusionIterm_singlePhase_DarcyFlow_water_TPFA_BoundaryFace(
	MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
	const GravUpwind& gu, double k_default, const PressureBCAdapter& bc,
	const std::string& a_name = "a_f_Diff_p_w",
	const std::string& s_name = "s_f_Diff_p_w",
	const std::string& x_name = "p_w",
	bool enable_buoy = true, int gradSmoothIters = 0,
	bool rho_in_matrix = false      // <<< 新增：默认不把ρ进矩阵
) {
	DarcyWaterMobility_singlePhase mob_w(k_default);
	UpwindDensityByPotential_water rhoPol_w(true);
	Diffusion_TPFA_BoundaryFace_SinglePhase(
		mgr, reg, freg, gu, mob_w, rhoPol_w, bc,
		/*a*/a_name, /*s*/s_name, /*x*/x_name,
		/*buoy*/enable_buoy, /*gradSmooth*/gradSmoothIters,
		/*rho_in_matrix*/ rho_in_matrix
	);
	debugPrintBoundaryFaces_TPFA(
		mgr, reg, freg, gu, mob_w, rhoPol_w, bc,
		a_name.c_str(), s_name.c_str(), x_name.c_str(),
		enable_buoy, gradSmoothIters, /*max_to_print*/100
	);
}

inline void DiffusionIterm_TPFA_water_singlePhase_DarcyFlow
(
	MeshManager& mgr,
	const FieldRegistry& reg,
	FaceFieldRegistry& freg,
	const GravUpwind& gu,
	double k_default, // 各向同性渗透率
	const PressureBCAdapter& bc, // BC 提供器
	const std::string& a_name = "a_f_Diff_p_w",
	const std::string& s_name = "s_f_Diff_p_w",
	const std::string& x_name = "p_w",
	bool enable_buoy = true,
	int gradSmoothIters = 0
)
{
	DiffusionIterm_singlePhase_DarcyFlow_water_TPFA_Innerface (mgr, reg, freg, gu, k_default, a_name, s_name, x_name, enable_buoy, gradSmoothIters);
	DiffusionIterm_singlePhase_DarcyFlow_water_TPFA_BoundaryFace (mgr, reg, freg, gu, k_default, bc, a_name, s_name, x_name, enable_buoy, gradSmoothIters);
}
/************************************************************************************************/

/************************************CO2单相达西流扩散项离散封装*********************************/

////Step1: CO2单相达西流内部面扩散项离散
inline void DiffusionIterm_singlePhase_DarcyFlow_CO2_TPFA_Innerface
(
	MeshManager& mgr,
	const FieldRegistry& reg,
	FaceFieldRegistry& freg,
	const GravUpwind& gu,
	double k_default, // 各向同性渗透率
	const std::string& a_name = "a_f_Diff_p_g",
	const std::string& s_name = "s_f_Diff_p_g",
	const std::string& x_name = "p_g",
	bool enable_buoy = true,
	int gradSmoothIters = 0)
{
	DarcyCO2Mobility_singlePhase mob_g (k_default);
	UpwindDensityByPotential_CO2 rhoPol_g (true);
	Diffusion_TPFA_InnerFace_SinglePhase(
		mgr, reg, freg, gu, mob_g, rhoPol_g,
		a_name, s_name, x_name,
		enable_buoy, gradSmoothIters,
		/*rho_in_matrix*/ false);
	debugPrintInnerFaces_TPFA(mgr, reg, freg, gu, mob_g, rhoPol_g,
		a_name.c_str(), s_name.c_str(), x_name.c_str(),
		enable_buoy, gradSmoothIters, 100);

}
////Step2: CO2单相达西流边界面扩散项离散
inline void DiffusionIterm_singlePhase_DarcyFlow_CO2_TPFA_BoundaryFace(
	MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
	const GravUpwind& gu, double k_default, const PressureBCAdapter& bc,
	const std::string& a_name = "a_f_Diff_p_g",
	const std::string& s_name = "s_f_Diff_p_g",
	const std::string& x_name = "p_g",
	bool enable_buoy = true, int gradSmoothIters = 0,
	bool rho_in_matrix = false
) {
	DarcyCO2Mobility_singlePhase mob_g(k_default);
	UpwindDensityByPotential_CO2 rhoPol_g(true);
	Diffusion_TPFA_BoundaryFace_SinglePhase(
		mgr, reg, freg, gu, mob_g, rhoPol_g, bc,
		a_name, s_name, x_name,
		enable_buoy, gradSmoothIters,
		rho_in_matrix
	);
	debugPrintBoundaryFaces_TPFA(
		mgr, reg, freg, gu, mob_g, rhoPol_g, bc,
		a_name.c_str(), s_name.c_str(), x_name.c_str(),
		enable_buoy, gradSmoothIters, 100
	);
}

inline void DiffusionIterm_TPFA_CO2_singlePhase_DarcyFlow
(
	MeshManager& mgr,
	const FieldRegistry& reg,
	FaceFieldRegistry& freg,
	const GravUpwind& gu,
	double k_default, // 各向同性渗透率
	const PressureBCAdapter& bc, // BC 提供器
	const std::string& a_name = "a_f_Diff_p_g",
	const std::string& s_name = "s_f_Diff_p_g",
	const std::string& x_name = "p_g",
	bool enable_buoy = true,
	int gradSmoothIters = 0
)
{
	DiffusionIterm_singlePhase_DarcyFlow_CO2_TPFA_Innerface(mgr, reg, freg, gu, k_default, a_name, s_name, x_name, enable_buoy, gradSmoothIters);
	DiffusionIterm_singlePhase_DarcyFlow_CO2_TPFA_BoundaryFace (mgr, reg, freg, gu, k_default, bc, a_name, s_name, x_name, enable_buoy, gradSmoothIters);
}
/************************************************************************************************/

/*************************************单相能量守恒方程扩散项离散封装*****************************/

// ───────────────────────────────
//  纯扩散（例如温度扩散项）
//    ρ_up=1，ρ̄=0，浮力不计（enable_buoy=false）
// ───────────────────────────────

inline void DiffusionIterm_singlePhase_Temperature_TPFA_InnerFace
(
	MeshManager& mgr,
	const FieldRegistry& reg,
	FaceFieldRegistry& freg,
	const GravUpwind& gu,
	const std::string& lambda_name, // e.g. "lambda_eff"
	const std::string& a_name = "a_f_Diff_T",  //这里设置好了方便后续调用
	const std::string& s_name = "s_f_Diff_T",
	const std::string& x_name = "T",
	int gradSmoothIters = 0
)
{
	IsotropicLambdaFromField mob(lambda_name);
	NoDensity rhoPol;
	Diffusion_TPFA_InnerFace_SinglePhase (mgr, reg, freg, gu, mob, rhoPol, a_name, s_name, x_name, false, gradSmoothIters);
	debugPrintInnerFaces_TPFA(mgr, reg, freg, gu, mob, rhoPol,
		/*a_name*/"a_f_Diff_T", /*s_name*/"s_f_Diff_T", /*x_name*/"T",
		/*enable_buoy*/false, /*gradSmoothIters*/0,
		/*max_to_print*/100);

}

inline void DiffusionIterm_singlePhase_Temperature_TPFA_BoundaryFace
(
	MeshManager& mgr,
	const FieldRegistry& reg,
	FaceFieldRegistry& freg,
	const GravUpwind& gu,
	const std::string& lambda_name, // e.g. "lambda_eff"
	const TemperatureBCAdapter& Tbc,
	const std::string& a_name = "a_f_Diff_T",  //这里设置好了方便后续调用
	const std::string& s_name = "s_f_Diff_T",
	const std::string& x_name = "T",
	int gradSmoothIters = 0
)
{
	IsotropicLambdaFromField mob(lambda_name);
	NoDensity rhoPol;
	Diffusion_TPFA_BoundaryFace_SinglePhase(mgr, reg, freg, gu, mob, rhoPol, Tbc, a_name, s_name, x_name, false, gradSmoothIters);
	debugPrintBoundaryFaces_TPFA(mgr, reg, freg, gu, mob, rhoPol, Tbc,
		"a_f_Diff_T", "s_f_Diff_T", "T",
		/*enable_buoy*/false, /*gradSmoothIters*/0,
		/*max_to_print*/100);
}

inline void DiffusionIterm_TPFA_Temperature_singlePhase
(
	MeshManager& mgr,
	const FieldRegistry& reg,
	FaceFieldRegistry& freg,
	const GravUpwind& gu,
	const std::string& lambda_name, // e.g. "lambda_eff"
	const TemperatureBCAdapter& Tbc,
	const std::string& a_name = "a_f_Diff_T",  //这里设置好了方便后续调用
	const std::string& s_name = "s_f_Diff_T",
	const std::string& x_name = "T",
	int gradSmoothIters = 0
)
{
	DiffusionIterm_singlePhase_Temperature_TPFA_InnerFace(mgr, reg, freg, gu, lambda_name, a_name, s_name, x_name, gradSmoothIters);
	DiffusionIterm_singlePhase_Temperature_TPFA_BoundaryFace(mgr, reg, freg, gu, lambda_name, Tbc, a_name, s_name, x_name, gradSmoothIters);
}
/************************************************************************************************/

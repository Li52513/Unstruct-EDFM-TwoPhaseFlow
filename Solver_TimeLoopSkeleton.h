#pragma once
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "PhysicalPropertiesManager.h"
#include "Diff_TPFA_Operators.h"
#include "Conv_FirstOrder_Operators.h"
#include "TimeIterm_Euler_SinglePhase_PressureEq.h"     // TimeTerm_Euler_SinglePhase_Flow(...)
#include "TimeIterm_Euler_SinglePhase_TemperatureEQ.h"  // TimeTerm_Euler_SinglePhase_Temperature(...)
#include "BCAdapter.h"
#include "TemperatureBCAdapter.h"
#include "Solver_AssemblerCOO.h"   // 你的装配与 SparseSystemCOO / buildUnknownMap
#include "Solver_TimeLoopUtils.h" // 你可能需要的迭代器/求解器等工具
#include "Solver_PostChecks.h"  


//================================小工具：向量-场 gather/scatter,规避ghost单元================================//
inline std::vector<int> buildUnknownMapChecked(Mesh& mesh, int& N)
{
	auto m = buildUnknownMap(mesh, N);
	if (N <= 0) { std::cerr << "[TimeLoop] empty unknown set.\n"; }
	return m;
}

inline std::vector<double> gatherFieldToVec(const FieldRegistry& reg, Mesh& mesh, const std::string& fld, const std::vector<int>& lid_of_cell, int N)
{
	std::vector<double> x(N, 0.0);
	auto f = reg.get<volScalarField>(fld);
	const auto& cells = mesh.getCells();
	const auto& id2idx = mesh.getCellId2Index();
	if (!f) return x;
	for (const auto& c : cells) {
		if (c.id < 0)continue;
		size_t i = id2idx.at(c.id);
		int r = lid_of_cell[i];
		if (r >= 0 && r < N) x[r] = (*f)[i]; // 只赋值自由度
	}
	return x;
}

inline void scatterVecToField(FieldRegistry& reg, Mesh& mesh,const std::string& fld, const std::vector<int>& lid_of_cell, const std::vector<double>& x)
{
	auto f = reg.get<volScalarField>(fld);
	const auto& cells = mesh.getCells();
	const auto& id2idx = mesh.getCellId2Index();
	if (!f) return;
	for (const auto& c : cells) 
	{
		if (c.id < 0) continue;
		size_t i = id2idx.at(c.id);
		int r = lid_of_cell[i];
		if (r >= 0) (*f)[i] = x[r];
	}
}

//==================== Jacobi 迭代（占位） ====================//
// 备注：便于最小可跑；后续你可换成 CG/AMG 等更强求解器。
struct JacobiOpts { int maxIt = 200; double omega = 0.8; double atol = 1e-8; };

inline void extractDiagonal(const SparseSystemCOO& sys, std::vector<double>& D) //提取对角线
{
	D.assign(sys.n, 0.0);
	for (const auto& t : sys.A) if (t.r == t.c) D[t.r] += t.v;
}

inline void spmv(const SparseSystemCOO& sys, const std::vector<double>& x, std::vector<double>& y) //r = 稀疏矩阵-向量乘积
{
	y.assign(sys.n, 0.0);
	for (const auto& t : sys.A) y[t.r] += t.v * x[t.c];
}

inline double infNorm(const std::vector<double>& v) 
{
	double m = 0.0; for (double a : v) m = std::max(m, std::abs(a)); return m;
}

inline void jacobiSolve
(
	const SparseSystemCOO& sys,		// in: 系统矩阵
	std::vector<double>& x,			// in: 初值; out: 解
	const JacobiOpts& opt,			// in: 迭代控制参数
	double* out_resid = nullptr,    // out: 最终残差
	int* out_iters = nullptr        // out: 实际迭代次数
)
{
	std::vector<double> D, r(sys.n), Dx(sys.n), Ax(sys.n);
	extractDiagonal(sys, D);
	if (x.size() != (size_t)sys.n) x.assign(sys.n, 0.0);

	double lastRes = 1e300;  //最后一次更新的无穷范数
	int k = 0;
	for (; k < opt.maxIt; ++k) {
		spmv(sys, x, Ax);
		for (int i = 0; i < sys.n; ++i) r[i] = sys.b[i] - Ax[i];
		// 更新
		double updInf = 0.0;
		for (int i = 0; i < sys.n; ++i) {
			double di = (std::abs(D[i]) > 0 ? r[i] / D[i] : 0.0);
			double dx = opt.omega * di;
			x[i] += dx;
			updInf = std::max(updInf, std::abs(dx));
		}
		lastRes = updInf;
		if (updInf < opt.atol)
		{
			cout << "内迭代次数" << k << endl;
			break;
		}
	}
	if (out_resid) *out_resid = lastRes;
	if (out_iters) *out_iters = k;
}


//==================== 欠松弛应用 ====================//
inline void underRelax(std::vector<double>& dest, const std::vector<double>& src, double urf)
{
	const size_t n = dest.size();
	for (size_t i = 0; i < n; ++i) dest[i] = dest[i] + urf * (src[i] - dest[i]);
}



///===============================外迭代骨架-控制参数==========================================//
struct SolverControls 
{
	int maxOuter = 50;	// 最大外迭代次数
	double tol_p_abs = 1e-6; // 压力方程绝对残差收敛
	double tol_T_abs = 1e-6; // 温度方程绝对残差收敛
	double urf_p = 0.7;	// 压力欠松弛因子 应用于压力场更新
	double urf_T = 0.7; // 温度欠松弛因子 应用于温度度场更新
	double c_phi_const = 1e-12; // 孔隙度可压缩性
	JacobiOpts jac_p; //压力线性解
	JacobiOpts jac_T; //温度线性解
	// 可选：是否每次迭代打印装配报告（默认 false）
	bool reportPerIter = true;
	// 可选：是否导出 MatrixMarket（仅在 reportPerIter 为 true 时，且只导出最后一次迭代）
	bool dumpMMOnLastIter = true;

};

// ============ 一次外迭代：CO2相 ============
/**
 * @brief CO2：做“一次外迭代”
 * @param ...（见参数名：网格/场/面场/物性/边界/重力/岩石/步长/控制）
 * @param dp_inf  输出：本次外迭代 |Δp|_inf
 * @param dT_inf  输出：本次外迭代 |ΔT|_inf
 * @param lastSysP 可选输出：本次迭代“压力方程”压缩后的系统，用于进一步导出/检查
 * @param lastSysT 可选输出：本次迭代“温度方程”压缩后的系统
 * @return true/false
 */
inline bool doOneOuter_CO2
(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const PressureBCAdapter& Pbc,
	const TemperatureBCAdapter& Tbc,
	const GravUpwind& gu,
	const RockDefaults& rock,
	double dt,
	const SolverControls& ctrl,
	// out
	double& dp_inf,
	double& dT_inf,
	// optional out
	SparseSystemCOO* lastSysP = nullptr,
	SparseSystemCOO* lastSysT = nullptr
)
{
	Mesh& mesh = mgr.mesh();
	// 1) 把上一迭代层拷到“工作场”（p_g,T）
	if (!startOuterIteration(reg, /*p*/"p_g", /*T*/"T",/*p_prev*/"p_g_prev", /*T_prev*/"T_prev")) return false;

	// 2) 物性（全隐：在工作场 p_g/T 处评估）
	ppm.UpdateMatrixRockAt(mgr, reg, "p_g", "T");
	ppm.UpdateMatrixFluidAt(mgr, reg, "p_g", "T", "CO2");
	ppm.ComputeMatrixEffectiveThermalsAt(mgr, reg, "p_g", "T", "CO2", 1e-12);

	// 3) 压力：面账本 + 时间项 + 组装 + 解
	SparseSystemCOO sysP;
	{
		// 扩散项
		DiffusionIterm_TPFA_CO2_singlePhase_DarcyFlow(mgr, reg, freg, gu, rock.k_iso, Pbc, "a_f_Diff_p_g", "s_f_Diff_p_g", "p_g", /*enable_buoy=*/true, /*gradSmoothIters=*/0);

		//时间项
		computeRhoAndDrhoDpAt(mgr, reg,/*p_eval*/"p_g", /*T_eval*/"T",/*phase*/"CO2",/*rhoOut*/"rho_time", /*drhoOut*/"drho_dp_time");
		TimeTerm_Euler_SinglePhase_Flow(mgr, reg, dt, ctrl.c_phi_const, "phi", "p_g_old", "rho_time", "drho_dp_time", "aC_time_p", "bC_time_p");

		//装配
		assemblePressure_CO2_singlePhase_COO(mgr, reg, freg, &sysP);
		sysP.compressInPlace(0.0);

		if (ctrl.reportPerIter) 
		{
			auto R = PostChecks::reportAssembly(sysP, /*force_compress*/false);
			PostChecks::printAssemblyReport(R, "P");
		}
		int N = 0; auto lid = buildUnknownMap(mesh, N);
		auto p_vec = gatherFieldToVec(reg, mesh, "p_g", lid, N);
		double resP = 0; int itP = 0;
		jacobiSolve(sysP, p_vec, ctrl.jac_p, &resP, &itP);
		scatterVecToField(reg, mesh, "p_g", lid, p_vec);
	}

	// 4) 温度：扩散+对流（基于 CO2 质量通量）+ 时间项 + 解

	SparseSystemCOO sysT;
	{
		DiffusionIterm_TPFA_Temperature_singlePhase( mgr, reg, freg, gu,/*lambda_eff*/"lambda_eff", Tbc,/*a_f*/"a_f_Diff_T", /*s_f*/"s_f_Diff_T",/*T_field*/"T");

		Convective_FirstOrder_SinglePhase_Temperature( mgr, reg, freg, Tbc,/*cp*/"cp_g",/*p*/ "p_g",/*T*/ "T",/*a_f^p*/"a_f_Diff_p_g", /*s_f^p*/"s_f_Diff_p_g",/*aPP*/"aPP_conv", /*aPN*/"aPN_conv", /*bP*/"bP_conv");

		TimeTerm_Euler_SinglePhase_Temperature(mgr, reg, dt, /*rock comp*/1e-12,/*phi*/"phi",/*rho_r*/"rho_r", /*cp_r*/"cp_r",/*rho_f*/"rho_g", /*cp_f*/"cp_g",/*T_now*/"T",/*T_old*/"T_old", /*aC*/"aC_time_T", /*bC*/"bC_time_T");

		assembleTemperature_singlePhase_COO(mgr, reg, freg,"a_f_Diff_T", "s_f_Diff_T","aPP_conv", "aPN_conv", "bP_conv","aC_time_T", "bC_time_T",&sysT);
		
		sysT.compressInPlace(0.0);

		if (ctrl.reportPerIter)
		{
			auto R = PostChecks::reportAssembly(sysT, /*force_compress*/false);
			PostChecks::printAssemblyReport(R, "T");
		}

		int N = 0; auto lid = buildUnknownMap(mesh, N);
		auto T_vec = gatherFieldToVec(reg, mesh, "T", lid, N);
		double resT = 0; int itT = 0;
		jacobiSolve(sysT, T_vec, ctrl.jac_T, &resT, &itT);
		scatterVecToField(reg, mesh, "T", lid, T_vec);
	}

	// 5) 欠松弛 + 收敛衡量（相对上一迭代层 prev）
	underRelaxInPlace(reg, "p_g", "p_g_prev", ctrl.urf_p);
	underRelaxInPlace(reg, "T", "T_prev", ctrl.urf_T);

	dp_inf = maxAbsDiff(reg, "p_g", "p_g_prev");
	dT_inf = maxAbsDiff(reg, "T", "T_prev");

	// 6) 写回 prev
	updatePrevIterates(reg, /*p*/"p_g", /*T*/"T", /*p_prev*/"p_g_prev", /*T_prev*/"T_prev");

	// 可选：把最后一次迭代的系统拷出去供上层导出
	if (lastSysP) *lastSysP = sysP;
	if (lastSysT) *lastSysT = sysT;
	return true;
}

// ============ 一次外迭代：WATER ============
inline bool doOneOuter_WATER(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const PressureBCAdapter& Pbc,
	const TemperatureBCAdapter& Tbc,
	const GravUpwind& gu,
	const RockDefaults& rock,
	double dt,
	const SolverControls& ctrl,
	// out
	double& dp_inf,
	double& dT_inf,
	// optional out
	SparseSystemCOO* lastSysP = nullptr,
	SparseSystemCOO* lastSysT = nullptr
) {
	Mesh& mesh = mgr.mesh();

	if (!startOuterIteration(reg, /*p*/"p_w", /*T*/"T",
		/*p_prev*/"p_w_prev", /*T_prev*/"T_prev"))
		return false;

	ppm.UpdateMatrixRockAt(mgr, reg, "p_w", "T");
	ppm.UpdateMatrixFluidAt(mgr, reg, "p_w", "T", "water");
	ppm.ComputeMatrixEffectiveThermalsAt(mgr, reg, "p_w", "T", "water", 1e-12);

	SparseSystemCOO sysP;
	{
		DiffusionIterm_TPFA_water_singlePhase_DarcyFlow(
			mgr, reg, freg, gu, rock.k_iso, Pbc,
			/*a_f*/"a_f_Diff_p_w", /*s_f*/"s_f_Diff_p_w",
			/*p_field*/"p_w", /*enable_buoy*/true, /*gradSmoothIters*/0);

		computeRhoAndDrhoDpAt(mgr, reg, "p_w", "T", "water", "rho_time", "drho_dp_time");

		TimeTerm_Euler_SinglePhase_Flow(
			mgr, reg, dt, ctrl.c_phi_const,
			/*phi*/"phi", /*p_old*/"p_w_old",
			/*rho_time*/"rho_time", /*drho_dp_time*/"drho_dp_time",
			/*aC*/"aC_time_p", /*bC*/"bC_time_p");

		assemblePressure_water_singlePhase_COO(mgr, reg, freg, &sysP);
		sysP.compressInPlace(0.0);

		if (ctrl.reportPerIter) {
			auto R = PostChecks::reportAssembly(sysP, false);
			PostChecks::printAssemblyReport(R, "P(w)");
		}

		int N = 0; auto lid = buildUnknownMap(mesh, N);
		auto p_vec = gatherFieldToVec(reg, mesh, "p_w", lid, N);
		double resP = 0; int itP = 0;
		jacobiSolve(sysP, p_vec, ctrl.jac_p, &resP, &itP);
		scatterVecToField(reg, mesh, "p_w", lid, p_vec);
	}

	SparseSystemCOO sysT;
	{
		DiffusionIterm_TPFA_Temperature_singlePhase(
			mgr, reg, freg, gu,
			"lambda_eff", Tbc,
			"a_f_Diff_T", "s_f_Diff_T",
			"T");

		Convective_FirstOrder_SinglePhase_Temperature(
			mgr, reg, freg, Tbc,
			"cp_w", "p_w", "T",
			"a_f_Diff_p_w", "s_f_Diff_p_w",
			"aPP_conv", "aPN_conv", "bP_conv");

		TimeTerm_Euler_SinglePhase_Temperature(
			mgr, reg, dt, 1e-12,
			"phi", "rho_r", "cp_r",
			"rho_w", "cp_w",
			"T", "T_old",
			"aC_time_T", "bC_time_T");

		assembleTemperature_singlePhase_COO(
			mgr, reg, freg,
			"a_f_Diff_T", "s_f_Diff_T",
			"aPP_conv", "aPN_conv", "bP_conv",
			"aC_time_T", "bC_time_T",
			&sysT);
		sysT.compressInPlace(0.0);

		if (ctrl.reportPerIter) {
			auto R = PostChecks::reportAssembly(sysT, false);
			PostChecks::printAssemblyReport(R, "T(w)");
		}

		int N = 0; auto lid = buildUnknownMap(mesh, N);
		auto T_vec = gatherFieldToVec(reg, mesh, "T", lid, N);
		double resT = 0; int itT = 0;
		jacobiSolve(sysT, T_vec, ctrl.jac_T, &resT, &itT);
		scatterVecToField(reg, mesh, "T", lid, T_vec);
	}

	underRelaxInPlace(reg, "p_w", "p_w_prev", ctrl.urf_p);
	underRelaxInPlace(reg, "T", "T_prev", ctrl.urf_T);

	dp_inf = maxAbsDiff(reg, "p_w", "p_w_prev");
	dT_inf = maxAbsDiff(reg, "T", "T_prev");

	updatePrevIterates(reg, "p_w", "T", "p_w_prev", "T_prev");

	if (lastSysP) *lastSysP = sysP;
	if (lastSysT) *lastSysT = sysT;
	return true;
}


// ============ 外迭代驱动（单步） ============
/**
 * @brief 单步外迭代驱动（根据 phase 选择 CO2/WATER 的一次外迭代）
 * @param phase 字符串 "CO2" 或 "water"
 * 其余参数同上。
 */

inline bool outerIter_OneStep_singlePhase(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const PressureBCAdapter& Pbc,
	const TemperatureBCAdapter& Tbc,
	const GravUpwind& gu,
	const RockDefaults& rock,
	double dt,
	const SolverControls& ctrl,
	const std::string& phase = "CO2"
) {
	// 时间步开始（冻结 *_old，初始化 *_prev）
	if (phase == "CO2" || phase == "co2") {
		if (!startTimeStep(mgr.mesh(), reg, "p_g", "T", "p_g_old", "T_old", "p_g_prev", "T_prev")) return false;
	}
	else if (phase == "water" || phase == "WATER") {
		if (!startTimeStep(mgr.mesh(), reg, "p_w", "T", "p_w_old", "T_old", "p_w_prev", "T_prev")) return false;
	}
	else {
		std::cerr << "[outer] unknown phase: " << phase << "\n";
		return false;
	}

	// 外迭代
	for (int it = 0; it < ctrl.maxOuter; ++it) 
	{
		double dp = 0, dT = 0;
		SparseSystemCOO lastP, lastT;

		bool ok = (phase == "CO2" || phase == "co2")
			? doOneOuter_CO2(mgr, reg, freg, ppm, Pbc, Tbc, gu, rock, dt, ctrl, dp, dT,
				ctrl.dumpMMOnLastIter ? &lastP : nullptr,
				ctrl.dumpMMOnLastIter ? &lastT : nullptr)
			: doOneOuter_WATER(mgr, reg, freg, ppm, Pbc, Tbc, gu, rock, dt, ctrl, dp, dT,
				ctrl.dumpMMOnLastIter ? &lastP : nullptr,
				ctrl.dumpMMOnLastIter ? &lastT : nullptr);
		if (!ok) return false;

		std::cout << "[Outer " << it << "]  |Δp|_inf=" << dp << "  |ΔT|_inf=" << dT << "\n";

		if (dp < ctrl.tol_p_abs && dT < ctrl.tol_T_abs) {
			std::cout << "Converged at outer iter " << it << "\n";

			// 如果需要，把最后一次迭代的系统导出（只导出最后一轮）
			if (ctrl.dumpMMOnLastIter) {
#if __cplusplus >= 201703L
				try { std::filesystem::create_directories("mm"); }
				catch (...) {}
#endif
				const char* tag = (phase == "CO2" || phase == "co2") ? "CO2" : "WAT";
				PostChecks::dumpCOO_to_matrix_market(lastP, std::string("mm/A_P_") + tag + ".mtx",
					std::string("mm/b_P_") + tag + ".txt", false);
				PostChecks::dumpCOO_to_matrix_market(lastT, std::string("mm/A_T_") + tag + ".mtx",
					std::string("mm/b_T_") + tag + ".txt", false);
			}
			break;
		}
		if (it == ctrl.maxOuter - 1) {
			std::cout << "Reached maxOuter without meeting tolerances.\n";
		}
	}
	return true;
}
























// 返回：是否收敛
//inline bool outerIter_OneStep_singlePhase
//(
//	MeshManager& mgr,
//	FieldRegistry& reg,
//	FaceFieldRegistry& freg,
//	PhysicalPropertiesManager& ppm,
//	const PressureBCAdapter& bcAdapter,
//	const TemperatureBCAdapter& Tbc,
//	const GravUpwind& gu,
//	const RockDefaults& rock,
//	double dt,
//	const SolverControls& ctrl,
//	const std::string& phase = "water"
//)
//{
//	Mesh& mesh = mgr.mesh();
//	
//	if (phase == "CO2")
//	{
//		// 1) 时间步开始：把当前时层的值复制给 *_old；利用*_old初始化 *_prev = *_old
//		{
//			auto p = reg.get<volScalarField>("p_g");
//			auto T = reg.get<volScalarField>("T");
//			auto p0 = reg.getOrCreate<volScalarField>("p_g_old", p ? p->data.size() : 0, 0.0); //上一时间步步压力
//			auto T0 = reg.getOrCreate<volScalarField>("T_old", T ? T->data.size() : 0, 0.0);   //上一时间步步温度
//			auto pk = reg.getOrCreate<volScalarField>("p_g_prev", p ? p->data.size() : 0, 0.0);
//			auto Tk = reg.getOrCreate<volScalarField>("T_prev", T ? T->data.size() : 0, 0.0);
//			if (!p || !T) { std::cerr << "[advanceOneStep] missing p_g/T.\n"; return false; }
//			*p0 = *p; *T0 = *T; // 记录上一时间步
//			*pk = *p0; *Tk = *T0; // 初始化迭代初值  
//		}
//
//		// 2）外迭代
//		int N = 0;
//		auto lid_of_cell = buildUnknownMapChecked(mesh, N); //自由度映射
//		std::vector<double>p_prev = gatherFieldToVec(reg, mesh, "p_g_prev", lid_of_cell, N);
//		std::vector<double>T_prev = gatherFieldToVec(reg, mesh, "T_prev", lid_of_cell, N);
//
//
//
//		for (int iter = 0; iter < ctrl.maxOuter; ++iter)
//		{
//			const auto p_prev_before = p_prev;  // 迭代前的向量快照
//			const auto T_prev_before = T_prev;  // 迭代前的向量快照
//
//			//2.1 创建迭代层和上一时层的时间中点  （看情况选取）
//			//bulidTimeLinearizationPoint(mgr, reg, "p_g_old", "T_old", "p_g_prev", "T_prev", "p_time_lin", "T_time_lin");
//
//			//2.2 物性更新（用当前迭代层）
//			ppm.UpdateMatrixRockAt(mgr, reg, "p_g_prev", "T_prev");
//			ppm.UpdateMatrixFluidAt(mgr, reg, "p_g_prev", "T_prev", "CO2");
//			ppm.ComputeMatrixEffectiveThermalsAt(mgr, reg, "p_g_prev", "T_prev", "CO2", 1e-12);
//
//			//2.3 连续性方程+单相Darcy定律 方程:面账本 + 时间项 + 装配 + 线性解 + 欠松弛
//			{
//				//扩散项
//				DiffusionIterm_TPFA_CO2_singlePhase_DarcyFlow(mgr, reg, freg, gu, rock.k_iso, bcAdapter, "a_f_Diff_p_g", "s_f_Diff_p_g", "p_g_prev", /*enable_buoy=*/true, /*gradSmoothIters=*/0);
//
//				//时间项
//				TimeTermStats stP;
//				computeRhoAndDrhoDpAt(
//					mgr, reg,
//					/*p_eval*/ "p_g_prev",
//					/*T_eval*/ "T_prev",
//					/*phase */ "CO2",
//					/*rhoOut*/ "rho_time",
//					/*drhoOut*/"drho_dp_time"
//				);
//				TimeTerm_Euler_SinglePhase_Flow(mgr, reg, dt, ctrl.c_phi_const, "phi", "p_g_old", "rho_time", "drho_dp_time", "aC_time_p", "bC_time_p");
//
//				// 装配
//				SparseSystemCOO sysP;
//				assemblePressure_CO2_singlePhase_COO(mgr, reg, freg, &sysP);
//				sysP.compressInPlace(0.0);
//
//				// 初值 = 当前迭代层
//				std::vector<double> p_sol = p_prev;
//
//				// 线性解
//				double resP = 0; int itP = 0;
//				jacobiSolve(sysP, p_sol, ctrl.jac_p, &resP, &itP);
//
//				std::cout << "[sizes] N(map)=" << N
//					<< "  sysP.n=" << sysP.n
//					<< "  p_prev=" << p_prev.size()
//					<< "  p_sol=" << p_sol.size() << "\n";
//
//				// 欠松弛并写回迭代层
//				underRelax(p_prev, p_sol, ctrl.urf_p);  //当前p_prev为k+1层返回值
//				scatterVecToField(reg, mesh, "p_g_prev", lid_of_cell, p_prev);
//			}
//
//			//2.4 能量守恒方程 ：扩散项+对流（基于上一步 pressure 质量通量）+ 时间项 + 装配 + 解 + 欠松弛
//			{
//				//扩散项
//				DiffusionIterm_TPFA_Temperature_singlePhase( mgr, reg, freg, gu, "lambda_eff", Tbc, "a_f_Diff_T", "s_f_Diff_T", "T_prev");
//				//对流项
//				Convective_FirstOrder_SinglePhase_Temperature(mgr, reg, freg, Tbc, "cp_g", "p_g_prev", "T_prev", "a_f_Diff_p_g", "s_f_Diff_p_g", "aPP_conv", "aPN_conv", "bP_conv");
//				//时间项
//				TimeTerm_Euler_SinglePhase_Temperature(mgr, reg, dt, 1e-12,"phi", "rho_r", "cp_r","rho_g", "cp_g","T_prev","T_old","aC_time_T", "bC_time_T");
//				// 装配
//				SparseSystemCOO sysT;
//				assembleTemperature_singlePhase_COO(mgr, reg, freg, "a_f_Diff_T", "s_f_Diff_T","aPP_conv", "aPN_conv", "bP_conv","aC_time_T", "bC_time_T", &sysT);
//				sysT.compressInPlace(0.0);
//				// 初值 = 当前迭代层
//				std::vector<double> T_sol = T_prev;
//				// 线性解
//				double resT = 0; int itT = 0;
//				jacobiSolve(sysT, T_sol, ctrl.jac_T, &resT, &itT);
//				// 欠松弛并写回迭代层
//				underRelax(T_prev, T_sol, ctrl.urf_T);
//				scatterVecToField(reg, mesh, "T_prev", lid_of_cell, T_prev);
//			}
//
//			// 2.5 收敛检查（∞-范数）
//			{
//				// 以“本次解与上次迭代层”的差值来衡量。这里直接用欠松弛前后的差（上一节已更新）
//				// 再做一次 gather 当前场与上次迭代前备份对比
//
//				auto p_prev_now = gatherFieldToVec(reg, mesh, "p_g_prev", lid_of_cell, N);
//				auto T_prev_now = gatherFieldToVec(reg, mesh, "T_prev", lid_of_cell, N);
//				double dp_inf = 0.0, dT_inf = 0.0;
//				for (int i = 0; i < N; ++i) dp_inf = std::max(dp_inf, std::abs(p_prev_now[i] - p_prev_before[i]));
//				for (int i = 0; i < N; ++i) dT_inf = std::max(dT_inf, std::abs(T_prev_now[i] - T_prev_before[i]));
//
//				// 更新备份
//				p_prev = std::move(p_prev_now);
//				T_prev = std::move(T_prev_now);
//				std::cout << "[Outer " << iter << "]  |Δp|_inf=" << dp_inf << "  |ΔT|_inf=" << dT_inf << "\n";
//				if (dp_inf < ctrl.tol_p_abs && dT_inf < ctrl.tol_T_abs)
//				{
//					std::cout << "Converged at outer iter " << iter << "\n";
//					break;
//				}
//				if (iter == ctrl.maxOuter - 1) {
//					std::cout << "Reached maxOuter without meeting tolerances.\n";
//				}
//
//			}
//
//		}
//
//		// 2.6 收敛/或到上限：把迭代层提交为当前层 (n+1)
//		{
//			auto p = reg.get<volScalarField>("p_g");  auto pk = reg.get<volScalarField>("p_g_prev");
//			auto T = reg.get<volScalarField>("T");    auto Tk = reg.get<volScalarField>("T_prev");
//			if (p && pk) *p = *pk;
//			if (T && Tk) *T = *Tk;
//		}
//		return true;
//
//	};
//
//
//}
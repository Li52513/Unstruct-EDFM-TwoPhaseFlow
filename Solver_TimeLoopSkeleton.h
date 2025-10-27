#pragma once
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
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
#include "LinearSolver_Eigen.h"


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
	double tol_p_rel = 1e-6;   // 新增：压力相对容差
	double tol_T_rel = 1e-6;   // 新增：温度相对容差
	double urf_p = 0.7;	// 压力欠松弛因子 应用于压力场更新
	double urf_T = 0.7; // 温度欠松弛因子 应用于温度度场更新
	double c_phi_const = 1e-12; // 孔隙度可压缩性
	double theta_p = 0.5;	// 时间项theta_p
	double theta_T = 0.5;	// 时间项theta_T
	JacobiOpts jac_p; //压力线性解
	JacobiOpts jac_T; //温度线性解
	LinearSolverOptions lin_p; // 压力求解设置
	LinearSolverOptions lin_T; // 温度求解设置
	bool useJacobi = false;    // 默认用 Krylov，必要时退回 Jacobi


	// 可选：是否每次迭代打印装配报告（默认 false）
	bool reportPerIter = false;
	// 可选：是否导出 MatrixMarket（仅在 reportPerIter 为 true 时，且只导出最后一次迭代）
	bool dumpMMOnLastIter = false;

};


// ============ 一次外迭代：CO2相 ============
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

	// 1) k层初值 -> 工作场
	if (!startOuterIteration(reg, "p_g", "T", "p_g_prev", "T_prev")) return false;

	// 2) 构建 θ-评估点（只控制一次，干净）
	const double theta_p = std::min(1.0, std::max(0.0, ctrl.theta_p));
	const double theta_T = std::min(1.0, std::max(0.0, ctrl.theta_T));
	buildEvalFields(
		reg,
		/*p_old*/ "p_g_old", /*p_iter*/ "p_g",
		/*T_old*/ "T_old",   /*T_iter*/ "T",
		/*p_eval*/ "p_eval", /*T_eval*/ "T_eval",
		theta_p, theta_T
	);

	// 3) 物性在“评估场”处评估（一次性）
	ppm.UpdateMatrixRockAt(mgr, reg, "p_eval", "T_eval");
	ppm.UpdateMatrixFluidAt(mgr, reg, "p_eval", "T_eval", "CO2");
	ppm.ComputeMatrixEffectiveThermalsAt(mgr, reg, "p_eval", "T_eval", "CO2", 1e-12);

	// —— 压力时间项需要：rho_eval 与 drho_dp_eval（评估点），以及 rho_old（旧层） ——
	// 3.1 评估点的 rho 和 drho/dp
	computeRhoAndDrhoDpAt(
		mgr, reg,
		/*p_name*/ "p_eval",
		/*T_name*/ "T_eval",
		/*phase*/  "CO2",
		/*rho_out*/     "rho_eval",
		/*drhodp_out*/  "drho_dp_eval"
	);
	// 3.2 旧层密度（若你已有 computeRhoAt 就用它；没有的话用 AndDrhoDp 退而求其次）
	// computeRhoAt(mgr, reg, "p_g_old", "T_old", "CO2", "rho_old");
	{
		// 临时：再算一次，导数丢掉即可
		computeRhoAndDrhoDpAt(
			mgr, reg,
			/*p_name*/ "p_g_old",
			/*T_name*/ "T_old",
			/*phase*/  "CO2",
			/*rho_out*/     "rho_old",
			/*drhodp_out*/  "_junk_drdp_old"
		);
	}

	// 4) 压力：扩散 + （θ-版）时间项 + 组装 + 解
	SparseSystemCOO sysP;
	{
		DiffusionIterm_TPFA_CO2_singlePhase_DarcyFlow
		(
			mgr, reg, freg, gu, rock.k_iso, Pbc,
			/*a_name*/ "a_f_Diff_p_g",
			/*s_name*/ "s_f_Diff_p_g",
			/*x_name*/ "p_g",
			/*enable_buoy*/ true,
			/*gradSmoothIters*/ 1
		);

		// θ-评估统一版时间项（与温度风格一致）
		TimeTerm_Theta_SinglePhase_Flow
		(
			mgr, reg, dt, ctrl.c_phi_const,
			/*phi_name*/        "phi",
			/*p_old_name*/      "p_g_old",
			/*rho_old_name*/    "rho_old",
			/*p_eval_name*/     "p_eval",
			/*T_eval_name*/     "T_eval",
			/*rho_eval_name*/   "rho_eval",
			/*drdp_eval_name*/  "drho_dp_eval",
			/*aC_name*/         "aC_time_p",
			/*bC_name*/         "bC_time_p"
		);

		assemblePressure_CO2_singlePhase_COO(mgr, reg, freg, &sysP);
		
		const bool needCompressP =
			(ctrl.lin_p.type == LinearSolverOptions::Type::SparseLU) ||
			(ctrl.lin_p.type == LinearSolverOptions::Type::LDLT) ||
			(lastSysP != nullptr);
		
		if (needCompressP) sysP.compressInPlace(0.0);

		if (ctrl.reportPerIter)
		{
			auto R = PostChecks::reportAssembly(sysP, false);
			PostChecks::printAssemblyReport(R, "P(CO2)");
		}

		int N = 0; auto lid = buildUnknownMap(mesh, N);
		auto p_vec = gatherFieldToVec(reg, mesh, "p_g", lid, N);
		double resP = 0; int itP = 0; bool okP = false;

		if (!ctrl.useJacobi) {
			auto opt = ctrl.lin_p; if (opt.tol <= 0.0) opt.tol = ctrl.tol_p_abs;
			okP = solveCOO_Eigen(sysP, p_vec, opt, &itP, &resP);
			if (!okP) std::cerr << "[LinearSolver] pressure failed, fallback Jacobi.\n";
		}
		if (ctrl.useJacobi || !okP) {
			jacobiSolve(sysP, p_vec, ctrl.jac_p, &resP, &itP);
		}
		scatterVecToField(reg, mesh, "p_g", lid, p_vec);
	}

	// 5) 温度：扩散 + 对流 + （θ-版）时间项 + 解
	SparseSystemCOO sysT;
	{
		DiffusionIterm_TPFA_Temperature_singlePhase(
			mgr, reg, freg, gu,
			/*lambda_eff*/ "lambda_eff",
			Tbc,
			/*a_f*/ "a_f_Diff_T",
			/*s_f*/ "s_f_Diff_T",
			/*T_field*/ "T"
		);

		Convective_FirstOrder_SinglePhase_Temperature(
			mgr, reg, freg, Tbc,
			/*cp*/ "cp_g",
			/*p*/  "p_g",
			/*T*/  "T",
			/*a_f^p*/ "a_f_Diff_p_g",
			/*s_f^p*/ "s_f_Diff_p_g",
			/*aPP*/ "aPP_conv",
			/*aPN*/ "aPN_conv",
			/*bP*/  "bP_conv"
		);

		// θ-评估统一版温度时间项
		// 此时 rho_r/cp_r/rho_g/cp_g 已由 ppm 在 (p_eval,T_eval) 处更新到字段里
		TimeTerm_Theta_SinglePhase_Temperature
		(
			mgr, reg, dt,
			/*Ceff_floor*/ 1e-12,
			/*phi_name*/   "phi",
			/*rho_r_name*/ "rho_r",
			/*cp_r_name*/  "cp_r",
			/*rho_f_name*/ "rho_g",
			/*cp_f_name*/  "cp_g",
			/*T_old_name*/ "T_old",
			/*aC_name*/    "aC_time_T",
			/*bC_name*/    "bC_time_T",
			/*Ceff_out*/   "Ceff_T"  // 可选诊断输出
		);

		assembleTemperature_singlePhase_COO(
			mgr, reg, freg,
			"a_f_Diff_T", "s_f_Diff_T",
			"aPP_conv", "aPN_conv", "bP_conv",
			"aC_time_T", "bC_time_T",
			&sysT
		);
		const bool needCompressT =
			(ctrl.lin_T.type == LinearSolverOptions::Type::SparseLU) ||
			(ctrl.lin_T.type == LinearSolverOptions::Type::LDLT) ||
			(lastSysT != nullptr);
		if (needCompressT) sysT.compressInPlace(0.0);

		if (ctrl.reportPerIter) {
			auto R = PostChecks::reportAssembly(sysT, false);
			PostChecks::printAssemblyReport(R, "T(CO2)");
		}

		int N = 0; auto lid = buildUnknownMap(mesh, N);
		auto T_vec = gatherFieldToVec(reg, mesh, "T", lid, N);
		double resT = 0; int itT = 0; bool okT = false;

		if (!ctrl.useJacobi) {
			auto opt = ctrl.lin_T; if (opt.tol <= 0.0) opt.tol = ctrl.tol_T_abs;
			okT = solveCOO_Eigen(sysT, T_vec, opt, &itT, &resT);
			if (!okT) std::cerr << "[LinearSolver] temperature failed, fallback Jacobi.\n";
		}
		if (ctrl.useJacobi || !okT) {
			jacobiSolve(sysT, T_vec, ctrl.jac_T, &resT, &itT);
		}
		scatterVecToField(reg, mesh, "T", lid, T_vec);

		// 可选：温度限幅（若仍想保留）
		//if (auto Tfield = reg.get<volScalarField>("T")) {
		//	for (double& val : Tfield->data) {
		//		if (val < 373.15)      val = 373.15;
		//		else if (val > 450.0)  val = 450.0;
		//	}
		//}
	}

	// 6) 欠松弛 + 收敛衡量
	underRelaxInPlace(reg, "p_g", "p_g_prev", ctrl.urf_p);
	underRelaxInPlace(reg, "T", "T_prev", ctrl.urf_T);

	dp_inf = maxAbsDiff(reg, "p_g", "p_g_prev");
	dT_inf = maxAbsDiff(reg, "T", "T_prev");

	updatePrevIterates(reg, "p_g", "T", "p_g_prev", "T_prev");

	if (lastSysP) *lastSysP = sysP;
	if (lastSysT) *lastSysT = sysT;
	return true;
}

// ============ 一次外迭代：WATER ============
inline bool doOneOuter_WATER
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

	// 1) k层初值 -> 工作场
	if (!startOuterIteration(reg, "p_w", "T", "p_w_prev", "T_prev")) return false;

	// 2) 构建 θ-评估点（只控制一次，干净）
	const double theta_p = std::min(1.0, std::max(0.0, ctrl.theta_p));
	const double theta_T = std::min(1.0, std::max(0.0, ctrl.theta_T));
	buildEvalFields(
		reg,
		/*p_old*/ "p_w_old", /*p_iter*/ "p_w",
		/*T_old*/ "T_old",   /*T_iter*/ "T",
		/*p_eval*/ "p_eval", /*T_eval*/ "T_eval",
		theta_p, theta_T
	);

	// 3) 物性在“评估场”处评估（一次性）
	ppm.UpdateMatrixRockAt(mgr, reg, "p_eval", "T_eval");
	ppm.UpdateMatrixFluidAt(mgr, reg, "p_eval", "T_eval", "water");
	ppm.ComputeMatrixEffectiveThermalsAt(mgr, reg, "p_eval", "T_eval", "water", 1e-12);

	// —— 压力时间项需要：rho_eval 与 drho_dp_eval（评估点），以及 rho_old（旧层） ——
	// 3.1 评估点的 rho 和 drho/dp
	computeRhoAndDrhoDpAt(
		mgr, reg,
		/*p_name*/ "p_eval",
		/*T_name*/ "T_eval",
		/*phase*/  "water",
		/*rho_out*/     "rho_eval",
		/*drhodp_out*/  "drho_dp_eval"
	);
	// 3.2 旧层密度
	{
		// 临时：再算一次，导数丢掉即可
		computeRhoAndDrhoDpAt(
			mgr, reg,
			/*p_name*/ "p_w_old",
			/*T_name*/ "T_old",
			/*phase*/  "water",
			/*rho_out*/     "rho_old",
			/*drhodp_out*/  "_junk_drdp_old"
		);
	}

	// 4) 压力：扩散 + （θ-版）时间项 + 组装 + 解
	SparseSystemCOO sysP;
	{
		DiffusionIterm_TPFA_water_singlePhase_DarcyFlow(
			mgr, reg, freg, gu, rock.k_iso, Pbc,
			/*a_name*/ "a_f_Diff_p_w",
			/*s_name*/ "s_f_Diff_p_w",
			/*x_name*/ "p_w",
			/*enable_buoy*/ true,
			/*gradSmoothIters*/ 1
		);

		// θ-评估统一版时间项（与温度风格一致）
		TimeTerm_Theta_SinglePhase_Flow(
			mgr, reg, dt, ctrl.c_phi_const,
			/*phi_name*/        "phi",
			/*p_old_name*/      "p_w_old",
			/*rho_old_name*/    "rho_old",
			/*p_eval_name*/     "p_eval",
			/*T_eval_name*/     "T_eval",
			/*rho_eval_name*/   "rho_eval",
			/*drdp_eval_name*/  "drho_dp_eval",
			/*aC_name*/         "aC_time_p",
			/*bC_name*/         "bC_time_p"
		);

		assemblePressure_water_singlePhase_COO(mgr, reg, freg, &sysP);
		const bool needCompressP =
			(ctrl.lin_p.type == LinearSolverOptions::Type::SparseLU) ||
			(ctrl.lin_p.type == LinearSolverOptions::Type::LDLT) ||
			(lastSysP != nullptr);
		if (needCompressP) sysP.compressInPlace(0.0);

		if (ctrl.reportPerIter) {
			auto R = PostChecks::reportAssembly(sysP, false);
			PostChecks::printAssemblyReport(R, "P(water)");
		}

		int N = 0; auto lid = buildUnknownMap(mesh, N);
		auto p_vec = gatherFieldToVec(reg, mesh, "p_w", lid, N);
		double resP = 0; int itP = 0; bool okP = false;

		if (!ctrl.useJacobi) {
			auto opt = ctrl.lin_p; if (opt.tol <= 0.0) opt.tol = ctrl.tol_p_abs;
			okP = solveCOO_Eigen(sysP, p_vec, opt, &itP, &resP);
			if (!okP) std::cerr << "[LinearSolver] pressure failed, fallback Jacobi.\n";
		}
		if (ctrl.useJacobi || !okP) {
			jacobiSolve(sysP, p_vec, ctrl.jac_p, &resP, &itP);
		}
		scatterVecToField(reg, mesh, "p_w", lid, p_vec);
	}

	// 5) 温度：扩散 + 对流 + （θ-版）时间项 + 解
	SparseSystemCOO sysT;
	{
		DiffusionIterm_TPFA_Temperature_singlePhase(
			mgr, reg, freg, gu,
			/*lambda_eff*/ "lambda_eff",
			Tbc,
			/*a_f*/ "a_f_Diff_T",
			/*s_f*/ "s_f_Diff_T",
			/*T_field*/ "T"
		);

		Convective_FirstOrder_SinglePhase_Temperature(
			mgr, reg, freg, Tbc,
			/*cp*/ "cp_w",
			/*p*/  "p_w",
			/*T*/  "T",
			/*a_f^p*/ "a_f_Diff_p_w",
			/*s_f^p*/ "s_f_Diff_p_w",
			/*aPP*/ "aPP_conv",
			/*aPN*/ "aPN_conv",
			/*bP*/  "bP_conv"
		);

		// θ-评估统一版温度时间项
		// 此时 rho_r/cp_r/rho_g/cp_g 已由 ppm 在 (p_eval,T_eval) 处更新到字段里
		TimeTerm_Theta_SinglePhase_Temperature
		(
			mgr, reg, dt,
			/*Ceff_floor*/ 1e-12,
			/*phi_name*/   "phi",
			/*rho_r_name*/ "rho_r",
			/*cp_r_name*/  "cp_r",
			/*rho_f_name*/ "rho_w",
			/*cp_f_name*/  "cp_w",
			/*T_old_name*/ "T_old",
			/*aC_name*/    "aC_time_T",
			/*bC_name*/    "bC_time_T",
			/*Ceff_out*/   "Ceff_T"  // 可选诊断输出
		);

		assembleTemperature_singlePhase_COO(
			mgr, reg, freg,
			"a_f_Diff_T", "s_f_Diff_T",
			"aPP_conv", "aPN_conv", "bP_conv",
			"aC_time_T", "bC_time_T",
			&sysT
		);
		const bool needCompressT =
			(ctrl.lin_T.type == LinearSolverOptions::Type::SparseLU) ||
			(ctrl.lin_T.type == LinearSolverOptions::Type::LDLT) ||
			(lastSysT != nullptr);
		if (needCompressT) sysT.compressInPlace(0.0);

		if (ctrl.reportPerIter) {
			auto R = PostChecks::reportAssembly(sysT, false);
			PostChecks::printAssemblyReport(R, "T(water)");
		}

		int N = 0; auto lid = buildUnknownMap(mesh, N);
		auto T_vec = gatherFieldToVec(reg, mesh, "T", lid, N);
		double resT = 0; int itT = 0; bool okT = false;

		if (!ctrl.useJacobi) {
			auto opt = ctrl.lin_T; if (opt.tol <= 0.0) opt.tol = ctrl.tol_T_abs;
			okT = solveCOO_Eigen(sysT, T_vec, opt, &itT, &resT);
			if (!okT) std::cerr << "[LinearSolver] temperature failed, fallback Jacobi.\n";
		}
		if (ctrl.useJacobi || !okT) {
			jacobiSolve(sysT, T_vec, ctrl.jac_T, &resT, &itT);
		}
		scatterVecToField(reg, mesh, "T", lid, T_vec);

		//// 可选：温度限幅（若仍想保留）
		//if (auto Tfield = reg.get<volScalarField>("T")) {
		//	for (double& val : Tfield->data) {
		//		if (val < 373.15)      val = 373.15;
		//		else if (val > 450.0)  val = 450.0;
		//	}
		//}
	}

	// 6) 欠松弛 + 收敛衡量
	underRelaxInPlace(reg, "p_w", "p_g_prev", ctrl.urf_p);
	underRelaxInPlace(reg, "T", "T_prev", ctrl.urf_T);

	dp_inf = maxAbsDiff(reg, "p_w", "p_g_prev");
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

		static double prev_dp = 1e300;
		static double prev_dT = 1e300;

		if (it > 0) {
			double rp = dp / std::max(prev_dp, 1e-30);
			double rT = dT / std::max(prev_dT, 1e-30);

			double up = (rp < 0.7) ? +0.05 : (rp > 0.95 ? -0.05 : 0.0);
			double uT = (rT < 0.7) ? +0.05 : (rT > 0.95 ? -0.05 : 0.0);

			auto& ctrl_mut = const_cast<SolverControls&>(ctrl);
			ctrl_mut.urf_p = std::min(0.7, std::max(0.15, ctrl_mut.urf_p + up));
			ctrl_mut.urf_T = std::min(0.7, std::max(0.15, ctrl_mut.urf_T + uT));
		}

		prev_dp = dp;
		prev_dT = dT;

		auto maxAbsField = [&](const std::string& fld) -> double 
			{
			double m = 0.0;
			auto f = reg.get<volScalarField>(fld);
			if (!f) return 1.0;
			const auto& cells = mgr.mesh().getCells();
			const auto& id2idx = mgr.mesh().getCellId2Index();
			for (const auto& c : cells) {
				if (c.id < 0) continue;
				size_t i = id2idx.at(c.id);
				m = std::max(m, std::abs((*f)[i]));
			}
			return std::max(1.0, m);
			};

		const std::string pName = (phase == "CO2" || phase == "co2") ? "p_g" : "p_w";
		double pScale = maxAbsField(pName);
		double TScale = maxAbsField("T");

		bool convP = dp < std::max(ctrl.tol_p_abs, ctrl.tol_p_rel * pScale);
		bool convT = dT < std::max(ctrl.tol_T_abs, ctrl.tol_T_rel * TScale);
		if (convP && convT) {
			std::cout << "Converged at outer iter " << it << "\n";
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



// 测试函数：2D-常物性-单相-CO2-热扩散问题
inline bool doOneOuters_test_singlePhase_CO2_T_diffusion
(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const TemperatureBCAdapter& Tbc,
	const GravUpwind& gu,
	double dt,
	const SolverControls& ctrl,
	// out
	double& dT_inf,
	// optional out
	SparseSystemCOO* lastSysP = nullptr,
	SparseSystemCOO* lastSysT = nullptr

)
{
	Mesh& mesh = mgr.mesh();
	if (!startOuterIteration_T(reg, "T", "T_prev")) return false;

	// 利用当前外迭代层的温度场计算物性参数（常物性） 
	ppm.RockProperties_test_constProperties_singlePhase_CO2_T_diffusion(mgr, reg);
	ppm.CO2Properties_test_constProperties_singlePhase_CO2_T_diffusion(mgr, reg);
	ppm.ComputeEffectiveThermalProperties_test_constProperties_singlePhase_CO2_T_diffusion(mgr, reg);

	// 温度扩散方程离散与组装
	// 离散
	{
		// 扩散项
		DiffusionIterm_TPFA_Temperature_singlePhase
		(
			mgr, reg, freg, gu,
			/*lambda_eff*/ "lambda_eff",
			Tbc,
			/*a_f*/ "a_f_Diff_T",
			/*s_f*/ "s_f_Diff_T",
			/*T_field*/ "T"
		);

		//时间项
		TimeTerm_Theta_SinglePhase_Temperature
		(
			mgr, reg, dt,
			/*Ceff_floor*/ 1e-12,
			/*phi_name*/   "phi",
			/*rho_r_name*/ "rho_r",
			/*cp_r_name*/  "cp_r",
			/*rho_f_name*/ "rho_g",
			/*cp_f_name*/  "cp_g",
			/*T_old_name*/ "T_old",
			/*aC_name*/    "aC_time_T",
			/*bC_name*/    "bC_time_T",
			/*Ceff_out*/   "Ceff_T"  // 可选诊断输出
		);
	}
	//组装
	SparseSystemCOO sysT_test;
	{
		const OperatorFieldNames nmT = makeNames("T");  //	取出温度场相关的算子场名
		assemble_COO(mgr, reg, freg, "ddt+laplacian", nmT, &sysT_test);  //
		if (lastSysT) *lastSysT = sysT_test; //

		//打印装配报告
		if (ctrl.reportPerIter)
		{
			auto R = PostChecks::reportAssembly(sysT_test, false);
			PostChecks::printAssemblyReport(R, "T(diffusion, implicit)");
		}
	}

	//求解
	{
		int N = 0; auto lid = buildUnknownMap(mesh, N); //生成待求解变量向量

		auto T_vec = gatherFieldToVec(reg, mesh, "T", lid, N); //从场中提取待求解变量向量

		//求解线性系统

		double resT = 0.0; int itT = 0; bool okT = false;
		auto opt = ctrl.lin_T; //传入线性求解器设置，位于  LinearSolverOptions 结构体内
		if (opt.tol <= 0.0) opt.tol = ctrl.tol_T_abs;

		//调用Eigen求解器,实现方程组的求解
		okT = solveCOO_Eigen(sysT_test, T_vec, opt, &itT, &resT); //输入参数为系数矩阵，待求解向量，线性求解器设置，输出迭代次数和最终残差
		if (!okT) {
			std::cerr << "[LinearSolver] temperature solve failed.\n";
			return false;
		}

		//散射回场
		scatterVecToField(reg, mesh, "T", lid, T_vec);

	}

	//外迭代松弛
	underRelaxInPlace(reg, "T", "T_prev", ctrl.urf_T);
	//收敛性检查
	dT_inf = maxAbsDiff(reg, "T", "T_prev");
	updatePrevIterates(reg, { {"T","T_prev"} });

	return true;

}
// 仅温度的外迭代驱动器：2D / 单相 CO2 / 非稳态 / 纯隐式热扩散
inline bool outerIter_test_singlePhase_CO2_T_diffusion(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const TemperatureBCAdapter& Tbc,
	const GravUpwind& gu,
	double dt,
	const SolverControls& ctrl
)
{

	// 外迭代
	double prev_dT = 1e300;

	for (int it = 0; it < ctrl.maxOuter; ++it) {

		double dT = 0.0;
		SparseSystemCOO lastT; //储存上一时层的线性系统矩阵
		bool ok = doOneOuters_test_singlePhase_CO2_T_diffusion
		(
			mgr, reg, freg, ppm, Tbc, gu, dt, ctrl, dT, (ctrl.dumpMMOnLastIter ? &lastT : nullptr)
		);

		if (!ok) return false;

		// —— 进度输出 —— 
		std::cout << "[Outer " << it << "]  |ΔT|_inf=" << dT << "\n";

		// —— 松弛因子自适应调整 ——
		if (it > 0) 
		{
			double rT = dT / std::max(prev_dT, 1e-30);
			double uT = (rT < 0.7) ? +0.05 : (rT > 0.95 ? -0.05 : 0.0);
			auto& ctrl_mut = const_cast<SolverControls&>(ctrl);
			ctrl_mut.urf_T = std::min(0.7, std::max(0.15, ctrl_mut.urf_T + uT));
		}
		prev_dT = dT;

		auto maxAbsField = [&](const std::string& fld) -> double 
			{
			double m = 0.0;
			auto f = reg.get<volScalarField>(fld);
			if (!f) return 1.0;
			const auto& cells = mgr.mesh().getCells();
			const auto& id2idx = mgr.mesh().getCellId2Index();
			for (const auto& c : cells) {
				if (c.id < 0) continue;
				size_t i = id2idx.at(c.id);
				m = std::max(m, std::abs((*f)[i]));
			}
			return std::max(1.0, m);
			};

		double TScale = maxAbsField("T");
		bool convT = dT < std::max(ctrl.tol_T_abs, ctrl.tol_T_rel * TScale);

		if (convT) {
			std::cout << "Converged at outer iter " << it << "\n";

			if (ctrl.dumpMMOnLastIter) {
#if __cplusplus >= 201703L
				try { std::filesystem::create_directories("mm"); }
				catch (...) {}
#endif
				PostChecks::dumpCOO_to_matrix_market(
					lastT,
					"mm/A_T_CO2_diff.mtx",
					"mm/b_T_CO2_diff.txt",
					/*sym=*/false
				);
			}
			break;
		}

		if (it == ctrl.maxOuter - 1) {
			std::cout << "Reached maxOuter without meeting T tolerances.\n";
		}


	}

	return true;

}
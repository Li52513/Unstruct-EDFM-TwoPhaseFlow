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
#include"DiffusionCentral.h"
#include"ConvectionUpwind_Flux.h"
#include"ConvectionUpwind.h"
#include"Timeterm_BDF.h"
#include "FVM_SourceTerm_Sources.h"
#include "FVM_SourceTerm_StrongDirichlet.h"
#include"FVM_SourceTerm_ProdWellOps.h"
#include "FVM_WellCoupling.h"
#include "FVM_WellDOF.h"
#include "FVM_SourceTerm_WellHeat.h"
#include "WellConfig.h"
#include "InflowPatch.h"


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

	struct InjectorPin {
		bool enable = false;
		double relativeWeight = 0.0;
		double minWeight = 0.0;
	} injectorPin;

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
inline bool doOneOuters_test_constProperties_singlePhase_CO2_T_diffusion
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
inline bool outerIter_test_constProperties_singlePhase_CO2_T_diffusion
(
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
		bool ok = doOneOuters_test_constProperties_singlePhase_CO2_T_diffusion
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


//测试函数：2D-变物性-单相-CO2-导热问题

inline bool doOneOuters_test_varyProperties_singlePhase_CO2_T_diffusion(
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

	ppm.RockProperties_test_constProperties_singlePhase_CO2_T_diffusion(mgr, reg);
	ppm.CO2Properties_test_varProperties_singlePhase_CO2_T_diffusion(mgr, reg, "T");
	ppm.ComputeEffectiveThermalProperties_test_varProperties_singlePhase_CO2_T_diffusion(mgr, reg, "T");

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
//输入参数 ：网格信息 场信息
inline bool outerIter_test_varyProperties_singlePhase_CO2_T_diffusion
(
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
		bool ok = doOneOuters_test_varyProperties_singlePhase_CO2_T_diffusion
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


//测试函数： 2D-常物性-单相-CO2-达西渗流问题
inline bool doOneOuters_test_constProperties_singlePhase_CO2_p_flow
(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const PressureBCAdapter& Pbc,
	const Vector& g,
	double dt,
	const SolverControls& ctrl,
	// out
	double& dp_inf,
	// optional out
	SparseSystemCOO* lastSysP = nullptr
)
{
	Mesh& mesh = mgr.mesh();

	// 1) 外迭代起步：保存上一迭代值 p_g_prev（同一时间层内）
	if (!startOuterIteration_p(reg, "p_g", "p_g_prev")) return false;

	// 2) 常物性装填（若已装过也没关系）
	ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
	ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);

	using FVM::Diffusion::RhoFaceMethod;

	// 3) 面系数（质量形）：a_f,s_f 直接包含“密度1”
	//    关键：mobility_tokens 里带 "rho:rho_g" → MASS_FORM 打开
	const OperatorFieldNames nmP = makeNames("p_g");
	bool ok = FVM::Diffusion::build_FaceCoeffs_Central(
		/*mgr,reg,freg*/ mgr, reg, freg,
		/*a_name*/ nmP.a_f_diff,             // "a_f_Diff_p_g"
		/*s_name*/ nmP.s_f_diff,             // "s_f_Diff_p_g"
		/*x_name*/ "p_g",
		/*mobility_tokens*/
		std::vector<std::string>{ "kxx:kxx", "kyy:kyy", "kzz:kzz", "/mu_g", "rho:rho_g" },
		/*rho_field_for_buoy*/ "rho_g",      // 浮力线性化用ρ（与 MASS_FORM 一致即可）
		/*rho_method*/ RhoFaceMethod::Linear,
		/*g*/ g,
		/*bc*/ Pbc,
		/*enable_buoy*/ true,
		/*gradSmoothIters*/ 0
	);
	if (!ok) return false;

	// 4) 时间项：全隐（用上一时间层 p_g_old 与当前线性化点 p_g_prev）
	ok = TimeTerm_FullyImplicit_SinglePhase_Flow(
		mgr, reg, dt,
		/*c_phi_name   */ "c_phi",
		/*phi_name     */ "phi",
		/*p_old_name   */ "p_g_old",     // 上一时间层 n
		/*rho_old_name */ "rho_g",
		/*p_lin_name   */ "p_g_prev",    // 本层迭代线性化点
		/*rho_lin_name */ "rho_g",
		/*drdp_lin_name*/ "Drho_Dp_g",
		/*aC_name      */ nmP.a_time,    // "aC_time_p_g"
		/*bC_name      */ nmP.b_time     // "bC_time_p_g"
	);
	if (!ok) return false;

	// 5) 组装并求解： TIME + DIFFUSION
	SparseSystemCOO sysP;
	{
		// 你已有两套装配入口，这里沿用“字符串表达式”版本
		assemble_COO(mgr, reg, freg, "ddt+diffusion", nmP, &sysP);
		if (lastSysP) *lastSysP = sysP;

		if (ctrl.reportPerIter) {
			auto R = PostChecks::reportAssembly(sysP, false);
			PostChecks::printAssemblyReport(R, "P(Time-implicit + Diffusion (mass-form Darcy))");
		}
	}

	// 6) 线性求解
	int N = 0; auto lid = buildUnknownMap(mesh, N);
	auto pvec = gatherFieldToVec(reg, mesh, "p_g", lid, N);

	auto opt = ctrl.lin_p; if (opt.tol <= 0.0) opt.tol = ctrl.tol_p_abs;
	double res = 0.0; int itL = 0;
	bool okLin = solveCOO_Eigen(sysP, pvec, opt, &itL, &res);
	if (!okLin) { std::cerr << "[LinearSolver] pressure solve failed.\n"; return false; }
	scatterVecToField(reg, mesh, "p_g", lid, pvec);

	// 7) 松弛与收敛度量（同一时间层内）
	underRelaxInPlace(reg, "p_g", "p_g_prev", ctrl.urf_p);
	dp_inf = maxAbsDiff(reg, "p_g", "p_g_prev");
	updatePrevIterates(reg, { {"p_g","p_g_prev"} });

	return true;

}

// 仅压力的外迭代驱动器：2D / 单相 CO2 / 非稳态 / 达西渗流
//输入参数 ：网格信息 场信息
// ================== 外迭代驱动：单相·常物性·CO2·达西渗流（压力） ==================
inline bool outerIter_test_constProperties_singlePhase_CO2_p_flow
(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const PressureBCAdapter& Pbc,
	const Vector& g,
	double dt,
	const SolverControls& ctrl
)
{
	double prev_dp = 1e300;

	for (int it = 0; it < ctrl.maxOuter; ++it) {

		double dp = 0.0;
		SparseSystemCOO lastP;

		// —— 单步外迭代：装配 TIME+CONVECTION 并求解 —— //
		bool ok = doOneOuters_test_constProperties_singlePhase_CO2_p_flow(
			mgr, reg, freg, ppm, Pbc, g, dt, ctrl,
			/*out*/ dp,
			/*lastSysP*/ (ctrl.dumpMMOnLastIter ? &lastP : nullptr)
		);
		if (!ok) return false;

		// —— 进度输出 —— //
		std::cout << "[Outer " << it << "]  |Δp|_inf=" << dp << "\n";

		// —— 自适应松弛（与温度版同策略）—— //
		if (it > 0) {
			double rp = dp / std::max(prev_dp, 1e-30);
			double up = (rp < 0.7) ? +0.05 : (rp > 0.95 ? -0.05 : 0.0);
			auto& ctrl_mut = const_cast<SolverControls&>(ctrl);
			ctrl_mut.urf_p = std::min(0.7, std::max(0.15, ctrl_mut.urf_p + up));
		}
		prev_dp = dp;

		// —— 归一化尺度：用当前 p_g 的最大幅值 —— //
		auto maxAbsField = [&](const std::string& fld)->double {
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

		double Pscale = maxAbsField("p_g");
		bool convP = dp < std::max(ctrl.tol_p_abs, ctrl.tol_p_rel * Pscale);

		if (convP) {
			std::cout << "Converged at outer iter " << it << "\n";
			if (ctrl.dumpMMOnLastIter) {
#if __cplusplus >= 201703L
				try { std::filesystem::create_directories("mm"); }
				catch (...) {}
#endif
				PostChecks::dumpCOO_to_matrix_market(
					lastP,
					"mm/A_P_CO2_flow.mtx",
					"mm/b_P_CO2_flow.txt",
					/*sym=*/false
				);
			}
			break;
		}

		if (it == ctrl.maxOuter - 1) {
			std::cout << "Reached maxOuter without meeting P tolerances.\n";
		}
	}
	return true;
}





//测试函数 导热问题
//测试函数： 2D-常物性-单相-CO2-达西渗流问题
inline bool doOneOuters_test_constProperties_singlePhase_CO2_T_newT
(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const TemperatureBCAdapter& Tbc,
	const Vector& g,
	double dt,
	const SolverControls& ctrl,
	// out
	double& dT_inf,
	// optional out
	SparseSystemCOO* lastSysP = nullptr
)
{
	Mesh& mesh = mgr.mesh();

	// 1) 外迭代起步：保存上一迭代值 p_g_prev（同一时间层内）
	if (!startOuterIteration_T(reg, "T", "T_prev")) return false;

	// 2) 常物性装填（若已装过也没关系）
	ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
	ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);
	ppm.ComputeEffectiveThermalProperties_test_constProperties_singlePhase_CO2_T_diffusion(mgr, reg);

	using FVM::Diffusion::RhoFaceMethod;

	// 3) 面系数（质量形）：a_f,s_f 直接包含“密度1”
	//    关键：mobility_tokens 里带 "rho:rho_g" → MASS_FORM 打开
	const OperatorFieldNames nmT = makeNames("T");
	bool ok = FVM::Diffusion::build_FaceCoeffs_Central(
		/*mgr,reg,freg*/ mgr, reg, freg,
		/*a_name*/ nmT.a_f_diff,             // "a_f_Diff_p_g"
		/*s_name*/ nmT.s_f_diff,             // "s_f_Diff_p_g"
		/*x_name*/ "T",
		/*mobility_tokens*/
		std::vector<std::string>{ "iso:lambda_eff" },
		/*rho_field_for_buoy*/ "",      // 浮力线性化用ρ（与 MASS_FORM 一致即可）
		/*rho_method*/ RhoFaceMethod::Linear,
		/*g*/ Vector{0,0,0},
		/*bc*/ Tbc,
		/*enable_buoy*/ false,
		/*gradSmoothIters*/ 0
	);
	if (!ok) return false;

	// 4) 时间项：全隐（用上一时间层 p_g_old 与当前线性化点 p_g_prev）
	ok = TimeTerm_Theta_SinglePhase_Temperature
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
	if (!ok) return false;

	// 5) 组装并求解： TIME + DIFFUSION
	SparseSystemCOO sysT;
	{
		// 你已有两套装配入口，这里沿用“字符串表达式”版本
		assemble_COO(mgr, reg, freg, "ddt+diffusion", nmT, &sysT);
		if (lastSysP) *lastSysP = sysT;

		if (ctrl.reportPerIter) {
			auto R = PostChecks::reportAssembly(sysT, false);
			PostChecks::printAssemblyReport(R, "P(Time-implicit + Diffusion (mass-form Darcy))");
		}
	}

	// 6) 线性求解
	int N = 0; auto lid = buildUnknownMap(mesh, N);
	auto pvec = gatherFieldToVec(reg, mesh, "T", lid, N);

	auto opt = ctrl.lin_T; if (opt.tol <= 0.0) opt.tol = ctrl.tol_T_abs;
	double res = 0.0; int itL = 0;
	bool okLin = solveCOO_Eigen(sysT, pvec, opt, &itL, &res);
	if (!okLin) { std::cerr << "[LinearSolver] pressure solve failed.\n"; return false; }
	scatterVecToField(reg, mesh, "T", lid, pvec);

	// 7) 松弛与收敛度量（同一时间层内）
	underRelaxInPlace(reg, "T", "T_prev", ctrl.urf_T);
	dT_inf = maxAbsDiff(reg, "T", "T_prev");
	updatePrevIterates(reg, { {"T","T_prev"} });

	return true;

}

//












// 仅压力的外迭代驱动器：2D / 单相 CO2 / 非稳态 / 达西渗流
//输入参数 ：网格信息 场信息
// ================== 外迭代驱动：单相·常物性·CO2·达西渗流（压力） ==================
inline bool outerIter_test_constProperties_singlePhase_CO2_T_newT
(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const TemperatureBCAdapter& Tbc,
	const Vector& g,
	double dt,
	const SolverControls& ctrl
)
{
	double prev_dT = 1e300;

	for (int it = 0; it < ctrl.maxOuter; ++it) {

		double dT = 0.0;
		SparseSystemCOO lastT;

		// —— 单步外迭代：装配 TIME+CONVECTION 并求解 —— //
		bool ok = doOneOuters_test_constProperties_singlePhase_CO2_T_newT(
			mgr, reg, freg, ppm, Tbc, {0,0,0}, dt, ctrl,
			/*out*/ dT,
			/*lastSysP*/ (ctrl.dumpMMOnLastIter ? &lastT : nullptr)
		);
		if (!ok) return false;

		// —— 进度输出 —— //
		std::cout << "[Outer " << it << "]  |ΔT|_inf=" << dT << "\n";

		// —— 自适应松弛（与温度版同策略）—— //
		if (it > 0) {
			double rT = dT / std::max(prev_dT, 1e-30);
			double uT = (rT < 0.7) ? +0.05 : (rT > 0.95 ? -0.05 : 0.0);
			auto& ctrl_mut = const_cast<SolverControls&>(ctrl);
			ctrl_mut.urf_T = std::min(0.7, std::max(0.15, ctrl_mut.urf_T + uT));
		}
		prev_dT = dT;

		// —— 归一化尺度：用当前 p_g 的最大幅值 —— //
		auto maxAbsField = [&](const std::string& fld)->double {
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

		double Tscale = maxAbsField("T");
		bool convT = dT < std::max(ctrl.tol_T_abs, ctrl.tol_T_rel * Tscale);

		if (convT) {
			std::cout << "Converged at outer iter " << it << "\n";
			if (ctrl.dumpMMOnLastIter) {
#if __cplusplus >= 201703L
				try { std::filesystem::create_directories("mm"); }
				catch (...) {}
#endif
				PostChecks::dumpCOO_to_matrix_market(
					lastT,
					"mm/A_T_CO2_flow.mtx",
					"mm/b_T_CO2_flow.txt",
					/*sym=*/false
				);
			}
			break;
		}

		if (it == ctrl.maxOuter - 1) {
			std::cout << "Reached maxOuter without meeting P tolerances.\n";
		}
	}
	return true;
}



inline bool doOneOutert_constProperties_singlePhase_CO2_pressure
(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const PressureBCAdapter& Pbc,
	const Vector& g,
	double dt,
	const SolverControls& ctrl,
	// out
	double& dp_inf,
	// optional out
	SparseSystemCOO* lastSysP = nullptr
)
{
	Mesh& mesh = mgr.mesh();

	if (!startOuterIteration_scatter(reg, "p_g", "p_g_prev")) return false;

	ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
	ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);

	using FVM::Diffusion::RhoFaceMethod;
	const OperatorFieldNames nm = makeNames("p_g");

	bool ok = FVM::Diffusion::build_FaceCoeffs_Central(mgr, reg, freg, nm.a_f_diff, nm.s_f_diff, "p_g", {"kxx:kxx","kyy:kyy","kzz:kzz","/mu_g","rho:rho_g"}, "rho_g", RhoFaceMethod::Linear, g, Pbc, true, 0);
	if (!ok) return false;
	ok = TimeTerm_FullyImplicit_SinglePhase_Flow(mgr, reg, dt, "c_phi", "phi", "p_g_old", "rho_g", "p_g_prev", "rho_g", "Drho_Dp_g", nm.a_time, nm.b_time);
	if (!ok) return false;

	SparseSystemCOO sysp;
	{
		assemble_COO(mgr, reg, freg, "ddt+diffusion", nm, &sysp);
		if (lastSysP) *lastSysP = sysp;
		if (ctrl.reportPerIter) {
			auto R = PostChecks::reportAssembly(sysp, false);
			PostChecks::printAssemblyReport(R, "P(Time-implicit + Diffusion (mass-form Darcy))");
		}
	}
	int N = 0; auto lid = buildUnknownMap(mesh, N);
	auto pvec = gatherFieldToVec(reg, mesh, "p_g", lid, N);
	auto opt = ctrl.lin_p; if (opt.tol <= 0.0) opt.tol = ctrl.tol_p_abs;
	double res = 0.0; int itL = 0;
	bool okLin = solveCOO_Eigen(sysp, pvec, opt, &itL, &res);
	if (!okLin) { std::cerr << "[LinearSolver] pressure solve failed.\n"; return false; }
	scatterVecToField(reg, mesh, "p_g", lid, pvec);

	underRelaxInPlace(reg, "p_g", "p_g_prev", ctrl.urf_p);
	dp_inf = maxAbsDiff(reg, "p_g", "p_g_prev");
	updatePrevIterates(reg, { {"p_g","p_g_prev"} });
	return true;



}

inline bool outerIter_constProperties_singlePhase_CO2_Pressure
(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const PressureBCAdapter& Pbc,
	const Vector& g,
	double dt,
	const SolverControls& ctrl
)
{
	double prev_dp_g = 1e300;
	for (int it = 0; it < ctrl.maxOuter; ++it) 
	{
		double dp  = 0.0;
		SparseSystemCOO lastT;

		//一次外迭代推进
		bool ok = doOneOutert_constProperties_singlePhase_CO2_pressure(mgr, reg, freg, ppm, Pbc, g, dt, ctrl, dp, (ctrl.dumpMMOnLastIter ? &lastT : nullptr));
		if (!ok) return false;

		// —— 进度输出 —— //
		std::cout << "[Outer " << it << "]  |Δp|_inf=" << dp << "\n";

		// —— 自适应松弛—— //
		if (it > 0)
		{
			double rT = dp / std::max(prev_dp_g, 1e-30);
			double uT = (rT < 0.7) ? +0.05 : (rT > 0.95 ? -0.05 : 0.0);
			auto& ctrl_mut = const_cast<SolverControls&>(ctrl);
			ctrl_mut.urf_p = std::min(0.7, std::max(0.15, ctrl_mut.urf_p + uT));
		}
		prev_dp_g = dp;

		auto maxAbsField = [&](const std::string& fld)->double {
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
		double PScale = maxAbsField("p_g");

		bool convP = dp < std::max(ctrl.tol_p_abs, ctrl.tol_p_rel * PScale);

		if (convP) {
			std::cout << "Converged at outer iter " << it << "\n";
			if (ctrl.dumpMMOnLastIter) {
#if __cplusplus >= 201703L
				try { std::filesystem::create_directories("mm"); }
				catch (...) {}
#endif
				PostChecks::dumpCOO_to_matrix_market(
					lastT,
					"mm/A_P_CO2_flow.mtx",
					"mm/b_P_CO2_flow.txt",
					/*sym=*/false
				);
			}
			break;
		}

		if (it == ctrl.maxOuter - 1) {
			std::cout << "Reached maxOuter without meeting P tolerances.\n";
		}
	}
	return true;
}


inline bool doOneOutert_constProperties_singlePhase_CO2_T_H
(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const TemperatureBCAdapter& Tbc,
	const PressureBCAdapter& Pbc,
	const Vector& g,
	double dt,
	const SolverControls& ctrl,
	// out
	double& dT_inf,
	double& dp_inf,
	// optional out
	SparseSystemCOO* lastSysP = nullptr,
	SparseSystemCOO* lastSysT = nullptr
)
{
	Mesh& mesh = mgr.mesh();
	if (!startOuterIteration_scatter(reg, "p_g", "p_g_prev")) return false;
	if (!startOuterIteration_T(reg, "T", "T_prev")) return false;

	ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
	ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);
	ppm.ComputeEffectiveThermalProperties_constProperties_singlePhase_CO2_T_H(mgr, reg);

	using FVM::Diffusion::RhoFaceMethod;
	//压力方程面系数
	const OperatorFieldNames nmP = makeNames("p_g");
	bool ok = FVM::Diffusion::build_FaceCoeffs_Central(mgr, reg, freg, nmP.a_f_diff, nmP.s_f_diff, "p_g", { "kxx:kxx","kyy:kyy","kzz:kzz","/mu_g","rho:rho_g" }, "rho_g", RhoFaceMethod::Linear, g, Pbc, true, 0);
	if (!ok) return false;
	ok = TimeTerm_FullyImplicit_SinglePhase_Flow(mgr, reg, dt, "c_phi", "phi", "p_g_old", "rho_g", "p_g_prev", "rho_g", "Drho_Dp_g", nmP.a_time, nmP.b_time,
		/*strongBCmask=*/nullptr);
	if (!ok) return false;

	//压力方程组装求解
	SparseSystemCOO sysp;
	{
		assemble_COO(mgr, reg, freg, "ddt+diffusion", nmP, &sysp);
		if (lastSysP) *lastSysP = sysp;
		if (ctrl.reportPerIter) {
			auto R = PostChecks::reportAssembly(sysp, false);
			PostChecks::printAssemblyReport(R, "P(Time-implicit + Diffusion (mass-form Darcy))");
		}
	}
	int N = 0; auto lid = buildUnknownMap(mesh, N);
	auto pvec = gatherFieldToVec(reg, mesh, "p_g", lid, N);
	auto opt = ctrl.lin_p; if (opt.tol <= 0.0) opt.tol = ctrl.tol_p_abs;
	double res = 0.0; int itL = 0;
	bool okLin = solveCOO_Eigen(sysp, pvec, opt, &itL, &res);
	if (!okLin) { std::cerr << "[LinearSolver] pressure solve failed.\n"; return false; }
	scatterVecToField(reg, mesh, "p_g", lid, pvec);

	underRelaxInPlace(reg, "p_g", "p_g_prev", ctrl.urf_p);
	dp_inf = maxAbsDiff(reg, "p_g", "p_g_prev");
	updatePrevIterates(reg, { {"p_g","p_g_prev"} });

	//温度方程面系数
	std::vector<char>   maskT = mark_strong_BC_cells(mgr, Tbc);
	std::vector<double> Ttar; build_dirichlet_T_targets(mgr, Tbc, maskT, Ttar);

	const double PIN_W = 1e7;

	const OperatorFieldNames nmT = makeNames("T");
	ok = TimeTerm_FullyImplicit_SinglePhase_Temperature(mgr, reg, dt, "C_eff", "T_old", nmT.a_time, nmT.b_time, &maskT, &Ttar, PIN_W);
	if (!ok) return false;
	ok = FVM::Diffusion::build_FaceCoeffs_Central(mgr, reg, freg, nmT.a_f_diff, nmT.s_f_diff, "T", { "iso:lambda_eff" }, "", RhoFaceMethod::Linear, { 0,0,0 }, Tbc, false, 0);
	if (!ok) return false;
	ok = FVM::Convection::buildFlux_Darcy_Mass(mgr, reg, freg, nmP.a_f_diff, nmP.s_f_diff, "p_g", "rho_g", "mf_g", "Qf_g", "ufn_g", &Pbc, true);
	if (!ok) return false;
	ok = FVM::Convection::build_FaceCoeffs_Upwind(mgr, reg, freg, "T", "mf_g", { "cp_g" }, nmT, Tbc);
	if (!ok) return false;

	//温度方程组装求解
	SparseSystemCOO sysT;
	{
		assemble_COO(mgr, reg, freg, "ddt+diffusion+convection", nmT, &sysT);
		if (lastSysT) *lastSysT = sysT;
		if (ctrl.reportPerIter) {
			auto R = PostChecks::reportAssembly(sysT, false);
			PostChecks::printAssemblyReport(R, "T(Time-implicit + Diffusion + Convection)");
		}
	}
	{
		int N = 0; auto lid = buildUnknownMap(mesh, N);
		auto Tvec = gatherFieldToVec(reg, mesh, "T", lid, N);
		auto optT = ctrl.lin_T; if (optT.tol <= 0.0) optT.tol = ctrl.tol_T_abs;
		double resT = 0.0; int itLT = 0;
		bool okLinT = solveCOO_Eigen(sysT, Tvec, optT, &itLT, &resT);
		if (!okLinT) { std::cerr << "[LinearSolver] temperature solve failed.\n"; return false; }
		scatterVecToField(reg, mesh, "T", lid, Tvec);
	}

	underRelaxInPlace(reg, "T", "T_prev", ctrl.urf_T);
	dT_inf = maxAbsDiff(reg, "T", "T_prev");
	updatePrevIterates(reg, { {"T","T_prev"} });
	return true;

}
//================================================================================================================//

inline bool doOneOutert_constProperties_singlePhase_CO2_T_H_withouPIN
(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const TemperatureBCAdapter& Tbc,
	const PressureBCAdapter& Pbc,
	const Vector& g,
	double dt,
	const SolverControls& ctrl,
	// out
	double& dT_inf,
	double& dp_inf,
	// optional out
	SparseSystemCOO* lastSysP = nullptr,
	SparseSystemCOO* lastSysT = nullptr
)
{
	Mesh& mesh = mgr.mesh();

	// —— 外迭代起步：把上一轮解拷贝到 *_prev，便于欠松弛与收敛评估 —— //
	if (!startOuterIteration_scatter(reg, "p_g", "p_g_prev")) return false;
	if (!startOuterIteration_T(reg, "T", "T_prev")) return false;

	// —— 常物性：更新一次（保持接口一致，后续也可切换变物性） —— //
	ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
	ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);
	ppm.ComputeEffectiveThermalProperties_constProperties_singlePhase_CO2_T_H(mgr, reg);

	using FVM::Diffusion::RhoFaceMethod;

	// =========================
	// 1) 压力：时间项 + 达西扩散（质量式）
	// =========================
	const OperatorFieldNames nmP = makeNames("p_g");

	bool ok = FVM::Diffusion::build_FaceCoeffs_Central
	(
		mgr, reg, freg,
		nmP.a_f_diff, nmP.s_f_diff,
		"p_g",
		{ "kxx:kxx","kyy:kyy","kzz:kzz","/mu_g","rho:rho_g" },  // K/mu ⋅ rho_g
		"rho_g", RhoFaceMethod::Linear,
		g, Pbc,
		/*massForm=*/true, /*alpha_anisotropy=*/0
	);
	if (!ok) return false;

	ok = TimeTerm_FullyImplicit_SinglePhase_Flow(
		mgr, reg, dt,
		"c_phi",             // 孔隙度压缩性
		"phi",               // φ^n
		"p_g_old",           // p^n
		"rho_g",             // ρ^n
		"p_g_prev",          // p^⋆
		"rho_g",             // ρ^⋆
		"Drho_Dp_g",         // (∂ρ/∂p)^⋆
		nmP.a_time, nmP.b_time
		/* strongBCmask = */ // 这里不做行覆盖；采用标准边界处理
	);
	if (!ok) return false;

	SparseSystemCOO sysp;
	assemble_COO(mgr, reg, freg, "ddt+diffusion", nmP, &sysp);
	if (lastSysP) *lastSysP = sysp;
	if (ctrl.reportPerIter) {
		auto R = PostChecks::reportAssembly(sysp, /*brief=*/false);
		PostChecks::printAssemblyReport(R, "P(Time-implicit + Diffusion (mass-form Darcy))");
	}

	int Np = 0; auto lid_p = buildUnknownMap(mesh, Np);
	auto pvec = gatherFieldToVec(reg, mesh, "p_g", lid_p, Np);

	auto optP = ctrl.lin_p;
	if (optP.tol <= 0.0) optP.tol = ctrl.tol_p_abs;

	double resP = 0.0; int itP = 0;
	bool okLinP = solveCOO_Eigen(sysp, pvec, optP, &itP, &resP);
	if (!okLinP) return false;

	scatterVecToField(reg, mesh, "p_g", lid_p, pvec);

	underRelaxInPlace(reg, "p_g", "p_g_prev", ctrl.urf_p);
	dp_inf = maxAbsDiff(reg, "p_g", "p_g_prev");
	updatePrevIterates(reg, { {"p_g","p_g_prev"} });

	// =========================
	// 2) 温度：时间项 + 传导 + 对流（无 PIN，强 Dirichlet 行覆盖）
	// =========================

	// 2.1 边界对应的“需要强制”的单元，以及目标 T_b（面积加权）
	std::vector<char>   maskT = mark_strong_BC_cells(mgr, Tbc);
	std::vector<double> Ttar;build_dirichlet_T_targets(mgr, Tbc, maskT, Ttar);

	const OperatorFieldNames nmT = makeNames("T");

	// 时间项（不做 PIN），保持统一装配接口
	ok = TimeTerm_FullyImplicit_SinglePhase_Temperature(
		mgr, reg, dt,
		"C_eff", "T_old",
		nmT.a_time, nmT.b_time
		/* 无需传 mask/pin */
	);
	if (!ok) return false;

	// 传导项
	ok = FVM::Diffusion::build_FaceCoeffs_Central(
		mgr, reg, freg,
		nmT.a_f_diff, nmT.s_f_diff,
		"T",
		{ "iso:lambda_eff" },  // 等效导热系数
		"",
		RhoFaceMethod::Linear,
		{ 0,0,0 }, Tbc,
		/*massForm=*/false, /*alpha_anisotropy=*/0
	);
	if (!ok) return false;

	// 先算达西质量通量/体积通量/法向速度
	ok = FVM::Convection::buildFlux_Darcy_Mass(
		mgr, reg, freg,
		nmP.a_f_diff, nmP.s_f_diff,
		"p_g", "rho_g",
		"mf_g", "Qf_g", "ufn_g", &Pbc, true
	);
	if (!ok) return false;

	// 对流项（一阶迎风；mf × cp）
	ok = FVM::Convection::build_FaceCoeffs_Upwind(
		mgr, reg, freg,
		"T",
		"mf_g",            // 质量通量：内部已避免再次乘 rho
		{ "cp_g" },        // 上风携带物性
		nmT, Tbc
	);
	if (!ok) return false;

	// 组装线性系统
	SparseSystemCOO sysT;
	assemble_COO(mgr, reg, freg, "ddt+diffusion+convection", nmT, &sysT);
	if (lastSysT) *lastSysT = sysT;
	if (ctrl.reportPerIter) {
		auto R = PostChecks::reportAssembly(sysT, /*brief=*/false);
		PostChecks::printAssemblyReport(R, "T(Time-implicit + Diffusion + Convection)");
	}

	// 强 Dirichlet 行覆盖（无 PIN）
	int Nt = 0; auto lid_t = buildUnknownMap(mesh, Nt);
	apply_strong_dirichlet_rows_T(lid_t, maskT, Ttar, sysT);

	// 求解
	auto Tvec = gatherFieldToVec(reg, mesh, "T", lid_t, Nt);

	auto optT = ctrl.lin_T;
	if (optT.tol <= 0.0) optT.tol = ctrl.tol_T_abs;

	double resT = 0.0; int itT = 0;
	bool okLinT = solveCOO_Eigen(sysT, Tvec, optT, &itT, &resT);
	if (!okLinT) return false;

	scatterVecToField(reg, mesh, "T", lid_t, Tvec);

	// 欠松弛 + 收敛
	underRelaxInPlace(reg, "T", "T_prev", ctrl.urf_T);
	dT_inf = maxAbsDiff(reg, "T", "T_prev");
	updatePrevIterates(reg, { {"T","T_prev"} });

	return true;
}
//================================================================================================================//

inline bool doOneOutert_constProperties_singlePhase_CO2_T_H_closed
(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const TemperatureBCAdapter& Tbc,
	const PressureBCAdapter& Pbc,
	const Vector& g,
	double dt,
	const SolverControls& ctrl,
	// out
	double& dT_inf,
	double& dp_inf,
	// optional out
	SparseSystemCOO* lastSysP = nullptr,
	SparseSystemCOO* lastSysT = nullptr
)
{
	Mesh& mesh = mgr.mesh();

	// —— 外迭代起步：把上一轮解拷贝到 *_prev，便于欠松弛与收敛评估 —— //
	if (!startOuterIteration_scatter(reg, "p_g", "p_g_prev")) return false;
	if (!startOuterIteration_T(reg, "T", "T_prev")) return false;

	// —— 常物性：更新一次（保持接口一致，后续也可切换变物性） —— //
	ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
	ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);
	ppm.ComputeEffectiveThermalProperties_constProperties_singlePhase_CO2_T_H(mgr, reg);

	using FVM::Diffusion::RhoFaceMethod;

	// =========================
	// 1) 压力：时间项 + 达西扩散（质量式）
	// =========================
	const OperatorFieldNames nmP = makeNames("p_g");

	bool ok = FVM::Diffusion::build_FaceCoeffs_Central
	(
		mgr, reg, freg,
		nmP.a_f_diff, nmP.s_f_diff,
		"p_g",
		{ "kxx:kxx","kyy:kyy","kzz:kzz","/mu_g","rho:rho_g" },  // K/mu ⋅ rho_g
		"rho_g", RhoFaceMethod::Linear,
		g, Pbc,
		/*massForm=*/false, /*alpha_anisotropy=*/0
	);
	if (!ok) return false;

	ok = FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Flow
	(
		mgr, reg, dt,
		"c_phi",             // 孔隙度压缩性
		"phi",               // φ^n
		"p_g_old",           // p^n
		"rho_g",             // ρ^n
		"p_g_prev",          // p^⋆
		"rho_g",             // ρ^⋆
		"Drho_Dp_g",         // (∂ρ/∂p)^⋆
		nmP.a_time, nmP.b_time
	);
	if (!ok) return false;

	SparseSystemCOO sysp;
	assemble_COO(mgr, reg, freg, "ddt+diffusion", nmP, &sysp);
	if (lastSysP) *lastSysP = sysp;
	if (ctrl.reportPerIter) {
		auto R = PostChecks::reportAssembly(sysp, /*brief=*/false);
		PostChecks::printAssemblyReport(R, "P(Time-implicit + Diffusion (mass-form Darcy))");
	}

	int Np = 0; auto lid_p = buildUnknownMap(mesh, Np);
	auto pvec = gatherFieldToVec(reg, mesh, "p_g", lid_p, Np);

	auto optP = ctrl.lin_p;
	if (optP.tol <= 0.0) optP.tol = ctrl.tol_p_abs;

	double resP = 0.0; int itP = 0;
	bool okLinP = solveCOO_Eigen(sysp, pvec, optP, &itP, &resP);
	if (!okLinP) return false;

	scatterVecToField(reg, mesh, "p_g", lid_p, pvec);

	underRelaxInPlace(reg, "p_g", "p_g_prev", ctrl.urf_p);
	dp_inf = maxAbsDiff(reg, "p_g", "p_g_prev");
	updatePrevIterates(reg, { {"p_g","p_g_prev"} });

	// =========================
	// 2) 温度：时间项 + 传导 + 对流
	// =========================

	const OperatorFieldNames nmT = makeNames("T");

	// 时间项，保持统一装配接口
	ok = FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Temperature
	(
		mgr, reg, dt,
		"C_eff", "T_old",
		nmT.a_time, nmT.b_time
	);
	if (!ok) return false;

	// 传导项
	ok = FVM::Diffusion::build_FaceCoeffs_Central
	(
		mgr, reg, freg,
		nmT.a_f_diff, nmT.s_f_diff,
		"T",
		{ "iso:lambda_eff" },  // 等效导热系数
		"",
		RhoFaceMethod::Linear,
		{ 0,0,0 }, Tbc,
		/*massForm=*/false, /*alpha_anisotropy=*/0
	);
	if (!ok) return false;

	// 先算达西质量通量/体积通量/法向速度
	ok = FVM::Convection::buildFlux_Darcy_Mass
	(
		mgr, reg, freg,
		nmP.a_f_diff, nmP.s_f_diff,
		"p_g", "rho_g",
		"mf_g", "Qf_g", "ufn_g", &Pbc, true
	);

	if (!ok) return false;

	debugCheckMassFlux(mgr, freg, "mf_g", 1e-20);

	

	//// 1) patch 查询
	//const auto& bfaces = mgr.boundaryFaces();
	//auto lookup = BoundaryUtils::makePatchFaceLookup(bfaces);

	//// 2) 条件型入边温度：仅在【温度绝热 && 压力非零通量 && flux<0 && 属于右边界】时给 T_prev(owner)
	//auto provT_cond_right = BoundaryUtils::makeConditionalPrevTimeValueProviderForPatches(
	//	mgr, reg, freg,
	//	lookup, { "xL" },
	//	/*fieldName_prev*/ "T_prev",
	//	/*T*/ Tbc,
	//	/*P*/ Pbc,
	//	/*flux for temperature*/ "mf_g",   // 或 "Qf_g"，与温度方程使用的通量一致
	//	/*epsFlux*/ 1e-18
	//);

	//// 5) 注册到 scalarBC（关键：键名必须与 var_key 一致，这里 var_key="T"）
	//std::unordered_map<std::string, std::function<bool(int, double&)>> scalarBC;
	//scalarBC["T"] = provT_cond_right;


	// 对流项（一阶迎风；mf × cp）
	ok = FVM::Convection::build_FaceCoeffs_Upwind
	(
		mgr, reg, freg,
		"T",
		"mf_g",            // 质量通量：内部已避免再次乘 rho
		{ "cp_g" },        // 上风携带物性
		nmT, Tbc
	);
	if (!ok) return false;

	// 组装线性系统
	SparseSystemCOO sysT;
	assemble_COO(mgr, reg, freg, "ddt+diffusion+convection", nmT, &sysT);
	if (lastSysT) *lastSysT = sysT;
	if (ctrl.reportPerIter) {
		auto R = PostChecks::reportAssembly(sysT, /*brief=*/false);
		PostChecks::printAssemblyReport(R, "T(Time-implicit + Diffusion + Convection)");
	}

	int Nt = 0; auto lid_t = buildUnknownMap(mesh, Nt);

	// 求解
	auto Tvec = gatherFieldToVec(reg, mesh, "T", lid_t, Nt);

	auto optT = ctrl.lin_T;
	if (optT.tol <= 0.0) optT.tol = ctrl.tol_T_abs;

	double resT = 0.0; int itT = 0;
	bool okLinT = solveCOO_Eigen(sysT, Tvec, optT, &itT, &resT);
	if (!okLinT) return false;

	scatterVecToField(reg, mesh, "T", lid_t, Tvec);

	// 欠松弛 + 收敛
	underRelaxInPlace(reg, "T", "T_prev", ctrl.urf_T);
	dT_inf = maxAbsDiff(reg, "T", "T_prev");
	updatePrevIterates(reg, { {"T","T_prev"} });

	return true;
}





inline bool outerIter_constProperties_singlePhase_CO2_T_H
(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const TemperatureBCAdapter& Tbc,
	const PressureBCAdapter& Pbc,
	const Vector& g,
	double dt,
	const SolverControls& ctrl
)
{
	double prev_dp_g = 1e300;
	double prev_dT = 1e300;

	for (int it = 0; it < ctrl.maxOuter; ++it)
	{
		double dp = 0.0;
		double dT = 0.0;
		SparseSystemCOO lastP;
		SparseSystemCOO lastT;
		//一次外迭代推进
		bool ok = doOneOutert_constProperties_singlePhase_CO2_T_H_closed(mgr, reg, freg, ppm, Tbc, Pbc, g, dt, ctrl, dT, dp, (ctrl.dumpMMOnLastIter ? &lastP : nullptr),(ctrl.dumpMMOnLastIter ? &lastT : nullptr));
		if (!ok) return false;

		// —— 进度输出 —— //
		std::cout << "[Outer " << it << "]  |Δp|_inf=" << dp << "  |ΔT|_inf=" << dT << "\n";

		// —— 自适应松弛—— //
		if (it > 0)
		{
			double rT = dT / std::max(prev_dT, 1e-30);
			double uT = (rT < 0.7) ? +0.05 : (rT > 0.95 ? -0.05 : 0.0);
			auto& ctrl_mut = const_cast<SolverControls&>(ctrl);
			ctrl_mut.urf_T = std::min(0.7, std::max(0.15, ctrl_mut.urf_T + uT));
			double rp = dp / std::max(prev_dp_g, 1e-30);
			double up = (rp < 0.7) ? +0.05 : (rp > 0.95 ? -0.05 : 0.0);
			ctrl_mut.urf_p = std::min(0.7, std::max(0.15, ctrl_mut.urf_p + up));
		}

		prev_dp_g = dp;
		prev_dT = dT;
		// 收敛判据：绝对/相对（相对采用场量标度）
		auto maxAbsField = [&](const std::string& fld)->double {
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

		const double PScale = maxAbsField("p_g");
		const double TScale = maxAbsField("T");

		const bool convP = (dp < std::max(ctrl.tol_p_abs, ctrl.tol_p_rel * PScale));
		const bool convT = (dT < std::max(ctrl.tol_T_abs, ctrl.tol_T_rel * TScale));

		if (convP && convT) {
			std::cout << "Converged at outer iter " << it << "\n";

			if (ctrl.dumpMMOnLastIter) {
#if __cplusplus >= 201703L
				try { std::filesystem::create_directories("mm"); }
				catch (...) {}
#endif
				// 同时导出 P/T 的矩阵账本（便于核对边界行与对流行）
				PostChecks::dumpCOO_to_matrix_market(
					lastP,
					"mm/A_P_CO2_pTH.mtx",
					"mm/b_P_CO2_pTH.txt",
					/*sym=*/false
				);
				PostChecks::dumpCOO_to_matrix_market(
					lastT,
					"mm/A_T_CO2_pTH.mtx",
					"mm/b_T_CO2_pTH.txt",
					/*sym=*/false
				);
			}
			break;
		}

		if (it == ctrl.maxOuter - 1) {
			std::cout << "Reached maxOuter without meeting P/T tolerances.\n";
		}
	}

	return true;


}



inline bool doOneOutert_constProperties_singlePhase_CO2_T_H_closed_withWell
(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const TemperatureBCAdapter& Tbc,
	const PressureBCAdapter& Pbc,
	const Vector& g,
	double dt,
	const SolverControls& ctrl,
	const std::vector<WellConfig>& wellsCfg_in,
	// out
	double& dT_inf,
	double& dp_inf,
	// optional out
	SparseSystemCOO* lastSysP = nullptr,
	SparseSystemCOO* lastSysT = nullptr
)
{

	Mesh& mesh = mgr.mesh();
	const int Nc = (int)mesh.getCells().size();



	// —— 外迭代起步（保持原有） —— //
	if (!startOuterIteration_scatter(reg, "p_g", "p_g_prev")) return false;
	if (!startOuterIteration_T(reg, "T", "T_prev")) return false;

	// —— 常物性刷新（你工程内已有三件套，名字按你项目） —— //
	ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
	ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);
	ppm.ComputeEffectiveThermalProperties_constProperties_singlePhase_CO2_T_H(mgr, reg);

	//离散系数和源项场名称后缀  可以从外部输入 在main函数定义
	const OperatorFieldNames nmP = makeNames("p_g");
	const OperatorFieldNames nmT = makeNames("T");


	//case1 不包含注/采井
	if (wellsCfg_in.size() == 0)
	{
		// ======================================================
		// 1.1) 压力方程：离散（ddt+Darcy扩散）
		// ======================================================

		// a 扩散项
		if (!FVM::Diffusion::build_FaceCoeffs_Central(mgr,reg,freg,nmP.a_f_diff,nmP.s_f_diff,"p_g" ,{ "kxx:kxx","kyy:kyy","kzz:kzz","/mu_g","rho:rho_g" },"rho_g", FVM::Diffusion::RhoFaceMethod::Linear, g, Pbc,false,3 )) return false;

		// b 时间项 
		if (!FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Flow(mgr, reg, dt, "c_phi", "phi", "p_g_old", "rho_g", "p_g_prev", "rho_g", "Drho_Dp_g", nmP.a_time, nmP.b_time)) return false;

		// ======================================================
		// 1.2) 压力方程：装配与求解
		// ======================================================

		SparseSystemCOO sysp;
		if (!assemble_COO(mgr, reg, freg, "ddt+diffusion", nmP, &sysp)) return false;

		int Np = 0; auto lid_p = buildUnknownMap(mesh, Np);
		auto pvec = gatherFieldToVec(reg, mesh, "p_g", lid_p, Np);

		auto optP = ctrl.lin_p;
		if (optP.tol <= 0.0) optP.tol = ctrl.tol_p_abs;

		double resP = 0.0; int itP = 0;
		bool okLinP = solveCOO_Eigen(sysp, pvec, optP, &itP, &resP);
		if (!okLinP) return false;

		scatterVecToField(reg, mesh, "p_g", lid_p, pvec);

		underRelaxInPlace(reg, "p_g", "p_g_prev", ctrl.urf_p);
		dp_inf = maxAbsDiff(reg, "p_g", "p_g_prev");
		updatePrevIterates(reg, { {"p_g","p_g_prev"} });

		// ======================================================
		// 2.1) 温度方程离散：时间项 + 传导 + 对流
		// ======================================================

		std::vector<char>   maskT = mark_strong_BC_cells(mgr, Tbc);
		std::vector<double> Ttar; build_dirichlet_T_targets(mgr, Tbc, maskT, Ttar);

		//a 时间项
		if (! FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Temperature(mgr, reg, dt, "C_eff", "T_old", nmT.a_time, nmT.b_time)) return false;

		//b 扩散项
		if (!FVM::Diffusion::build_FaceCoeffs_Central(mgr, reg, freg, nmT.a_f_diff, nmT.s_f_diff, "T", { "iso:lambda_eff" }, "", FVM::Diffusion::RhoFaceMethod::Linear, g, Tbc, false, 3)) return true;

		//c.1 计算质量通量
		if (!FVM::Convection::buildFlux_Darcy_Mass(mgr, reg, freg,nmP.a_f_diff, nmP.s_f_diff,"p_g", "rho_g","mf_g", "Qf_g", "ufn_g", &Pbc, true)) return false;
		
		//c.2 对流项，利用质量通量，采用一阶迎风格式计算面系数
		if (!FVM::Convection::build_FaceCoeffs_Upwind(mgr, reg, freg,"T","mf_g",{ "cp_g" },nmT, Tbc)) return false;

		// ======================================================
		// 2.2) 温度方程：装配与求解
		// ======================================================

		SparseSystemCOO sysT;
		assemble_COO(mgr, reg, freg, "ddt+diffusion+convection", nmT, &sysT);
		if (lastSysT) *lastSysT = sysT;
		if (ctrl.reportPerIter) {
			auto R = PostChecks::reportAssembly(sysT, /*brief=*/false);
			PostChecks::printAssemblyReport(R, "T(Time-implicit + Diffusion + Convection)");
		}

		// 强 Dirichlet 行覆盖
		int Nt = 0; auto lid_t = buildUnknownMap(mesh, Nt);
		apply_strong_dirichlet_rows_T(lid_t, maskT, Ttar, sysT);

		// 求解
		auto Tvec = gatherFieldToVec(reg, mesh, "T", lid_t, Nt);

		auto optT = ctrl.lin_T;
		if (optT.tol <= 0.0) optT.tol = ctrl.tol_T_abs;

		double resT = 0.0; int itT = 0;
		bool okLinT = solveCOO_Eigen(sysT, Tvec, optT, &itT, &resT);
		if (!okLinT) return false;

		scatterVecToField(reg, mesh, "T", lid_t, Tvec);

		// 欠松弛 + 收敛
		underRelaxInPlace(reg, "T", "T_prev", ctrl.urf_T);
		dT_inf = maxAbsDiff(reg, "T", "T_prev");
		updatePrevIterates(reg, { {"T","T_prev"} });

		return true;
	}
	else {
		// 1.1 扩散（Darcy mass-form）
		if (!FVM::Diffusion::build_FaceCoeffs_Central(
			mgr, reg, freg,
			nmP.a_f_diff, nmP.s_f_diff,
			"p_g",
			{ "kxx:kxx","kyy:kyy","kzz:kzz","/mu_g","rho:rho_g" }, // (K/mu)*rho
			"rho_g", FVM::Diffusion::RhoFaceMethod::Linear,
			g, Pbc,
			/*massForm=*/true, /*alpha_aniso=*/0)) return false;

		// 1.2 时间项
		if (!FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Flow(
			mgr, reg, dt,
			"c_phi", "phi", "p_g_old",
			"rho_g", "p_g_prev",
			"rho_g", "Drho_Dp_g",
			nmP.a_time, nmP.b_time)) return false;

		// 1.3 组装基础系统（不含井）
		SparseSystemCOO sysp;
		if (!assemble_COO(mgr, reg, freg, "ddt+diffusion", nmP, &sysp)) return false;

		// 1.4 多井：构建/刷新掩码与 PI、注册井 DoF、扩维、列耦合与井行
		std::vector<WellConfig> wellsCfg = wellsCfg_in;   // 本地拷贝，便于补齐字段名
		build_masks_and_PI_for_all(mgr, reg, wellsCfg);   // 每口井：mask_<name>/PI_<name>

		std::vector<WellDOF> wells;
		const int Ntot = register_well_dofs_for_all(Nc, wellsCfg, wells); // 回填 wells[k].lid
		extend_linear_system_size(sysp, Ntot);                             // !!! 必须在 addA/addb 前

		std::vector<int> lid_cell(Nc); for (int i = 0; i < Nc; ++i) lid_cell[i] = i;

		for (const auto& w : wells) {
			// 列耦合：质量导纳 PI，故 scaleByCellMeasure=false
			add_peaceman_coupling_cell_rows(
				sysp, mesh, reg, w.PI_field, w.mask_field, lid_cell, w.lid, /*false*/false);
			// 井行：Pressure/Rate 通吃
			add_well_row(
				sysp, mesh, reg, w.PI_field, w.mask_field, lid_cell, w.lid, w.mode, w.target, /*false*/false);
		}
		if (lastSysP) *lastSysP = sysp;

		// 1.5 —— 压力求解段：严格保留你的原有结构 —— //
		int Np = 0; auto lid_p = buildUnknownMap(mesh, Np);
		auto pvec = gatherFieldToVec(reg, mesh, "p_g", lid_p, Np);

		// 追加井未知到同一解向量尾部（保持原结构，仅加两步：resize + 初猜）
		pvec.resize(Ntot, 0.0);
		for (const auto& w : wells) {
			if (w.mode == WellDOF::Mode::Pressure) pvec[w.lid] = w.target; // BHP 作为井压初猜
			else                                   pvec[w.lid] = pvec[0];   // 定流：随便给个初猜
		}

		auto optP = ctrl.lin_p; if (optP.tol <= 0.0) optP.tol = ctrl.tol_p_abs;
		double resP = 0.0; int itP = 0;
		bool okLinP = solveCOO_Eigen(sysp, pvec, optP, &itP, &resP);
		if (!okLinP) return false;

		// 只把“前 Np 分量（单元压力）”散回 p_g（与你原逻辑一致）
		{
			std::vector<double> p_cells(Np);
			std::copy(pvec.begin(), pvec.begin() + Np, p_cells.begin());
			scatterVecToField(reg, mesh, "p_g", lid_p, p_cells);
		}

		// 写回每口井的井底压 p_w_<name>（单值场）
		writeback_pw_fields_for_all(reg, wells, pvec);

		// 欠松弛与收敛评估（原结构不变）
		underRelaxInPlace(reg, "p_g", "p_g_prev", ctrl.urf_p);
		dp_inf = maxAbsDiff(reg, "p_g", "p_g_prev");
		updatePrevIterates(reg, { {"p_g","p_g_prev"} });

		// ======================================================
		// 2) 温度方程：ddt + 传导 + 对流 + 基于 q_m 的井热源
		// ======================================================
		std::vector<char>   maskT = mark_strong_BC_cells(mgr, Tbc);
		std::vector<double> Ttar; build_dirichlet_T_targets(mgr, Tbc, maskT, Ttar);



		if (!FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Temperature(
			mgr, reg, dt, "C_eff", "T_old", nmT.a_time, nmT.b_time)) return false;

		// 1) 纯正交离散

		if (!FVM::Diffusion::build_FaceCoeffs_Central(
			mgr, reg, freg,
			nmT.a_f_diff, nmT.s_f_diff,
			"T",
			{ "iso:lambda_eff" },
			"", FVM::Diffusion::RhoFaceMethod::Linear,
			Vector(0, 0, 0), Tbc,
			/*enable_buoy=*/false,
			/*gradSmoothIters=*/1)) return false;


		if (!FVM::Convection::buildFlux_Darcy_Mass(
			mgr, reg, freg,
			nmP.a_f_diff, nmP.s_f_diff,
			"p_g", "rho_g",
			"mf_g", "Qf_g", "ufn_g", &Pbc, true)) return false;

		if (!FVM::Convection::build_FaceCoeffs_Upwind(
			mgr, reg, freg,
			"T",
			"mf_g",
			{ "cp_g" },
			nmT, Tbc)) return false;

		// —— 多井 q_m 热源：注井加 b（焓输入），采井加 a（抽能） —— //
		const FVM::SourceTerm::InjectorDirichletPinOptions pinOpts
		{
		ctrl.injectorPin.enable,
		ctrl.injectorPin.relativeWeight,
		ctrl.injectorPin.minWeight
		};

		bool first = true;
		for (const auto& w : wells) {
			const double Tin = (w.role == WellDOF::Role::Injector) ? w.Tin : 0.0;
			if (!FVM::SourceTerm::add_temperature_source_from_qm_single(
				mgr, reg,
				w.role,
				w.PI_field.c_str(), w.mask_field.c_str(), w.name.c_str(),
				Tin,
				/*p*/ "p_g", /*cp*/ "cp_g",
				/*thickness*/ 1.0,
				/*accumulate*/ !first,
				/*verbose*/ false,
				pinOpts))                                  // <- new argument
				return false;
			first = false;
		}

		// 2.1 —— 温度求解段：严格保留你的原有结构 —— //
		SparseSystemCOO sysT;
		if (!assemble_COO(mgr, reg, freg, "ddt+diffusion+convection+source", nmT, &sysT)) return false;
		if (lastSysT) *lastSysT = sysT;

		int Nt = 0; auto lid_t = buildUnknownMap(mesh, Nt);
		apply_strong_dirichlet_rows_T(lid_t, maskT, Ttar, sysT);

		auto Tvec = gatherFieldToVec(reg, mesh, "T", lid_t, Nt);

		auto optT = ctrl.lin_T; if (optT.tol <= 0.0) optT.tol = ctrl.tol_T_abs;
		double resT = 0.0; int itT = 0;
		bool okLinT = solveCOO_Eigen(sysT, Tvec, optT, &itT, &resT);
		if (!okLinT) return false;

		scatterVecToField(reg, mesh, "T", lid_t, Tvec);

		underRelaxInPlace(reg, "T", "T_prev", ctrl.urf_T);
		dT_inf = maxAbsDiff(reg, "T", "T_prev");
		updatePrevIterates(reg, { {"T","T_prev"} });


		return true;
	}
	

	

}


inline bool outerIter_constProperties_singlePhase_CO2_T_H_withWell
(
	MeshManager& mgr,
	FieldRegistry& reg,
	FaceFieldRegistry& freg,
	PhysicalPropertiesManager& ppm,
	const TemperatureBCAdapter& Tbc,
	const PressureBCAdapter& Pbc,
	const Vector& g,
	const std::vector<WellConfig>& wellsCfg_in,
	double dt,
	const SolverControls& ctrl
)
{
	double prev_dp_g = 1e300;
	double prev_dT = 1e300;

	for (int it = 0; it < ctrl.maxOuter; ++it)
	{
		double dp = 0.0;
		double dT = 0.0;
		SparseSystemCOO lastP;
		SparseSystemCOO lastT;
		//一次外迭代推进
		bool ok = doOneOutert_constProperties_singlePhase_CO2_T_H_closed_withWell(mgr, reg, freg, ppm, Tbc, Pbc, g, dt, ctrl, wellsCfg_in ,dT, dp, (ctrl.dumpMMOnLastIter ? &lastP : nullptr), (ctrl.dumpMMOnLastIter ? &lastT : nullptr));
		if (!ok) return false;

		// —— 进度输出 —— //
		std::cout << "[Outer " << it << "]  |Δp|_inf=" << dp << "  |ΔT|_inf=" << dT << "\n";

		// —— 自适应松弛—— //
		if (it > 0)
		{
			double rT = dT / std::max(prev_dT, 1e-30);
			double uT = (rT < 0.7) ? +0.05 : (rT > 0.95 ? -0.05 : 0.0);
			auto& ctrl_mut = const_cast<SolverControls&>(ctrl);
			ctrl_mut.urf_T = std::min(0.7, std::max(0.15, ctrl_mut.urf_T + uT));
			double rp = dp / std::max(prev_dp_g, 1e-30);
			double up = (rp < 0.7) ? +0.05 : (rp > 0.95 ? -0.05 : 0.0);
			ctrl_mut.urf_p = std::min(0.7, std::max(0.15, ctrl_mut.urf_p + up));
		}

		prev_dp_g = dp;
		prev_dT = dT;
		// 收敛判据：绝对/相对（相对采用场量标度）
		auto maxAbsField = [&](const std::string& fld)->double {
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

		const double PScale = maxAbsField("p_g");
		const double TScale = maxAbsField("T");

		const bool convP = (dp < std::max(ctrl.tol_p_abs, ctrl.tol_p_rel * PScale));
		const bool convT = (dT < std::max(ctrl.tol_T_abs, ctrl.tol_T_rel * TScale));

		if (convP && convT) {
			std::cout << "Converged at outer iter " << it << "\n";

			if (ctrl.dumpMMOnLastIter) {
#if __cplusplus >= 201703L
				try { std::filesystem::create_directories("mm"); }
				catch (...) {}
#endif
				// 同时导出 P/T 的矩阵账本（便于核对边界行与对流行）
				PostChecks::dumpCOO_to_matrix_market(
					lastP,
					"mm/A_P_CO2_pTH.mtx",
					"mm/b_P_CO2_pTH.txt",
					/*sym=*/false
				);
				PostChecks::dumpCOO_to_matrix_market(
					lastT,
					"mm/A_T_CO2_pTH.mtx",
					"mm/b_T_CO2_pTH.txt",
					/*sym=*/false
				);
			}
			break;
		}

		if (it == ctrl.maxOuter - 1) {
			std::cout << "Reached maxOuter without meeting P/T tolerances.\n";
		}
	}

	return true;


}

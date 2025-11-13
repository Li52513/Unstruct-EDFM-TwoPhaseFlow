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
#include"FVM_SourceTerm_ProdWellOps.h"
#include "FVM_WellCoupling.h"
#include "FVM_WellDOF.h"
#include "FVM_SourceTerm_WellHeat.h"
#include "WellConfig.h"
#include "InflowPatch.h"


//==================== 欠松弛应用 ====================//
inline void underRelax(std::vector<double>& dest, const std::vector<double>& src, double urf)
{
	const size_t n = dest.size();
	for (size_t i = 0; i < n; ++i) dest[i] = dest[i] + urf * (src[i] - dest[i]);
}

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

inline void scatterVecToField(FieldRegistry& reg, Mesh& mesh, const std::string& fld, const std::vector<int>& lid_of_cell, const std::vector<double>& x)
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
	LinearSolverOptions lin_p; // 压力求解设置
	LinearSolverOptions lin_T; // 温度求解设置
	bool useJacobi = false;    // 默认用 Krylov，必要时退回 Jacobi


	// 可选：是否每次迭代打印装配报告（默认 false）
	bool reportPerIter = false;
	// 可选：是否导出 MatrixMarket（仅在 reportPerIter 为 true 时，且只导出最后一次迭代）
	bool dumpMMOnLastIter = false;

	//**修改
	// Pressure inner sweeps & flux refresh controls
	int    NsweepP_max = 3;          // 每次外迭代允许的压力子扫次数
	double tol_p_inner = 5e-7;       // 结束压力子扫的 L∞ 阈值
	double kappa_p_for_flux_update = 0.3;   // 用于判定“需全量重建通量”的放宽系数
	bool   enable_incremental_convection = true;
	double incremental_flip_ratio = 0.01;   // 若翻转面占比低于该值，则按需重装对流面
	double flux_sign_epsilon = 1e-14;       // 面通量符号判定的零阈值

	// Linear solver reuse / preconditioner刷新
	bool reuse_linear_pattern = true;
	int  rebuild_precond_every = 10;

	// Aitken 欠松弛
	bool   enable_Aitken_p = true;
	bool   enable_Aitken_T = true;
	double aitken_omega_min = 0.05;
	double aitken_omega_max = 1.5;

	// 温度 CFL 守门配置
	bool   enable_CFL_guard = true;
	double CFL_T_threshold = 0.7;
	double CFL_dt_scale_min = 0.2;

	// 时间步自适应 (基于 CFL)
	bool   enable_dt_adapt = true;
	double dt_adapt_CFL_target = 0.6;
	double dt_adapt_CFL_hysteresis = 0.3;
	double dt_adapt_grow = 1;
	double dt_adapt_shrink = 0.5;
	double dt_min = 1e-3;
	double dt_max = 1e6;
};


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
		if (!FVM::Diffusion::build_FaceCoeffs_Central(mgr,reg,freg,nmP.a_f_diff,nmP.s_f_diff,"p_g" ,{ "kxx:kxx","kyy:kyy","kzz:kzz","/mu_g","rho:rho_g" },"rho_g", FVM::Diffusion::RhoFaceMethod::Linear, g, Pbc,false,0 )) return false;

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

		//a 时间项
		if (! FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Temperature(mgr, reg, dt, "C_eff", "T_old", nmT.a_time, nmT.b_time)) return false;

		//b 扩散项
		if (!FVM::Diffusion::build_FaceCoeffs_Central(mgr, reg, freg, nmT.a_f_diff, nmT.s_f_diff, "T", { "iso:lambda_eff" }, "", FVM::Diffusion::RhoFaceMethod::Linear, g, Tbc, false, 0)) return false;

		//c.1 计算质量通量
		if (!FVM::Convection::buildFlux_Darcy_Mass(mgr, reg, freg,nmP.a_f_diff, nmP.s_f_diff,"p_g", "rho_g","mf_g", "Qf_g", "ufn_g", &Pbc, true)) return false;

		//检查通量正确性
		debugCheckMassFlux(mgr, freg, "mf_g", 1e-20);
		
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
	else {
		// 1.1 扩散（Darcy mass-form）
		if (!FVM::Diffusion::build_FaceCoeffs_Central(mgr, reg, freg, nmP.a_f_diff, nmP.s_f_diff, "p_g", { "kxx:kxx","kyy:kyy","kzz:kzz","/mu_g","rho:rho_g" }, "rho_g", FVM::Diffusion::RhoFaceMethod::Linear, g, Pbc, false, 0)) return false;

		// 1.2 时间项
		if (!FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Flow(mgr, reg, dt, "c_phi", "phi", "p_g_old", "rho_g", "p_g_prev", "rho_g", "Drho_Dp_g", nmP.a_time, nmP.b_time)) return false;

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

		for (const auto& w : wells) 
		{
			// 列耦合：质量导纳 PI，故 scaleByCellMeasure=false
			add_peaceman_coupling_cell_rows(
				sysp, mesh, reg, w.PI_field, w.mask_field, lid_cell, w.lid, /*false*/false);
			// 井行：Pressure/Rate 通吃
			add_well_row(
				sysp, mesh, reg, w.PI_field, w.mask_field, lid_cell, w.lid, w.mode, w.target, /*false*/false);
		}
		if (lastSysP) *lastSysP = sysp;

		// 1.5 —— 压力求解
		int Np = 0; auto lid_p = buildUnknownMap(mesh, Np);
		auto pvec = gatherFieldToVec(reg, mesh, "p_g", lid_p, Np);

		// 追加井未知到同一解向量尾部（保持原结构，仅加两步：resize + 初猜）
		pvec.resize(Ntot, 0.0);
		for (const auto& w : wells) 
		{
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



		//a 时间项
		if (!FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Temperature(mgr, reg, dt, "C_eff", "T_old", nmT.a_time, nmT.b_time)) return false;

		//b 扩散项
		if (!FVM::Diffusion::build_FaceCoeffs_Central(mgr, reg, freg, nmT.a_f_diff, nmT.s_f_diff, "T", { "iso:lambda_eff" }, "", FVM::Diffusion::RhoFaceMethod::Linear, g, Tbc, false, 0)) return false;

		//c.1 计算质量通量
		if (!FVM::Convection::buildFlux_Darcy_Mass(mgr, reg, freg, nmP.a_f_diff, nmP.s_f_diff, "p_g", "rho_g", "mf_g", "Qf_g", "ufn_g", &Pbc, true)) return false;

		//检查通量正确性
		debugCheckMassFlux(mgr, freg, "mf_g", 1e-20);

		//c.2 对流项，利用质量通量，采用一阶迎风格式计算面系数
		if (!FVM::Convection::build_FaceCoeffs_Upwind(mgr, reg, freg, "T", "mf_g", { "cp_g" }, nmT, Tbc)) return false;

		bool first = true;
		for (const auto& w : wells) 
		{
			const double Tin = (w.role == WellDOF::Role::Injector) ? w.Tin : 0.0;
			if (!FVM::SourceTerm::add_temperature_source_from_qm_single(
				mgr, reg,
				w.role,
				w.PI_field.c_str(), w.mask_field.c_str(), w.name.c_str(),
				Tin,
				/*p*/ "p_g", /*cp*/ "cp_g",
				/*thickness*/ 1.0,
				/*accumulate*/ !first,
				/*verbose*/ false))                                  // <- new argument
				return false;
			first = false;
		}

		// 2.1 —— 温度求解段：严格保留你的原有结构 —— //
		SparseSystemCOO sysT;
		if (!assemble_COO(mgr, reg, freg, "ddt+diffusion+convection+source", nmT, &sysT)) return false;
		if (lastSysT) *lastSysT = sysT;

		int Nt = 0; auto lid_t = buildUnknownMap(mesh, Nt);

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







#pragma once
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"

#include "SolverContrlStrName.h"
#include "Solver_AssemblerCOO.h"
#include "Solver_Tools.h"
#include "BCAdapter.h"

#include "DiffusionCentral.h"
#include "Timeterm_BDF.h"

#include "LinearSolver_Eigen.h"
#include "PhysicalPropertiesManager.h"
#include "WellConfig.h"
#include "FVM_WellCoupling.h"

namespace SinglePhase {

	struct PressureAssemblyConfig 
	{
		PressureEquation_String						P_Eq_str;
		std::string operator_tag =					P_Eq_str.operator_tag;                   // pressure operator tag for nm
		std::string pressure_field =				P_Eq_str.pressure_field;                 // current eval pressure field
		std::string pressure_old_field =			P_Eq_str.pressure_old_field;
		std::string pressure_prev_field =			P_Eq_str.pressure_prev_field;
		Vector gravity = { 0.0, 0.0, 0.0 };
		bool enable_buoyancy = false;
		int gradient_smoothing = 0;
	};

	struct PressureSolveControls
	{
		PressureAssemblyConfig assembly;
		LinearSolverOptions    linear;
		double under_relax = 1.0;		///< 欠松弛系数 (0<urf<=1, 建议 0.3~0.7)
		double tol_abs = 1e-6;			///绝对残差
		double tol_rel = 1e-6;			///相对残差
		int max_outer = 10;				///最大外迭代次数
		bool   verbose = false;			///< 是否打印每轮外迭代的 dp_inf / linRes
	};
	/**
	 * \brief 单次压力组装 + 线性求解 的报告
	 */
	struct PressureStepReport
	{
		double lin_residual = 0.0;      ///< 线性求解器返回残差（如果可用）
		int    lin_iterations = 0;      ///< 线性求解器迭代步数
		double dp_inf = 0.0;            ///< 本次 outer 中压力场的无穷范数变化
		int    outer_iterations = 0;   ///< 实际使用的外迭代次数
	};

	struct PressureAssemblyResult
	{
		SparseSystemCOO system;
		std::vector<int> cell_lid;
	};
	/**
	*  brief: 该函数用于单相压力方程系数矩阵组装
	*/
	inline bool AssembleandSolver_PressureEq_SinglePhase(
		MeshManager& mgr,
		FieldRegistry& reg,
		FaceFieldRegistry& freg,
		PhysicalPropertiesManager& ppm,
		PressureAssemblyConfig& ctrl,
		PressureAssemblyResult& result,
		PressureSolveControls& sol_ctrl,
		PressureStepReport& rep,
		const PressureBCAdapter& Pbc,
		double dt,
		const std::vector<WellConfig>& wellsCfg_in,
		// per-outer 输出量：
		double& dp_inf,
		double& lin_residual,
		int& lin_iterations
	)
	{
		Mesh& mesh = mgr.mesh();
		const int Nc = (int)mesh.getCells().size();
		const OperatorFieldNames nmP = makeNames(ctrl.operator_tag);
		if (!GeneralTools::startOuterIteration_scatter(reg,ctrl.pressure_field,ctrl.pressure_prev_field)) return false;
		if (wellsCfg_in.size() == 0)
		{
			// a 扩散项离散
			std::vector<std::string> mobility_tokens = {"kxx:kxx", "kyy:kyy", "kzz:kzz" ,"/mu_g","rho:rho_g"};
			if (!FVM::Diffusion::build_FaceCoeffs_Central(mgr, reg, freg, nmP.a_f_diff, nmP.s_f_diff, ctrl.pressure_field, mobility_tokens, "rho_g", FVM::Diffusion::RhoFaceMethod::Linear, ctrl.gravity, Pbc, ctrl.enable_buoyancy, ctrl.gradient_smoothing))
			{
				std::cerr << "[SinglePhase][Pressure] diffusion operator build failed.\n";
				return false;
			}

			// b 时间项离散
			if (!FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Flow(mgr, reg, dt, "c_phi", "phi_r", ctrl.pressure_old_field, "rho_g", ctrl.pressure_prev_field, "rho_g", "Drho_Dp_g", nmP.a_time, nmP.b_time))
			{
				std::cerr << "[SinglePhase][Pressure] Timeterm operator build failed.\n";
				return false;
			}

			// c 系数矩阵组装
			SparseSystemCOO sys;
			if (!assemble_COO(mgr, reg, freg, "ddt+diffusion", nmP, &sys))
			{
				std::cerr << "[SinglePhase][Pressure] Assemble process build failed.\n";
				return false;
			}
			result.system = sys;	//储存进PressureAssemblyResult
			int N = 0;
			auto lid_cell = buildUnknownMap(mgr.mesh(), N);
			result.cell_lid = lid_cell;	//储存进PressureAssemblyResult

			// d 备份上一外迭代层的值并完成求解
			auto pvec = GeneralTools::gatherFieldToVec(reg, mesh, ctrl.pressure_field, lid_cell, N);

			double linRes = 0.0;
			int    linIters = 0;
			if (!solveCOO_Eigen(sys, pvec, sol_ctrl.linear, &linIters, &linRes))
			{
				std::cerr << "[SinglePhase][Pressure] Pressure equation solution failed.\n";
				return false;
			}

			GeneralTools::scatterVecToField(reg, mesh, ctrl.pressure_field, lid_cell, pvec);
			double dpInf = 0.0;
			dpInf = GeneralTools::maxAbsDiff(reg, ctrl.pressure_field, ctrl.pressure_prev_field);
			GeneralTools::underRelaxInPlace(reg, ctrl.pressure_field, ctrl.pressure_prev_field, sol_ctrl.under_relax);

			// e 输出结果
			dp_inf = dpInf;
			lin_residual = linRes;
			lin_iterations = linIters;
			return true;
		}
		else 
		{
			// a 扩散项离散
			std::vector<std::string> mobility_tokens = { "kxx:kxx", "kyy:kyy", "kzz:kzz" ,"/mu_g","rho:rho_g" };
			if (!FVM::Diffusion::build_FaceCoeffs_Central(mgr, reg, freg, nmP.a_f_diff, nmP.s_f_diff, ctrl.pressure_field, mobility_tokens, "rho_g", FVM::Diffusion::RhoFaceMethod::Linear, ctrl.gravity, Pbc, ctrl.enable_buoyancy, ctrl.gradient_smoothing))
			{
				std::cerr << "[SinglePhase][Pressure] diffusion operator build failed.\n";
				return false;
			}

			// b 时间项离散
			if (!FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Flow(mgr, reg, dt, "c_phi", "phi_r", ctrl.pressure_old_field, "rho_g", ctrl.pressure_prev_field, "rho_g", "Drho_Dp_g", nmP.a_time, nmP.b_time))
			{
				std::cerr << "[SinglePhase][Pressure] Timeterm operator build failed.\n";
				return false;
			}

			// c 系数矩阵组装
			SparseSystemCOO sys;
			if (!assemble_COO(mgr, reg, freg, "ddt+diffusion", nmP, &sys))
			{
				std::cerr << "[SinglePhase][Pressure] Assemble process build failed.\n";
				return false;
			}
			int N = 0;
			auto lid_cell = buildUnknownMap(mgr.mesh(), N);
			result.cell_lid = lid_cell;

			// d 耦合井系统,对系数矩阵
			std::vector<WellConfig> wellsCfg = wellsCfg_in;   // 本地拷贝，便于补齐字段名
			build_masks_and_PI_for_all(mgr, reg, wellsCfg);		//WellConfig.h
			std::vector<WellDOF> wells;
			const int Ntot = register_well_dofs_for_all(Nc, wellsCfg, wells);//WellConfig.h
			GeneralTools::extend_linear_system_size(sys, Ntot);

			for (const auto& w : wells)
			{ 
				add_peaceman_coupling_cell_rows(sys, mesh, reg, w.PI_field, w.mask_field, lid_cell, w.lid);		//#include "FVM_WellCoupling.h"
				add_well_row(sys, mesh, reg, w.PI_field, w.mask_field, lid_cell, w.lid, w.mode, w.target);		//#include "FVM_WellCoupling.h"
			}
			result.system = sys;

			// d 备份上一外迭代层的值并完成求解
			auto pvec = GeneralTools::gatherFieldToVec(reg, mesh, ctrl.pressure_field, lid_cell, N);
			pvec.resize(Ntot, 0.0);
			for (const auto& w : wells)
			{
				if (w.mode == WellDOF::Mode::Pressure) pvec[w.lid] = w.target; // BHP 作为井压初猜
				else                                   pvec[w.lid] = pvec[0];   // 定流：随便给个初猜
			}

			double linRes = 0.0;
			int    linIters = 0;
			if (!solveCOO_Eigen(sys, pvec, sol_ctrl.linear, &linIters, &linRes))
			{
				std::cerr << "[SinglePhase][Pressure] Pressure equation solution failed.\n";
				return false;
			}

			std::vector<double> p_cells(N);
			std::copy(pvec.begin(), pvec.begin() + N, p_cells.begin());
			GeneralTools::scatterVecToField(reg, mesh, ctrl.pressure_field, lid_cell, p_cells);
			writeback_pw_fields_for_all(reg, wells, pvec);   ///#include "WellConfig.h"

			double dpInf = 0.0;
			dpInf = GeneralTools::maxAbsDiff(reg, ctrl.pressure_field, ctrl.pressure_prev_field);
			GeneralTools::underRelaxInPlace(reg, ctrl.pressure_field, ctrl.pressure_prev_field, sol_ctrl.under_relax);

			// e 输出结果
			dp_inf = dpInf;
			lin_residual = linRes;
			lin_iterations = linIters;
			return true;
		}	
	}
}
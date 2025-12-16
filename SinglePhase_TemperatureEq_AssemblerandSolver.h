#pragma once
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"

#include "SolverContrlStrName.h"
#include "Solver_AssemblerCOO.h"
#include "Solver_Tools.h"

#include "TemperatureBCAdapter.h"

#include "DiffusionCentral.h"
#include "Timeterm_BDF.h"
#include "ConvectionUpwind.h"
#include "SinglePhase_FaceMassRateCalculate.h"

#include "LinearSolver_Eigen.h"
#include "WellConfig.h"
//#include "FVM_WellCoupling.h"
#include "FVM_SourceTerm_WellHeat.h"

namespace SinglePhase 
{
	struct TemperatureAssemblyConfig
	{
		PhysicalParameters_String				pmm_str;
		TemperatureEquation_String				T_eq_string;
		std::string operator_tag =				T_eq_string.operator_tag;
		std::string temperature_field =			T_eq_string.temperatue_field;
		std::string temperature_old_field =		T_eq_string.temperatue_old_field;
		std::string temperature_prev_field =	T_eq_string.temperatue_prev_field;
		Vector gravity = { 0.0, 0.0, 0.0 };
		int gradient_smoothing = 0;
	};

	struct TemperatureSolveControls
	{
		TemperatureAssemblyConfig assembly;
		LinearSolverOptions    linear;
		double	under_relax = 1.0;		///< 欠松弛系数 (0<urf<=1, 建议 0.3~0.7)
		double	tol_abs = 1e-6;			///绝对残差
		double	tol_rel = 1e-6;			///相对残差
		int		max_outer = 10;			///最大外迭代次数
		bool	verbose = false;			///< 是否打印每轮外迭代的 dp_inf / linRes
		
	};

	struct TemperatureStepReport
	{
		double lin_residual = 0.0;      ///< 线性求解器返回残差（如果可用）
		int    lin_iterations = 0;      ///< 线性求解器迭代步数
		double dT_inf = 0.0;            ///< 本次 outer 中压力场的无穷范数变化
		int    outer_iterations = 0;   ///< 实际使用的外迭代次数
	};

	struct TemperatureAssemblyResult
	{
		SparseSystemCOO system;
		std::vector<int> cell_lid;
	};

	/**
	*  brief: 该函数用于单相温度方程系数矩阵组装及求解，并输出迭代残差和求解残差及次数
	*/
	inline bool AssembleandSolver_TemperatureEq_SinglePhase
	(
		MeshManager& mgr,
		FieldRegistry& reg,
		FaceFieldRegistry& freg,
		PressureAssemblyConfig& p_ctrl,
		TemperatureAssemblyConfig& ctrl,
		const FaceMassRateConfig mf,
		TemperatureAssemblyResult& result,
		TemperatureSolveControls& sol_ctrl,
		TemperatureStepReport& rep,
		const TemperatureBCAdapter& Tbc,
		double dt,
		const std::vector<WellDOF>& wells_in
	)
	{
		Mesh& mesh = mgr.mesh();
		const int Nc = (int)mesh.getCells().size();
		const OperatorFieldNames nmT = makeNames(ctrl.operator_tag);
		if (!GeneralTools::startOuterIteration_scatter(reg, ctrl.temperature_field, ctrl.temperature_prev_field)) return false;
		if (wells_in.size() == 0)
		{
			// a 扩散项离散
			std::vector<std::string> mobility_tokens = { "iso:lambda_eff" };
			if (!FVM::Diffusion::build_FaceCoeffs_Central(
				mgr, reg, freg,
				nmT.a_f_diff, nmT.s_f_diff,
				ctrl.temperature_field,
				mobility_tokens,
				"", FVM::Diffusion::RhoFaceMethod::Linear,
				ctrl.gravity,
				Tbc,
				false,
				ctrl.gradient_smoothing))
			{
				std::cerr << "[SinglePhase][Temperature] diffusion operator build failed.\n";
				return false;
			}

			// b 时间项离散  其中C_eff 还是没有剥离出来，并且强烈依赖Timeterm_BDF.h
			if (!FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Temperature_revised(mgr, reg, dt, ctrl.pmm_str.C_eff_field, ctrl.pmm_str.C_eff_old_field, ctrl.temperature_old_field, nmT.a_time, nmT.b_time))
			{
				std::cerr << "[SinglePhase][Temperature] diffusion operator build failed.\n";
				return false;
			}

			// c 对流项离散 
			if (!FVM::Convection::build_FaceCoeffs_Upwind(mgr, reg, freg, ctrl.temperature_field, mf.total_mass_flux, { ctrl.pmm_str.cp_fluid_field }, nmT, Tbc))
			{
				std::cerr << "[SinglePhase][Temperature] convection operator build failed.\n";
				return false;
			}

			// d 系数矩阵组装
			SparseSystemCOO sys;
			if (!assemble_COO(mgr, reg, freg, "ddt+diffusion+convection", nmT, &sys))
			{
				std::cerr << "[SinglePhase][Temperature] Assemble process build failed.\n";
				return false;
			}
			result.system = sys;
			int N = 0;
			auto lid_cell = buildUnknownMap(mgr.mesh(), N);
			result.cell_lid = lid_cell;

			// e 上一外迭代层的值作为初值
			auto pvec = GeneralTools::gatherFieldToVec(reg, mesh, ctrl.temperature_field, lid_cell, N);

			double linRes = 0.0;
			int    linIters = 0;
			if (!solveCOO_Eigen(sys, pvec, sol_ctrl.linear, &linIters, &linRes))
			{
				std::cerr << "[SinglePhase][Temperature] Temperature equation solution failed.\n";
				return false;
			}
			GeneralTools::scatterVecToField(reg, mesh, ctrl.temperature_field, lid_cell, pvec);
			double dTInf = 0.0;
			dTInf = GeneralTools::maxAbsDiff(reg, ctrl.temperature_field, ctrl.temperature_prev_field);
			GeneralTools::underRelaxInPlace(reg, ctrl.temperature_field, ctrl.temperature_prev_field, sol_ctrl.under_relax);

			// f 输出结果
			rep.dT_inf = dTInf;
			rep.lin_residual = linRes;
			rep.lin_iterations = linIters;
			return true;
		}
		else
		{
			// a 扩散项离散
			std::vector<std::string> mobility_tokens = { "iso:lambda_eff" };
			if (!FVM::Diffusion::build_FaceCoeffs_Central(
				mgr, reg, freg,
				nmT.a_f_diff, nmT.s_f_diff,
				ctrl.temperature_field,
				mobility_tokens,
				"", FVM::Diffusion::RhoFaceMethod::Linear,
				ctrl.gravity,
				Tbc,
				false,
				ctrl.gradient_smoothing))
			{
				std::cerr << "[SinglePhase][Temperature] diffusion operator build failed.\n";
				return false;
			}

			// b 时间项离散  其中C_eff 还是没有剥离出来，并且强烈依赖Timeterm_BDF.h
			if (!FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Temperature_revised(mgr, reg, dt, ctrl.pmm_str.C_eff_field, ctrl.pmm_str.C_eff_old_field,ctrl.temperature_old_field, nmT.a_time, nmT.b_time))
			{
				std::cerr << "[SinglePhase][Temperature] diffusion operator build failed.\n";
				return false;
			}

			// c 对流项离散 
			if (!FVM::Convection::build_FaceCoeffs_Upwind(mgr, reg, freg, ctrl.temperature_field, mf.total_mass_flux, { ctrl.pmm_str.cp_fluid_field }, nmT, Tbc))
			{
				std::cerr << "[SinglePhase][Temperature] convection operator build failed.\n";
				return false;
			}

			bool first = true;
			for (const auto& w : wells_in)
			{
				const double Tin = (w.role == WellDOF::Role::Injector) ? w.Tin : 0.0;
				if (!FVM::SourceTerm::add_temperature_source_from_qm_single(
					mgr, reg,
					w.role,
					w.PI_field.c_str(), w.mask_field.c_str(), w.name.c_str(),
					Tin,
					/*p*/ p_ctrl.pressure_field.c_str(), /*cp*/ ctrl.pmm_str.cp_fluid_field.c_str(),
					/*thickness*/ 1.0,
					nmT.a_src.c_str(),
					nmT.b_src.c_str(),
					/*accumulate*/ !first,
					/*verbose*/ sol_ctrl.verbose))                                  // <- new argument
					return false;
				first = false;
			}

			// d 系数矩阵组装
			SparseSystemCOO sys;
			if (!assemble_COO(mgr, reg, freg, "ddt+diffusion+convection+source", nmT, &sys))
			{
				std::cerr << "[SinglePhase][Temperature] Assemble process build failed.\n";
				return false;
			}
			result.system = sys;
			int N = 0;
			auto lid_cell = buildUnknownMap(mgr.mesh(), N);
			result.cell_lid = lid_cell;

			// f 备份上一外迭代层的值并完成求解
			auto pvec = GeneralTools::gatherFieldToVec(reg, mesh, ctrl.temperature_field, lid_cell, N);

			double linRes = 0.0;
			int    linIters = 0;
			if (!solveCOO_Eigen(sys, pvec, sol_ctrl.linear, &linIters, &linRes))
			{
				std::cerr << "[SinglePhase][Temperature] Temperature equation solution failed.\n";
				return false;
			}
			GeneralTools::scatterVecToField(reg, mesh, ctrl.temperature_field, lid_cell, pvec);
			double dTInf = 0.0;
			dTInf = GeneralTools::maxAbsDiff(reg, ctrl.temperature_field, ctrl.temperature_prev_field);
			GeneralTools::underRelaxInPlace(reg, ctrl.temperature_field, ctrl.temperature_prev_field, sol_ctrl.under_relax);

			// f 输出结果
			rep.dT_inf = dTInf;
			rep.lin_residual = linRes;
			rep.lin_iterations = linIters;
			return true;
		}
	}

}
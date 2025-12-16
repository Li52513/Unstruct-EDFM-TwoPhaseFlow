#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#if __cplusplus >= 201703L
#include <filesystem>
#endif

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "PhysicalPropertiesManager.h"
#include "BCAdapter.h"
#include "TemperatureBCAdapter.h"

#include "SinglePhase_FaceMassRateCalculate.h"
#include "SinglePhase_PressureEq_AssemblerandSolver.h"
#include "SinglePhase_TemperatureEq_AssemblerandSolver.h"

#include "PostProcess_.h"
#include "PostProcesstoCSV.h"

namespace SinglePhase {


	inline bool runTransient_SinglePhase_HT_Iteration(
		MeshManager& mgr,
		FieldRegistry& reg,
		FaceFieldRegistry& freg,
		PhysicalPropertiesManager& ppm,
		const TemperatureBCAdapter& Tbc,
		const PressureBCAdapter& Pbc,
		const std::vector<WellConfig>& wellsCfg_in,
		PressureSolveControls& P_sol_ctrl,
		FaceMassRateConfig& mf_sol_ctrl,
		TemperatureSolveControls& T_sol_ctrl,
		double dt,
		int nSteps,
		const Vector& g,

		int                        writeEveryP = 0,
		int                        writeEveryT = 0,
		const std::string& outPrefixP = "",
		const std::string& outPrefixT = "",
		int                        snapshotEveryCsv = 0,
		const std::string& snapshotPrefix = "")
	{
#if __cplusplus >= 201703L
		try {
			if (!outPrefixP.empty()) {
				auto dirP = std::filesystem::path(outPrefixP).parent_path();
				if (!dirP.empty()) std::filesystem::create_directories(dirP);
			}
			if (!outPrefixT.empty()) {
				auto dirT = std::filesystem::path(outPrefixT).parent_path();
				if (!dirT.empty()) std::filesystem::create_directories(dirT);
			}
			if (!snapshotPrefix.empty()) {
				auto dirS = std::filesystem::path(snapshotPrefix).parent_path();
				if (!dirS.empty()) std::filesystem::create_directories(dirS);
			}
		}
		catch (...) {
			std::cerr << "[SinglePhase-HT] cannot create output directories.\n";
			return false;
		}
#endif

		if (dt <= 0.0 || nSteps <= 0) {
			std::cerr << "[SinglePhase-HT] invalid dt or nSteps.\n";
			return false;
		}
		auto& mesh = mgr.mesh();
		const size_t nCells = mesh.getCells().size();
		// ---------- string shortcuts (canonical names live in solve-controls) ----------
		const auto& Pasm = P_sol_ctrl.assembly;
		const auto& Tasm = T_sol_ctrl.assembly;
		const std::string p_name = Pasm.pressure_field;
		const std::string p_old_name = Pasm.pressure_old_field;
		const std::string p_prev_name = Pasm.pressure_prev_field;
		const std::string T_name = Tasm.temperature_field;      
		const std::string T_old_name = Tasm.temperature_old_field;
		const std::string T_prev_name = Tasm.temperature_prev_field;
		// fully-implicit time terms require freezing these "old" property fields:
		const std::string rho_name = Pasm.pmm_str.rho_fluid_field;
		const std::string rho_old_name2 = Pasm.pmm_str.rho_fluid_old_field;     

		const std::string Ceff_name = Tasm.pmm_str.C_eff_field;
		const std::string Ceff_old_name2 = Tasm.pmm_str.C_eff_old_field;

		auto makePltName = [](const std::string& prefix, const std::string& tag, int step) -> std::string {
			std::ostringstream oss;
			oss << prefix << tag << "_"
				<< std::setw(6) << std::setfill('0') << step << ".plt";
			return oss.str();
			};

		// ---------- a lightweight property-update wrapper ----------
		auto updateProps_all = [&]() -> bool
		{
				// 你当前用的是“常物性模块”，这里保持接口不变；
				// 后续切换变物性，只替换这三行即可。
				ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
				ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);
				ppm.ComputeEffectiveThermalProperties_constProperties_singlePhase_CO2_T_H(mgr, reg);
				return true;
		};
		// 初始物性更新一次（确保 rho / Ceff 等场已存在）
		if (!updateProps_all()) return false;

		double t = 0.0;

		// ============== time loop ==============
		for (int step = 0; step < nSteps; ++step)
		{
			// ---------- Step 1) start time step: freeze p_old/T_old and set p_prev/T_prev ----------
			if (!GeneralTools::startTimeStep_scalar(mesh, reg, p_name, p_old_name, p_prev_name)) return false;
			if (!GeneralTools::startTimeStep_scalar(mesh, reg, T_name, T_old_name, T_prev_name)) return false;

			// ---------- Step 2) freeze rho_old and Ceff_old (strict fully-implicit needs these) ----------
					// Ensure they exist, then copy current → old (only once per time step!)
			reg.getOrCreate<volScalarField>(rho_old_name2, nCells, 0.0);
			reg.getOrCreate<volScalarField>(Ceff_old_name2, nCells, 0.0);

			// make sure current rho/Ceff are consistent with the just-frozen p_prev/T_prev
			if (!GeneralTools::copyField1(reg, p_prev_name, p_name)) return false;
			if (!GeneralTools::copyField1(reg, T_prev_name, T_name)) return false;
			if (!updateProps_all()) return false;

			if (!GeneralTools::copyField1(reg, rho_name, rho_old_name2)) return false;
			if (!GeneralTools::copyField1(reg, Ceff_name, Ceff_old_name2)) return false;

			// ---------- Step 3) outer iteration inside time step ----------
			bool accept_step_P = false;
			bool accept_step_T = false;
			PressureStepReport    prep{};
			TemperatureStepReport trep{};
			PressureAssemblyResult p_assemblyresult{};
			TemperatureAssemblyResult   T_assemblyresult{};
			std::vector<WellDOF>   wells_output;

			for (int outer = 0; outer < P_sol_ctrl.max_outer; ++outer)
			{
				// --- (a) Pressure Picard step ---
				
				if (!updateProps_all()) return false;

				if (!AssembleandSolver_PressureEq_SinglePhase(
					mgr, reg, freg,
					P_sol_ctrl.assembly,
					p_assemblyresult,
					P_sol_ctrl,
					prep,
					Pbc,
					dt,
					wellsCfg_in,
					wells_output
				))
				{
					std::cerr << "[SinglePhase-HT] pressure solve failed.\n";
					return false;
				}

				if (P_sol_ctrl.verbose)
				{
					std::cout << "[SinglePhase-HT] step " << step
						<< " outer " << outer
						<< " dp_inf=" << prep.dp_inf << "\n";
				}
				if (prep.dp_inf <= P_sol_ctrl.tol_abs)
				{
					if (!updateProps_all()) return false;
					std:: cout << "Pressure Equation convergence at Step " << outer << endl;
					if (!FinalizeFaceCoeffsForFlux(mgr, reg, freg, P_sol_ctrl.assembly, Pbc)) return false;
					accept_step_P = true;
					break;
				}
				GeneralTools::updatePrevIterates(reg, { { Pasm.pressure_field,Pasm.pressure_prev_field } });
				if (outer == P_sol_ctrl.max_outer - 1)
				{
					std::cout << "Reached maxOuter at" << outer << "without meeting P tolerances.\n";
					return false;
				}			
			}
			if (!accept_step_P) return false;

			// ---(b)Face mass rates (must be AFTER pressure convergence) ---
			if (!buildFaceMassRates(mgr, reg, freg, Pasm, Pbc, mf_sol_ctrl))
			{
				std::cerr << "[SinglePhase-HT] buildFaceMassRates failed.\n";
				return false;
			}

			// ---(c)Temperature Equation solve----
			bool temperature_converged = false;
			for (int outer = 0; outer < T_sol_ctrl.max_outer; ++outer)
			{
				if (!updateProps_all()) return false;

				if (!AssembleandSolver_TemperatureEq_SinglePhase(
					mgr, reg, freg,
					P_sol_ctrl.assembly, T_sol_ctrl.assembly,
					mf_sol_ctrl,
					T_assemblyresult,
					T_sol_ctrl, trep,
					Tbc, dt,
					wells_output))
				{
					std::cerr << "[SinglePhase-HT] temperature solve failed.\n";
					return false;
				}
				if (T_sol_ctrl.verbose)
				{
					std::cout << "[SinglePhase-HT] step " << step
						<< " outer " << outer
						<< " dT_inf=" << trep.dT_inf << "\n";
				}
				if (trep.dT_inf <= T_sol_ctrl.tol_abs)
				{
					temperature_converged = true;
					accept_step_T = true;
					cout << "Temperature Equation convergence at Step " << outer << endl;
					break;
				}
				GeneralTools::updatePrevIterates(reg, { { Tasm.temperature_field,Tasm.temperature_prev_field } });
				if (outer == T_sol_ctrl.max_outer-1)
				{
					std::cout << "Reached maxOuter at" << outer << "without meeting T tolerances.\n";
					return false;
				}
			}

			if (accept_step_P && accept_step_T)
			{
				if (writeEveryP > 0 && !outPrefixP.empty() && (step % writeEveryP == 0)) {
					const std::string fn = outPrefixP + std::to_string(step) + ".plt";
					outputTecplot_cellToFaceToNode_BC_P(
						mgr, reg, freg,
						/*Tbc*/ nullptr,
						/*Pbc*/ &Pbc,
						P_sol_ctrl.assembly.pressure_field,
						P_sol_ctrl.assembly.pressure_field + "_face",
						/*grad*/ nullptr,
						fn
					);
				}
				if (writeEveryT > 0 && !outPrefixT.empty() && (step % writeEveryT == 0)) {
					const std::string fn = outPrefixT + std::to_string(step) + ".plt";
					outputTecplot_cellToFaceToNode_BC(
						mgr, reg, freg,
						/*Tbc*/ &Tbc,
						/*Pbc*/ nullptr,
						T_sol_ctrl.assembly.temperature_field,
						T_sol_ctrl.assembly.temperature_field + "_face",
						/*grad*/ nullptr,
						fn
					);
				}
				if (snapshotEveryCsv > 0 && !snapshotPrefix.empty() && (step % snapshotEveryCsv == 0))
				{
					// 1) 单元场：cell_id + center + p + T
					if (!SinglePhase::OutputCSV::writeCellPTSnapshotCSV(
						snapshotPrefix, step, /*simTime*/ t,
						mgr, reg,
						P_sol_ctrl.assembly.pressure_field,
						T_sol_ctrl.assembly.temperature_field))
						return false;

					// 2) 井：pw + T_bh + Qm
					if (!SinglePhase::OutputCSV::writeWellBHP_T_Qm_SnapshotCSV(
						snapshotPrefix, step, /*simTime*/ t,
						mgr, reg,
						wells_output,
						P_sol_ctrl.assembly.pressure_field,
						T_sol_ctrl.assembly.temperature_field))
						return false;
				}
			}
			t += dt;
			std::cout << "[SinglePhase-HT] step " << step << " done. t=" << t << "\n";

		}
		return true;
	}
}
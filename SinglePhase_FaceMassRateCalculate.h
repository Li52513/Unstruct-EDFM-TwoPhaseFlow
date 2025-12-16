#pragma once
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"

#include "SolverContrlStrName.h"
#include "Solver_AssemblerCOO.h"
#include "BCAdapter.h"
#include "ConvectionUpwind_Flux.h"

#include "SinglePhase_PressureEq_AssemblerandSolver.h"

namespace SinglePhase
{
	struct FaceMassRateConfig 
	{
		FaceMassRate_String ctrl;
		std::string total_mass_flux = ctrl.total_mass_flux_name;
		std::string total_vol_flux = ctrl.total_vol_flux_name;						/// 对应的体积通量 Q_total [m^3/s]
		std::string total_velocity = ctrl.total_velocity_name;						/// 对应的法向速度 u_n_total [m/s]
		bool clamp_dirichlet_backflow = true;										/// 是否对 Dirichlet 出流边界做“回流截断”（交给 buildFlux_Darcy_Mass）
		double dirichlet_zero_flux_tol = 1e-10;										/// |p_owner - p_bc| <= tol 时将边界通量视为 0（交给 buildFlux_Darcy_Mass）
		bool enable_conservation_check = true;										/// 是否启用通量守恒检查（内部+边界）
		double flux_check_tol = 1e-10;												/// 用于检查“理论零通量”的容差（如 no-flow 边界、cap/gravity 边界） 若 <=0，则仅打印最大值，不做阈值判断
	};

	//==================check函数===================//
	// /// \brief 对给定面质量通量 mf，做简单的守恒性诊断。
	///
	/// - 对每个 cell 计算 netFlux[cell] = Σ_faces sign * mf_face，
	///   其中 sign = +1(作为 owner) 或 -1(作为 neighbor)；
	/// - 统计 max |netFlux|，打印出来；
	///   * 对稳态、无源、无井、无流边界的算例，理论上应趋于 0；
	///   * 一般非稳态/有井时，这个值 ≈ time-term+源项，不一定为 0，仅作参考。
	/// - 若提供 flux_check_tol > 0，可对“理论零通量”边界做阈值检查。
	
	inline void diagnoseFaceFluxConservation_internal
	(
		MeshManager& mgr,
		const faceScalarField& mf,
		const std::string& label
	)
	{
		const Mesh& mesh = mgr.mesh();
		const auto& faces = mesh.getFaces();
		const auto& cells = mesh.getCells();
		const auto& id2idx = mesh.getCellId2Index();
		std::vector<double> netFlux(cells.size(), 0.0);

		// 汇总：owner 加，neighbor 减
		for (const auto& F : faces)
		{
			const int iF = F.id - 1;
			if (iF < 0 || iF >= static_cast<int>(mf.data.size())) continue;
			const double flux = mf.data[static_cast<size_t>(iF)];
			const int Pid = F.ownerCell;
			const int Nid = F.neighborCell;
			if (Pid >= 0)
			{
				auto itP = id2idx.find(Pid);
				if (itP != id2idx.end())
				{
					netFlux[static_cast<size_t>(itP->second)] += flux;
				}
			}
			if (Nid >= 0)
			{
				auto itN = id2idx.find(Nid);
				if (itN != id2idx.end())
				{
					netFlux[static_cast<size_t>(itN->second)] -= flux;
				}
			}
		}
		double maxAbsNet = 0.0;
		double sumAbsNet = 0.0;
		for (double v : netFlux)
		{
			const double av = std::abs(v);
			maxAbsNet = std::max(maxAbsNet, av);
			sumAbsNet += av;
		}
		std::cout << "[FluxCheck][" << label << "] "
			<< "max |net flux per cell| = " << maxAbsNet
			<< ", sum |net flux per cell| = " << sumAbsNet << "\n";
	}

	/// \brief 检查总压力通量在 no-flow Neumann 边界上的大小。
	inline void diagnoseNoFlowBoundaryFlux
	(
		MeshManager& mgr,
		const faceScalarField& mf,
		const PressureBCAdapter& pbc,
		double a_eps,
		double c_eps,
		double flux_tol,
		const std::string& label
	)
	{
		const Mesh& mesh = mgr.mesh();
		const auto& faces = mesh.getFaces();

		const double aTol = std::max(a_eps, 0.0);
		const double cTol = std::max(c_eps, 0.0);

		double maxAbsFluxNoFlow = 0.0;

		for (const auto& F : faces)
		{
			if (!F.isBoundary()) continue;
			const int iF = F.id - 1;
			if (iF < 0 || iF >= static_cast<int>(mf.data.size())) continue;

			double a = 0.0, b = 0.0, c = 0.0;
			if (!pbc.getABC(F.id, a, b, c)) continue;

			const bool isZeroGrad = (std::abs(a) <= aTol && std::abs(b) > aTol);
			const bool isZeroFlux = (std::abs(c) <= cTol);
			if (!isZeroGrad || !isZeroFlux) continue;

			const double flux = mf.data[static_cast<size_t>(iF)];
			maxAbsFluxNoFlow = std::max(maxAbsFluxNoFlow, std::abs(flux));
		}

		std::cout << "[FluxCheck][" << label << "] no-flow Neumann boundary: "
			<< "max |mf| = " << maxAbsFluxNoFlow << "\n";

		if (flux_tol > 0.0 && maxAbsFluxNoFlow > flux_tol)
		{
			std::cout << "  WARNING: max |mf| on no-flow boundary > flux_check_tol = "
				<< flux_tol << "\n";
		}
	}

	/**
	 * @brief 压力收敛后用于通量计算的“最终面系数刷新”。
	 *
	 * 只做：
	 *  - 重新 build_FaceCoeffs_Central(...) 写入 nm.a_f_diff / nm.s_f_diff
	 * 不做：
	 *  - 不组装矩阵、不求解
	 *
	 * 用途：压力 Picard 收敛后，你通常会再 updateProps_all() 一次，
	 *       这会改变 rho/mu 等场；为了保证通量用的是“收敛 p + 最终物性”，
	 *       就用本函数刷新面系数，然后 buildFaceMassRates() 只读系数算通量。
	 */
	inline bool FinalizeFaceCoeffsForFlux(
		MeshManager& mgr,
		FieldRegistry& reg,
		FaceFieldRegistry& freg,
		const PressureAssemblyConfig& pcfg,
		const PressureBCAdapter& pbc)
	{
		const OperatorFieldNames nm = makeNames(pcfg.operator_tag);

		// 单相 CO2：kxx/kyy/kzz 与 mu_g、rho_g（或你 pmm_str 中定义的 rho_fluid_field）
		std::vector<std::string> mobility_tokens = {
			"kxx:kxx", "kyy:kyy", "kzz:kzz",
			"/mu_g",
			"rho:rho_g"
		};

		if (!FVM::Diffusion::build_FaceCoeffs_Central(
			mgr, reg, freg,
			nm.a_f_diff, nm.s_f_diff,
			pcfg.pressure_field,
			mobility_tokens,
			pcfg.pmm_str.rho_fluid_field,                // 你当前用的是 rho_g
			FVM::Diffusion::RhoFaceMethod::Linear,
			pcfg.gravity,
			pbc,
			pcfg.enable_buoyancy,
			pcfg.gradient_smoothing))
		{
			std::cerr << "[SinglePhase][FinalizeFaceCoeffsForFlux] build_FaceCoeffs_Central failed.\n";
			return false;
		}
		return true;
	}

	//================================================================//

	
	/**
	 * @brief 基于主扩散算子的 a_f / s_f，计算总达西质量通量。
	 *
	 * - 使用 PressureAssemblyConfig::operator_tag 生成 nm；
	 * - 读取 nm.a_f_diff / nm.s_f_diff；
	 * - 以 pressure_field / rho_mix_field 为输入，调用 buildFlux_Darcy_Mass；
	 * - 输出 total_mass_flux (+ Q_total, u_n)。
	 */
	inline bool buildFaceMassRates
	(
		MeshManager& mgr,
		FieldRegistry& reg,
		FaceFieldRegistry& freg,
		const PressureAssemblyConfig& pcfg,
		const PressureBCAdapter& pbc,
		const FaceMassRateConfig& cfg
	)
	{
		if (cfg.total_mass_flux.empty())
		{
			std::cerr << "[FaceMassRate] total_mass_flux name is empty.\n";
			return false;
		}
		const OperatorFieldNames nm = makeNames(pcfg.operator_tag);

		if (!FVM::Convection::buildFlux_Darcy_Mass(
			mgr, reg, freg,
			nm.a_f_diff, nm.s_f_diff,
			pcfg.pressure_field, pcfg.pmm_str.rho_fluid_field,
			cfg.total_mass_flux, cfg.total_vol_flux, cfg.total_velocity,
			&pbc,
			cfg.clamp_dirichlet_backflow,
			cfg.dirichlet_zero_flux_tol))
		{
			std::cerr << "[FaceMassRate] buildFlux_Darcy_Mass failed.\n";
			return false;
		}

		// ======= 守恒性检查（可选） ======= //
		if (cfg.enable_conservation_check)
		{
			auto mf = freg.get<faceScalarField>(cfg.total_mass_flux.c_str());
			if (mf)
			{
				diagnoseFaceFluxConservation_internal(mgr, *mf, "total_mass");
				diagnoseNoFlowBoundaryFlux(
					mgr, *mf, pbc,
					1e-30,   // 这里你可以把 no-flow 判据从别处挪进来
					1e-30,
					cfg.flux_check_tol,
					"total_mass");
			}
		}
		return true;
	}
}
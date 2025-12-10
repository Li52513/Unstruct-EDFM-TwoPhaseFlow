#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <memory>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "FaceSignMask.hpp"
#include "BCAdapter.h"
#include "IMPES_CommonUtils.h"

#include "SolverContrlStrName.h"
#include "ConvectionUpwind_Flux.h"
#include "PressureEqAssemblerandSolver.h"
namespace IMPES_Iteration
{
	struct FaceMassRateConfig
	{
		FaceMassRate_String ctrl;
		std::string total_mass_flux = ctrl.total_mass_flux_name;					// 总压力驱动质量通量 m_total [kg/s]
		std::string total_vol_flux = ctrl.total_vol_flux_name;						/// 对应的体积通量 Q_total [m^3/s]
		std::string total_velocity = ctrl.total_velocity_name;						/// 对应的法向速度 u_n_total [m/s]
		std::string capillary_mass_flux = ctrl.capillary_correction_flux_name;		/// 毛细校正质量通量 m_cap [kg/s]
		std::string gravity_mass_flux = ctrl.gravity_correction_flux_name;			/// 重力校正质量通量 m_grav [kg/s]
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
	/// \brief 检查“理论应为零”的边界通量（用于 cap/gravity 通量）。
	inline void diagnoseZeroBoundaryFlux
	(
		MeshManager& mgr,
		const faceScalarField& mf,
		double flux_tol,
		const std::string& label
	)
	{
		const Mesh& mesh = mgr.mesh();
		const auto& faces = mesh.getFaces();

		double maxAbsFluxBnd = 0.0;

		for (const auto& F : faces)
		{
			if (!F.isBoundary()) continue;
			const int iF = F.id - 1;
			if (iF < 0 || iF >= static_cast<int>(mf.data.size())) continue;

			const double flux = mf.data[static_cast<size_t>(iF)];
			maxAbsFluxBnd = std::max(maxAbsFluxBnd, std::abs(flux));
		}

		std::cout << "[FluxCheck][" << label << "] boundary faces: "
			<< "max |mf| = " << maxAbsFluxBnd << "\n";

		if (flux_tol > 0.0 && maxAbsFluxBnd > flux_tol)
		{
			std::cout << "  WARNING: max |mf| on boundary > flux_check_tol = "
				<< flux_tol << "\n";
		}
	}

	//===============================================//
	/**
	 * @brief 基于主扩散算子的 a_f / s_f，计算总达西质量通量。
	 *
	 * - 使用 PressureAssemblyConfig::operator_tag 生成 nm；
	 * - 读取 nm.a_f_diff / nm.s_f_diff；
	 * - 以 pressure_field / rho_mix_field 为输入，调用 buildFlux_Darcy_Mass；
	 * - 输出 total_mass_flux (+ Q_total, u_n)。
	 */
	inline bool buildTotalDarcyMassFlux(
		MeshManager& mgr,
		const FieldRegistry& reg,
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
		
		const std::string Q_name = cfg.total_vol_flux.empty()
			? (cfg.total_mass_flux + "_Q")
			: cfg.total_vol_flux;
		
		const std::string ufn_name = cfg.total_velocity.empty()
			? (cfg.total_mass_flux + "_ufn")
			: cfg.total_velocity;
		bool ok = FVM::Convection::buildFlux_Darcy_Mass(
			mgr, reg, freg,
			nm.a_f_diff,
			nm.s_f_diff,
			pcfg.pressure_field,
			pcfg.rho_mix_field,
			cfg.total_mass_flux,
			Q_name,
			ufn_name,
			&pbc,
			cfg.clamp_dirichlet_backflow,
			cfg.dirichlet_zero_flux_tol);

		if (!ok) {
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
	/**
	 * @brief 基于毛细算子的 a_f / s_f 和 Pc 场，计算毛细质量通量。
	 *
	 * - 内部面：m_cap = a_f_cap (Pc_P - Pc_N) - s_f_cap
	 * - 边界面：m_cap = 0（因为 buildAndApplyOperator 时使用 nullPbc）
	 */
	inline bool buildCapillaryMassFlux(
		MeshManager& mgr,
		const FieldRegistry& reg,
		FaceFieldRegistry& freg,
		const PressureAssemblyConfig& pcfg,
		const FaceMassRateConfig& cfg)
	{
		if (cfg.capillary_mass_flux.empty()) 
		{
			// 不需要毛细质量通量，直接跳过
			return true;
		}
		const OperatorFieldNames nmCap = makeNames(pcfg.operator_tag + "_cap");
		auto aF = freg.get<faceScalarField>(nmCap.a_f_diff.c_str());
		auto sF = freg.get<faceScalarField>(nmCap.s_f_diff.c_str());
		if (!aF || !sF) 
		{
			std::cerr << "[FaceMassRate] capillary operator face fields not found: "
				<< nmCap.a_f_diff << ", " << nmCap.s_f_diff << "\n";
			return false;
		}
		auto mCap = freg.getOrCreate<faceScalarField>(
			cfg.capillary_mass_flux.c_str(), aF->data.size(), 0.0);
		
		const Mesh& mesh = mgr.mesh();
		const auto& faces = mesh.getFaces();
		const auto& id2idx = mesh.getCellId2Index();

		for (const auto& F : faces) {
			const int iF = F.id - 1;
			double flux = 0.0;

			if (!F.isBoundary() && F.neighborCell >= 0)
			{
				const double PcP = cellScalar(reg, mesh, pcfg.Pc_field.c_str(), F.ownerCell, 0.0);
				const double PcN = cellScalar(reg, mesh, pcfg.Pc_field.c_str(), F.neighborCell, 0.0);
				flux = (*aF)[iF] * (PcP - PcN) - (*sF)[iF];
			}
			else
			{
				// 与 assemblePressureTwoPhase 中 nullPbc 的用法保持一致，
				// 边界毛细通量统一视为 0。
				flux = 0.0;
			}
			(*mCap)[iF] = flux;
		}
		if (cfg.enable_conservation_check)
		{
			diagnoseFaceFluxConservation_internal(mgr, *mCap, "capillary_mass");
			diagnoseZeroBoundaryFlux(mgr, *mCap, cfg.flux_check_tol, "capillary_mass");
		}
		return true;
	}
	/**
	 * @brief 基于重力算子的 s_f，计算重力质量通量。
	 *
	 * 理论上：m_grav = - s_f_grav（owner -> neighbor 为正方向）。
	 * 由于装配时采用 nullPbc，重力项仅作为内部体力校正，这里对边界直接设为 0。
	 */
	inline bool buildGravityMassFlux(
		MeshManager& mgr,
		FaceFieldRegistry& freg,
		const PressureAssemblyConfig& pcfg,
		const FaceMassRateConfig& cfg)
	{
		if (cfg.gravity_mass_flux.empty()) {
			return true; // 不需要重力通量
		}

		const OperatorFieldNames nmGrav = makeNames(pcfg.operator_tag + "_grav");
		auto aF = freg.get<faceScalarField>(nmGrav.a_f_diff.c_str());
		auto sF = freg.get<faceScalarField>(nmGrav.s_f_diff.c_str());
		if (!aF || !sF) {
			std::cerr << "[FaceMassRate] gravity operator face fields not found: "
				<< nmGrav.a_f_diff << ", " << nmGrav.s_f_diff << "\n";
			return false;
		}

		auto mGrav = freg.getOrCreate<faceScalarField>(
			cfg.gravity_mass_flux.c_str(), aF->data.size(), 0.0);

		const Mesh& mesh = mgr.mesh();
		const auto& faces = mesh.getFaces();

		for (const auto& F : faces) {
			const int iF = F.id - 1;
			double flux = 0.0;

			if (!F.isBoundary() && F.neighborCell >= 0)
			{
				// 与 buildAndApplyOperator / assemble_diffusion_faces 的符号保持一致：
				// 对应 cell 方程中，b += ± s_f；对面通量而言就是 F_grav = - s_f。
				flux = -(*sF)[iF];
			}
			else
			{
				// 重力项仅作为内部体力校正，边界不施加额外重力通量
				flux = 0.0;
			}
			(*mGrav)[iF] = flux;
		}
		if (cfg.enable_conservation_check)
		{
			diagnoseFaceFluxConservation_internal(mgr, *mGrav, "gravity_mass");
			diagnoseZeroBoundaryFlux(mgr, *mGrav, cfg.flux_check_tol, "gravity_mass");
		}
		return true;
	}
	/**
	 * @brief 在压力收敛之后构建全部面质量通量：
	 *        总达西质量通量 + 毛细质量通量 + 重力质量通量。
	 *
	 * 使用前置条件：
	 *  - assemblePressureTwoPhase 已经被调用一次，且：
	 *    * 主扩散算子的 a_f/s_f:  nm.a_f_diff, nm.s_f_diff 已存在；
	 *    * 毛细算子的   a_f/s_f:  nmCap.a_f_diff, nmCap.s_f_diff 已存在；
	 *    * 重力算子的   a_f/s_f:  nmGrav.a_f_diff, nmGrav.s_f_diff 已存在；
	 *  - 压力场 / Pc 场 / rho_mix 场已更新为当前时间步/外迭代的解。
	 */
	inline bool buildFaceMassRates(
		MeshManager& mgr,
		FieldRegistry& reg,
		FaceFieldRegistry& freg,
		const PressureAssemblyConfig& pcfg,
		const PressureBCAdapter& pbc,
		const FaceMassRateConfig& cfg
	)
	{
		// 1) 主压力驱动通量（走 buildFlux_Darcy_Mass）
		if (!buildTotalDarcyMassFlux(mgr, reg, freg, pcfg, pbc, cfg)) {
			std::cerr << "[FaceMassRate] buildTotalDarcyMassFlux failed.\n";
			return false;
		}

		// 2) 毛细质量通量
		if (!buildCapillaryMassFlux(mgr, reg, freg, pcfg, cfg)) {
			std::cerr << "[FaceMassRate] buildCapillaryMassFlux failed.\n";
			return false;
		}

		// 3) 重力质量通量
		if (!buildGravityMassFlux(mgr, freg, pcfg, cfg)) {
			std::cerr << "[FaceMassRate] buildGravityMassFlux failed.\n";
			return false;
		}

		return true;
	}

}// namespace IMPES_Iteration
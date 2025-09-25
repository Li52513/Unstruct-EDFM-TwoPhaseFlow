#pragma once
#include <iostream>
#include <iomanip>
#include <string>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "FieldAcessForDiscre.h"               // cellScalar(...)
#include "TPFA_PermeabilityOperation.h"        // kEffAlong/getKdiag（若 Mobility 里用到）
#include "TPFA_GradientsOperation.h"           // computeCellGradients_LSQ_with_GG(...)
#include "TPFA_Mobility.h"                     // DarcyWaterMobility_singlePhase 等
#include "TPFA_UpwindforGravityandDensity.h"   // UpwindDensityByPotential_* / GravUpwind
#include "BCAdapter.h"   

// 小工具 ：统一几何（内部/边界）

inline void faceGeomPack( Mesh& mesh, const Face& F, double& dperp, Vector& ehat, double& Eabs, double& Aabs) {
	const Vector Aj = F.vectorE + F.vectorT;
	if (!F.isBoundary())
	{
		dperp = std::max(F.ownerToNeighbor.Mag(), 1e-14);
		ehat = (dperp > 1e-14 ? F.ownerToNeighbor / dperp : F.normal);
		
	}
	else
	{
		const auto& id2idx = mesh.getCellId2Index();
		const Vector r = F.midpoint - mesh.getCells()[id2idx.at(F.ownerCell)].center;
		dperp = std::max(r.Mag(), 1e-14);
		ehat = (dperp > 1e-14) ? (r / dperp) : F.normal;
	}
	Eabs = F.vectorE.Mag();
	Aabs = Aj.Mag(); // = length(2D) or area(3D)

}

// ============ 内部面打印/校核 ============
// MobilityProvider: 需实现 mobilityAlong(mesh, reg, cellId, ehat)
// DensityPolicy   : 需实现 rhoUp(...) 与 rhoBar(...)

template<class MobilityProvider, class DensityPolicy>
inline void debugPrintInnerFaces_TPFA(MeshManager& mgr,
	const FieldRegistry& reg,
	const FaceFieldRegistry& freg,
	const GravUpwind& gu,
	const MobilityProvider& mob,
	const DensityPolicy& rhoPol,
	const std::string& a_name = "a_f_Diff",
	const std::string& s_name = "s_f_Diff",
	const std::string& x_name = "p_w",
	bool enable_buoy = true,
	int  gradSmoothIters = 0,
	int  max_to_print = -1)
{
	Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
	const auto& faces = mesh.getFaces();
	const auto& cells = mesh.getCells();
	const auto& id2idx = mesh.getCellId2Index();

	auto a_f_Diff = freg.get<faceScalarField>(a_name);
	auto s_f_Diff = freg.get<faceScalarField>(s_name);
	if (!a_f_Diff || !s_f_Diff) {
		std::cout << "[DebugInner] 缺少面场 '" << a_name << "' 或 '" << s_name << "'\n";
		return;
	}

	// 预计算每个单元的 ∇p（LSQ 优先，GG 兜底，可选平滑）,用于交叉项计算
	const std::vector<Vector> grad = computeCellGradients_LSQ_with_GG(mesh, reg, x_name.c_str(), gradSmoothIters);

	int printed = 0;
	std::cout << std::scientific << std::setprecision(6);
	for (const auto& F : faces)
	{
		if (F.isBoundary()) continue;

		const int P = F.ownerCell, N = F.neighborCell;
		const size_t iF = F.id - 1;
		const size_t iP = id2idx.at(P), iN = id2idx.at(N);

		double dperp, Eabs, Aabs; Vector ehat;
		faceGeomPack(mesh, F, dperp, ehat, Eabs, Aabs);
		// gamma 兜底
		double gamma = F.f_linearInterpolationCoef;
		if (!(gamma > 0.0 && gamma < 1.0)) {
			const Vector& CP = cells[iP].center;
			const Vector& CN = cells[iN].center;
			const double dPf = (F.midpoint - CP).Mag();
			const double dFn = (CN - F.midpoint).Mag();
			gamma = dPf / std::max(dPf + dFn, 1e-14);
		}

		// λ 调和到面
		const double lamP = std::max(mob.mobilityAlong(mesh, reg, P, ehat), 1e-30);
		const double lamN = std::max(mob.mobilityAlong(mesh, reg, N, ehat), 1e-30);
		const double lam_f = 1.0 / std::max(gamma / lamP + (1.0 - gamma) / lamN, 1e-30);

		// 变量与密度
		const double xP = cellScalar(reg, mesh, x_name.c_str(), P, 0.0);
		const double xN = cellScalar(reg, mesh, x_name.c_str(), N, 0.0);
		const Vector CPc = cells[iP].center, CNc = cells[iN].center;

		const double rho_bar = rhoPol.rhoBar(mesh, reg, P, N,gamma);
		const double rho_up = rhoPol.rhoUp(mesh, reg, P, N, xP, xN, CPc, CNc, gu);

		const double beta_f = std::max(rho_up * lam_f, 0.0);
		const double a_calc = beta_f * (Eabs / dperp);

		// 源项组成
		const Vector gF = grad[iP] * (1.0 - gamma) + grad[iN] * gamma;
		const double s_cross = beta_f * (gF * F.vectorT);
		double s_buoy = 0.0;
		if (enable_buoy) {
			s_buoy = -beta_f * rho_bar * (gu.g * ehat) * Eabs;
		}
		const double s_calc = s_cross + s_buoy;

		const double a_store = (*a_f_Diff)[iF];
		const double s_store = (*s_f_Diff)[iF];

		std::cout << "[IFace " << F.id << "] "
			<< "P=" << P << " N=" << N
			<< "  a_store=" << a_store << " a_calc=" << a_calc
			<< "  |Δa|=" << std::abs(a_store - a_calc) << "\n"
			<< "           s_store=" << s_store
			<< "  s_cross=" << s_cross
			<< "  s_buoy=" << s_buoy
			<< "  s_calc=" << s_calc
			<< "  |Δs|=" << std::abs(s_store - s_calc) << "\n";

		if (max_to_print > 0 && ++printed >= max_to_print) break;


	}
}

// ============ 边界面打印/校核 ============
// 需要 BCProvider（PressureBCAdapter 等）来拿原始 a,b,c
template<class MobilityProvider, class DensityPolicy, class BCProvider>
inline void debugPrintBoundaryFaces_TPFA(MeshManager& mgr,
	const FieldRegistry& reg,
	const FaceFieldRegistry& freg,
	const GravUpwind& gu,
	const MobilityProvider& mob,
	const DensityPolicy& rhoPol,
	const BCProvider& bc,
	const std::string& a_name = "a_f_Diff",
	const std::string& s_name = "s_f_Diff",
	const std::string& x_name = "p_w",
	bool enable_buoy = true,
	int  gradSmoothIters = 0,
	int  max_to_print = -1)
{
	Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
	const auto& faces = mesh.getFaces();
	const auto& cells = mesh.getCells();
	const auto& id2idx = mesh.getCellId2Index();

	auto a_f_Diff = freg.get<faceScalarField>(a_name);
	auto s_f_Diff = freg.get<faceScalarField>(s_name);
	if (!a_f_Diff || !s_f_Diff) {
		std::cout << "[DebugBFace] 缺少面场 '" << a_name << "' 或 '" << s_name << "'\n";
		return;
	}

	// 交叉项使用 cell 梯度（注意：边界 T=0 => 交叉应为 0，但保留同一套写法）
	const auto grad = computeCellGradients_LSQ_with_GG(mesh, reg, x_name.c_str(), gradSmoothIters);

	int printed = 0;
	std::cout << std::scientific << std::setprecision(6);
	for (const auto& F : faces) {
		if (!F.isBoundary()) continue;

		const int P = F.ownerCell;
		const size_t iF = F.id - 1;
		const size_t iP = id2idx.at(P);

		double dperp, Eabs, Aabs; Vector ehat;
		faceGeomPack(mesh, F, dperp, ehat, Eabs, Aabs);

		// 面移动率（取 P）
		const double lamP = std::max(mob.mobilityAlong(mesh, reg, P, ehat), 1e-30);
		const double xP = cellScalar(reg, mesh, x_name.c_str(), P, 0.0);
		const Vector CPc = cells[iP].center;

		const double rho_up = rhoPol.rhoUp(mesh, reg, P, P, xP, xP, CPc, CPc, gu);
		const double beta_f = std::max(rho_up * lamP, 0.0);
		const double a_face = beta_f * (Eabs / dperp); // 修正前的面系数

		// Robin 原始 a,b,c
		double a = 0.0, b = 0.0, c = 0.0;
		const bool has = bc.getABC(F.id, a, b, c);

		double alpha_f = 0.0, betaB_f = 0.0, a_PB = a_face;
		if (has) {
			const double num = b * (Eabs / dperp);
			const double den = a * std::max(Aabs, 1e-30) + num;
			alpha_f = (den > 1e-30) ? (num / den) : 0.0;

			if (std::abs(a) > 1e-30) betaB_f = (c * Aabs) / den;
			else                      betaB_f = 0.0;      // 纯 Neumann

			a_PB = a_face * (1.0 - alpha_f);
		}

		// 源项：交叉 + 浮力 + BC
		const double s_cross = beta_f * (grad[iP] * F.vectorT);        // T=0 => 应该是 0
		double s_buoy = 0.0;
		if (enable_buoy) {
			const double rhoP = rhoPol.rhoUp(mesh, reg, P, P, xP, xP, CPc, CPc, gu);
			s_buoy = -beta_f * rhoP * (gu.g * ehat) * Eabs;
		}
		const double s_BC = a_face * betaB_f;
		const double a_calc = a_PB;
		const double s_calc = s_cross + s_buoy + s_BC;

		const double a_store = (*a_f_Diff)[iF];
		const double s_store = (*s_f_Diff)[iF];

		std::cout << "[BFace " << F.id << "] P=" << P
			<< "  a_store=" << a_store << " a_calc=" << a_calc
			<< "  |Δa|=" << std::abs(a_store - a_calc) << "\n"
			<< "           s_store=" << s_store
			<< "  s_cross=" << s_cross
			<< "  s_buoy=" << s_buoy
			<< "  s_BC=" << s_BC
			<< "  s_calc=" << s_calc
			<< "  |Δs|=" << std::abs(s_store - s_calc) << "\n";

		if (max_to_print > 0 && ++printed >= max_to_print) break;
	}
}

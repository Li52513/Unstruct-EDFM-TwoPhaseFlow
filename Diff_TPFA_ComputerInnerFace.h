#pragma once
#include <string> // std::string
#include <memory> // std::shared_ptr
#include <tuple> // std::tie
#include <utility> // std::forward
#include <algorithm> 
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "FieldAcessForDiscre.h"
#include "Diff_TPFA_PermeabilityOperation.h"
#include "Diff_TPFA_GradientsOperation.h"
#include "Diff_TPFA_UpwindforGravityandDensity.h"
#include "Diff_TPFA_Mobility.h"

// ───────────────────────────────
//  内部面：TPFA（含可选“达西/纯扩散”两类）
//    号约定与现有实现保持一致：
//    β_f > 0； s_cross = +β_f (∇φ)_f·T_f； s_buoy = −β_f ρ̄ (g·ê) |E|-β_f ρg·T
//  并将计算得到的网格面离散系数和源项，存入面场 a_f_Diff, s_f_Diff
// ───────────────────────────────
template<class MobilityProvider, class DensityPolicy>
inline void Diffusion_TPFA_InnerFace_SinglePhase
(
	MeshManager& mgr,
	const FieldRegistry& reg,
	FaceFieldRegistry& freg,
	const GravUpwind& gu,
	const MobilityProvider& mob,
	const DensityPolicy& rhoPol,
	/*存放字段名*/ const std::string& a_name = "a_f_Diff",
	const std::string& s_name = "s_f_Diff",
	const std::string& x_name = "p",
	/*是否计浮力*/ bool enable_buoy = true,
	/*梯度平滑*/ int gradSmoothIters = 0
)
{
	Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
	auto faces = const_cast<std::vector<Face>&>(mesh.getFaces());
	const auto& id2idx = mesh.getCellId2Index();

	//输出字段（面场）
	auto a_f_Diff = freg.getOrCreate<faceScalarField>(a_name, faces.size(), 0.0);
	auto s_f_Diff = freg.getOrCreate<faceScalarField>(s_name, faces.size(), 0.0);

	// 预计算每个单元的 ∇p（LSQ 优先，GG 兜底，可选平滑）,用于交叉项计算
	const std::vector<Vector> grad = computeCellGradients_LSQ_with_GG(mesh, reg, x_name.c_str(), gradSmoothIters);

	const double eps_d = 1e-14, eps_mu = 1e-30, eps_l = 1e-30;
#pragma omp parallel for schedule(static)
	for (int f = 0; f < static_cast<int>(faces.size()); ++f)
	{
		const Face& F = faces[f];
		if (F.isBoundary()) { (*a_f_Diff)[F.id - 1] = 0; (*s_f_Diff)[F.id - 1] = 0.0; continue; }

		const int P = F.ownerCell;
		const int N = F.neighborCell;

		const double dperp = std::max(F.ownerToNeighbor.Mag(), eps_d);
		const Vector ehat = (dperp > 1e-14 ? F.ownerToNeighbor / dperp : F.normal);
		const Vector Aj = F.vectorE + F.vectorT;
		const double Eabs = F.vectorE.Mag();

		//插值系数gamma
		double gamma = F.f_linearInterpolationCoef;

		if (!(gamma > 0.0 && gamma < 1.0))
		{
			const Vector CO = mesh.getCells()[P].center;
			const Vector CN = mesh.getCells()[N].center;
			const Vector dON = CN - CO;
			const double D = dON.Mag();
			if (D > eps_d)
			{
				const Vector eON = dON / D;
				const double s = (F.midpoint - CO) * eON; // 面心到 owner 的投影长度
				gamma = std::min(1.0, std::max(0.0, s / D));
			}
			else
			{
				gamma = 0.5;
			}
		}

		// λ（沿 e 的等效）：调和到面
		const double lamP = std::max(mob.mobilityAlong(mesh, reg, P, ehat), eps_l);  //涉及到的物性参数为黏度mu
		const double lamN = std::max(mob.mobilityAlong(mesh, reg, N, ehat), eps_l);
		const double lam_f = 1.0 / std::max(gamma / lamP + (1.0 - gamma) / lamN, eps_l);

		//用于达西类引入的扩散项，需要计算rho_up和rho_inter（来源于浮力）；对于纯扩散类，rho_up=1，rho_inter=0
		double p_p = cellScalar(reg, mesh, x_name.c_str(), P, 0.0);
		double p_n = cellScalar(reg, mesh, x_name.c_str(), N, 0.0);

		const Vector& CP = mesh.getCells()[id2idx.at(P)].center;
		const Vector& CN = mesh.getCells()[id2idx.at(N)].center;

		const double rho_inter = rhoPol.rhoBar(mesh, reg, P, N, gamma);
		const double rho_up = rhoPol.rhoUp(mesh, reg, P, N, p_p, p_n, CP, CN, gu);  //涉及到的物性参数为密度rho

		// β_f = λ_e * |E| / d_perp
		const double beta_f = std::max(rho_up * lam_f, 0.0);
		const double a_face = beta_f * Eabs / dperp;
		(*a_f_Diff)[F.id - 1] = a_face;

		// s_cross = +β_f (∇φ)_f·T_f
		const Vector& gP = grad[id2idx.at(P)];
		const Vector& gN = grad[id2idx.at(N)];
		const Vector grad_f = (1.0 - gamma) * gP + gamma * gN; //线性插值

		const double s_cross = beta_f * (grad_f * F.vectorT);
		double s_buoy = 0.0;

		//如果包含浮力项
		if (enable_buoy)
		{
			// s_buoy = −β_f ρ̄ (g·ê) |E|-β_f ρg·T
			const double gdot = gu.g * ehat; // g·ê
			auto s_buoy_e = -beta_f * rho_inter * gdot * Eabs;
			auto s_buoy_t = -beta_f * rho_inter * gu.g * F.vectorT;
			s_buoy = s_buoy_e + s_buoy_t;
		}
		const double s_face_raw = s_cross + s_buoy;
		const double s_face_lim =
			std::max(-0.5 * a_face, std::min(0.5 * a_face, s_face_raw));
		(*s_f_Diff)[F.id - 1] = s_face_lim;

	}

}
#include <tuple> // std::tie
#include <utility> // std::forward
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "BCAdapter.h"
#include "FieldAcessForDiscre.h"
#include "Diff_TPFA_PermeabilityOperation.h"
#include "Diff_TPFA_GradientsOperation.h"
#include "Diff_TPFA_UpwindforGravityandDensity.h"
#include "Diff_TPFA_Mobility.h"

// ───────────────────────────────
//  边界面：TPFA + 通用 ABC（面积化）
// a_PB = a_f (1 − α_f)，其中
//   a_f   = β_f |E| / d⊥ （β_f=ρ_up λ_P）
//   α_f   = (b |E|/d⊥) / ( a |A| + b |E|/d⊥ )
// s_BC   = + a_f * βB_f
// 交叉/浮力与内部约定：s_cross = β_f (∇φ_P)·T_f；s_buoy = −β_f ρ̄ (g·ê) |E|-β_f ρg·T 
// 
// 并将计算得到的网格面离散系数和源项，存入面场 a_f_Diff, s_f_Diff
// ───────────────────────────────

template<class MobilityProvider, class DensityPolicy, class BCProvider>
inline void Diffusion_TPFA_BoundaryFace_SinglePhase(
    MeshManager& mgr,
    const FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const GravUpwind& gu,
    const MobilityProvider& mob,
    const DensityPolicy& rhoPol,   // 纯扩散可用 NoDensity
    const BCProvider& bc,
    const std::string& a_name = "a_f_Diff",
    const std::string& s_name = "s_f_Diff",
    const std::string& x_name = "p",
    bool enable_buoy = true,
    int  gradSmoothIters = 0,
    bool rho_in_matrix = false      // <<< 新增：与内部面策略一致，默认不把ρ放进矩阵
)
{
    Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
    auto faces = const_cast<std::vector<Face>&>(mesh.getFaces());
    const auto& id2idx = mesh.getCellId2Index();

    auto a_f_Diff = freg.getOrCreate<faceScalarField>(a_name, faces.size(), 0.0);
    auto s_f_Diff = freg.getOrCreate<faceScalarField>(s_name, faces.size(), 0.0);

    const std::vector<Vector> grad =
        computeCellGradients_LSQ_with_GG(mesh, reg, x_name.c_str(), gradSmoothIters);

    const double eps_d = 1e-14, eps_l = 1e-30;

    //取消并行计算
    for (size_t f = 0; f < faces.size(); ++f)
    {
        const Face& F = faces[f];
        if (!F.isBoundary()) continue;

        const int P = F.ownerCell;
        const size_t iP = id2idx.at(P);

        const double dperp = std::max(F.ownerToNeighbor.Mag(), eps_d);
        const Vector ehat = (dperp > 1e-14 ? F.ownerToNeighbor / dperp : F.normal);
        const double Eabs = F.vectorE.Mag();
        const double Aabs = (F.vectorE + F.vectorT).Mag();

        // λ_P（沿 ehat 的等效导通率）
        const double lamP = std::max(mob.mobilityAlong(mesh, reg, P, ehat), eps_l);

        // 迎风密度（边界：单边，取 P）
        const double pP = cellScalar(reg, mesh, x_name.c_str(), P, 0.0);
        const Vector CPc = mesh.getCells()[iP].center;
        const double rho_up = rhoPol.rhoUp(mesh, reg, P, /*N*/P, pP, pP, CPc, CPc, gu);

        // ——关键分离：矩阵用 beta_face，源项用 beta_src——
        const double beta_face = (rho_in_matrix ? rho_up * lamP : lamP);
        const double beta_src = (rho_up * lamP);

        // a_f：只用 beta_face
        const double a_face = beta_face * (Eabs / dperp);

        // Robin 统一（面积化）
        double a = 0.0, b = 0.0, c = 0.0;
        const bool has = bc.getABC(F.id, a, b, c);

        double alpha_f = 0.0, betaB_f = 0.0;
        double a_PB = a_face; // 默认 Dirichlet(b=0) => alpha=0

        if (has)
        {
            const double num = b * (Eabs / dperp);
            const double den = a * std::max(Aabs, 1e-30) + num;

            alpha_f = (den > 1e-30) ? (num / den) : 0.0;
            if (std::abs(a) > 1e-30) {
                const double C = c * Aabs;
                betaB_f = C / den;
            }
            else {
                // 纯 Neumann：没有“等效壁值”，只通过 s_f 注入；这里 betaB_f=0
                betaB_f = 0.0;
            }

            a_PB = a_face * (1.0 - alpha_f);
        }

        // ——源项：交叉 + 浮力 + BC——
        // 交叉项用单元梯度（边界：只用P侧）
        const Vector& gP = grad[iP];
        const double s_cross = beta_src * (gP * F.vectorT);

        double s_buoy = 0.0;
        if (enable_buoy)
        {
            // 用 beta_src × ρ 与 g 的投影；NoDensity 时 rho_up=1，自然退化
            const double g_dot_e = gu.g * ehat;
            s_buoy = -beta_src * g_dot_e * Eabs
                + -beta_src * (gu.g * F.vectorT);
        }

        const double s_BC = a_face * betaB_f;

        (*a_f_Diff)[F.id - 1] = a_PB;
        const double s_aux = s_cross + s_buoy;
        const double s_aux_limited = std::max(-0.5 * a_face, std::min(0.5 * a_face, s_aux));
        (*s_f_Diff)[F.id - 1] = s_aux_limited + s_BC;
        //(*s_f_Diff)[F.id - 1] = s_cross + s_buoy + s_BC;
    }

//#pragma omp parallel for schedule(static)
//    for (int f = 0; f < static_cast<int>(faces.size()); ++f)
//    {
//        const Face& F = faces[f];
//        if (!F.isBoundary()) continue;
//
//        const int P = F.ownerCell;
//        const size_t iP = id2idx.at(P);
//
//        const double dperp = std::max(F.ownerToNeighbor.Mag(), eps_d);
//        const Vector ehat = (dperp > 1e-14 ? F.ownerToNeighbor / dperp : F.normal);
//        const double Eabs = F.vectorE.Mag();
//        const double Aabs = (F.vectorE + F.vectorT).Mag();
//
//        // λ_P（沿 ehat 的等效导通率）
//        const double lamP = std::max(mob.mobilityAlong(mesh, reg, P, ehat), eps_l);
//
//        // 迎风密度（边界：单边，取 P）
//        const double pP = cellScalar(reg, mesh, x_name.c_str(), P, 0.0);
//        const Vector CPc = mesh.getCells()[iP].center;
//        const double rho_up = rhoPol.rhoUp(mesh, reg, P, /*N*/P, pP, pP, CPc, CPc, gu);
//
//        // ——关键分离：矩阵用 beta_face，源项用 beta_src——
//        const double beta_face = (rho_in_matrix ? rho_up * lamP : lamP);
//        const double beta_src = (rho_up * lamP);
//
//        // a_f：只用 beta_face
//        const double a_face = beta_face * (Eabs / dperp);
//
//        // Robin 统一（面积化）
//        double a = 0.0, b = 0.0, c = 0.0;
//        const bool has = bc.getABC(F.id, a, b, c);
//
//        double alpha_f = 0.0, betaB_f = 0.0;
//        double a_PB = a_face; // 默认 Dirichlet(b=0) => alpha=0
//
//        if (has)
//        {
//            const double num = b * (Eabs / dperp);
//            const double den = a * std::max(Aabs, 1e-30) + num;
//
//            alpha_f = (den > 1e-30) ? (num / den) : 0.0;
//            if (std::abs(a) > 1e-30) {
//                const double C = c * Aabs;
//                betaB_f = C / den;
//            }
//            else {
//                // 纯 Neumann：没有“等效壁值”，只通过 s_f 注入；这里 betaB_f=0
//                betaB_f = 0.0;
//            }
//
//            a_PB = a_face * (1.0 - alpha_f);
//        }
//
//        // ——源项：交叉 + 浮力 + BC——
//        // 交叉项用单元梯度（边界：只用P侧）
//        const Vector& gP = grad[iP];
//        const double s_cross = beta_src * (gP * F.vectorT);
//
//        double s_buoy = 0.0;
//        if (enable_buoy)
//        {
//            // 用 beta_src × ρ 与 g 的投影；NoDensity 时 rho_up=1，自然退化
//            const double g_dot_e = gu.g * ehat;
//            s_buoy = -beta_src * g_dot_e * Eabs
//                + -beta_src * (gu.g * F.vectorT);
//        }
//
//        const double s_BC = a_face * betaB_f;
//
//        (*a_f_Diff)[F.id - 1] = a_PB;
//        const double s_aux = s_cross + s_buoy;
//        const double s_aux_limited = std::max(-0.5 * a_face, std::min(0.5 * a_face, s_aux));
//        (*s_f_Diff)[F.id - 1] = s_aux_limited + s_BC;
//        //(*s_f_Diff)[F.id - 1] = s_cross + s_buoy + s_BC;
//    }
}



//template<class MobilityProvider, class DensityPolicy, class BCProvider>
//inline void Diffusion_TPFA_BoundaryFace_SinglePhase(MeshManager& mgr,
//	const FieldRegistry& reg, 
//	FaceFieldRegistry& freg,
//	const GravUpwind& gu,
//	const MobilityProvider& mob,
//	const DensityPolicy& rhoPol,   // 纯扩散可用 NoDensity（其 ρ_P=1，ρ̄无关）
//	const BCProvider& bc,
//	const std::string& a_name = "a_f_Diff",
//	const std::string& s_name = "s_f_Diff",
//	const std::string& x_name = "p",
//	bool enable_buoy = true,
//	int gradSmoothIters = 0)
//{
//
//	//获取网格和网格面拓扑信息
//	Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
//	auto faces = const_cast<std::vector<Face>&>(mesh.getFaces());
//	const auto& id2idx = mesh.getCellId2Index();
//
//	//建立（获得）存储离散系数和源项的面场
//	auto a_f_Diff = freg.getOrCreate<faceScalarField>(a_name, faces.size(), 0.0);
//	auto s_f_Diff = freg.getOrCreate<faceScalarField>(s_name, faces.size(), 0.0);
//
//	// 预计算每个单元的 ∇p（LSQ 优先，GG 兜底，可选平滑）,用于交叉项计算
//	const std::vector<Vector> grad = computeCellGradients_LSQ_with_GG(mesh, reg, x_name.c_str(), gradSmoothIters);
//
//	const double eps_d = 1e-14, eps_mu = 1e-30, eps_l = 1e-30;
//
//#pragma omp parallel for schedule(static)
//	for (int f = 0; f < static_cast<int>(faces.size()); ++f)
//	{
//		const Face& F = faces[f];
//		if (!F.isBoundary()) continue;
//		const int P = F.ownerCell;
//		const size_t iP = id2idx.at(P);
//
//		const double dperp = std::max(F.ownerToNeighbor.Mag(), eps_d);
//		const Vector ehat = (dperp > 1e-14 ? F.ownerToNeighbor / dperp : F.normal);
//		const Vector Aj = F.vectorE + F.vectorT;
//		const double Eabs = F.vectorE.Mag();
//		const double Aabs = Aj.Mag(); // = length(2D) or area(3D)
//
//		//对于边界系数，直接取E，插值系数gamma 因为T=0，所以不影响，不需要计算 s_cross =0
//		const double lamP = std::max(mob.mobilityAlong(mesh, reg, P, ehat), eps_l);
//		const double lam_f = lamP; //边界面只用 owner 单元的渗透率、黏度
//		// ρ_up：单边，取 P（或达西策略下按定义返回 ρ_P）
//		const double pP = cellScalar(reg, mesh, x_name.c_str(), P, 0.0);
//		const Vector CPc = mesh.getCells()[iP].center;
//		// N 无，下面给一个影子参数
//		const double rho_up = rhoPol.rhoUp(mesh, reg, P, /*N*/P, pP, pP, CPc, CPc, gu);  // 单边取 P
//		const double beta_f = std::max(rho_up * lamP, 0.0);
//		const double a_face = beta_f * (Eabs / dperp); // a_f
//
//		// Robin 统一（面积化）：
//		double a = 0.0, b = 0.0, c = 0.0;
//		bool has = bc.getABC(F.id, a, b, c); // a,b,c 已是 ∮(...)dA 的系数
//		double alpha_f = 0.0, betaB_f = 0.0, a_PB = a_face; // 默认 Dirichlet(b=0)=>alpha=0
//
//		if (has)
//		{
//			// α_f
//			const double num = b * (Eabs / dperp);
//			const double den = a * std::max(Aabs, 1e-30) + num;
//			alpha_f = (den > 1e-30) ? (num / den) : 0.0;
//
//			// βB_f
//			const double C = c * Aabs;
//
//			if (std::abs(a) > 1e-30) betaB_f = C / den;
//			else betaB_f = 0.0; // pure Neumann
//
//			a_PB = a_face * (1.0 - alpha_f); // a_PB
//		}
//
//		//源项=交叉+浮力+BC
//		const Vector& gP = grad[iP];
//		double s_cross = beta_f * (gP * F.vectorT); // s_cross
//		double s_buoy = 0.0;
//
//		if (enable_buoy)
//		{
//			// 单边取 ρ_P（达西类）；纯扩散则 ρ_P=1，s_buoy 不会被使用
//			const double rhoP = rhoPol.rhoUp(mesh, reg, P, P, pP, pP, CPc, CPc, gu);
//			const double g_dot_e = gu.g * ehat; // g·ê
//			s_buoy = -beta_f * rhoP * g_dot_e * Eabs + -beta_f * rhoP * gu.g * F.vectorT; // s_buoy
//		}
//
//		const double s_BC = a_face * betaB_f; // s_BC
//		(*a_f_Diff)[F.id - 1] = a_PB;
//		(*s_f_Diff)[F.id - 1] = s_cross + s_buoy + s_BC;
//	}
//
//}
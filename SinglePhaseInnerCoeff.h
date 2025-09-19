#pragma once
#include <algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"

//===========小工具==================//
//开关：重力是否用迎风格式，计算重力势能
struct GravUpwind {
	Vector g;
	bool use_potential = true;
	inline double dot_pos(const Vector& p)const { return g.m_x * p.m_x + g.m_y * p.m_y + g.m_z * p.m_z; }
};

//读取网格单元中的标量场
inline double cellScalar(const FieldRegistry& reg, const Mesh& mesh, const char* name, int cellId, double fallback = 0.0)
{
	auto fld = reg.get<volScalarField>(name);
	if (!fld)
	{
		return fallback;
		cout << "field of " << name << "cannot get" << endl;
	}
	auto it = mesh.getCellId2Index().find(cellId);
	if (it == mesh.getCellId2Index().end())
	{
		return fallback;
		cout << "Id of "<< cellId<< "cannot get" << endl;

	}
	return fld->data[it->second];

}


//取各向异性K的三个主值（若未写入则用 k_default 兜底）
inline void getKdiag(const FieldRegistry& reg, const Mesh& mesh,int cellId, double k_default,double& kxx, double& kyy, double& kzz)
{
	kxx = cellScalar(reg, mesh, "kxx", cellId, k_default);
	kyy = cellScalar(reg, mesh, "kyy", cellId, k_default);
	kzz = cellScalar(reg, mesh, "kzz", cellId, k_default);

}

//================= 工具 1：Green–Gauss（GG）梯度 =================//
// Green–Gauss：∇p ≈ (1/V) Σ pf * A_out
// - 仅依赖 Face.vectorE / vectorT（方法无关），退化时用 normal*length
// - A_out 逐面对“当前单元 CP”外指
// - 内部面 pf 线性插值；边界面 pf = p_owner（先不考虑边界BC）
inline Vector greenGaussGrad_P(Mesh& mesh,
	const FieldRegistry& reg,
	const Cell& CP)
{
	const auto& faces = mesh.getFaces();
	const auto& cells = mesh.getCells();
	const auto& id2idx = mesh.getCellId2Index();
	auto p = reg.get<volScalarField>("p_w");
	const double eps = 1e-14;

	auto itP = id2idx.find(CP.id);
	if (!p || itP == id2idx.end() || CP.volume <= 0.0) return Vector(0, 0, 0);
	const size_t iP = itP->second;
	const double pP = (*p)[iP];

	Vector sumA(0, 0, 0);   // 自检：封闭曲面应≈0
	Vector flux(0, 0, 0);   // Σ pf * A_out

	for (int fid_ext : CP.CellFaceIDs) {
		const int idx = fid_ext - 1;                // 你当前拓扑是 1 基
		if (idx < 0 || idx >= (int)faces.size()) continue;
		const Face& F = faces[idx];

		// —— 方法无关的面积矢量 —— //
		Vector AjAligned = F.vectorE + F.vectorT;
		if (AjAligned.Mag() < eps) {                // 极端退化兜底
			Vector Aj0 = F.normal * F.length;
			if (F.neighborCell >= 0) {
				const Vector d = cells[id2idx.at(F.neighborCell)].center
					- cells[id2idx.at(F.ownerCell)].center;
				if ((Aj0 * d) < 0.0) Aj0 = -Aj0;   // 与 owner->neighbor 同向
			}
			AjAligned = Aj0;
		}
		// 对“当前单元 CP”外指
		const Vector Aout = (F.ownerCell == CP.id ? AjAligned : (AjAligned * -1.0));

		// —— 面上标量 pf —— //
		double pf = pP; // 默认先用 owner = CP
		if (!F.isBoundary()) {
			const size_t iO = id2idx.at(F.ownerCell);
			const size_t iN = id2idx.at(F.neighborCell);
			const double pO = (*p)[iO], pN = (*p)[iN];

			// gamma 若未预置，用“owner→neighbor 的距离投影”兜底
			double gamma = F.f_linearInterpolationCoef;
			if (!(gamma > 0.0 && gamma < 1.0)) {
				const Vector CO = cells[iO].center;
				const Vector CN = cells[iN].center;
				Vector dON = CN - CO;
				const double D = dON.Mag();
				if (D > eps) {
					const Vector ehat = dON / D;
					const double s = (F.midpoint - CO) * ehat;
					gamma = std::min(1.0, std::max(0.0, s / std::max(D, eps)));
				}
				else gamma = 0.5;
			}
			pf = (1.0 - gamma) * pO + gamma * pN;  // 内部面线性插值
		}
		// 边界面：先不考虑 BC，pf=pP 即可（不会破坏常量场）

		sumA = sumA + Aout;
		flux = flux + Aout * pf;
	}

	// 自检（可选）：封闭曲面 ΣA_out≈0
	// if (sumA.Mag() / std::max(1.0, CP.volume) > 1e-10) { ... }

	return flux * (1.0 / CP.volume);
}




//============= 工具 2：加权最小二乘（LSQ）梯度，失败则回退GG ===========//
struct Mat3 {
	double a[3][3];
	void zero() { for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) a[i][j] = 0.0; }
	Mat3(){ zero(); }
	void addOuter(const Vector& r, double w) {
		const double x = r.m_x, y = r.m_y, z = r.m_z;
		a[0][0] += w * x * x; a[0][1] += w * x * y; a[0][2] += w * x * z;
		a[1][0] += w * y * x; a[1][1] += w * y * y; a[1][2] += w * y * z;
		a[2][0] += w * z * x; a[2][1] += w * z * y; a[2][2] += w * z * z;
	}
	//求解Ag=b(3x3):返回是否成功
	bool solve(const Vector& b, Vector& g)const {
		const double A00 = a[0][0], A01 = a[0][1], A02 = a[0][2];
		const double A10 = a[1][0], A11 = a[1][1], A12 = a[1][2];
		const double A20 = a[2][0], A21 = a[2][1], A22 = a[2][2];

		// 行列式
		const double det =
			A00 * (A11 * A22 - A12 * A21) -
			A01 * (A10 * A22 - A12 * A20) +
			A02 * (A10 * A21 - A11 * A20);
		if (std::abs(det) < 1e-18) return false;

		// 伴随矩阵 * b / det
		const double i00 = (A11 * A22 - A12 * A21) / det;
		const double i01 = -(A01 * A22 - A02 * A21) / det;
		const double i02 = (A01 * A12 - A02 * A11) / det;
		const double i10 = -(A10 * A22 - A12 * A20) / det;
		const double i11 = (A00 * A22 - A02 * A20) / det;
		const double i12 = -(A00 * A12 - A02 * A10) / det;
		const double i20 = (A10 * A21 - A11 * A20) / det;
		const double i21 = -(A00 * A21 - A01 * A20) / det;
		const double i22 = (A00 * A11 - A01 * A10) / det;
		g.m_x = i00 * b.m_x + i01 * b.m_y + i02 * b.m_z;
		g.m_y = i10 * b.m_x + i11 * b.m_y + i12 * b.m_z;
		g.m_z = i20 * b.m_x + i21 * b.m_y + i22 * b.m_z;
		return true;
	}

};

inline std::vector<Vector>
computeCellGradients_LSQ_with_GG(Mesh& mesh, const FieldRegistry& reg, int smoothIters = 0)
{
	const auto& cells = mesh.getCells();
	const auto& faces = mesh.getFaces();
	const auto& id2idx = mesh.getCellId2Index();
	auto p = reg.get<volScalarField>("p_w");
	std::vector<Vector> grad(cells.size(), Vector(0, 0, 0));
	if (!p) return grad;

	const double eps_r = 1e-12;

	// 判断维度，决定邻居最少数（2D≥2，3D≥3）
	bool is3D = false;
	for (auto& F : faces) { if (F.FaceNodeCoords.size() == 3) { is3D = true; break; } }
	const int minNb = is3D ? 3 : 2;

	struct Mat3 
	{
		double a[3][3];
		Mat3() { for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) a[i][j] = 0.0; }
		void addOuter(const Vector& r, double w)
		{
			const double x = r.m_x, y = r.m_y, z = r.m_z;
			a[0][0] += w * x * x; a[0][1] += w * x * y; a[0][2] += w * x * z;
			a[1][0] += w * y * x; a[1][1] += w * y * y; a[1][2] += w * y * z;
			a[2][0] += w * z * x; a[2][1] += w * z * y; a[2][2] += w * z * z;
		}
		bool solve(const Vector& b, Vector& g) const {
			const double A00 = a[0][0], A01 = a[0][1], A02 = a[0][2];
			const double A10 = a[1][0], A11 = a[1][1], A12 = a[1][2];
			const double A20 = a[2][0], A21 = a[2][1], A22 = a[2][2];
			const double det = A00 * (A11 * A22 - A12 * A21) - A01 * (A10 * A22 - A12 * A20) + A02 * (A10 * A21 - A11 * A20);
			if (std::abs(det) < 1e-18) return false;
			const double i00 = (A11 * A22 - A12 * A21) / det;
			const double i01 = -(A01 * A22 - A02 * A21) / det;
			const double i02 = (A01 * A12 - A02 * A11) / det;
			const double i10 = -(A10 * A22 - A12 * A20) / det;
			const double i11 = (A00 * A22 - A02 * A20) / det;
			const double i12 = -(A00 * A12 - A02 * A10) / det;
			const double i20 = (A10 * A21 - A11 * A20) / det;
			const double i21 = -(A00 * A21 - A01 * A20) / det;
			const double i22 = (A00 * A11 - A01 * A10) / det;
			g.m_x = i00 * b.m_x + i01 * b.m_y + i02 * b.m_z;
			g.m_y = i10 * b.m_x + i11 * b.m_y + i12 * b.m_z;
			g.m_z = i20 * b.m_x + i21 * b.m_y + i22 * b.m_z;
			return true;
		}
	};

	// 第一次：LSQ，失败则 GG
	for (const auto& CP : cells) 
	{
		auto itP = id2idx.find(CP.id);
		if (itP == id2idx.end()) continue;
		const size_t iP = itP->second;
		const double phiP = (*p)[iP];

		Mat3 G; Vector b(0, 0, 0);
		int nbCount = 0;

		for (int fid_ext : CP.CellFaceIDs) {
			const int idx = fid_ext - 1;                          // 1基→0基
			if (idx < 0 || idx >= (int)faces.size()) continue;
			const Face& F = faces[idx];

			const int N = (F.ownerCell == CP.id ? F.neighborCell : F.ownerCell);
			if (N < 0) continue;                                  // 只用内部邻居做 LSQ
			const size_t iN = id2idx.at(N);
			const Vector r = cells[iN].center - CP.center;
			const double rmag = std::max(r.Mag(), eps_r);
			const double w = 1.0 / (rmag * rmag);
			const double dphi = (*p)[iN] - phiP;

			G.addOuter(r, w);
			b = b + r * (w * dphi);
			++nbCount;
		}

		bool ok = false;
		if (nbCount >= minNb) 
		{              // 2D≥2, 3D≥3
			// 轻微对角正则，避免近奇异
			G.a[0][0] += 1e-18; G.a[1][1] += 1e-18; G.a[2][2] += 1e-18;
			ok = G.solve(b, grad[iP]);
			cout << "对单元" << CP.id << "采用了LSQ算法计算梯度" << endl;
		}
		if (!ok) {
			grad[iP] = greenGaussGrad_P(mesh, reg, CP);  // GG 兜底
			cout << "对单元" << CP.id << "采用了GG算法计算梯度" << endl;;
		}
	}

	// 可选：邻域平滑（提升鲁棒性）
	for (int it = 0; it < smoothIters; ++it)
	{
		std::vector<Vector> newg = grad;
		for (const auto& CP : cells) {
			auto itP = id2idx.find(CP.id);
			if (itP == id2idx.end()) continue;
			const size_t iP = itP->second;
			Vector acc = grad[iP]; int cnt = 1;
			for (int fid_ext : CP.CellFaceIDs) {
				const int idx = fid_ext - 1;
				if (idx < 0 || idx >= (int)faces.size()) continue;
				const Face& F = faces[idx];
				const int N = (F.ownerCell == CP.id ? F.neighborCell : F.ownerCell);
				if (N < 0) continue;
				const size_t iN = id2idx.at(N);
				acc = acc + grad[iN]; ++cnt;
			}
			newg[iP] = acc * (1.0 / double(cnt));
		}
		grad.swap(newg);
	}
	return grad;
}


// 将各相渗透率等效到e方向上 k_e = e^T K e（K=diag）
inline double kEffAlong(const Vector& e, double kxx, double kyy, double kzz)
{
	// 假设 Vector 有成员 m_x/m_y/m_z
	const double ex = e.m_x, ey = e.m_y, ez = e.m_z;
	return kxx * ex * ex + kyy * ey * ey + kzz * ez * ez;
}


//=============内部面a_f(faceDiscreCoef)和离散源项的计算===========================///

inline void computerInnerFaceCoeffAndSources_PW(MeshManager& mgr, const FieldRegistry& reg, const GravUpwind& gu, double k_default, bool enable_buoy, int gradSmoothIters = 0)
{
	Mesh& mesh = const_cast<Mesh&>(mgr.mesh());  //获取网格mesh
	auto& faces = const_cast<std::vector<Face>&>(mesh.getFaces()); //获取网格面face
	const auto& id2idx = mesh.getCellId2Index();

	//根据当前reg中存储的压力场 预计算每个单元的∇p（LSQ 优先，GG 兜底）
	const std::vector<Vector> gradP = computeCellGradients_LSQ_with_GG(mesh, reg, gradSmoothIters);

	const double eps_mu = 1e-30;
	const double eps_l = 1e-30;
	const double eps_d = 1e-14;

	for (auto& F : faces)
	{
		if (F.isBoundary()) continue;

		const int P = F.ownerCell;
		const int N = F.neighborCell;

		// --- 几何：d⊥ 与 ê ，以及 |E|
		const double dperp = std::max(F.ownerToNeighbor.Mag(), eps_d);
		const Vector ehat = (dperp > 1e-14 ? F.ownerToNeighbor / dperp : F.normal);
		const double Eabs = F.vectorE.Mag();


		//--- gamma：若未赋值，按“到面心距离”兜底
		double gamma = F.f_linearInterpolationCoef;
		if (!(gamma > 0.0 && gamma < 1.0))
		{
			const Vector& CP = mesh.getCells()[id2idx.at(P)].center;
			const Vector& CN = mesh.getCells()[id2idx.at(N)].center;
			const double dPf = (F.midpoint - CP).Mag();
			const double dfn = (CN - F.midpoint).Mag();
			const double denom = std::max(dPf + dfn, eps_d);
			gamma = dPf / denom;
		}

		// ---两侧物性
		const double pP = cellScalar(reg, mesh, "p_w", P);
		cout << "第" << F.id << "个面Owner单元的水压力为" << pP << endl;
		const double pN = cellScalar(reg, mesh, "p_w", N);
		cout << "第" << F.id << "个面Neighbor单元的水压力为" << pN << endl;
		const double rhoP = cellScalar(reg, mesh, "rho_w", P);
		cout << "第" << F.id << "个面Owner单元的水的密度为" << rhoP << endl;
		const double rhoN = cellScalar(reg, mesh, "rho_w", N);
		cout << "第" << F.id << "个面Neighbor单元的水的密度为" << rhoN << endl;
		const double muP = std::max(cellScalar(reg, mesh, "mu_w", P), eps_mu);
		cout << "第" << F.id << "个面Owner单元的水的黏度为" << muP << endl;
		const double muN = std::max(cellScalar(reg, mesh, "mu_w", N), eps_mu);
		cout << "第" << F.id << "个面Neighbor单元的水的黏度为" << muN << endl;

		double kxxP, kyyP, kzzP, kxxN, kyyN, kzzN;
		getKdiag(reg, mesh, P, k_default, kxxP, kyyP, kzzP);  //取出P网格内各个方向的渗透率
		getKdiag(reg, mesh, N, k_default, kxxN, kyyN, kzzN); //取出N网格内各个方向的渗透率

		const double k_eP = std::max(kEffAlong(ehat, kxxP, kyyP, kzzP), eps_l); //计算等效渗透率
		cout << "第" << F.id << "个面Owner单元的等效渗透系数为" << k_eP << endl;
		const double k_eN = std::max(kEffAlong(ehat, kxxN, kyyN, kzzN), eps_l);
		cout << "第" << F.id << "个面Neighbor单元的等效渗透系数为" << k_eN << endl;
		const double lamP = std::max(k_eP / muP, eps_l); //计算p网格上的流度
		cout << "第" << F.id << "个面Owner单元的流度lamda为" << lamP << endl;
		const double lamN = std::max(k_eN / muN, eps_l); //计算N网格上的流度
		cout << "第" << F.id << "个面Neighbor单元的流度lamda为" << lamN << endl;

		// --- 面上流度（调和平均）
		const double lam_f = 1.0 / std::max(gamma / lamP + (1.0 - gamma) / lamN, eps_l);

		// --- 上游密度（相势差或压差）
		double rho_up;
		if (gu.use_potential) {
			const Vector CPc = mesh.getCells()[id2idx.at(P)].center;
			const Vector CNc = mesh.getCells()[id2idx.at(N)].center;
			const double rhoBar = 0.5 * (rhoP + rhoN);  //计算重力势能的时候密度采用平均
			const double phiP = pP - rhoBar * gu.dot_pos(CPc);
			const double phiN = pN - rhoBar * gu.dot_pos(CNc);
			rho_up = (phiP - phiN > 0.0 ? rhoP : rhoN); // P->N 流则取 P 侧
		}
		else {
			rho_up = (pP - pN > 0.0 ? rhoP : rhoN);
		}
		const double beta_f = rho_up * lam_f;

		// ======内部面系数计算=============//
		F.faceDiscreCoef = beta_f * (Eabs / dperp);
		cout << "第" << F.id << "个面的离散系数为：" << F.faceDiscreCoef << endl;

		// ======源项的计算： 交叉项 + 浮力项 =======//

		//首先要计算面上梯度
		const Vector& gP = gradP[id2idx.at(P)];
		cout << "第" << F.id << "个面Owner单元的压力梯度为" << "("<<gP.m_x<<","<< gP.m_y << "," << gP.m_z << ")" << endl;
		const Vector& gN = gradP[id2idx.at(N)];
		cout << "第" << F.id << "个面Neighbor单元的压力梯度为" << "(" << gN.m_x << "," << gN.m_y << "," << gN.m_z << ")" << endl;
		const Vector gF = gP * (1 - gamma) + gN * gamma;
		cout << "第" << F.id << "个面的压力梯度为" << "(" << gF.m_x << "," << gF.m_y << "," << gF.m_z << ")" << endl;

		//计算交叉项
		double s_cross = 0.0;
		// 非正交交叉项：β_f * (∇p)_f · T_f
		s_cross = beta_f * (gF * F.vectorT);
		cout << "第" << F.id << "个面的非正交交叉源项为：" << s_cross << endl;

		//计算 浮力项：β_f * ρ̄ * (g·ê) * |E|  浮力项只取正交
		double s_buoy = 0.0;
		if (enable_buoy)
		{
			const double rhoBar = 0.5 * (rhoP + rhoN);
			const double g_dot_e = gu.g.m_x * ehat.m_x + gu.g.m_y * ehat.m_y + gu.g.m_z * ehat.m_z;
			s_buoy = -beta_f * rhoBar * g_dot_e * Eabs;
			cout << "第" << F.id << "个面的浮力源项为：" << s_buoy << endl;
		}
		F.faceSource = s_cross + s_buoy;
		cout << "第" << F.id << "个面的离散源项为：" << F.faceSource << endl;

	}

}


// ========================= 边界面 a_f与源项计算 ===================================//
//约定：
// beta_f = rho_up * lambda_f  (>0, 单边取 owner 值)
// a_f = beta_f  * |E| / d_perp
// alpha_f = (b * |E| / d_perp) / (a * |A| + b * |E| / d_perp)
//  betaB_f = (c * |A|)         / (a * |A| + b * |E| / d_perp)  这里由于T=0
//  a_PB = a_f * (1 - alpha_f)
//   s_cross = + beta_f * (gradP_P · T)                 // 与内部面同号约定
//   s_buoy  = - beta_f * rhoBar * (g·e_hat) * |E|      // 与内部面同号约定
//   s_BC    = + a_f * betaB_f   我认为这里应该是 ；s_BC    =  a_f * betaB_f
//   F.faceDiscreCoef = a_PB
//   F.faceSource     = s_cross + s_buoy + s_BC
//
// 说明：PressureBC::Registry 需能通过 F.id 查询到 {a,b,c}（面积积分意义下）。
//       若该面未设置 BC，则退化为“Neumann=0”（即 a=0,b=1,c=0）。

inline void computeBoundaryFaceCoeffAndSources_PW(
	MeshManager& mgr,
	const FieldRegistry& reg,
	const GravUpwind& gu,
	double k_default,
	const PressureBC::Registry& pbc,     // 按面积积分形式的 a,b,c
	bool enable_buoy,
	int  gradSmoothIters = 0
){
	Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
	auto& faces = const_cast<std::vector<Face>&>(mesh.getFaces());
	const auto& cells = mesh.getCells();
	const auto& id2idx = mesh.getCellId2Index();

	// 预计算每个单元的 ∇p（与内部面一致：LSQ 优先，GG 兜底，可选平滑）
	const std::vector<Vector> gradP = computeCellGradients_LSQ_with_GG(mesh, reg, gradSmoothIters);  //一定会用GG 因为边界面上不符合LSQ的计算逻辑

	const double eps_mu = 1e-30;
	const double eps_l = 1e-30;
	const double eps_d = 1e-14;
	const double eps_bc = 1e-30;

	for (auto& F : faces)
	{
		if (!F.isBoundary()) continue;            // 只处理边界面
		const int P = F.ownerCell;
		auto itP = id2idx.find(P);
		if (itP == id2idx.end()) continue;
		const size_t iP = itP->second;  //找到全局编号

		// ----计算几何量---

		//d_perp 取 P 到面心的距离；Eabs=|E|；Aabs=|A|（方法无关：E+T）
		const Vector CPc = cells[iP].center;
		const double dperp = std::max((F.midpoint - CPc).Mag(), eps_d);

		const Vector Aj = F.vectorE + F.vectorT;     // 方法无关面积矢量
		const double Aabs = std::max(Aj.Mag(), eps_d); 
		const double Eabs = std::max(F.vectorE.Mag(), eps_d);

		//e_hat：正交方向，优先用 E 的方向，退化用法向
		Vector ehat = (Eabs > 1e-14) ? (F.vectorE * (1.0 / Eabs)) : F.normal;

		// ---物性（对于边界面取Owner单元的值）
		const double pP = cellScalar(reg, mesh, "p_w", P);
		const double rhoP = cellScalar(reg, mesh, "rho_w", P);
		const double muP = std::max(cellScalar(reg, mesh, "mu_w", P), eps_mu);

		double kxxP, kyyP, kzzP;
		getKdiag(reg, mesh, P, k_default, kxxP, kyyP, kzzP);
		const double k_eP = std::max(kEffAlong(ehat, kxxP, kyyP, kzzP), eps_l);

		const double lamP = std::max(k_eP / muP, eps_l);
		const double rho_up = rhoP;                      // 边界面：用 owner 密度
		const double beta_f = rho_up * lamP;             // 单边流度 * 上游密度 > 0

		// --- a_f（不含 alpha_f） ---
		const double a_f = beta_f * (Eabs / dperp);

		// --- 读取面积积分形式的边界系数 (a,b,c) ---
		double a = 0.0, b = 1.0, c = 0.0;               // 默认 Neumann=0 的退化（稳健）
		if (auto bc = pbc.find(F.id))
		{
			a = bc->a; b = bc->b; c = bc->c;            // 注意：这里的 a,b,c 已是“按面积积分”的系数
		}
		// denom = a*|A| + b*|E|/d_perp
		const double denom = std::max(a * Aabs + b * (Eabs / dperp), eps_bc);
		double alpha_f = (b * (Eabs / dperp)) / denom;   // α_f ∈ [0,1]（数值上做个夹取更稳）
		if (alpha_f < 0.0) alpha_f = 0.0;
		if (alpha_f > 1.0) alpha_f = 1.0;
		
		cout << "第" << F.id << "面上的alpha等于" << alpha_f << endl;

		const double betaB_f = (c * Aabs) / denom;

		cout << "第" << F.id << "面上的betaB_f等于" << betaB_f << endl;

		//------计算面对角系数（装配到P的主对角）——
		const double a_PB = a_f * (1.0 - alpha_f);
		F.faceDiscreCoef = a_PB;
		cout << "第" << F.id << "个面的离散系数为：" << F.faceDiscreCoef << endl;

		// --- 源项三块：交叉 + 浮力 + BC 常数 ---
		// 交叉项（边界面常用 gP，若 T=0 则此项自动为 0）
		const Vector& gP = gradP[iP];
		double s_cross = beta_f * (gP * F.vectorT);

		// 浮力项（与内部面一致的号约定）
		double s_buoy = 0.0;
		if (enable_buoy)
		{
			const double rhoBar = rhoP;  // 单边，取 P
			const double g_dot_e = gu.g.m_x * ehat.m_x + gu.g.m_y * ehat.m_y + gu.g.m_z * ehat.m_z;
			s_buoy = -beta_f * rhoBar * g_dot_e * Eabs;    // 负号：与内部面保持一致
		}

		// BC 常数项
		const double s_BC = a_f * betaB_f;

		F.faceSource = s_cross + s_buoy + s_BC;


		// ---- 调试输出（可保留/可屏蔽） ----
		std::cout << "[BFace " << F.id << "] "
			<< "pP=" << pP
			<< "  rhoP=" << rhoP
			<< "  muP=" << muP
			<< "  k_eP=" << k_eP
			<< "  lamP=" << lamP
			<< "  beta_f=" << beta_f
			<< "  Eabs=" << Eabs
			<< "  dperp=" << dperp
			<< "  a_f=" << a_f
			<< "  (a,b,c)=(" << a << "," << b << "," << c << ")"
			<< "  alpha_f=" << alpha_f
			<< "  betaB_f=" << betaB_f
			<< "  a_PB=" << a_PB
			<< "  s_cross=" << s_cross
			<< "  s_buoy=" << s_buoy
			<< "  s_BC=" << s_BC
			<< "  faceSource=" << F.faceSource
			<< std::endl;



	}


}
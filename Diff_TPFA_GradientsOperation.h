#pragma once
#include<vector>
#include<string>
#include<algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FieldAcessForDiscre.h"

// ======工具：网格单元上梯度计算：LSQ + GG 兜底（声明即可；用你已有定义） ======
//================= 工具 1：Green–Gauss（GG）梯度 =================//
// Green–Gauss：∇p ≈ (1/V) Σ pf * A_out
// - 仅依赖 Face.vectorE / vectorT（方法无关），退化时用 normal*length
// - A_out 逐面对“当前单元 CP”外指
// - 内部面 pf 线性插值；边界面 pf = p_owner（先不考虑边界BC）
inline Vector greenGaussGrad(Mesh& mesh, const FieldRegistry& reg, const char* name, const Cell& CP)
{
	const auto& faces = mesh.getFaces();
	const auto& cells = mesh.getCells();
	const auto& id2idx = mesh.getCellId2Index();
	auto phi = reg.get<volScalarField>(name);
	const double eps = 1e-14;

	auto itP = id2idx.find(CP.id);
	if (!phi || itP == id2idx.end() || CP.volume <= 0.0) return Vector(0, 0, 0);
	const size_t iP = itP->second;
	const double phiP = (*phi)[iP];

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
		double pf = phiP; // 默认先用 owner = CP
		if (!F.isBoundary()) {
			const size_t iO = id2idx.at(F.ownerCell);
			const size_t iN = id2idx.at(F.neighborCell);
			const double pO = (*phi)[iO], pN = (*phi)[iN];

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
struct Mat3
{
	double a[3][3];
	void zero() { for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) a[i][j] = 0.0; }
	Mat3() { zero(); }
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
computeCellGradients_LSQ_with_GG(Mesh& mesh, const FieldRegistry& reg, const char* name, int smoothIters = 0)
{
	const auto& cells = mesh.getCells();
	const auto& faces = mesh.getFaces();
	const auto& id2idx = mesh.getCellId2Index();
	auto phi = reg.get<volScalarField>(name);
	std::vector<Vector> grad(cells.size(), Vector(0, 0, 0));
	if (!phi) return grad;

	const double eps_r = 1e-12;

	// 判断维度，决定邻居最少数（2D≥2，3D≥3）
	bool is3D = false;
	for (auto& F : faces) { if (F.FaceNodeCoords.size() > 2) { is3D = true; break; } }
	const int minNb = is3D ? 3 : 2;

	// 第一次：LSQ，失败则 GG
	for (const auto& CP : cells)
	{
		auto itP = id2idx.find(CP.id);
		if (itP == id2idx.end()) continue;
		const size_t iP = itP->second;
		const double phiP = (*phi)[iP];

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
			const double dphi = (*phi)[iN] - phiP;

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
			//cout << "对单元" << CP.id << "采用了LSQ算法计算梯度" << endl;
		}
		if (!ok) {
			grad[iP] = greenGaussGrad(mesh, reg, name, CP);  // GG 兜底
			//cout << "对单元" << CP.id << "采用了GG算法计算梯度" << endl;;
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
#pragma once
#include<vector>
#include<string>
#include<algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FieldAcessForDiscre.h"

// ======工具：网格单元上梯度计算：LSQ + GG 兜底 ======
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

inline Vector greenGaussGrad_linearConsistent(
    Mesh& mesh, const FieldRegistry& reg, const char* name,
    const Cell& CP, bool is3D,
    const std::vector<int>& nbCellsCached) // 传入已收集的邻居，便于在边界做临时LSQ
{
    const auto& cells = mesh.getCells();
    const auto& faces = mesh.getFaces();
    const auto& id2idx = mesh.getCellId2Index();
    auto phi = reg.get<volScalarField>(name);
    if (!phi) return Vector(0, 0, 0);

    // 1) 若需要，用一次局部LSQ估计 gLSQ 供边界面线性外推 φ_f
    auto lsq_local = [&](Vector& g)->bool {
        if (!is3D) {
            double Gxx = 0, Gxy = 0, Gyy = 0, bx = 0, by = 0;
            auto itP = id2idx.find(CP.id); if (itP == id2idx.end()) return false;
            const size_t iP = itP->second; const double phiP = (*phi)[iP];
            for (int N : nbCellsCached) {
                auto itN = id2idx.find(N); if (itN == id2idx.end()) continue;
                const size_t iN = itN->second;
                const Vector r3 = cells[iN].center - CP.center;
                const double rx = r3.m_x, ry = r3.m_y;
                const double r2 = std::max(rx * rx + ry * ry, 1e-20);
                const double w = 1.0 / r2;
                const double dphi = (*phi)[iN] - phiP;
                Gxx += w * rx * rx; Gxy += w * rx * ry; Gyy += w * ry * ry;
                bx += w * rx * dphi; by += w * ry * dphi;
            }
            const double tiny = 1e-20;
            Gxx += tiny; Gyy += tiny;
            const double det = Gxx * Gyy - Gxy * Gxy;
            if (std::abs(det) < 1e-30) return false;
            const double inv = 1.0 / det;
            g.m_x = (Gyy * bx - Gxy * by) * inv;
            g.m_y = (-Gxy * bx + Gxx * by) * inv;
            g.m_z = 0.0;
            return true;
        }
        else {
            // 3D 简化：复用一次 3×3 法方程
            Mat3 G; Vector b(0, 0, 0);
            auto itP = id2idx.find(CP.id); if (itP == id2idx.end()) return false;
            const size_t iP = itP->second; const double phiP = (*phi)[iP];
            for (int N : nbCellsCached) {
                auto itN = id2idx.find(N); if (itN == id2idx.end()) continue;
                const size_t iN = itN->second;
                const Vector r = cells[iN].center - CP.center;
                const double rmag2 = std::max(r.m_x * r.m_x + r.m_y * r.m_y + r.m_z * r.m_z, 1e-24);
                const double w = 1.0 / rmag2;
                const double dphi = (*phi)[iN] - phiP;
                G.addOuter(r, w); b = b + r * (w * dphi);
            }
            // 微正则
            G.a[0][0] += 1e-18; G.a[1][1] += 1e-18; G.a[2][2] += 1e-18;
            return G.solve(b, g);
        }
        };

    Vector gLSQ(0, 0, 0);
    bool haveLSQ = lsq_local(gLSQ);

    // 2) Green–Gauss：φ_f 线性一致重构，grad ≈ (1/V) Σ φ_f * S_f
    auto itP = id2idx.find(CP.id); if (itP == id2idx.end()) return Vector(0, 0, 0);
    const size_t iP = itP->second; const double phiP = (*phi)[iP];
    const double V = CP.volume; // 2D: area; 3D: volume
    if (V <= 0) return Vector(0, 0, 0);

    Vector fluxSum(0, 0, 0);
    for (int fid_ext : CP.CellFaceIDs) {
        const int idx = fid_ext - 1;
        if (idx < 0 || idx >= (int)faces.size()) continue;
        const Face& F = faces[idx];

        // 面矢量 S_f = n * A（2D 用 length 替代）
        const double Af = is3D ? /*F.area*/ F.length : /*F.length*/ F.length;   // <-- 若字段名不同请替换
        const Vector Sf = F.normal * Af;                                      // <-- F.normal 为单位法向；若非单位需改为 F.areaVector

        // 面心 φ_f
        double phi_f = phiP;
        if (F.ownerCell == CP.id && F.neighborCell >= 0) {
            // 内部面：简单平均（对线性场在对称几何/均匀网格线性一致；非均匀几何仍然OK作为GG的面值）
            const size_t iN = id2idx.at(F.neighborCell);
            const double phiN = (*phi)[iN];
            phi_f = 0.5 * (phiP + phiN);
        }
        else if (F.neighborCell == CP.id && F.ownerCell >= 0) {
            const size_t iN = id2idx.at(F.ownerCell);
            const double phiN = (*phi)[iN];
            phi_f = 0.5 * (phiP + phiN);
        }
        else {
            // 边界面：用局部LSQ梯度做一次线性外推，保证线性一致
            if (haveLSQ) {
                const Vector dx = F.midpoint - CP.center;                       // <-- 若无 F.center，请替换为可用的面几何中心
                phi_f = phiP + gLSQ*dx;
            }
            else {
                // 兜底：仍用 φP，至少保证稳定（线性一致性会变差，但通常不至于炸）
                phi_f = phiP;
            }
        }

        fluxSum = fluxSum + Sf * phi_f;
    }
    return fluxSum * (1.0 / V);
}

inline std::vector<Vector>
computeCellGradients_LSQ_with_GG(Mesh& mesh, const FieldRegistry& reg, const char* name, int smoothIters = 0)
{
    const auto& cells = mesh.getCells();
    const auto& faces = mesh.getFaces();
    const auto& id2idx = mesh.getCellId2Index();
    auto phi = reg.get<volScalarField>(name);
    std::vector<Vector> grad(cells.size(), Vector(0, 0, 0));
    if (!phi) return grad;

    // 判维
    bool is3D = false;
    for (const auto& c : cells) { if (std::abs(c.center.m_z) > 1e-14) { is3D = true; break; } }
    const int minNb = is3D ? 3 : 2;
    const int targetNb = is3D ? 6 : 4;   // 经验值：2D 至少 4 个样本，3D 至少 6 个样本更稳

    // ---- 主循环：每个单元做 LSQ；失败则 GG 兜底 ----
    for (const auto& CP : cells)
    {
        auto itP = id2idx.find(CP.id);
        if (itP == id2idx.end()) continue;
        const size_t iP = itP->second;
        const double phiP = (*phi)[iP];

        // 1) 一环邻居
        std::vector<int> nbCells;
        nbCells.reserve(16);
        auto contains = [&](const std::vector<int>& v, int cid)->bool {
            for (int x : v) if (x == cid) return true;
            return false;
            };
        for (int fid_ext : CP.CellFaceIDs) {
            const int idx = fid_ext - 1;
            if (idx < 0 || idx >= (int)faces.size()) continue;
            const Face& F = faces[idx];
            const int N = (F.ownerCell == CP.id ? F.neighborCell : F.ownerCell);
            if (N >= 0 && !contains(nbCells, N)) nbCells.push_back(N);
        }

        // 2) 无条件扩展到 targetNb：二环邻居
        if ((int)nbCells.size() < targetNb) {
            std::vector<int> firstRing = nbCells;
            for (int Nid : firstRing) {
                auto itN = id2idx.find(Nid);
                if (itN == id2idx.end()) continue;
                const auto& CN = cells[itN->second];
                for (int fid2 : CN.CellFaceIDs) {
                    const int idx2 = fid2 - 1;
                    if (idx2 < 0 || idx2 >= (int)faces.size()) continue;
                    const Face& F2 = faces[idx2];
                    const int NN = (F2.ownerCell == Nid ? F2.neighborCell : F2.ownerCell);
                    if (NN < 0 || NN == CP.id) continue;
                    if (!contains(nbCells, NN)) nbCells.push_back(NN);
                    if ((int)nbCells.size() >= targetNb) break;
                }
                if ((int)nbCells.size() >= targetNb) break;
            }
        }

        // 3) LSQ 求解
        bool ok = false;
        if ((int)nbCells.size() >= minNb)
        {
            if (!is3D) {
                // ---- 2D 显式 2×2 LSQ ----
                double Gxx = 0, Gxy = 0, Gyy = 0, bx = 0, by = 0;
                for (int N : nbCells) {
                    auto itN = id2idx.find(N);
                    if (itN == id2idx.end()) continue;
                    const size_t iN = itN->second;
                    const Vector r3 = cells[iN].center - CP.center;
                    const double rx = r3.m_x, ry = r3.m_y;
                    const double r2 = std::max(rx * rx + ry * ry, 1e-20);
                    const double w = 1.0 / r2;               // r^-2 权
                    const double dphi = (*phi)[iN] - phiP;
                    Gxx += w * rx * rx; Gxy += w * rx * ry; Gyy += w * ry * ry;
                    bx += w * rx * dphi; by += w * ry * dphi;
                }
                const double tiny = 1e-20;
                Gxx += tiny; Gyy += tiny;
                const double det = Gxx * Gyy - Gxy * Gxy;
                if (std::abs(det) > 1e-30) {
                    const double inv = 1.0 / det;
                    grad[iP].m_x = (Gyy * bx - Gxy * by) * inv;
                    grad[iP].m_y = (-Gxy * bx + Gxx * by) * inv;
                    grad[iP].m_z = 0.0;
                    ok = true;
                }
            }
            else {
                // ---- 3D 仍用 3×3 法方程 ----
                Mat3 G; Vector b(0, 0, 0);
                for (int N : nbCells) {
                    auto itN = id2idx.find(N);
                    if (itN == id2idx.end()) continue;
                    const size_t iN = itN->second;
                    const Vector r = cells[iN].center - CP.center;
                    const double rmag2 = std::max(r.m_x * r.m_x + r.m_y * r.m_y + r.m_z * r.m_z, 1e-24);
                    const double w = 1.0 / rmag2;
                    const double dphi = (*phi)[iN] - phiP;
                    G.addOuter(r, w);
                    b = b + r * (w * dphi);
                }
                G.a[0][0] += 1e-18; G.a[1][1] += 1e-18; G.a[2][2] += 1e-18;  // 轻微对角正则
                ok = G.solve(b, grad[iP]);
            }
        }

        // 4) 兜底：线性一致 Green–Gauss
        if (!ok) {
            grad[iP] = greenGaussGrad_linearConsistent(mesh, reg, name, CP, is3D, nbCells);
        }
    }

    // ---- 可选平滑（单元测试时建议 smoothIters=0）----
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

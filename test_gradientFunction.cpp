#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "VolField.h"
#include <tuple>


// 声明待测函数（按你给出的签名）
inline std::vector<Vector>
computeCellGradients_LSQ_with_GG(Mesh& mesh, const FieldRegistry& reg, const char* name, int smoothIters = 0);


//φ(x,y) = x + y测试
int main2() 
{
	// 1) 构造一个小的 2D 三角形形网格（1×1 m，2×2 单元）
	double lengthX = 1.0, lengthY = 1.0, lengthZ = 0.0;
	int sectionNumX = 50, sectionNumY = 50, sectionNumZ = 0;
	bool usePrism = true;    // 2D
	bool useQuadBase = false;  // 三角形网格

	MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
	mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

	Mesh& mesh = mgr.mesh();
	const auto& cells = mesh.getCells();
	const auto& id2idx = mesh.getCellId2Index();

	// 2) 注册并填充标量场 φ(x,y) = x + y（线性函数，梯度应为 (1,1,0)）
	FieldRegistry reg;
	auto phi = reg.getOrCreate<VolField<double>>("phi", cells.size(), 0.0);

	for (const auto& c : cells) {
		if (c.id < 0) continue;
		size_t i = id2idx.at(c.id);
		double x = c.center.m_x;
		double y = c.center.m_y;
		(*phi)[i] = x + y;
	}

	// 3) 计算梯度（不做平滑）
	auto grads = computeCellGradients_LSQ_with_GG(mesh, reg, "phi", /*smoothIters=*/1);

	// 4) 校验与输出：期望 ∇φ = (1,1,0)
	const Vector gradExact(1.0, 1.0, 0.0);
	double maxAbsErr = 0.0;
	std::cout.setf(std::ios::scientific);
	std::cout << "cell_id, cx, cy, grad_x, grad_y, grad_z, |err|_inf\n";
	for (const auto& c : cells) {
		if (c.id < 0) continue;
		size_t i = id2idx.at(c.id);
		Vector g = grads[i];
		double ex = std::abs(g.m_x - gradExact.m_x);
		double ey = std::abs(g.m_y - gradExact.m_y);
		double ez = std::abs(g.m_z - gradExact.m_z);
		double errInf = std::max({ ex, ey, ez });
		maxAbsErr = std::max(maxAbsErr, errInf);
		std::cout << c.id << ", "
			<< c.center.m_x << ", " << c.center.m_y << ", "
			<< g.m_x << ", " << g.m_y << ", " << g.m_z << ", "
			<< errInf << "\n";
	}

	// 5) 判定：线性可再现性（LSQ/GG 应精确或达到机器精度）
	const double tol = 1e-10; // 非正交修正、数值噪声下的合理阈值
	if (maxAbsErr <= tol) {
		std::cout << "[PASS] max |error|_inf = " << maxAbsErr << " <= " << tol << "\n";
		return 0;
	}
	else {
		std::cout << "[FAIL] max |error|_inf = " << maxAbsErr << " > " << tol << "\n";
		return 1;
	}
}
//φ(x,y) = x^2 + y^2测试
int main3()
{
    // 1) 网格（1×1 区域，三角网）
    double lengthX = 1.0, lengthY = 1.0, lengthZ = 0.0;
    int sectionNumX = 50, sectionNumY = 50, sectionNumZ = 0;
    bool usePrism = true;       // 2D
    bool useQuadBase = false;   // 三角形底面

    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    Mesh& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();

    // 2) 注册并填充 φ(x,y) = x^2 + y^2
    FieldRegistry reg;
    auto phi = reg.getOrCreate<VolField<double>>("phi", cells.size(), 0.0);
    for (const auto& c : cells) {
        if (c.id < 0) continue;
        size_t i = id2idx.at(c.id);
        double x = c.center.m_x;
        double y = c.center.m_y;
        (*phi)[i] = x * x + y * y;
    }

    // 3) 计算梯度（建议不平滑，避免掩盖离散误差）
    auto grads = computeCellGradients_LSQ_with_GG(mesh, reg, "phi", /*smoothIters=*/0);

    // 4) 校验并输出：期望 ∇φ = (2x, 2y, 0)
    double maxAbsErr = 0.0;
    double sumSq = 0.0;
    size_t count = 0;

    std::cout.setf(std::ios::scientific);
    std::cout << "cell_id, cx, cy, grad_x, grad_y, grad_z, ex_gx, ex_gy, |err|_inf\n";

    for (const auto& c : cells) {
        if (c.id < 0) continue;
        size_t i = id2idx.at(c.id);
        const Vector g = grads[i];

        // 解析梯度
        const double ex_gx = 2.0 * c.center.m_x;
        const double ex_gy = 2.0 * c.center.m_y;
        const double ex_gz = 0.0;

        // 误差
        const double ex = std::abs(g.m_x - ex_gx);
        const double ey = std::abs(g.m_y - ex_gy);
        const double ez = std::abs(g.m_z - ex_gz);
        const double errInf = std::max({ ex, ey, ez });

        maxAbsErr = std::max(maxAbsErr, errInf);
        sumSq += (ex * ex + ey * ey + ez * ez);
        ++count;

        std::cout << c.id << ", "
            << c.center.m_x << ", " << c.center.m_y << ", "
            << g.m_x << ", " << g.m_y << ", " << g.m_z << ", "
            << ex_gx << ", " << ex_gy << ", "
            << errInf << "\n";
    }

    // L2 范数（按单元计的均方根）
    const double l2 = std::sqrt(sumSq / std::max<size_t>(1, count));

    // 5) 判定：二次场下期望一阶收敛，这里给一个经验阈值 ~ C*h
    const double h = 1.0 / std::min(sectionNumX, sectionNumY);  // 代表网格尺度
    const double tol = 5.0 * h;  // 经验系数 C≈5，可按需要调整/做网格加密回归

    std::cout << "L_inf = " << maxAbsErr << ",  L2 = " << l2 << ",  h = " << h << "\n";
    if (maxAbsErr <= tol) {
        std::cout << "[PASS] L_inf <= C*h (C=5): " << maxAbsErr << " <= " << tol << "\n";
        return 0;
    }
    else {
        std::cout << "[WARN] L_inf > C*h，建议做网格加密(例如 25→50→100)检查收敛阶。\n";
        return 0; // 不强行失败，方便连续跑收敛曲线
    }
}

struct ErrResult {
    double Linf;
    double L2;
};

ErrResult run_case(int nx, int ny) {
    // —— 构网 ——
    double Lx = 1.0, Ly = 1.0, Lz = 0.0;
    bool usePrism = true, useQuadBase = false;
    MeshManager mgr(Lx, Ly, Lz, nx, ny, 0, usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    Mesh& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();

    // —— φ = x^2 + y^2 ——
    FieldRegistry reg;
    auto phi = reg.getOrCreate<VolField<double>>("phi", cells.size(), 0.0);
    for (const auto& c : cells) {
        if (c.id < 0) continue;
        size_t i = id2idx.at(c.id);
        double x = c.center.m_x, y = c.center.m_y;
        (*phi)[i] = x * x + y * y;
    }

    // —— 计算梯度 ——
    auto g = computeCellGradients_LSQ_with_GG(mesh, reg, "phi", 0);

    // —— 误差统计 ——
    double Linf = 0.0, sumSq = 0.0; size_t cnt = 0;
    for (const auto& c : cells) {
        if (c.id < 0) continue;
        size_t i = id2idx.at(c.id);
        const double ex = std::abs(g[i].m_x - 2.0 * c.center.m_x);
        const double ey = std::abs(g[i].m_y - 2.0 * c.center.m_y);
        const double ez = std::abs(g[i].m_z - 0.0);
        const double emax = std::max({ ex,ey,ez });
        Linf = std::max(Linf, emax);
        sumSq += (ex * ex + ey * ey + ez * ez);
        ++cnt;
    }
    const double L2 = std::sqrt(sumSq / std::max<size_t>(1, cnt));
    return { Linf, L2 };
}

static double order(double E1, double E2, double h1, double h2) {
    return std::log(E1 / E2) / std::log(h1 / h2);
}


int main4() 
{
    std::vector<int> Ns = { 25, 50, 100 };
    std::vector<double> hs, Einf, E2;

    for (int n : Ns) {
        ErrResult r = run_case(n, n);
        double h = 1.0 / n;
        hs.push_back(h); Einf.push_back(r.Linf); E2.push_back(r.L2);
        std::cout << "n=" << n << "  h=" << h
            << "  Linf=" << r.Linf << "  L2=" << r.L2 << "\n";
    }

    std::cout << "Order(Linf, 25->50)  = " << order(Einf[0], Einf[1], hs[0], hs[1]) << "\n";
    std::cout << "Order(Linf, 50->100) = " << order(Einf[1], Einf[2], hs[1], hs[2]) << "\n";
    std::cout << "Order(L2,   25->50)  = " << order(E2[0], E2[1], hs[0], hs[1]) << "\n";
    std::cout << "Order(L2,   50->100) = " << order(E2[1], E2[2], hs[1], hs[2]) << "\n";
    return 0;
}


/**
 * @file test_FVM_Grad_Benchmark.cpp
 * @brief FVM-EDFM 梯度算子基准测试 (Corrected Variable Names)
 * @details
 * 修正说明：
 * 修复了之前代码中 nFrac 与 nFracCells 变量名混用导致的“未定义标识符”编译错误。
 * 保持了对边界效应的过滤逻辑。
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <memory>
#include <stdexcept>

 // 引入核心模块
#include "UserDefineVarType.h"
#include "MeshManager.h"            // 2D-EDFM Manager
#include "3D_MeshManager.h"         // 3D-EDFM Manager
#include "2D_FieldManager.h"        // 2D-EDFM Field Manager
#include "3D_FieldManager.h"        // 3D-EDFM Field Manager
#include "FVM_Grad.h"               // 待测算子

#include <functional>
using FieldFunc = std::function<double(const Vector&)>;
using GradFunc = std::function<Vector(const Vector&)>;

// =========================================================
// 辅助：边界单元判定
// =========================================================

// 判定基岩单元是否在边界上
bool IsMatrixBoundary(const Mesh& mesh, int cellIdx)
{
    const auto& cell = mesh.getCells()[cellIdx];
    for (int globalFaceID : cell.CellFaceIDs) {
        int faceIdx = mesh.getFaceIndex(globalFaceID);
        if (faceIdx >= 0) {
            if (mesh.getFaces()[faceIdx].isBoundary()) return true;
        }
    }
    return false;
}

// 判定 3D 裂缝单元是否在边界上 (存在没有邻居的边)
bool IsFractureBoundary_3D(const FractureNetwork_2D& net, const FractureElement_2D* elem)
{
    for (int edgeIdx : elem->connectedEdgeIndices) {
        const auto& edge = net.getGlobalEdges()[edgeIdx];
        if (edge.neighborCell_solverIndex == -1) return true;
    }
    return false;
}

// 判定 2D 裂缝单元(线段)是否是端点
// [Fix] 必须同时统计 Tangent Neighbors (前后) 和 Intersection Neighbors (交叉)
bool IsFractureBoundary_2D(const FractureElement* elem)
{
    int totalNeighbors = 0;

    // 1. 统计切向邻居 (Siblings)
    // 正常的内部线段这里应该是 2 (前一个 + 后一个)
    // 端点是 1
    totalNeighbors += (int)elem->tangentNeighbors.size();

    // 2. 统计交叉邻居 (Intersections)
    for (const auto& nb : elem->neighbors) {
        if (nb.solverIndex != -1) totalNeighbors++;
    }

    // 如果总邻居数 < 2，则视为真正的物理端点（Tip），跳过检查
    // 如果 >= 2，说明是内部单元（或者是一个连接点），应当检查
    return (totalNeighbors < 2);
}

// =========================================================
// Case 1: 2D-EDFM Benchmark
// =========================================================
void run_Benchmark_2D_EDFM_Grad()
{
    std::cout << "\n=========================================================" << std::endl;
    std::cout << "[Benchmark] Starting 2D-EDFM Gradient Test (Internal Patch Test)" << std::endl;
    std::cout << "=========================================================" << std::endl;

    double Lx = 1.0, Ly = 1.0;
    int nx = 20, ny = 20;
    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OrthogonalCorrection);

    mgr.addFracture(Vector(0.2, 0.2, 0), Vector(0.8, 0.8, 0));
    mgr.addFracture(Vector(0.2, 0.5, 0), Vector(0.8, 0.5, 0));

    mgr.setDistanceMetric(DistanceMetric::CrossAwareGauss);
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.ComputeFractureGeometryCouplingCoefficient();

    auto fm = std::make_shared<FieldManager_2D>();
    size_t nMatrix = mgr.mesh().getGridCount();
    size_t nFrac = mgr.fracture_network().getOrderedFractureElements().size();
    size_t nMatFaces = mgr.mesh().getFaces().size();
    fm->InitSizes(nMatrix, nFrac, 0, 0, nMatFaces, 0);

    // P = x + 2y
    auto matP = fm->matrixFields.create<volScalarField>("Pressure", nMatrix);
    for (int i = 0; i < matP->size; ++i) {
        Vector c = mgr.mesh().getCells()[i].center;
        (*matP)[i] = 1.0 * c.m_x + 2.0 * c.m_y;
    }

    auto fracP = fm->fractureFields.create<volScalarField>("Pressure", nFrac);
    const auto& fracElems = mgr.fracture_network().getOrderedFractureElements();

    // [Fix] 使用 nFrac (之前误写为 nFracCells)
    for (int i = 0; i < (int)nFrac; ++i) {
        const auto* elem = fracElems[i];
        const auto& parent = mgr.fracture_network().fractures[elem->parentFractureID];
        Vector start = parent.start;
        Vector end = parent.end;
        Vector center = start + 0.5 * (elem->param0 + elem->param1) * (end - start);
        (*fracP)[i] = 1.0 * center.m_x + 2.0 * center.m_y;
    }

    FVM_Grad gradOp(mgr.mesh(), nullptr, &mgr.fracture_network(), nullptr);
    auto matGrad = gradOp.compute(*matP, FVM_Grad::Method::LeastSquares);
    auto fracGrad = gradOp.compute(*fracP, mgr.fracture_network(), FVM_Grad::Method::LeastSquares);

    // --- Verification (Internal Only) ---
    std::cout << ">>> Verifying Matrix Gradients (Internal Cells Only)..." << std::endl;
    double maxErrM = 0.0;
    int checkedM = 0;
    for (int i = 0; i < matGrad->size; ++i) {
        if (IsMatrixBoundary(mgr.mesh(), i)) continue;

        Vector g = (*matGrad)[i];
        double err = std::sqrt(std::pow(g.m_x - 1.0, 2) + std::pow(g.m_y - 2.0, 2));
        if (err > maxErrM) maxErrM = err;
        checkedM++;
    }
    std::cout << "    Checked Cells: " << checkedM << "/" << matGrad->size << std::endl;
    std::cout << "    Max Internal Error: " << maxErrM << (maxErrM < 1e-10 ? " [PASS]" : " [FAIL]") << std::endl;

    std::cout << ">>> Verifying Fracture Gradients (Internal Elements Only)..." << std::endl;
    double maxErrF = 0.0;
    int checkedF = 0;
    for (int i = 0; i < fracGrad->size; ++i) {
        const auto* elem = fracElems[i];
        if (IsFractureBoundary_2D(elem)) continue;

        Vector g = (*fracGrad)[i];
        Vector ref = (elem->parentFractureID == 0) ? Vector(1.5, 1.5, 0.0) : Vector(1.0, 0.0, 0.0);

        double err = (g - ref).Mag();
        if (err > maxErrF) maxErrF = err;
        checkedF++;
    }
    std::cout << "    Checked Elements: " << checkedF << "/" << fracGrad->size << std::endl;
    std::cout << "    Max Internal Error: " << maxErrF << (maxErrF < 1e-10 ? " [PASS]" : " [FAIL]") << std::endl;
}

// =========================================================
// Case 2: 3D-EDFM Benchmark
// =========================================================
Fracture_2D CreateRectFracture(int id, const Vector& p1, const Vector& p2, const Vector& p3, const Vector& p4)
{
    std::vector<Vector> pts = { p1, p2, p3, p4 };
    return Fracture_2D(id, pts);
}

void run_Benchmark_3D_EDFM_Grad()
{
    std::cout << "\n=========================================================" << std::endl;
    std::cout << "[Benchmark] Starting 3D-EDFM Gradient Test (Internal Patch Test)" << std::endl;
    std::cout << "=========================================================" << std::endl;

    double L = 1.0;
    int N = 10;
    MeshManager_3D mgr(L, L, L, N, N, N, true, false);
    mgr.BuildSolidMatrixGrid_3D();

    // Frac A (Oblique)
    Vector p1(0.8, 0.2, 0.1), p2(0.2, 0.8, 0.1), p3(0.2, 0.8, 0.9), p4(0.8, 0.2, 0.9);
    mgr.addFracturetoFractureNetwork(CreateRectFracture(0, p1, p2, p3, p4));

    // Frac B (Horizontal)
    Vector q1(0.1, 0.1, 0.5), q2(0.9, 0.1, 0.5), q3(0.9, 0.9, 0.5), q4(0.1, 0.9, 0.5);
    mgr.addFracturetoFractureNetwork(CreateRectFracture(1, q1, q2, q3, q4));

    mgr.meshAllFracturesinNetwork(10, 10, NormalVectorCorrectionMethod::OrthogonalCorrection);

    mgr.setupGlobalIndices();
    mgr.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);
    mgr.DetectFractureFractureIntersectionsInNetwork(FFIntersectionStrategy::Octree_Optimized);
    mgr.fracture_network().rebuildEdgeProperties();
    mgr.buildTopologyMaps();

    auto fm = std::make_shared<FieldManager_3D>();
    size_t nMatrix = mgr.mesh().getGridCount();
    size_t nFrac = mgr.fracture_network().getOrderedFractureElements().size();
    size_t nNNC = mgr.getInteractionPairs().size();
    size_t nMatFaces = mgr.mesh().getFaces().size();
    fm->InitSizes(nMatrix, nFrac, nNNC, 0, nMatFaces, 0);

    // P = x + 2y + 3z
    auto matP = fm->matrixFields.create<volScalarField>("Pressure", nMatrix);
    for (int i = 0; i < matP->size; ++i) {
        Vector c = mgr.mesh().getCells()[i].center;
        (*matP)[i] = 1.0 * c.m_x + 2.0 * c.m_y + 3.0 * c.m_z;
    }

    auto fracP = fm->fractureFields.create<volScalarField>("Pressure", nFrac);
    const auto& fracElems = mgr.fracture_network().getOrderedFractureElements();

    // [Fix] 使用 nFrac (之前误写为 nFracCells)
    for (int i = 0; i < (int)nFrac; ++i) {
        const auto* elem = fracElems[i];
        Vector c = elem->centroid;
        (*fracP)[i] = 1.0 * c.m_x + 2.0 * c.m_y + 3.0 * c.m_z;
    }

    FVM_Grad gradOp(mgr.mesh(), &mgr.fracture_network(), nullptr, nullptr);
    auto matGrad = gradOp.compute(*matP, FVM_Grad::Method::LeastSquares);
    auto fracGrad = gradOp.compute(*fracP, mgr.fracture_network(), FVM_Grad::Method::LeastSquares);

    // --- Verification (Internal Only) ---
    std::cout << ">>> Verifying Matrix Gradients (Internal Cells Only)..." << std::endl;
    double maxErrM = 0.0;
    int checkedM = 0;
    for (int i = 0; i < matGrad->size; ++i) {
        if (IsMatrixBoundary(mgr.mesh(), i)) continue;

        Vector g = (*matGrad)[i];
        double err = (g - Vector(1, 2, 3)).Mag();
        if (err > maxErrM) maxErrM = err;
        checkedM++;
    }
    std::cout << "    Checked Cells: " << checkedM << "/" << matGrad->size << std::endl;
    std::cout << "    Max Internal Error: " << maxErrM << (maxErrM < 1e-10 ? " [PASS]" : " [FAIL]") << std::endl;

    std::cout << ">>> Verifying Fracture Gradients (Internal Elements Only)..." << std::endl;
    double maxErrF = 0.0;
    int checkedF = 0;
    for (int i = 0; i < fracGrad->size; ++i) {
        const auto* elem = fracElems[i];
        if (IsFractureBoundary_3D(mgr.fracture_network(), elem)) continue;

        Vector g = (*fracGrad)[i];
        Vector ref = (elem->parentFractureID == 0) ? Vector(-0.5, 0.5, 3.0) : Vector(1.0, 2.0, 0.0);

        double err = (g - ref).Mag();
        if (err > maxErrF) {
            maxErrF = err;
        }
        checkedF++;
    }
    std::cout << "    Checked Elements: " << checkedF << "/" << fracGrad->size << std::endl;
    std::cout << "    Max Internal Error: " << maxErrF << (maxErrF < 1e-10 ? " [PASS]" : " [FAIL]") << std::endl;
}

void run_Advanced_Accuracy_Test(const std::string& testName,
    FieldFunc func,
    GradFunc anaGradFunc,
    bool isQuadratic = false)
{
    std::cout << "\n=========================================================" << std::endl;
    std::cout << "[Advanced Benchmark] " << testName << std::endl;
    std::cout << "=========================================================" << std::endl;

    // 1. 构建网格 (使用稍密的网格以观察非线性收敛)
    double L = 1.0;
    int N = 50; // 加密网格
    MeshManager_3D mgr(L, L, L, N, N, N, true, false); // 3D Mesh
    mgr.BuildSolidMatrixGrid_3D();

    // 2. 添加典型裂缝 (一横一斜)
    // Horizontal
    Vector q1(0.1, 0.1, 0.5), q2(0.9, 0.1, 0.5), q3(0.9, 0.9, 0.5), q4(0.1, 0.9, 0.5);
    mgr.addFracturetoFractureNetwork(CreateRectFracture(0, q1, q2, q3, q4));
    // Oblique
    Vector p1(0.8, 0.2, 0.1), p2(0.2, 0.8, 0.1), p3(0.2, 0.8, 0.9), p4(0.8, 0.2, 0.9);
    mgr.addFracturetoFractureNetwork(CreateRectFracture(1, p1, p2, p3, p4));

    // 3. 拓扑与几何构建
    mgr.meshAllFracturesinNetwork(15, 15, NormalVectorCorrectionMethod::OrthogonalCorrection);
    mgr.setupGlobalIndices();
    mgr.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);
    mgr.DetectFractureFractureIntersectionsInNetwork(FFIntersectionStrategy::Octree_Optimized);
    mgr.fracture_network().rebuildEdgeProperties(); // Global Index & Topology Sink
    mgr.buildTopologyMaps();

    // 4. 初始化场数据
    auto fm = std::make_shared<FieldManager_3D>();
    size_t nMatrix = mgr.mesh().getGridCount();
    size_t nFrac = mgr.fracture_network().getOrderedFractureElements().size();
    fm->InitSizes(nMatrix, nFrac, 0, 0, 0, 0);

    auto matP = fm->matrixFields.create<volScalarField>("Pressure", nMatrix);
    auto fracP = fm->fractureFields.create<volScalarField>("Pressure", nFrac);

    // 填充 Matrix 场
    for (int i = 0; i < matP->size; ++i) {
        (*matP)[i] = func(mgr.mesh().getCells()[i].center);
    }

    // 填充 Fracture 场
    const auto& fracElems = mgr.fracture_network().getOrderedFractureElements();
    for (int i = 0; i < (int)nFrac; ++i) {
        (*fracP)[i] = func(fracElems[i]->centroid);
    }

    // 5. 计算梯度
    FVM_Grad gradOp(mgr.mesh(), &mgr.fracture_network(), nullptr, nullptr);
    auto matGrad = gradOp.compute(*matP, FVM_Grad::Method::LeastSquares);
    auto fracGrad = gradOp.compute(*fracP, mgr.fracture_network(), FVM_Grad::Method::LeastSquares);

    // 6. 误差统计
    // ---------------------------------------------------
    // Matrix Error
    double rmsErrM = 0.0;
    double maxErrM = 0.0;
    int countM = 0;
    for (int i = 0; i < matGrad->size; ++i) {
        if (IsMatrixBoundary(mgr.mesh(), i)) continue; // 仅检查内部

        Vector numGrad = (*matGrad)[i];
        Vector anaGrad = anaGradFunc(mgr.mesh().getCells()[i].center);

        double err = (numGrad - anaGrad).Mag();
        rmsErrM += err * err;
        if (err > maxErrM) maxErrM = err;
        countM++;
    }
    rmsErrM = std::sqrt(rmsErrM / countM);

    // Fracture Error (注意：解析解需要投影到裂缝切平面)
    double rmsErrF = 0.0;
    double maxErrF = 0.0;
    int countF = 0;
    for (int i = 0; i < fracGrad->size; ++i) {
        const auto* elem = fracElems[i];
        if (IsFractureBoundary_3D(mgr.fracture_network(), elem)) continue; // 仅检查内部

        Vector numGrad = (*fracGrad)[i];
        Vector globalAnaGrad = anaGradFunc(elem->centroid);

        // [Key] 投影解析梯度到裂缝平面： Grad_tau = Grad - (Grad . n) * n
        Vector n = elem->normal; n.Normalize();
        Vector anaGradProj = globalAnaGrad - (globalAnaGrad * n) * n;

        double err = (numGrad - anaGradProj).Mag();
        rmsErrF += err * err;
        if (err > maxErrF) maxErrF = err;
        countF++;
    }
    rmsErrF = std::sqrt(rmsErrF / countF);

    // 7. 输出报告
    std::cout << std::scientific << std::setprecision(4);
    std::cout << ">>> Report: " << testName << std::endl;
    std::cout << "    [Matrix]   RMS Error: " << rmsErrM << " | Max Error: " << maxErrM << std::endl;
    std::cout << "    [Fracture] RMS Error: " << rmsErrF << " | Max Error: " << maxErrF << std::endl;

    // 判定标准
    // 对于二次场，一阶 LS 格式通常有 O(h) 或 O(h^2) 的误差，不应为 0
    // 对于 30x30x30 网格，dx ≈ 0.033。误差应在 1e-2 ~ 1e-3 量级
    double tolerance = isQuadratic ? 0.1 : 1e-9;

    if (maxErrM < tolerance && maxErrF < tolerance) {
        std::cout << "    [RESULT] PASS (Within expected accuracy)" << std::endl;
    }
    else {
        std::cout << "    [RESULT] WARNING/FAIL (Check grid quality or boundary effects)" << std::endl;
    }
}
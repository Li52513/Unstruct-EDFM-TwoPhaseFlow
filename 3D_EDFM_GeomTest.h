#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <functional>

// 引入核心头文件
#include "3D_MeshManager.h"

// ---------------------------------------------------------
// 辅助函数：快速创建四边形裂缝
// ---------------------------------------------------------
Fracture_2D CreateFracture_GeoTest(int id, const Vector& p1, const Vector& p2, const Vector& p3, const Vector& p4)
{
    std::vector<Vector> boundaryVertices = { p1, p2, p3, p4 };
    return Fracture_2D(id, boundaryVertices);
}

// ---------------------------------------------------------
// 辅助类：统计误差信息的累加器
// ---------------------------------------------------------
struct ErrorStats {
    double maxError = 0.0;
    double sumError = 0.0;
    int count = 0;
    int failCount = 0;

    void update(double computed, double expected, double tol = 1e-5) {
        double err = std::abs(computed - expected);
        maxError = std::max(maxError, err);
        sumError += err;
        count++;
        if (err > tol) failCount++;
    }

    void printReport(const std::string& prefix) const {
        std::cout << "[" << prefix << "] Verified " << count << " pairs." << std::endl;
        std::cout << "    -> Max Error: " << std::scientific << maxError << std::defaultfloat << std::endl;
        std::cout << "    -> Avg Error: " << (count > 0 ? sumError / count : 0.0) << std::endl;
        if (failCount == 0) {
            std::cout << "    -> Result: \033[32mPASS\033[0m (All within tolerance)" << std::endl;
        }
        else {
            std::cout << "    -> Result: \033[31mFAIL\033[0m (" << failCount << " pairs exceeded tolerance)" << std::endl;
        }
    }
};

/**
 * @brief 验证 Step 1 几何参数准确性的 Benchmark 测试函数
 * @details 包含两个子测试用例：
 * Case A: 垂直裂缝 (验证常数距离)
 * Case B: 45度斜裂缝 (验证投影公式)
 */
int RunBenchmark_Step1_Distance_Accuracy()
{
    std::cout << "\n=========================================================" << std::endl;
    std::cout << "   Benchmark: Step 1 - NNC Geometry Accuracy Verification" << std::endl;
    std::cout << "=========================================================" << std::endl;

    // -------------------------------------------------------------------------
    // 基础网格配置 (100x100x100, 10x10x10 Cells)
    // -------------------------------------------------------------------------
    // 基岩网格尺寸: 10m x 10m x 10m
    // 单元中心通式: Center(i,j,k) = (10*i+5, 10*j+5, 10*k+5)
    double L = 100.0;
    int N = 10;

    // 初始化管理器
    std::cout << "-> Initializing Matrix Mesh (100x100x100)..." << std::endl;
    MeshManager_3D manager(L, L, L, N, N, N, true, false); // usePrism=true (Hex), QuadBase
    manager.BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod::OrthogonalCorrection, "Bench_Matrix");

    // -------------------------------------------------------------------------
    // Case A: 垂直裂缝测试 (Vertical Fracture at x = 52.0)
    // -------------------------------------------------------------------------
    // 理论推导:
    // 裂缝平面: x = 52
    // 相交基岩单元: i=5 (x范围 [50, 60])
    // 基岩中心 X坐标: 10*5 + 5 = 55.0
    // 理论 d_NNC = |55.0 - 52.0| = 3.0
    // -------------------------------------------------------------------------
    std::cout << "\n[Test Case A] Vertical Fracture at x=52.0" << std::endl;

    Fracture_2D fracVert = CreateFracture_GeoTest(1,
        Vector(52, 0, 0), Vector(52, 100, 0),
        Vector(52, 100, 100), Vector(52, 0, 100));

    manager.addFracturetoFractureNetwork(fracVert);
    manager.meshAllFracturesinNetwork(5, 5); // 5x5 sub-grids

    // 执行求交 (使用目前最高效的策略)
    manager.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);

    // [修正] 验证逻辑
    const auto& pairsA = manager.getInteractionPairs();
    ErrorStats statsA;

    // 获取网格引用以便查询真实中心
    const auto& cells_A = manager.mesh().getCells();

    for (const auto& pair : pairsA) {
        // 1. 获取基岩单元的真实几何中心
        int cellIndex = manager.mesh().getCellIndex(pair.matrixCellGlobalID);
        Vector center = cells_A[cellIndex].center;

        // 2. 动态计算理论距离
        // 裂缝面 x = 52.0，距离 d = |x_center - 52.0|
        double theoryDist = std::abs(center.m_x - 52.0);

        // 3. 比对代码计算值与理论值
        statsA.update(pair.distMatrixToFracPlane, theoryDist);
    }
    statsA.printReport("Case A");

    // 导出以便人工复核
    manager.exportInteractionPairsToCSV("./3D_EDFM_PostTest/Bench_Step1_CaseA.csv");


    // -------------------------------------------------------------------------
    // 重置环境 (为了 Case B，重新构建一个新的 Manager 实例)
    // -------------------------------------------------------------------------
    // 注意：实际工程中通常会清空列表，这里为了代码清晰直接新建对象
    std::cout << "\n[Test Case B] Diagonal Fracture (Plane x - z = 0)" << std::endl;
    MeshManager_3D managerB(L, L, L, N, N, N, true, false);
    managerB.BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod::OrthogonalCorrection, "Bench_Matrix_B");

    // -------------------------------------------------------------------------
    // Case B: 斜裂缝测试 (Diagonal Fracture, 45 degree)
    // -------------------------------------------------------------------------
    // 理论推导:
    // 裂缝平面方程: x - z = 0 (法向 n = (1, 0, -1)/sqrt(2))
    // 基岩中心: C = (cx, cy, cz)
    // 理论 d_NNC = |(C - P_on_plane) . n| = |(cx - cz) / sqrt(2)|
    // -------------------------------------------------------------------------
    Fracture_2D fracDiag = CreateFracture_GeoTest(2,
        Vector(0, 0, 0), Vector(100, 0, 100),   // Bottom-Left to Top-Right
        Vector(100, 100, 100), Vector(0, 100, 0) // Extruded along Y
    );

    managerB.addFracturetoFractureNetwork(fracDiag);
    managerB.meshAllFracturesinNetwork(5, 5);

    // 执行求交
    managerB.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);

    // 验证逻辑
    const auto& pairsB = managerB.getInteractionPairs();
    ErrorStats statsB;
    const double sqrt2 = std::sqrt(2.0);

    // 获取基岩单元列表以查询中心坐标 (需要 Mesh 提供访问接口)
    // 假设 Mesh 有 getCells() 方法
    const auto& cells = managerB.mesh().getCells();

    for (const auto& pair : pairsB) {
        // 1. 获取基岩中心
        // 注意：interactionPair 中存储的是 matrixCellGlobalID
        // 需通过 ID 获取 Cell 对象。这里假设 cellId2Index 映射或直接索引
        // 简化起见，直接利用 geometry 算出的距离进行反向验证，
        // 或者如果 Pair 里没存 Center，我们需要去 Mesh 里查。
        // *关键*：Step 1 修改中，我们在 Pair 里并没有存 Matrix Center，只存了结果 d_NNC。
        // 因此这里需要根据 ID 反算中心进行验证。

        // 假设 cells 是按 ID 顺序存储的，或者通过 getCellIndex 查询
        int cellIndex = managerB.mesh().getCellIndex(pair.matrixCellGlobalID);
        Vector center = cells[cellIndex].center;

        // 2. 计算理论值
        double theoryDist = std::abs(center.m_x - center.m_z) / sqrt2;

        // 3. 比对
        statsB.update(pair.distMatrixToFracPlane, theoryDist);
    }
    statsB.printReport("Case B");

    // 导出结果
    managerB.exportInteractionPairsToCSV("./3D_EDFM_PostTest/Bench_Step1_CaseB.csv");

    std::cout << "\n=========================================================" << std::endl;
    std::cout << "   Benchmark Completed." << std::endl;
    std::cout << "=========================================================" << std::endl;
	return 0;
}


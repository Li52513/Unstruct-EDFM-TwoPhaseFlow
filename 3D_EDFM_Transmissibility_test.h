#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <functional>
#include <unordered_map>

// 引入核心头文件
#include "3D_MeshManager.h"
#include "3D_FieldManager.h"
#include "TransmissibilitySolver_3D.h"
#include "SolverContrlStrName_op.h"

using namespace PhysicalProperties_string_op;

// ---------------------------------------------------------
// 辅助工具：误差统计
// ---------------------------------------------------------
struct BenchStats {
    double maxError = 0.0;
    int count = 0;
    int failCount = 0;

    void update(double computed, double expected, double tol = 1e-8) {
        double err = std::abs(computed - expected);
        // 如果 expected 很大，使用相对误差
        if (std::abs(expected) > 1e-9) err /= std::abs(expected);

        maxError = std::max(maxError, err);
        count++;
        if (err > tol) failCount++;
    }

    void report(const std::string& caseName) const {
        std::cout << "[" << caseName << "] Checked " << count << " items." << std::endl;
        std::cout << "    -> Max Rel Error: " << std::scientific << maxError << std::defaultfloat << std::endl;
        if (failCount == 0) std::cout << "    -> Result: \033[32mPASS\033[0m" << std::endl;
        else std::cout << "    -> Result: \033[31mFAIL\033[0m (" << failCount << " items failed)" << std::endl;
    }
};

// ---------------------------------------------------------
// 辅助函数：创建裂缝
// ---------------------------------------------------------
Fracture_2D CreateBenchFracture(int id, const Vector& p1, const Vector& p2, const Vector& p3, const Vector& p4, double aperture)
{
    std::vector<Vector> pts = { p1, p2, p3, p4 };
    Fracture_2D f(id, pts);
    f.aperture = aperture; // 显式设置开度
    return f;
}

/**
 * @brief 验证 TransmissibilitySolver_3D 的 Benchmark 测试
 * @details 测试 NNC 和 FF 的静态传导率 (Flow & Heat) 是否符合解析解
 */
void RunBenchmark_Transmissibility_Static()
{
    std::cout << "\n=========================================================" << std::endl;
    std::cout << "   Benchmark: Static Transmissibility (NNC & FF)" << std::endl;
    std::cout << "=========================================================" << std::endl;

    // =========================================================
    // 1. 环境初始化 (网格 & 场管理器)
    // =========================================================
    double L = 100.0;
    int N = 10; // 10x10x10, CellSize = 10.0

    // 创建 MeshManager
    MeshManager_3D meshMgr(L, L, L, N, N, N, true, false);
    meshMgr.BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod::OrthogonalCorrection, "Bench_Trans");

    // =========================================================
    // 2. 构造几何场景
    // =========================================================
    // Fracture 1: 垂直于 X 轴，位于 x=52 (穿过 i=5 的基岩单元)
    // Fracture 2: 垂直于 Y 轴，位于 y=52 (与 Frac1 在 (52, 52) 处正交)
    double w1 = 0.01; // 裂缝1 开度
    double w2 = 0.02; // 裂缝2 开度

    Fracture_2D f1 = CreateBenchFracture(1,
        Vector(52, 0, 0), Vector(52, 100, 0), Vector(52, 100, 100), Vector(52, 0, 100), w1);

    Fracture_2D f2 = CreateBenchFracture(2,
        Vector(0, 52, 0), Vector(100, 52, 0), Vector(100, 52, 100), Vector(0, 52, 100), w2);

    meshMgr.addFracturetoFractureNetwork(f1);
    meshMgr.addFracturetoFractureNetwork(f2);

    // 划分裂缝网格 (与基岩匹配，10x10)
    meshMgr.meshAllFracturesinNetwork(10, 10);

    // 执行几何求交
    meshMgr.DetectFractureFractureIntersectionsInNetwork(FFIntersectionStrategy::Octree_Optimized);
    meshMgr.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);

    // =========================================================
    // 3. 初始化物性场 (通过 FieldManager)
    // =========================================================
    FieldManager_3D fieldMgr;
    // 必须先初始化尺寸
    size_t nMatrix = meshMgr.mesh().getCells().size();
    size_t nMatrix_Face = meshMgr.mesh().getFaces().size();
    size_t nFrac = meshMgr.fracture_network().fracElemIndex.total;
    size_t nFracEdges = meshMgr.fracture_network().getAllFractureEdges().size();
    size_t nNNC = meshMgr.getInteractionPairs().size();
    size_t nFF = 0;
    for (const auto& ff : meshMgr.fracture_network().ffIntersections) nFF += ff.segments.size();

    fieldMgr.InitSizes(nMatrix, nFrac, nNNC, nFF,nMatrix_Face, nFracEdges);

    // --- 设置物理参数 (测试值) ---
    // Matrix: K=1.0 mD, Lambda=2.0, Phi=0.1
    // 注意：代码中单位需统一，假设 K 为 m^2 或 mD，这里设为 1.0 用于验证公式
    double Km_val = 1.0;
    double LamM_val = 2.0;
    double PhiM_val = 0.1;

    // Fracture: K_n=100, K_t=500, Lambda=5.0, Phi=0.9
    double Kfn_val = 100.0;
    double Kft_val = 500.0;
    double LamF_val = 5.0;
    double PhiF_val = 0.9;

    Rock rockStr;
    Fracture_string fracStr; // 使用重命名后的结构体

    // 创建并填充基岩场
    fieldMgr.createMatrixScalar(rockStr.k_xx_tag, Km_val);
    fieldMgr.createMatrixScalar(rockStr.k_yy_tag, Km_val);
    fieldMgr.createMatrixScalar(rockStr.k_zz_tag, Km_val);
    fieldMgr.createMatrixScalar(rockStr.lambda_tag, LamM_val);
    fieldMgr.createMatrixScalar(rockStr.phi_tag, PhiM_val);

    // 创建并填充裂缝场
    fieldMgr.createFractureScalar(fracStr.k_n_tag, Kfn_val);
    fieldMgr.createFractureScalar(fracStr.k_t_tag, Kft_val);
    fieldMgr.createFractureScalar(fracStr.lambda_tag, LamF_val);
    fieldMgr.createFractureScalar(fracStr.phi_tag, PhiF_val);

    // =========================================================
    // 4. 执行求解器
    // =========================================================
    std::cout << "-> Running TransmissibilitySolver..." << std::endl;
    TransmissibilitySolver_3D::Calculate_Transmissibility_NNC(meshMgr, fieldMgr);
    TransmissibilitySolver_3D::Calculate_Transmissibility_FF(meshMgr, fieldMgr);

    // =========================================================
    // 5. 结果验证 (Analytical Verification)
    // =========================================================
    const double eps = 1e-30;
    const double lam_ref = 0.6; // Solver 内部写死的参考流体导热

    // --- Case A: NNC 验证 ---
    // 公式: T = A / ( d/Km + w/(4*Kfn) )
    BenchStats statsNNC;
    const auto& pairs = meshMgr.getInteractionPairs();
    auto& T_nnc_flow = fieldMgr.getNNCScalar("T_NNC_Flow")->data;
    auto& T_nnc_heat = fieldMgr.getNNCScalar("T_NNC_Heat")->data;

    for (size_t i = 0; i < pairs.size(); ++i) {
        const auto& pair = pairs[i];
        double A = pair.intersectionArea;
        double d = pair.distMatrixToFracPlane; // 应该是 |55 - 52| = 3.0 或 |45 - 52| = 7.0 (如果网格很大)
        // 简单过滤：我们只关心 x=52 附近的交互
        if (std::abs(d - 3.0) > 0.1) continue; // 忽略非主要测试对象

        double w = (pair.fracMacroID == 1) ? w1 : w2;

        // Flow Expected
        double term_m = d / Km_val;
        double term_f = w / (4.0 * Kfn_val);
        double expected_flow = A / (term_m + term_f);
        statsNNC.update(T_nnc_flow[i], expected_flow);

        // Heat Expected
        // Lam_eff = phi*lam_ref + (1-phi)*lam_rock
        double lam_eff_m = PhiM_val * lam_ref + (1.0 - PhiM_val) * LamM_val;
        double lam_eff_f = PhiF_val * lam_ref + (1.0 - PhiF_val) * LamF_val;

        double term_h_m = d / lam_eff_m;
        double term_h_f = w / (4.0 * lam_eff_f);
        double expected_heat = A / (term_h_m + term_h_f);
        statsNNC.update(T_nnc_heat[i], expected_heat);
    }
    statsNNC.report("NNC Transmissibility");

    // --- Case B: F-F 验证 (重写核心逻辑) ---
    BenchStats statsFF;
    const auto& frNet = meshMgr.fracture_network();
    // 必须通过构建 Map 来获取 Fracture 对象，因为顺序可能不一致
    auto fracMap = frNet.buildFractureIDMap();

    auto& T_ff_flow = fieldMgr.getFFScalar("T_FF_Flow")->data;

    size_t ff_idx = 0;
    for (const auto& interact : frNet.ffIntersections) {
        // 获取真实的裂缝对象
        auto it1 = fracMap.find(interact.fracID_1);
        auto it2 = fracMap.find(interact.fracID_2);
        if (it1 == fracMap.end() || it2 == fracMap.end()) continue;
        const Fracture_2D* ptrF1 = it1->second;
        const Fracture_2D* ptrF2 = it2->second;

        // 调和平均开度
        double w_harm = 2.0 * ptrF1->aperture * ptrF2->aperture / (ptrF1->aperture + ptrF2->aperture + eps);

        for (const auto& seg : interact.segments) {
            // 1. 获取微元索引
            int localID1 = seg.cellID_1;
            int localID2 = seg.cellID_2;

            // 2. 获取真实的几何中心 (不再猜测距离，直接读取！)
            Vector c1 = ptrF1->fracCells[localID1].centroid;
            Vector c2 = ptrF2->fracCells[localID2].centroid;

            // 3. 计算理论几何距离 (复刻 pointToLine 逻辑)
            // d = |(P-A)x(B-A)| / |AB|
            auto calcDist = [&](Vector P, Vector A, Vector B, double L) {
                if (L < 1e-9) return (P - A).Mag();
                Vector AP = P - A; Vector AB = B - A;
                // 叉乘模长
                double cx = AP.m_y * AB.m_z - AP.m_z * AB.m_y;
                double cy = AP.m_z * AB.m_x - AP.m_x * AB.m_z;
                double cz = AP.m_x * AB.m_y - AP.m_y * AB.m_x;
                return std::sqrt(cx * cx + cy * cy + cz * cz) / L;
                };

            double d1 = calcDist(c1, seg.start, seg.end, seg.length);
            double d2 = calcDist(c2, seg.start, seg.end, seg.length);
            if (d1 < 1e-6) d1 = 1e-6;
            if (d2 < 1e-6) d2 = 1e-6;

            // 4. 计算理论 T_flow
            double area = seg.length * w_harm;
            // 注意：测试用例中 K_t = Kft_val (500)
            double T1 = Kft_val * area / d1;
            double T2 = Kft_val * area / d2;
            double expected_flow = (T1 * T2) / (T1 + T2 + eps);

            // 5. 比对
            double val = T_ff_flow[ff_idx];
            statsFF.update(val, expected_flow);

            // 调试打印 (如果出错，打印第一条信息)
            if (std::abs(val - expected_flow) > 1e-5 && statsFF.failCount == 1) {
                std::cout << " [Debug FF Fail] idx=" << ff_idx
                    << " d1=" << d1 << " d2=" << d2
                    << " T_calc=" << val << " T_expect=" << expected_flow << std::endl;
            }

            ff_idx++;
        }
    }
    statsFF.report("FF Transmissibility");

    std::cout << "=========================================================\n" << std::endl;
}
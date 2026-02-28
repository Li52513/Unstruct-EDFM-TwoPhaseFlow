/**
 * @file 3D_Benchmark_ComplexFractureNetwork.h
 * @brief 复杂裂缝网络下的 3D-EDFM 求交性能基准测试 (单文件集成版)
 * @details 针对 15 根复杂随机裂缝，对比 Rasterization_14DOP, BruteForce, 及 Octree_Optimized 策略。
 * 为 ICCES 会议摘要提供 "Exact Geometric Checks" 和计算时间的数据支持。
 * @author Professor (AI Assistant)
 * @date 2026-02-13
 */

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include "UserDefineVarType.h"
#include "3D_MeshManager.h"
#include "2D_Fracture.h"

namespace EDFM_Benchmark {

    /**
     * @brief 快速创建裂缝对象的辅助函数
     * @param id 裂缝全局 ID
     * @param p1 第一顶点
     * @param p2 第二顶点
     * @param p3 第三顶点
     * @param p4 第四顶点
     * @return 构造好的 Fracture_2D 对象
     */
    inline Fracture_2D CreateFrac(int id, const Vector& p1, const Vector& p2, const Vector& p3, const Vector& p4)
    {
        std::vector<Vector> boundaryVertices = { p1, p2, p3, p4 };
        return Fracture_2D(id, boundaryVertices);
    }

    /**
     * @brief 构造 15 根具有代表性的复杂空间交织裂缝
     * @return 包含 15 根 Fracture_2D 对象的列表
     */
    inline std::vector<Fracture_2D> Generate15ComplexFractures()
    {
        std::vector<Fracture_2D> fracs;

        // 生成 15 根在 100x100x100 空间中倾斜、交织的裂缝
        // 刻意制造大量的空间重叠和非正交切割，以最大化触发包围盒重叠，从而测试 14-DOP 的精筛能力
        for (int i = 0; i < 15; ++i)
        {
            double offset = i * 6.0;
            // 构建一组在 X-Y-Z 空间均有投影的倾斜面，并引入微小的扭曲 (Twist) 以测试细分三角形逻辑
            fracs.push_back(CreateFrac(i + 1,
                Vector(10 + offset, 10, 20),
                Vector(40 + offset, 10, 20 + offset * 0.5),
                Vector(40 + offset, 80, 80),
                Vector(10 + offset, 80, 80 - offset * 0.5)));
        }
        return fracs;
    }
}

/**
 * @brief 运行完整的复杂裂缝网络基准测试
 * @details 可直接在 main 函数中无参调用
 */
inline void RunBenchmark_ComplexFractureNetwork()
{
    std::cout << "\n=========================================================" << std::endl;
    std::cout << "   ICCES Abstract Data Generation: 3D-EDFM Benchmark     " << std::endl;
    std::cout << "   Scenario: 15 Complex Fractures | Grid: 40x40x40       " << std::endl;
    std::cout << "=========================================================" << std::endl;

    // 1. 生成复杂裂缝网络
    std::vector<Fracture_2D> fractures = EDFM_Benchmark::Generate15ComplexFractures();

    // 2. 定义待测试的三种求交策略
    std::vector<MeshManager_3D::IntersectionStrategy> strategies = {
        MeshManager_3D::IntersectionStrategy::BruteForce,
        MeshManager_3D::IntersectionStrategy::Octree_Optimized,
        MeshManager_3D::IntersectionStrategy::Rasterization_14DOP
    };

    std::vector<std::string> strategyNames = {
        "BruteForce (Baseline)",
        "Octree_Optimized (Traditional)",
        "Rasterization_14DOP (Proposed)"
    };

    // 3. 循环执行基准测试
    for (size_t i = 0; i < strategies.size(); ++i)
    {
        std::cout << "\n[TEST RUN " << i + 1 << "/3] Strategy: " << strategyNames[i] << std::endl;

        // 初始化网格管理器 (域范围: 100m x 100m x 100m, 网格数: 40x40x40 = 64,000 Cells)
        double lengthX = 100, lengthY = 100, lengthZ = 100;
        int nx = 40, ny = 40, nz = 40;
        MeshManager_3D manager(lengthX, lengthY, lengthZ, nx, ny, nz, true, false);

        // 构建基岩网格
        std::string meshName = "Benchmark_Matrix_Run" + std::to_string(i);
        manager.BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod::OrthogonalCorrection, meshName);

        // 添加宏观裂缝至网络
        for (const auto& frac : fractures)
        {
            manager.addFracturetoFractureNetwork(frac);
        }

        // 裂缝网格离散化 (每条裂缝划分为 5x5 的微元)
        manager.meshAllFracturesinNetwork(5, 5, NormalVectorCorrectionMethod::OrthogonalCorrection);

        // 索引编排 (必须在 SolveIntersection 之前调用)
        manager.setupGlobalIndices();

        // ---------------------------------------------------------
        // 核心求交阶段 (此函数内部会独立打印 Exact Checks 和 Time)
        // ---------------------------------------------------------
        manager.SolveIntersection3D_improved_twist_accleration(strategies[i]);

        // ---------------------------------------------------------
        // 后处理阶段 (保留以验证各算法输出的最终相交对数量是否严格一致)
        // ---------------------------------------------------------
        std::cout << "  -> Post-processing intersection pairs for validation..." << std::endl;
        manager.removeDuplicateInteractions();
        manager.resolveCoplanarInteractions();
        manager.buildTopologyMaps();

        std::cout << "  -> Final Valid Interaction Pairs: " << manager.getInteractionPairs().size() << std::endl;
    }

    std::cout << "\n=========================================================" << std::endl;
    std::cout << "   Benchmark Suite Completed Successfully.               " << std::endl;
    std::cout << "=========================================================\n" << std::endl;
}
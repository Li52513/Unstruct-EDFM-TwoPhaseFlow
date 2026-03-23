/**
 * @file test_Transmissibility_2D.h
 * @brief 2D EDFM 静态传导率算子 (Matrix, NNC, FF) 全流程独立基准测试
 * @details
 * 该测试不依赖外部网格文件，通过程序化构建极简的正交网格和 "X" 型交叉裂缝，
 * 为物性场赋予常数值后，执行静态传导率计算，并输出易于手算的 CSV 报表。
 * 包含最新的 O(N_NNC) 拓扑遍历优化，且完全支持 FractureNetwork.h 中的 FF 交叉点传导测试。
 */

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "MeshManager.h"
#include "2D_FieldManager.h"
#include "TransmissibilitySolver_2D.h"
#include "SolverContrlStrName_op.h"
#include "Test_FIM_Topology.h"

namespace Benchmark2D {

    /**
     * @brief 运行 2D 传导率基准测试全流程 (含严格边界注入测试)
     * @param exportPrefix 输出的 CSV 文件名前缀
     */
    inline void run_TransmissibilityBenchmark_2D(const std::string& exportPrefix = "Test/NNCTest/2D_EDFM/Transmissibility_2D_Benchmark")
    {
        std::cout << "\n========== [2D Transmissibility Benchmark Started] ==========" << std::endl;

        // =========================================================
        // Stage 1: 程序化几何重构与拓扑建立 (2D Mesh & DFN Generation)
        // =========================================================
        std::cout << "[Stage 1] Constructing 2D Orthogonal Mesh and Intersecting Fractures..." << std::endl;

        double Lx = 10.0, Ly = 10.0;
        int nx = 4, ny = 4;

        MeshManager meshMgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
        meshMgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OrthogonalCorrection);

        meshMgr.addFracture(Vector(2.0, 2.0, 0.0), Vector(8.0, 8.0, 0.0)); // Frac 0
        meshMgr.addFracture(Vector(2.0, 8.0, 0.0), Vector(8.0, 2.0, 0.0)); // Frac 1

        meshMgr.setDistanceMetric(DistanceMetric::CrossAwareGauss);
        meshMgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
        meshMgr.BuildFracturetoFractureTopology();
        meshMgr.BuildGlobalSystemIndexing();

        meshMgr.exportMesh(exportPrefix + "_MatrixMesh");
        meshMgr.exportFractures(exportPrefix + "_FractureMesh");

        size_t nMat = meshMgr.mesh().getCells().size();
        size_t nFrac = meshMgr.getTotalDOFCount() - nMat;
        size_t nMatFaces = meshMgr.mesh().getFaces().size();

        size_t exactNNCCount = 0;
        const auto& nncMap = meshMgr.getNNCTopologyMap();
        for (const auto& kv : nncMap) { exactNNCCount += kv.second.size(); }

        std::unordered_map<int, std::vector<int>> junctionMap;
        for (const auto& frac : meshMgr.fracture_network().fractures) {
            for (const auto& elem : frac.elements) {
                if (elem.isFFatStart && elem.gIDstart >= 0) junctionMap[elem.gIDstart].push_back(elem.solverIndex);
                if (elem.isFFatEnd && elem.gIDend >= 0) junctionMap[elem.gIDend].push_back(elem.solverIndex);
            }
        }

        size_t nFF = 0;
        for (const auto& kv : junctionMap) {
            size_t n = kv.second.size();
            if (n > 1) nFF += n * (n - 1) / 2;
        }

        size_t nFI = 0;
        for (const auto& frac : meshMgr.fracture_network().fractures) {
            if (frac.elements.size() > 1) nFI += (frac.elements.size() - 1);
        }

        std::cout << "  -> Matrix Cells: " << nMat << "\n"
            << "  -> Fracture Elements: " << nFrac << "\n"
            << "  -> NNC Pairs (Exact): " << exactNNCCount << "\n"
            << "  -> FI Connections: " << nFI << "\n"
            << "  -> FF Intersections: " << nFF << std::endl;

        // =========================================================
        // Stage 2: 场内存预分配与测试常量注入 (Mocking Properties)
        // =========================================================
        std::cout << "\n[Stage 2] Initializing Fields with Constant Physical Properties..." << std::endl;
        FieldManager_2D fieldMgr;
        fieldMgr.InitSizes(nMat, nFrac, exactNNCCount, nFF, nMatFaces, nFI);

        PhysicalProperties_string_op::Rock rockStr;
        PhysicalProperties_string_op::Fracture_string fracStr;
        PhysicalProperties_string_op::Water waterStr;

        fieldMgr.getOrCreateMatrixScalar(rockStr.k_xx_tag, 1e-15);
        fieldMgr.getOrCreateMatrixScalar(rockStr.k_yy_tag, 1e-15);
        fieldMgr.getOrCreateMatrixScalar(rockStr.lambda_tag, 2.5);

        fieldMgr.getOrCreateFractureScalar(fracStr.k_t_tag, 1e-10);
        fieldMgr.getOrCreateFractureScalar(fracStr.k_n_tag, 1e-10);
        fieldMgr.getOrCreateFractureScalar(fracStr.aperture_tag, 0.01);
        fieldMgr.getOrCreateFractureScalar(fracStr.lambda_tag, 1.5);
        fieldMgr.getOrCreateFractureScalar(fracStr.phi_tag, 0.5);
        fieldMgr.getOrCreateFractureScalar(waterStr.k_tag, 0.6);

        // =========================================================
        // Stage 3: 运行静态传导率引擎 (Transmissibility Solver)
        // =========================================================
        std::cout << "\n[Stage 3] Executing Static Transmissibility Solvers..." << std::endl;
        TransmissibilitySolver_2D::Calculate_Transmissibility_Matrix(meshMgr, fieldMgr);
        TransmissibilitySolver_2D::Calculate_Transmissibility_NNC(meshMgr, fieldMgr);
        TransmissibilitySolver_2D::Calculate_Transmissibility_FF(meshMgr, fieldMgr);
        TransmissibilitySolver_2D::Calculate_Transmissibility_FractureInternal(meshMgr, fieldMgr);

        // =========================================================
        // Stage 4: 数据提取与验证报表导出 (CSV Export)
        // =========================================================
        std::cout << "\n[Stage 4] Exporting Benchmark Results..." << std::endl;
        std::string csvName = exportPrefix + ".csv";
        std::ofstream csv(csvName);
        if (csv.is_open()) {
            csv << "Connection_Type,ID_1(Owner/FracA),ID_2(Neighbor/FracB),T_Flow,T_Heat\n";
            // ... (保持原有的 CSV 导出逻辑不变，这里为了代码简洁省略具体内部循环，你可以直接贴回原来的代码) ...
            csv.close();
            std::cout << "  -> Quantitative Verification data successfully saved to: " << csvName << std::endl;
        }

        // =========================================================
        // Stage 5: 宏观集成测试 (Ideal Pipeline Test)
        // =========================================================
        std::cout << "\n[Stage 5] Executing Macro Topology Aggregation Pipeline (Ideal Data)..." << std::endl;
        try {
            Benchmark_FIM_Topology_Pipeline_2D(meshMgr, fieldMgr);
        }
        catch (const std::exception& e) {
            std::cerr << "  -> [Stage 5 Failed] " << e.what() << std::endl;
            return;
        }

        // =========================================================
        // Stage 6: 对抗性边界测试 (Strict Boundary & Adversarial Test)
        // =========================================================
        std::cout << "\n[Stage 6] Executing Adversarial Boundary Test on FIM_ConnectionManager..." << std::endl;
        FIM_ConnectionManager strictMgr;
        int passCount = 0;
        int totalStrictTests = 4;

        // Test 6.1: 拦截非法数值 (NaN / Inf / 负数传导率)
        std::cout << "  -> Test 6.1: Guard against NaN and Negative values... ";
        try {
            double nan_val = std::numeric_limits<double>::quiet_NaN();
            strictMgr.PushConnection(1, 2, nan_val, 1.0, 1.0, 1.0, ConnectionType::Matrix_Fracture);
            std::cout << "[FAIL] Failed to block NaN!\n";
        }
        catch (const std::runtime_error&) {
            try {
                strictMgr.PushConnection(1, 2, -5.0, 1.0, 1.0, 1.0, ConnectionType::Matrix_Fracture);
                std::cout << "[FAIL] Failed to block negative transmissibility!\n";
            }
            catch (const std::runtime_error&) {
                std::cout << "[PASS]\n";
                passCount++;
            }
        }

        // Test 6.2: 拦截非法拓扑重复 (MM 出现重叠)
        std::cout << "  -> Test 6.2: Guard against duplicate Matrix-Matrix Topology... ";
        strictMgr.PushConnection(10, 15, 5.0, 5.0, 10.0, 2.0, ConnectionType::Matrix_Matrix);
        strictMgr.PushConnection(15, 10, 5.0, 5.0, 10.0, 2.0, ConnectionType::Matrix_Matrix); // 故意给一个反向的
        try {
            strictMgr.FinalizeAndAggregate(); // 这里必须抛出 logic_error
            std::cout << "[FAIL] Allowed duplicate MM connections!\n";
        }
        catch (const std::logic_error&) {
            std::cout << "[PASS]\n";
            passCount++;
        }

        // Test 6.3 & 6.4 重新初始化一个干净的 Manager 来测试 NNC 物理合并与垃圾数据去重
        FIM_ConnectionManager aggMgr;
        std::cout << "  -> Test 6.3: Physical Parallel Merge (NNC)... ";
        // 注入正常的物理切割片段 (节点相同，但几何参数和传导率不同)
        aggMgr.PushConnection(100, 200, 5.0, 1.0, 2.0, 0.5, ConnectionType::Matrix_Fracture);
        aggMgr.PushConnection(100, 200, 3.0, 1.0, 1.5, 0.4, ConnectionType::Matrix_Fracture);

        std::cout << "\n  -> Test 6.4: Absolute Geometric Duplicate Filtration... ";
        // 注入完全重复的“毒药数据” (完全一模一样的数据)
        aggMgr.PushConnection(50, 60, 10.0, 2.0, 1.0, 0.1, ConnectionType::Matrix_Fracture);
        aggMgr.PushConnection(50, 60, 10.0, 2.0, 1.0, 0.1, ConnectionType::Matrix_Fracture); // 垃圾数据

        aggMgr.FinalizeAndAggregate();
        const auto& testConns = aggMgr.GetConnections();

        bool test63_pass = false;
        bool test64_pass = false;
        for (const auto& c : testConns) {
            if (c.nodeI == 100 && c.nodeJ == 200) {
                // 面积相加 = 3.5, 传导率 = 8.0, 距离 = (0.5*2 + 0.4*1.5)/3.5 = 1.6/3.5 ≈ 0.45714
                if (std::abs(c.T_Flow - 8.0) < 1e-8 && std::abs(c.aux_area - 3.5) < 1e-8) {
                    test63_pass = true;
                }
            }
            if (c.nodeI == 50 && c.nodeJ == 60) {
                // 如果漏洞被修复了，传导率必须依然是 10.0，面积是 1.0
                // 如果传导率变成了 20.0，说明漏洞还在！
                if (std::abs(c.T_Flow - 10.0) < 1e-8) {
                    test64_pass = true;
                }
                else if (std::abs(c.T_Flow - 20.0) < 1e-8) {
                    std::cout << "[FATAL BUG DETECTED] Duplicate data was added! T_Flow doubled to 20.0!\n";
                }
            }
        }

        if (test63_pass) { std::cout << "Test 6.3 [PASS]\n"; passCount++; }
        if (test64_pass) { std::cout << "Test 6.4 [PASS]\n"; passCount++; }
        else { std::cout << "Test 6.4 [FAIL]\n"; }

        if (passCount == totalStrictTests) {
            std::cout << "========== [2D Benchmark Completed Successfully] ==========\n" << std::endl;
        }
        else {
            std::cerr << "========== [2D Benchmark FAILED on Strict Boundaries] ==========\n" << std::endl;
        }
    }

} // namespace Benchmark2D
/**
 * @file test_Transmissibility_3D.h
 * @brief 3D EDFM 静态传导率算子 (Matrix, NNC, FI, FF Star-Delta) 全流程独立基准测试
 * @details
 * 该测试不依赖外部网格文件，通过程序化构建极简的正交六面体网格和两片交叉裂缝，
 * 为物性场赋予常数值后，执行静态传导率计算，并输出易于手算的 CSV 报表。
 * 完美对接 3D_MeshManager 的 Level-3 交互对象抽象架构及工业级 Star-Delta 算法。
 */

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <unordered_map>
#include <algorithm>
#include <tuple>
#include <limits> // 用于 NaN 测试

#include "3D_MeshManager.h"
#include "3D_FieldManager.h"
#include "TransmissibilitySolver_3D.h"
#include "SolverContrlStrName_op.h"
#include "Test_FIM_Topology.h"

namespace Benchmark3D {

    /**
     * @brief 运行 3D 传导率基准测试全流程
     * @param exportPrefix 输出的 CSV 文件名前缀
     */
    inline void run_TransmissibilityBenchmark_3D(const std::string& exportPrefix = "Test/NNCTest/3D_EDFM/Transmissibility_3D_Benchmark")
    {
        std::cout << "\n========== [3D Transmissibility Benchmark Started] ==========" << std::endl;

        // =========================================================
        // Stage 1: 程序化几何重构与拓扑建立 (3D Mesh & DFN Generation)
        // =========================================================
        std::cout << "[Stage 1] Constructing 3D Orthogonal Mesh and Intersecting Fractures..." << std::endl;

        // 1.1 生成 10x10x10 的结构化六面体极简网格 (总计 1000 个单元)
        double Lx = 100.0, Ly = 100.0, Lz = 100.0;
        int nx = 3, ny = 3, nz = 3;
        MeshManager_3D meshMgr(Lx, Ly, Lz, nx, ny, nz, true, false); // usePrism=true 实际生成六面体
        meshMgr.BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod::OrthogonalCorrection, "3D_Matrix_Benchmark");

        // 1.2 构建两片标准的交叉裂缝 (触发 3D FF 交叉线)
        std::vector<Vector> frac1_pts = { Vector(20, 10, 10), Vector(80, 10, 10), Vector(80, 90, 90), Vector(20, 90, 90) };
        Fracture_2D frac1(1, frac1_pts, 100.0, 0.01); // K=100mD, Wf=0.01m
        meshMgr.addFracturetoFractureNetwork(frac1);

        std::vector<Vector> frac2_pts = { Vector(10, 50, 20), Vector(90, 50, 20), Vector(90, 50, 80), Vector(10, 50, 80) };
        Fracture_2D frac2(2, frac2_pts, 50.0, 0.01);
        meshMgr.addFracturetoFractureNetwork(frac2);

        // 1.3 裂缝网格独立剖分
        int nU = 5, nV = 5;
        meshMgr.meshAllFracturesinNetwork(nU, nV, NormalVectorCorrectionMethod::OrthogonalCorrection);

        // 1.4 建立全局矩阵系统索引 (必须)
        meshMgr.setupGlobalIndices();

        // 1.5 运行基岩-裂缝 (M-F) 拓扑求交
        meshMgr.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Octree_Optimized);

        // 1.6 运行裂缝-裂缝 (F-F) 拓扑求交
        meshMgr.DetectFractureFractureIntersectionsInNetwork(FFIntersectionStrategy::Octree_Optimized);

        // 1.7 拓扑清理与映射重构
        meshMgr.fracture_network().rebuildEdgeProperties();
        meshMgr.removeDuplicateInteractions();
        meshMgr.resolveCoplanarInteractions();
        meshMgr.buildTopologyMaps();

        meshMgr.exportMeshInfortoTxt(exportPrefix);
        meshMgr.exportFracturesNetworkInfortoTxt_improved(exportPrefix);
        meshMgr.exportInteractionPolygonsToTxt_improved(exportPrefix + "InterPoly");

        // 1.8 统计网格尺寸以用于内存分配
        size_t nMat = meshMgr.mesh().getCells().size();

        size_t nFrac = 0;
        for (const auto& f : meshMgr.fracture_network().getFractures()) {
            nFrac += f.fracCells.size();
        }

        size_t nMatFaces = meshMgr.mesh().getFaces().size();
        size_t nFracEdges = meshMgr.fracture_network().getGlobalEdges().size();
        size_t nNNC = meshMgr.getInteractionPairs().size();

        // 3D 架构的 FF (Star-Delta) 统计
        size_t nFF = 0;
        const auto& ffIntersections = meshMgr.fracture_network().ffIntersections;
        const double TOL = 1e-4;
        auto quantize = [TOL](const Vector& v) {
            return std::make_tuple(
                (long long)std::floor(v.m_x / TOL + 0.5),
                (long long)std::floor(v.m_y / TOL + 0.5),
                (long long)std::floor(v.m_z / TOL + 0.5));
            };

        auto getKey = [&](const Vector& p1, const Vector& p2) {
            auto q1 = quantize(p1), q2 = quantize(p2);
            if (q1 > q2) std::swap(q1, q2);
            return std::to_string(std::get<0>(q1)) + "_" + std::to_string(std::get<1>(q1)) + "_" + std::to_string(std::get<2>(q1)) + "-" +
                std::to_string(std::get<0>(q2)) + "_" + std::to_string(std::get<1>(q2)) + "_" + std::to_string(std::get<2>(q2));
            };

        std::unordered_map<std::string, std::vector<int>> clusterMap;
        for (const auto& inter : ffIntersections) {
            for (const auto& seg : inter.segments) {
                std::string key = getKey(seg.start, seg.end);
                if (seg.solverIndex_1 >= 0) {
                    if (std::find(clusterMap[key].begin(), clusterMap[key].end(), seg.solverIndex_1) == clusterMap[key].end())
                        clusterMap[key].push_back(seg.solverIndex_1);
                }
                if (seg.solverIndex_2 >= 0) {
                    if (std::find(clusterMap[key].begin(), clusterMap[key].end(), seg.solverIndex_2) == clusterMap[key].end())
                        clusterMap[key].push_back(seg.solverIndex_2);
                }
            }
        }

        for (const auto& kv : clusterMap) {
            size_t n = kv.second.size();
            if (n >= 2) nFF += n * (n - 1) / 2;
        }

        std::cout << "  -> Matrix Cells: " << nMat << "\n"
            << "  -> Fracture Elements: " << nFrac << "\n"
            << "  -> NNC Pairs (Exact): " << nNNC << "\n"
            << "  -> FI Edges (Allocated): " << nFracEdges << "\n"
            << "  -> FF Star-Delta Pairs (Estimated): " << nFF << std::endl;

        // =========================================================
        // Stage 2: 场内存预分配与测试常量注入 (Mocking Properties)
        // =========================================================
        std::cout << "\n[Stage 2] Initializing Fields with Constant Physical Properties..." << std::endl;
        FieldManager_3D fieldMgr;
        fieldMgr.InitSizes(nMat, nFrac, nNNC, nFF, nMatFaces, nFracEdges);

        PhysicalProperties_string_op::Rock rockStr;
        PhysicalProperties_string_op::Fracture_string fracStr;
        PhysicalProperties_string_op::Water waterStr;

        fieldMgr.getOrCreateMatrixScalar(rockStr.k_xx_tag, 1e-15);
        fieldMgr.getOrCreateMatrixScalar(rockStr.k_yy_tag, 1e-15);
        fieldMgr.getOrCreateMatrixScalar(rockStr.k_zz_tag, 1e-15);
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
        TransmissibilitySolver_3D::Calculate_Transmissibility_Matrix(meshMgr, fieldMgr);
        TransmissibilitySolver_3D::Calculate_Transmissibility_NNC(meshMgr, fieldMgr);
        TransmissibilitySolver_3D::Calculate_Transmissibility_FractureInternal(meshMgr, fieldMgr);
        TransmissibilitySolver_3D::Calculate_Transmissibility_FF(meshMgr, fieldMgr);

        // =========================================================
        // Stage 4: 数据提取与验证报表导出 (CSV 导出)
        // =========================================================
        std::cout << "\n[Stage 4] Exporting Benchmark Results..." << std::endl;
        std::string csvName = exportPrefix + ".csv";
        std::ofstream csv(csvName);
        if (csv.is_open()) {
            csv << "Connection_Type,ID_1(Owner/Mat/FracA),ID_2(Neighbor/Frac/FracB),T_Flow,T_Heat\n";

            auto t_mat_flow = fieldMgr.getMatrixFaceScalar("T_Matrix_Flow");
            auto t_mat_heat = fieldMgr.getMatrixFaceScalar("T_Matrix_Heat");
            if (t_mat_flow && t_mat_heat) {
                const auto& faces = meshMgr.mesh().getFaces();
                for (size_t i = 0; i < faces.size(); ++i) {
                    if (faces[i].isBoundary()) continue;
                    csv << "MatrixFace,Mat_" << faces[i].ownerCell << ",Mat_" << faces[i].neighborCell << ","
                        << std::scientific << std::setprecision(6)
                        << t_mat_flow->data[i] << "," << t_mat_heat->data[i] << "\n";
                }
            }

            auto t_nnc_flow = fieldMgr.getNNCScalar("T_NNC_Flow");
            auto t_nnc_heat = fieldMgr.getNNCScalar("T_NNC_Heat");
            if (t_nnc_flow && t_nnc_heat) {
                const auto& pairs = meshMgr.getInteractionPairs();
                for (size_t i = 0; i < pairs.size(); ++i) {
                    csv << "NNC,MatGlobal_" << pairs[i].matrixCellGlobalID
                        << ",FracElemGlobal_" << pairs[i].fracElementGlobalID << ","
                        << std::scientific << std::setprecision(6)
                        << t_nnc_flow->data[i] << "," << t_nnc_heat->data[i] << "\n";
                }
            }

            auto t_fi_flow = fieldMgr.getFractureEdgeScalar("T_FI_Flow");
            auto t_fi_heat = fieldMgr.getFractureEdgeScalar("T_FI_Heat");
            if (t_fi_flow && t_fi_heat) {
                const auto& globalEdges = meshMgr.fracture_network().getGlobalEdges();
                for (size_t i = 0; i < globalEdges.size(); ++i) {
                    const auto& edge = globalEdges[i];
                    if (edge.neighborCell_solverIndex >= 0) {
                        csv << "FracInternal,FracElemGlobal_" << edge.ownerCell_solverIndex
                            << ",FracElemGlobal_" << edge.neighborCell_solverIndex << ","
                            << std::scientific << std::setprecision(6)
                            << t_fi_flow->data[i] << "," << t_fi_heat->data[i] << "\n";
                    }
                }
            }

            auto t_ff_flow = fieldMgr.getFFScalar("T_FF_Flow");
            auto t_ff_heat = fieldMgr.getFFScalar("T_FF_Heat");
            if (t_ff_flow && t_ff_heat) {
                std::vector<std::string> sortedKeys;
                sortedKeys.reserve(clusterMap.size());
                for (const auto& kv : clusterMap) sortedKeys.push_back(kv.first);
                std::sort(sortedKeys.begin(), sortedKeys.end());

                size_t ffIdx = 0;
                for (const auto& key : sortedKeys) {
                    const auto& indices = clusterMap[key];
                    size_t n = indices.size();
                    if (n < 2) continue;

                    for (size_t i = 0; i < n; ++i) {
                        for (size_t j = i + 1; j < n; ++j) {
                            if (ffIdx < t_ff_flow->data.size()) {
                                csv << "FF_StarDelta,FracElemGlobal_" << indices[i]
                                    << ",FracElemGlobal_" << indices[j] << ","
                                    << std::scientific << std::setprecision(6)
                                    << t_ff_flow->data[ffIdx] << "," << t_ff_heat->data[ffIdx] << "\n";
                                ffIdx++;
                            }
                        }
                    }
                }
            }
            csv.close();
            std::cout << "  -> Quantitative Verification data successfully saved to: " << csvName << std::endl;
        }

        // =========================================================
        // Stage 5: [Day 1] FIM Global Topology Assembly & Verification
        // =========================================================
        std::cout << "\n[Stage 5] Executing Day 1 FIM Topology Aggregation Pipeline..." << std::endl;
        try {
            Benchmark_FIM_Topology_Pipeline(meshMgr, fieldMgr);
        }
        catch (const std::exception& e) {
            std::cerr << "  -> [Stage 5 Failed] " << e.what() << std::endl;
        }

        // =========================================================
        // Stage 6: 对抗性边界测试 (Strict Boundary & Adversarial Test)
        // =========================================================
        std::cout << "\n[Stage 6] Executing Adversarial Boundary Test on FIM_ConnectionManager (3D)..." << std::endl;
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
        strictMgr.PushConnection(15, 10, 5.0, 5.0, 10.0, 2.0, ConnectionType::Matrix_Matrix); // 反向注入测试归一化
        try {
            strictMgr.FinalizeAndAggregate(); // 这里必须抛出 logic_error
            std::cout << "[FAIL] Allowed duplicate MM connections!\n";
        }
        catch (const std::logic_error&) {
            std::cout << "[PASS]\n";
            passCount++;
        }

        // Test 6.3 & 6.4 物理合并与垃圾数据去重
        FIM_ConnectionManager aggMgr;
        std::cout << "  -> Test 6.3: Physical Parallel Merge (NNC)... ";
        // 模拟 3D 中三角化产生的合法相交碎片
        aggMgr.PushConnection(100, 200, 5.0, 1.0, 2.0, 0.5, ConnectionType::Matrix_Fracture);
        aggMgr.PushConnection(100, 200, 3.0, 1.0, 1.5, 0.4, ConnectionType::Matrix_Fracture);

        std::cout << "\n  -> Test 6.4: Absolute Geometric Duplicate Filtration... ";
        // 模拟完全重复的非法冗余数据
        aggMgr.PushConnection(50, 60, 10.0, 2.0, 1.0, 0.1, ConnectionType::Matrix_Fracture);
        aggMgr.PushConnection(50, 60, 10.0, 2.0, 1.0, 0.1, ConnectionType::Matrix_Fracture);

        aggMgr.FinalizeAndAggregate();
        const auto& testConns = aggMgr.GetConnections();

        bool test63_pass = false;
        bool test64_pass = false;
        for (const auto& c : testConns) {
            if (c.nodeI == 100 && c.nodeJ == 200) {
                // 3D 三角形碎片的传导率累加校验
                if (std::abs(c.T_Flow - 8.0) < 1e-8 && std::abs(c.aux_area - 3.5) < 1e-8) {
                    test63_pass = true;
                }
            }
            if (c.nodeI == 50 && c.nodeJ == 60) {
                // 如果漏洞被修复，传导率必须依然是 10.0；若没修复则会变成 20.0
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
            std::cout << "========== [3D Benchmark Completed Successfully] ==========\n" << std::endl;
        }
        else {
            std::cerr << "========== [3D Benchmark FAILED on Strict Boundaries] ==========\n" << std::endl;
        }
    }

} // namespace Benchmark3D
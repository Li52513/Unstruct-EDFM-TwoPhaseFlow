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

#include "3D_MeshManager.h"
#include "3D_FieldManager.h"
#include "TransmissibilitySolver_3D.h"
#include "SolverContrlStrName_op.h"

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

        // 3D 架构的 FI 直接来源于全局边的数量
        size_t nFracEdges = meshMgr.fracture_network().getGlobalEdges().size();

        // 3D 架构的 NNC 直接来源于扁平化的 interactionPairs
        size_t nNNC = meshMgr.getInteractionPairs().size();

        // 3D 架构的 FF (Star-Delta) 统计：采用端点量化无向键聚类
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

        // 聚类收集
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

        // 此处的统计保留合法网格对，去除越界数据（留给求解器最终裁决后会安全舍弃，预估只多不少）
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
        std::cout << "[Stage 2] Initializing Fields with Constant Physical Properties..." << std::endl;
        FieldManager_3D fieldMgr;

        // 基于精准的尺寸统计初始化场管理器 (第6个参数传 nFracEdges 给 FI)
        fieldMgr.InitSizes(nMat, nFrac, nNNC, nFF, nMatFaces, nFracEdges);

        PhysicalProperties_string_op::Rock rockStr;
        PhysicalProperties_string_op::Fracture_string fracStr;
        PhysicalProperties_string_op::Water waterStr;

        // 2.1 注入基岩物理量 (渗透率 1e-15 m2, 导热 2.5)
        fieldMgr.getOrCreateMatrixScalar(rockStr.k_xx_tag, 1e-15);
        fieldMgr.getOrCreateMatrixScalar(rockStr.k_yy_tag, 1e-15);
        fieldMgr.getOrCreateMatrixScalar(rockStr.k_zz_tag, 1e-15); // 3D 专属
        fieldMgr.getOrCreateMatrixScalar(rockStr.lambda_tag, 2.5);

        // 2.2 注入裂缝物理量 (渗透率 1e-10 m2, 缝宽 0.01m)
        fieldMgr.getOrCreateFractureScalar(fracStr.k_t_tag, 1e-10); // FF 使用
        fieldMgr.getOrCreateFractureScalar(fracStr.k_n_tag, 1e-10); // NNC 使用
        fieldMgr.getOrCreateFractureScalar(fracStr.aperture_tag, 0.01);
        fieldMgr.getOrCreateFractureScalar(fracStr.lambda_tag, 1.5);
        fieldMgr.getOrCreateFractureScalar(fracStr.phi_tag, 0.5);

        // 3D 求解器中明确使用了 Water::k_tag 作为流体热导率
        fieldMgr.getOrCreateFractureScalar(waterStr.k_tag, 0.6);

        // =========================================================
        // Stage 3: 运行静态传导率引擎 (Transmissibility Solver)
        // =========================================================
        std::cout << "[Stage 3] Executing Static Transmissibility Solvers..." << std::endl;
        TransmissibilitySolver_3D::Calculate_Transmissibility_Matrix(meshMgr, fieldMgr);
        TransmissibilitySolver_3D::Calculate_Transmissibility_NNC(meshMgr, fieldMgr);
        TransmissibilitySolver_3D::Calculate_Transmissibility_FractureInternal(meshMgr, fieldMgr); // 加入 FI 计算
        TransmissibilitySolver_3D::Calculate_Transmissibility_FF(meshMgr, fieldMgr);

        // =========================================================
        // Stage 4: 数据提取与验证报表导出 (CSV 导出)
        // =========================================================
        std::cout << "[Stage 4] Exporting Benchmark Results..." << std::endl;

        std::string csvName = exportPrefix + ".csv";
        std::ofstream csv(csvName);
        if (csv.is_open()) {
            csv << "Connection_Type,ID_1(Owner/Mat/FracA),ID_2(Neighbor/Frac/FracB),T_Flow,T_Heat\n";

            // 4.1 写入基岩内部面传导率 (Matrix-Matrix)
            auto t_mat_flow = fieldMgr.getMatrixFaceScalar("T_Matrix_Flow");
            auto t_mat_heat = fieldMgr.getMatrixFaceScalar("T_Matrix_Heat");
            if (t_mat_flow && t_mat_heat) {
                const auto& faces = meshMgr.mesh().getFaces();
                for (size_t i = 0; i < faces.size(); ++i) {
                    if (faces[i].isBoundary()) continue; // 跳过边界面
                    csv << "MatrixFace,Mat_" << faces[i].ownerCell << ",Mat_" << faces[i].neighborCell << ","
                        << std::scientific << std::setprecision(6)
                        << t_mat_flow->data[i] << "," << t_mat_heat->data[i] << "\n";
                }
            }

            // 4.2 写入基岩-裂缝窜流传导率 (NNC)
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

            // 4.3 写入裂缝内部传导率 (FI)
            auto t_fi_flow = fieldMgr.getFractureEdgeScalar("T_FI_Flow");
            auto t_fi_heat = fieldMgr.getFractureEdgeScalar("T_FI_Heat");
            if (t_fi_flow && t_fi_heat) {
                const auto& globalEdges = meshMgr.fracture_network().getGlobalEdges();
                for (size_t i = 0; i < globalEdges.size(); ++i) {
                    const auto& edge = globalEdges[i];
                    // 在对齐映射的数组中，仅剔除掉边界边的 0 值，保留有效内部流通对
                    if (edge.neighborCell_solverIndex >= 0) {
                        csv << "FracInternal,FracElemGlobal_" << edge.ownerCell_solverIndex
                            << ",FracElemGlobal_" << edge.neighborCell_solverIndex << ","
                            << std::scientific << std::setprecision(6)
                            << t_fi_flow->data[i] << "," << t_fi_heat->data[i] << "\n";
                    }
                }
            }

            // 4.4 写入裂缝-裂缝相交通导率 (FF - Deterministic Star-Delta)
            auto t_ff_flow = fieldMgr.getFFScalar("T_FF_Flow");
            auto t_ff_heat = fieldMgr.getFFScalar("T_FF_Heat");
            if (t_ff_flow && t_ff_heat) {
                // 保证导出顺序与求解器绝对一致的确定性字典序
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
        else {
            std::cerr << "  -> [Error] Failed to create benchmark CSV file: " << csvName << std::endl;
        }

        std::cout << "========== [3D Benchmark Completed Successfully] ==========\n" << std::endl;
    }

} // namespace Benchmark3D
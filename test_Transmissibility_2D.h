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

namespace Benchmark2D {

    /**
     * @brief 运行 2D 传导率基准测试全流程
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
        int nx = 3, ny = 3;

        // 1.1 生成 10x10 的结构化四边形极简网格 (用作完美对标)
        MeshManager meshMgr(Lx, Ly, 0.0, nx, ny, 0, false, true);
        meshMgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OrthogonalCorrection);

        // 1.2 构建一个标准的 "X" 型交叉裂缝，以完美触发 FF 交叉点计算
        meshMgr.addFracture(Vector(2.0, 2.0, 0.0), Vector(8.0, 8.0, 0.0)); // Frac 0
        meshMgr.addFracture(Vector(2.0, 8.0, 0.0), Vector(8.0, 2.0, 0.0)); // Frac 1

        // 1.3 裂缝求交与剖分
        meshMgr.setDistanceMetric(DistanceMetric::CrossAwareGauss);
        meshMgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);

        // 1.4 构建裂缝间交叉拓扑 (计算 FF 交点并填入 globalFFPts)
        meshMgr.BuildFracturetoFractureTopology();

        // 1.5 建立全局矩阵系统索引
        meshMgr.BuildGlobalSystemIndexing();

        meshMgr.exportMesh(exportPrefix + "_MatrixMesh");
		meshMgr.exportFractures(exportPrefix + "_FractureMesh");

        // 1.6 精准统计各类型连接规模，用于内存预分配
        size_t nMat = meshMgr.mesh().getCells().size();
        size_t nFrac = meshMgr.getTotalDOFCount() - nMat;
        size_t nMatFaces = meshMgr.mesh().getFaces().size();

        size_t exactNNCCount = 0;
        const auto& nncMap = meshMgr.getNNCTopologyMap();
        for (const auto& kv : nncMap) {
            exactNNCCount += kv.second.size();
        }

        // [修复核心] 使用 Star-Delta 模型重构 FF 数量统计
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
            if (n > 1) nFF += n * (n - 1) / 2; // 完全匹配星角变换对数
        }

        size_t nFI = 0;
        for (const auto& frac : meshMgr.fracture_network().fractures) {
            if (frac.elements.size() > 1) {
                nFI += (frac.elements.size() - 1);
            }
        }

        std::cout << "  -> Matrix Cells: " << nMat << "\n"
            << "  -> Fracture Elements: " << nFrac << "\n"
            << "  -> NNC Pairs (Exact): " << exactNNCCount << "\n"
            << "  -> FI Connections: " << nFI << "\n" // 新增打印
            << "  -> FF Intersections: " << nFF << std::endl;

        // =========================================================
        // Stage 2: 场内存预分配与测试常量注入 (Mocking Properties)
        // =========================================================
        std::cout << "[Stage 2] Initializing Fields with Constant Physical Properties..." << std::endl;
        FieldManager_2D fieldMgr;
        fieldMgr.InitSizes(nMat, nFrac, exactNNCCount, nFF, nMatFaces, nFI);

        PhysicalProperties_string_op::Rock rockStr;
        PhysicalProperties_string_op::Fracture_string fracStr;
        PhysicalProperties_string_op::Water waterStr;

        // 2.1 注入基岩常数物理量 (渗透率 1e-15 m2, 导热 2.5)
        fieldMgr.getOrCreateMatrixScalar(rockStr.k_xx_tag, 1e-15);
        fieldMgr.getOrCreateMatrixScalar(rockStr.k_yy_tag, 1e-15);
        fieldMgr.getOrCreateMatrixScalar(rockStr.lambda_tag, 2.5);

        // 2.2 注入裂缝常数物理量 (渗透率 1e-10 m2, 缝宽 0.01m)
        fieldMgr.getOrCreateFractureScalar(fracStr.k_t_tag, 1e-10); // FF 切向传导依赖此项
        fieldMgr.getOrCreateFractureScalar(fracStr.k_n_tag, 1e-10); // NNC 法向传导依赖此项
        fieldMgr.getOrCreateFractureScalar(fracStr.aperture_tag, 0.01);
        fieldMgr.getOrCreateFractureScalar(fracStr.lambda_tag, 1.5);
        fieldMgr.getOrCreateFractureScalar(fracStr.phi_tag, 0.5);
        fieldMgr.getOrCreateFractureScalar(waterStr.k_tag, 0.6);

        // =========================================================
        // Stage 3: 运行静态传导率引擎 (Transmissibility Solver)
        // =========================================================
        std::cout << "[Stage 3] Executing Static Transmissibility Solvers..." << std::endl;
        TransmissibilitySolver_2D::Calculate_Transmissibility_Matrix(meshMgr, fieldMgr);
        TransmissibilitySolver_2D::Calculate_Transmissibility_NNC(meshMgr, fieldMgr);
        TransmissibilitySolver_2D::Calculate_Transmissibility_FF(meshMgr, fieldMgr);
		TransmissibilitySolver_2D::Calculate_Transmissibility_FractureInternal(meshMgr, fieldMgr);

        // =========================================================
        // Stage 4: 数据提取与验证报表导出 (CSV 导出)
        // =========================================================
        std::cout << "[Stage 4] Exporting Benchmark Results..." << std::endl;

        std::string csvName = exportPrefix + ".csv";
        std::ofstream csv(csvName);
        if (csv.is_open()) {
            csv << "Connection_Type,ID_1(Owner/FracA),ID_2(Neighbor/FracB),T_Flow,T_Heat\n";

            // 4.1 Matrix-Matrix 内部面传导率导出
            auto t_mat_flow = fieldMgr.getMatrixFaceScalar("T_Matrix_Flow");
            auto t_mat_heat = fieldMgr.getMatrixFaceScalar("T_Matrix_Heat");
            if (t_mat_flow && t_mat_heat) {
                const auto& faces = meshMgr.mesh().getFaces();
                for (size_t i = 0; i < faces.size(); ++i) {
                    if (faces[i].isBoundary()) continue;
                    csv << "MatrixFace,Mat_" << faces[i].ownerCell_index << ",Mat_" << faces[i].neighborCell_index << ","
                        << std::scientific << std::setprecision(6)
                        << (*t_mat_flow)[i] << "," << (*t_mat_heat)[i] << "\n";
                }
            }

            // 4.2 NNC 窜流传导率导出
            auto t_nnc_flow = fieldMgr.getNNCScalar("T_NNC_Flow");
            auto t_nnc_heat = fieldMgr.getNNCScalar("T_NNC_Heat");
            if (t_nnc_flow && t_nnc_heat) {
                size_t nncIndex = 0;
                for (const auto& kv : nncMap) {
                    int matrixGlobalId = meshMgr.mesh().getCells()[kv.first].id;
					int matrixLocalId = kv.first;
                    for (int fSolverIdx : kv.second) {
                        csv << "NNC,MatGlobal_" << matrixLocalId << ",FracSolver_" << fSolverIdx << ","
                            << std::scientific << std::setprecision(6)
                            << (*t_nnc_flow)[nncIndex] << "," << (*t_nnc_heat)[nncIndex] << "\n";
                        nncIndex++;
                    }
                }
            }

            // 4.3 写入裂缝内部传导率 (FI)
            auto t_fi_flow = fieldMgr.getFractureFaceScalar("T_FI_Flow");
            auto t_fi_heat = fieldMgr.getFractureFaceScalar("T_FI_Heat");
            if (t_fi_flow && t_fi_heat) {
                size_t fiIdx = 0;
                for (const auto& frac : meshMgr.fracture_network().fractures) {
                    for (size_t i = 0; i < frac.elements.size() - 1; ++i) {
                        if (fiIdx < t_fi_flow->data.size()) {
                            csv << "FracInternal,Elem_" << frac.elements[i].solverIndex
                                << ",Elem_" << frac.elements[i + 1].solverIndex << ","
                                << std::scientific << std::setprecision(6)
                                << t_fi_flow->data[fiIdx] << "," << t_fi_heat->data[fiIdx] << "\n";
                            fiIdx++;
                        }
                    }
                }
            }

            // 4.3 FF 交叉传导率导出 (对齐 Star-Delta 展开后的微观单元对)
            auto t_ff_flow = fieldMgr.getFFScalar("T_FF_Flow");
            auto t_ff_heat = fieldMgr.getFFScalar("T_FF_Heat");
            if (t_ff_flow && t_ff_heat) {
                // 为了与计算顺序绝对一致，这里重新取一遍排序好的 JunctionID
                std::vector<int> junctionIDs;
                for (const auto& kv : junctionMap) junctionIDs.push_back(kv.first);
                std::sort(junctionIDs.begin(), junctionIDs.end());

                size_t ffIdx = 0;
                for (int jID : junctionIDs) {
                    const auto& elemIndices = junctionMap[jID];
                    size_t nElems = elemIndices.size();
                    if (nElems < 2) continue;

                    for (size_t i = 0; i < nElems; ++i) {
                        for (size_t j = i + 1; j < nElems; ++j) {
                            if (ffIdx < t_ff_flow->data.size()) {
                                // 现在 FF 输出的是微观段对段的真实连接，而非宏观裂缝！
                                csv << "FF_StarDelta,FracElem_" << elemIndices[i]
                                    << ",FracElem_" << elemIndices[j] << ","
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

        std::cout << "========== [2D Benchmark Completed Successfully] ==========\n" << std::endl;
    }

} // namespace Benchmark2D
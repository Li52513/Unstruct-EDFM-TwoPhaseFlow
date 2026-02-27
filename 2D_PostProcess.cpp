/**
 * @file 2D_PostProcess.cpp
 * @brief 2D-EDFM 后处理可视化模块 (Tecplot Exporter) 实现
 * @details 负责将 2D 基岩域和 1D 裂缝域的几何拓扑及所有物理场导出为 Tecplot BLOCK 格式，
 * 支持混合网格（三角形/四边形）自适应填充，并解决基岩与裂缝变量不对齐的问题。
 */

#include "2D_PostProcess.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <set>

 // ==============================================================================
 // 构造函数
 // ==============================================================================
PostProcess_2D::PostProcess_2D(const MeshManager& meshMgr, const FieldManager_2D& fieldMgr)
    : meshMgr_(meshMgr), fieldMgr_(fieldMgr)
{
}

// ==============================================================================
// 核心导出函数
// ==============================================================================
void PostProcess_2D::ExportTecplot(const std::string& filename, double time) const
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("[PostProcess_2D Error] Failed to open file for Tecplot export: " + filename);
    }

    std::cout << ">>> Exporting 2D Tecplot data to: " << filename << " (Time = " << time << ") ..." << std::endl;

    // 1. 获取全局唯一变量名
    std::vector<std::string> varNames = GetAllUniqueFieldNames();

    // 2. 写入全局 Header
    file << "TITLE = \"2D EDFM Multiphase Simulation Result\"\n";
    file << "VARIABLES = \"X\", \"Y\"";
    for (const auto& name : varNames) {
        file << ", \"" << name << "\"";
    }
    file << "\n";

    // 设置科学计数法及高精度输出
    file << std::scientific << std::setprecision(8);

    // =========================================================
    // Zone 1: Matrix Domain (基岩域 - 2D Cell)
    // =========================================================
    const auto& matrixNodes = meshMgr_.mesh().getNodes();
    const auto& matrixCells = meshMgr_.mesh().getCells();
    size_t numMatrixNodes = matrixNodes.size();
    size_t numMatrixCells = matrixCells.size();

    if (numMatrixCells > 0) {
        // 判定基岩网格的 Tecplot 单元类型 (Tri / Quad)
        int maxNodesPerCell = 0;
        for (const auto& cell : matrixCells) {
            maxNodesPerCell = std::max(maxNodesPerCell, static_cast<int>(cell.CellNodeIDs.size()));
        }
        std::string matrixZoneType = GetTecplotElementType(maxNodesPerCell);

        // 写入 Matrix Zone Header (BLOCK 格式)
        file << "ZONE T=\"Matrix\", SOLUTIONTIME=" << time << ", DATAPACKING=BLOCK, "
            << "NODES=" << numMatrixNodes << ", ELEMENTS=" << numMatrixCells << ", "
            << "ZONETYPE=" << matrixZoneType;

        // 声明哪些变量是在单元中心 (X,Y 为 Node，其余为 Cell-Centered)
        if (!varNames.empty()) {
            file << ", VARLOCATION=([3-" << 2 + varNames.size() << "]=CELLCENTERED)";
        }
        file << "\n";

        // (1.1) 输出基岩节点 X 坐标
        for (const auto& node : matrixNodes) {
            file << node.coord.m_x << "\n";
        }
        // (1.2) 输出基岩节点 Y 坐标
        for (const auto& node : matrixNodes) {
            file << node.coord.m_y << "\n";
        }

        // (1.3) 输出基岩变量场 (CELLCENTERED)
        for (const auto& name : varNames) {
            auto field = fieldMgr_.getMatrixScalar(name);
            if (field) {
                // 如果基岩有这个物理场，输出真实数据
                for (size_t i = 0; i < numMatrixCells; ++i) {
                    file << field->data[i] << "\n";
                }
            }
            else {
                // 如果不存在（例如只有裂缝特有的属性），补零对齐
                for (size_t i = 0; i < numMatrixCells; ++i) {
                    file << 0.0 << "\n";
                }
            }
        }

        // (1.4) 输出基岩网格连通性 (Connectivity)
        // 建立 Gmsh Global Node ID 到 局部按序输出的 1-based Index 映射
        std::unordered_map<int, int> gmshIdToLocalIdx;
        for (size_t i = 0; i < numMatrixNodes; ++i) {
            gmshIdToLocalIdx[matrixNodes[i].id] = static_cast<int>(i) + 1; // Tecplot 是 1-based
        }

        for (const auto& cell : matrixCells) {
            size_t nNodes = cell.CellNodeIDs.size();
            for (int i = 0; i < maxNodesPerCell; ++i) {
                // 如果当前单元是三角形（3个节点），但整体 Zone 是四边形（max=4），需重复最后一个节点
                int nodeGmshId = (i < nNodes) ? cell.CellNodeIDs[i] : cell.CellNodeIDs.back();

                // 安全获取并输出 1-based 索引
                auto it = gmshIdToLocalIdx.find(nodeGmshId);
                if (it != gmshIdToLocalIdx.end()) {
                    file << it->second << " ";
                }
                else {
                    file << "1 "; // 极端防崩后备方案
                }
            }
            file << "\n";
        }
    }

    // =========================================================
    // Zone 2: Fracture Domain (裂缝域 - 1D Segments)
    // =========================================================
    const auto& fractures = meshMgr_.fracture_network().fractures;

    // 统计裂缝微元数
    size_t totalFracCells = 0;
    for (const auto& frac : fractures) {
        totalFracCells += frac.elements.size();
    }

    // 采用独立节点策略：每个微元动态生成 2 个端点，总节点数为微元数 * 2
    size_t totalFracNodes = totalFracCells * 2;

    if (totalFracCells > 0) {
        // 写入 Fracture Zone Header
        file << "ZONE T=\"Fractures\", SOLUTIONTIME=" << time << ", DATAPACKING=BLOCK, "
            << "NODES=" << totalFracNodes << ", ELEMENTS=" << totalFracCells << ", "
            << "ZONETYPE=" << GetTecplotElementType(2);

        if (!varNames.empty()) {
            file << ", VARLOCATION=([3-" << 2 + varNames.size() << "]=CELLCENTERED)";
        }
        file << "\n";

        // 收集所有的物理节点坐标 (避免多次循环导致 I/O 错位)
        std::vector<double> fracX;
        std::vector<double> fracY;
        fracX.reserve(totalFracNodes);
        fracY.reserve(totalFracNodes);

        for (const auto& frac : fractures) {
            // 计算宏观裂缝的 X, Y 跨度
            double dx = frac.end.m_x - frac.start.m_x;
            double dy = frac.end.m_y - frac.start.m_y;

            for (const auto& cell : frac.elements) {
                // 基于归一化参数 [param0, param1] 线性插值计算微元端点坐标
                fracX.push_back(frac.start.m_x + dx * cell.param0); // 节点1 X
                fracX.push_back(frac.start.m_x + dx * cell.param1); // 节点2 X

                fracY.push_back(frac.start.m_y + dy * cell.param0); // 节点1 Y
                fracY.push_back(frac.start.m_y + dy * cell.param1); // 节点2 Y
            }
        }

        // (2.1) 输出裂缝节点 X 坐标
        for (double x : fracX) file << x << "\n";

        // (2.2) 输出裂缝节点 Y 坐标
        for (double y : fracY) file << y << "\n";

        // (2.3) 输出裂缝变量场 (CELLCENTERED)
        for (const auto& name : varNames) {
            auto field = fieldMgr_.getFractureScalar(name);
            if (field) {
                for (size_t i = 0; i < totalFracCells; ++i) {
                    file << field->data[i] << "\n";
                }
            }
            else {
                for (size_t i = 0; i < totalFracCells; ++i) {
                    file << 0.0 << "\n";
                }
            }
        }

        // (2.4) 输出裂缝连通性 (Connectivity)
        // 独立网格拓扑：第 i 个微元严格连接节点 (2*i + 1) 和 (2*i + 2) [Tecplot要求1-based]
        for (size_t i = 0; i < totalFracCells; ++i) {
            file << (2 * i + 1) << " " << (2 * i + 2) << "\n";
        }
    }

    file.close();
    std::cout << ">>> 2D Tecplot export completed successfully." << std::endl;
}

// ==============================================================================
// 内部辅助函数
// ==============================================================================

std::vector<std::string> PostProcess_2D::GetAllUniqueFieldNames() const
{
    std::set<std::string> uniqueNames;

    // 提取基岩域所有双精度标量场
    for (const auto& pair : fieldMgr_.matrixFields.fields) {
        // 利用 dynamic_pointer_cast 安全甄别 volScalarField (排除 Vector/ADVar 场)
        if (std::dynamic_pointer_cast<volScalarField>(pair.second)) {
            uniqueNames.insert(pair.first);
        }
    }

    // 提取裂缝域所有双精度标量场
    for (const auto& pair : fieldMgr_.fractureFields.fields) {
        if (std::dynamic_pointer_cast<volScalarField>(pair.second)) {
            uniqueNames.insert(pair.first);
        }
    }

    return std::vector<std::string>(uniqueNames.begin(), uniqueNames.end());
}

std::string PostProcess_2D::GetTecplotElementType(int numNodes) const
{
    if (numNodes == 2) {
        return "FELINESEG";         // 1D 线段 (裂缝)
    }
    else if (numNodes == 3) {
        return "FETRIANGLE";        // 纯三角形基岩
    }
    else if (numNodes == 4) {
        return "FEQUADRILATERAL";   // 四边形基岩 或 混合(三角形+四边形)
    }

    // 降级保护
    return "FEQUADRILATERAL";
}
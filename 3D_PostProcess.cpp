#include "3D_PostProcess.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <set>
#include <unordered_map> 
#include <algorithm> 

// =========================================================
// 辅助工具：获取 Tecplot 标准 ZoneType 及目标节点数
// =========================================================
struct TecplotElemInfo {
    std::string typeName;
    int targetNodeCount;
};

TecplotElemInfo GetTecplotInfo(const std::string& internalGridType) {
    if (internalGridType == "FETETRA")      return { "FETETRAHEDRON", 4 };
    if (internalGridType == "FEPYRAMID")    return { "FEBRICK", 8 };
    if (internalGridType == "FEPRISM")      return { "FEBRICK", 8 };
    if (internalGridType == "FEBRICK")      return { "FEBRICK", 8 };
    if (internalGridType == "FETRIANGLE")   return { "FETRIANGLE", 3 };
    if (internalGridType == "FEQUAD")       return { "FEQUADRILATERAL", 4 };
    return { "FETETRAHEDRON", 4 };
}

PostProcess_3D::PostProcess_3D(const MeshManager_3D& meshMgr, const FieldManager_3D& fieldMgr)
    : meshMgr_(meshMgr), fieldMgr_(fieldMgr)
{
}

std::vector<std::string> PostProcess_3D::GetAllUniqueFieldNames() const
{
    std::set<std::string> uniqueNames;
    for (const auto& pair : fieldMgr_.matrixFields.fields) {
        if (std::dynamic_pointer_cast<volScalarField>(pair.second)) uniqueNames.insert(pair.first);
    }
    for (const auto& pair : fieldMgr_.fractureFields.fields) {
        if (std::dynamic_pointer_cast<volScalarField>(pair.second)) uniqueNames.insert(pair.first);
    }
    return std::vector<std::string>(uniqueNames.begin(), uniqueNames.end());
}

std::string PostProcess_3D::GetTecplotElementType(int numNodes) const {
    switch (numNodes) {
    case 4: return "FETETRA";
    case 5: return "FEPYRAMID";
    case 6: return "FEPRISM";
    case 8: return "FEBRICK";
    default: return "FETETRA";
    }
}

// =========================================================
// 核心导出逻辑
// =========================================================
void PostProcess_3D::ExportTecplot(const std::string& filename, double time)
{
    std::cout << "[PostProcess] Exporting Ultra-Robust Tecplot file: " << filename << " (t=" << time << ")" << std::endl;

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "[Error] Cannot open file: " << filename << std::endl;
        return;
    }

    out << std::scientific << std::setprecision(6);
    std::vector<std::string> allVarNames = GetAllUniqueFieldNames();

    const auto& nodes = meshMgr_.mesh().getNodes();
    std::unordered_map<int, int> nodeID2Index;
    for (size_t i = 0; i < nodes.size(); ++i) {
        nodeID2Index[nodes[i].id] = static_cast<int>(i) + 1;
    }

    // Header
    out << "TITLE = \"3D-EDFM Simulation Results\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"Z\"";
    for (const auto& name : allVarNames) out << ", \"" << name << "\"";
    out << "\n";

    // ---------------------------------------------------------
    // Zone 1: Matrix
    // ---------------------------------------------------------
    const auto& cells = meshMgr_.mesh().getCells();
    size_t numCells = cells.size();

    if (numCells > 0) {
        std::string internalDesc = meshMgr_.mesh().getMatrixElementType();
        if (internalDesc == "Unknown" || internalDesc.empty()) {
            internalDesc = GetTecplotElementType(static_cast<int>(cells[0].CellNodeIDs.size()));
        }
        TecplotElemInfo tecInfo = GetTecplotInfo(internalDesc);

        // 强制写入 FEBRICK 类型
        out << "ZONE T=\"Matrix\"\n";
        out << " STRANDID=1, SOLUTIONTIME=" << time << "\n";
        out << " NODES=" << nodes.size() << ", ELEMENTS=" << numCells
            << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n";

        if (!allVarNames.empty()) {
            out << " VARLOCATION=([4-" << (3 + allVarNames.size()) << "]=CELLCENTERED)\n";
        }

        int col = 0;
        // Geometry Output 
        col = 0; for (const auto& node : nodes) { out << node.coord.m_x << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";
        col = 0; for (const auto& node : nodes) { out << node.coord.m_y << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";
        col = 0; for (const auto& node : nodes) { out << node.coord.m_z << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";

        // Fields Output (带 NaN 洗钱)
        for (const auto& name : allVarNames) {
            auto fieldPtr = fieldMgr_.getMatrixScalar(name);
            col = 0;
            if (fieldPtr) {
                for (double val : fieldPtr->data) {
                    out << (std::isfinite(val) ? val : 0.0) << " ";
                    if (++col % 10 == 0) out << "\n";
                }
            }
            else {
                for (size_t i = 0; i < numCells; ++i) {
                    out << "0.0 ";
                    if (++col % 10 == 0) out << "\n";
                }
            }
            if (col % 10 != 0) out << "\n";
        }

        // =========================================================
        // Connectivity (Mapping + Padding with Chirality Auto-Correction)
        // =========================================================
        for (const auto& cell : cells) {
            const auto& ids = cell.CellNodeIDs;
            int n = static_cast<int>(ids.size());

            std::vector<int> mappedIndices;
            mappedIndices.reserve(n);
            for (int id : ids) {
                if (nodeID2Index.find(id) != nodeID2Index.end()) mappedIndices.push_back(nodeID2Index[id]);
                else mappedIndices.push_back(1);
            }

            // --- 【核心引擎】：三维单元体积手性自动校验与修复 ---
            if (n >= 4) {
                int idx0 = mappedIndices[0] - 1;
                int idx1 = mappedIndices[1] - 1;
                int idx2 = mappedIndices[2] - 1;
                int idx3 = -1;

                if (n == 8) idx3 = mappedIndices[4] - 1;
                else if (n == 6) idx3 = mappedIndices[3] - 1;
                else if (n == 5) idx3 = mappedIndices[4] - 1;
                else if (n == 4) idx3 = mappedIndices[3] - 1;

                if (idx3 >= 0 && idx3 < nodes.size()) {
                    double ax = nodes[idx1].coord.m_x - nodes[idx0].coord.m_x;
                    double ay = nodes[idx1].coord.m_y - nodes[idx0].coord.m_y;
                    double az = nodes[idx1].coord.m_z - nodes[idx0].coord.m_z;

                    double bx = nodes[idx2].coord.m_x - nodes[idx0].coord.m_x;
                    double by = nodes[idx2].coord.m_y - nodes[idx0].coord.m_y;
                    double bz = nodes[idx2].coord.m_z - nodes[idx0].coord.m_z;

                    double cx = nodes[idx3].coord.m_x - nodes[idx0].coord.m_x;
                    double cy = nodes[idx3].coord.m_y - nodes[idx0].coord.m_y;
                    double cz = nodes[idx3].coord.m_z - nodes[idx0].coord.m_z;

                    double nx = ay * bz - az * by;
                    double ny = az * bx - ax * bz;
                    double nz = ax * by - ay * bx;

                    double vol = nx * cx + ny * cy + nz * cz;

                    // 如果体积小于0，强行翻转手性
                    if (vol < 0.0) {
                        if (n == 8) {
                            std::swap(mappedIndices[0], mappedIndices[4]);
                            std::swap(mappedIndices[1], mappedIndices[5]);
                            std::swap(mappedIndices[2], mappedIndices[6]);
                            std::swap(mappedIndices[3], mappedIndices[7]);
                        }
                        else if (n == 6) {
                            std::swap(mappedIndices[1], mappedIndices[2]);
                            std::swap(mappedIndices[4], mappedIndices[5]);
                        }
                        else if (n == 5) {
                            std::swap(mappedIndices[1], mappedIndices[3]);
                        }
                        else if (n == 4) {
                            std::swap(mappedIndices[1], mappedIndices[2]);
                        }
                    }
                }
            }

            // 执行退化映射
            if (n == 8) {
                for (int k = 0; k < 8; ++k) out << mappedIndices[k] << " ";
            }
            else if (n == 6) { // 三棱柱 6 -> 8
                out << mappedIndices[0] << " " << mappedIndices[1] << " " << mappedIndices[2] << " " << mappedIndices[2] << " "
                    << mappedIndices[3] << " " << mappedIndices[4] << " " << mappedIndices[5] << " " << mappedIndices[5] << " ";
            }
            else if (n == 5) { // 金字塔 5 -> 8
                out << mappedIndices[0] << " " << mappedIndices[1] << " " << mappedIndices[2] << " " << mappedIndices[3] << " "
                    << mappedIndices[4] << " " << mappedIndices[4] << " " << mappedIndices[4] << " " << mappedIndices[4] << " ";
            }
            else if (n == 4) { // 四面体 4 -> 8
                out << mappedIndices[0] << " " << mappedIndices[1] << " " << mappedIndices[2] << " " << mappedIndices[2] << " "
                    << mappedIndices[3] << " " << mappedIndices[3] << " " << mappedIndices[3] << " " << mappedIndices[3] << " ";
            }
            else {
                for (int k = 0; k < 8; ++k) out << mappedIndices[std::min(k, n - 1)] << " ";
            }
            out << "\n";
        }
    }

    // ---------------------------------------------------------
    // Zone 2: Fracture
    // ---------------------------------------------------------
    const auto& fractures = meshMgr_.fracture_network().getFractures();
    size_t totalFracElems = 0;
    size_t totalFracNodes = 0;

    for (const auto& frac : fractures) {
        totalFracNodes += frac.fracNodes.size();
        totalFracElems += frac.fracCells.size();
    }

    if (totalFracElems > 0) {
        // 强制写入 FEQUADRILATERAL 格式
        out << "ZONE T=\"Fracture\"\n";
        out << " STRANDID=2, SOLUTIONTIME=" << time << "\n";
        out << " NODES=" << totalFracNodes << ", ELEMENTS=" << totalFracElems
            << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n";

        if (!allVarNames.empty()) out << " VARLOCATION=([4-" << (3 + allVarNames.size()) << "]=CELLCENTERED)\n";

        int col = 0;
        // Fracture Geometry 
        col = 0; for (const auto& frac : fractures) for (const auto& n : frac.fracNodes) { out << n.coord.m_x << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";
        col = 0; for (const auto& frac : fractures) for (const auto& n : frac.fracNodes) { out << n.coord.m_y << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";
        col = 0; for (const auto& frac : fractures) for (const auto& n : frac.fracNodes) { out << n.coord.m_z << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";

        // Fracture Fields 
        for (const auto& name : allVarNames) {
            auto fieldPtr = fieldMgr_.getFractureScalar(name);
            col = 0;
            if (fieldPtr) {
                for (double val : fieldPtr->data) {
                    out << (std::isfinite(val) ? val : 0.0) << " ";
                    if (++col % 10 == 0) out << "\n";
                }
            }
            else {
                for (size_t i = 0; i < totalFracElems; ++i) {
                    out << "0.0 ";
                    if (++col % 10 == 0) out << "\n";
                }
            }
            if (col % 10 != 0) out << "\n";
        }

        // Connectivity
        size_t nodeOffset = 0;
        for (const auto& frac : fractures) {
            for (const auto& cell : frac.fracCells) {
                const auto& ids = cell.nodeIndices;
                int n = static_cast<int>(ids.size());

                for (int k = 0; k < 4; ++k) {
                    out << (nodeOffset + ids[std::min(k, n - 1)] + 1) << " ";
                }
                out << "\n";
            }
            nodeOffset += frac.fracNodes.size();
        }
    }

    out.close();
    std::cout << "[PostProcess] Export completed.\n";
}

// =========================================================
// 新一代渲染引擎：全量导出至 ParaView (VTK 格式)
// =========================================================
void PostProcess_3D::ExportVTK(const std::string& filename, double time)
{
    std::cout << "[PostProcess] Exporting ParaView VTK file: " << filename << " (t=" << time << ")" << std::endl;

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "[Error] Cannot open file: " << filename << std::endl;
        return;
    }

    out << std::scientific << std::setprecision(6);
    std::vector<std::string> allVarNames = GetAllUniqueFieldNames();

    // 1. 提取基岩数据
    const auto& nodes = meshMgr_.mesh().getNodes();
    const auto& cells = meshMgr_.mesh().getCells();
    size_t numMatrixNodes = nodes.size();
    size_t numMatrixCells = cells.size();

    // 建立基岩节点到 0-based 索引的映射 (VTK 必须从 0 开始)
    std::unordered_map<int, int> nodeID2Index;
    for (size_t i = 0; i < numMatrixNodes; ++i) {
        nodeID2Index[nodes[i].id] = static_cast<int>(i);
    }

    // 2. 提取裂缝数据
    const auto& fractures = meshMgr_.fracture_network().getFractures();
    size_t numFracNodes = 0;
    size_t numFracCells = 0;
    for (const auto& frac : fractures) {
        numFracNodes += frac.fracNodes.size();
        numFracCells += frac.fracCells.size();
    }

    // 总统计
    size_t totalNodes = numMatrixNodes + numFracNodes;
    size_t totalCells = numMatrixCells + numFracCells;

    // 计算 CELLS 列表的总长度 (每个单元记录：节点数 n + n个索引)
    size_t cellListSize = 0;
    for (const auto& cell : cells) cellListSize += (1 + cell.CellNodeIDs.size());
    for (const auto& frac : fractures) {
        for (const auto& cell : frac.fracCells) {
            cellListSize += (1 + cell.nodeIndices.size());
        }
    }

    // ==========================================
    // VTK Header
    // ==========================================
    out << "# vtk DataFile Version 3.0\n";
    out << "3D-EDFM Polyhedron Simulation (Time: " << time << ")\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // ==========================================
    // POINTS (节点坐标)
    // ==========================================
    out << "POINTS " << totalNodes << " double\n";
    int col = 0;
    // 基岩节点
    for (const auto& n : nodes) {
        out << n.coord.m_x << " " << n.coord.m_y << " " << n.coord.m_z << " ";
        if (++col % 3 == 0) out << "\n";
    }
    // 裂缝节点
    for (const auto& frac : fractures) {
        for (const auto& n : frac.fracNodes) {
            out << n.coord.m_x << " " << n.coord.m_y << " " << n.coord.m_z << " ";
            if (++col % 3 == 0) out << "\n";
        }
    }
    if (col % 3 != 0) out << "\n";

    // ==========================================
    // CELLS (单元连接表)
    // ==========================================
    out << "CELLS " << totalCells << " " << cellListSize << "\n";

    // 基岩单元
    for (const auto& cell : cells) {
        out << cell.CellNodeIDs.size() << " ";
        for (int id : cell.CellNodeIDs) {
            out << nodeID2Index[id] << " ";
        }
        out << "\n";
    }

    // 裂缝单元 (节点索引需加上基岩节点总数作为偏移量)
    size_t fracNodeOffset = numMatrixNodes;
    for (const auto& frac : fractures) {
        for (const auto& cell : frac.fracCells) {
            out << cell.nodeIndices.size() << " ";
            for (int localIdx : cell.nodeIndices) {
                out << (fracNodeOffset + localIdx) << " ";
            }
            out << "\n";
        }
        fracNodeOffset += frac.fracNodes.size();
    }

    // ==========================================
    // CELL_TYPES (VTK 单元类型定义)
    // ==========================================
    // 10=Tetra, 12=Hex, 13=Wedge(Prism), 14=Pyramid, 5=Triangle, 9=Quad
    out << "CELL_TYPES " << totalCells << "\n";

    // 基岩类型
    for (const auto& cell : cells) {
        size_t n = cell.CellNodeIDs.size();
        if (n == 8) out << "12\n";
        else if (n == 6) out << "13\n"; // 原生支持三棱柱，完美！
        else if (n == 5) out << "14\n";
        else if (n == 4) out << "10\n";
        else out << "10\n";
    }

    // 裂缝类型
    for (const auto& frac : fractures) {
        for (const auto& cell : frac.fracCells) {
            if (cell.nodeIndices.size() == 4) out << "9\n"; // Quad
            else out << "5\n"; // Triangle
        }
    }

    // ==========================================
    // CELL_DATA (物理场输出)
    // ==========================================
    out << "CELL_DATA " << totalCells << "\n";

    out << "SCALARS DomainID int 1\n";
    out << "LOOKUP_TABLE default\n";
    // 基岩单元全部标记为 0
    for (size_t i = 0; i < numMatrixCells; ++i) out << "0\n";
    // 裂缝单元全部标记为 1
    for (size_t i = 0; i < numFracCells; ++i) out << "1\n";

    for (const auto& name : allVarNames) {
        out << "SCALARS " << name << " double 1\n";
        out << "LOOKUP_TABLE default\n";

        // 输出基岩场
        auto matField = fieldMgr_.getMatrixScalar(name);
        if (matField) {
            for (double val : matField->data) out << (std::isfinite(val) ? val : 0.0) << "\n";
        }
        else {
            for (size_t i = 0; i < numMatrixCells; ++i) out << "0.0\n";
        }

        // 输出裂缝场
        auto fracField = fieldMgr_.getFractureScalar(name);
        if (fracField) {
            for (double val : fracField->data) out << (std::isfinite(val) ? val : 0.0) << "\n";
        }
        else {
            for (size_t i = 0; i < numFracCells; ++i) out << "0.0\n";
        }
    }

    out.close();
    std::cout << "[PostProcess] VTK Export completed successfully.\n";
}
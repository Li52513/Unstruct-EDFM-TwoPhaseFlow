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

std::vector<std::string> PostProcess_3D::GetAllUniqueFieldNames() const {
    std::set<std::string> uniqueNames;
    for (const auto& pair : fieldMgr_.matrixFields.fields) uniqueNames.insert(pair.first);
    for (const auto& pair : fieldMgr_.fractureFields.fields) uniqueNames.insert(pair.first);
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
// 核心导出逻辑 (含换行保护)
// =========================================================
void PostProcess_3D::ExportTecplot(const std::string& filename, double time)
{
    std::cout << "[PostProcess] Exporting Tecplot file: " << filename << " (t=" << time << ")" << std::endl;

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "[Error] Cannot open file: " << filename << std::endl;
        return;
    }

    out << std::scientific << std::setprecision(6);
    std::vector<std::string> allVarNames = GetAllUniqueFieldNames();

    // 构建映射表
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

        out << "ZONE T=\"Matrix\"\n";
        out << " STRANDID=1, SOLUTIONTIME=" << time << "\n";
        out << " NODES=" << nodes.size() << ", ELEMENTS=" << numCells
            << ", DATAPACKING=BLOCK, ZONETYPE=" << tecInfo.typeName << "\n";

        if (!allVarNames.empty()) {
            out << " VARLOCATION=([4-" << (3 + allVarNames.size()) << "]=CELLCENTERED)\n";
        }

        // =================================================
        // Geometry Output (修复：每10个换行)
        // =================================================
        int col = 0;
        // X
        col = 0; for (const auto& node : nodes) { out << node.coord.m_x << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";
        // Y
        col = 0; for (const auto& node : nodes) { out << node.coord.m_y << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";
        // Z
        col = 0; for (const auto& node : nodes) { out << node.coord.m_z << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";

        // =================================================
        // Fields Output (修复：每10个换行)
        // =================================================
        for (const auto& name : allVarNames) {
            auto fieldPtr = fieldMgr_.getMatrixScalar(name);
            col = 0;
            if (fieldPtr) {
                for (double val : fieldPtr->data) {
                    out << val << " ";
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

        // Connectivity (Mapping + Padding)
        // 这里的换行是天然的(每个Cell一行)，通常不会超长，除非是多面体。
        // 对于 Brick/Tetra，一行只有几个数，无需额外处理。
        for (const auto& cell : cells) {
            const auto& ids = cell.CellNodeIDs;
            int n = static_cast<int>(ids.size());

            std::vector<int> mappedIndices;
            mappedIndices.reserve(n);
            for (int id : ids) {
                if (nodeID2Index.find(id) != nodeID2Index.end()) mappedIndices.push_back(nodeID2Index[id]);
                else mappedIndices.push_back(1);
            }

            if (internalDesc == "FEPRISM" && tecInfo.typeName == "FEBRICK") {
                if (n == 6) {
                    out << mappedIndices[0] << " " << mappedIndices[1] << " " << mappedIndices[2] << " " << mappedIndices[2] << " "
                        << mappedIndices[3] << " " << mappedIndices[4] << " " << mappedIndices[5] << " " << mappedIndices[5] << " ";
                }
                else {
                    for (int k = 0; k < 8; ++k) out << mappedIndices[std::min(k, n - 1)] << " ";
                }
            }
            else {
                for (int k = 0; k < tecInfo.targetNodeCount; ++k) {
                    if (k < n) out << mappedIndices[k] << " ";
                    else out << mappedIndices[n - 1] << " ";
                }
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

    std::string internalFracType = "FETRIANGLE";
    for (const auto& frac : fractures) {
        if (!frac.fracCells.empty() && frac.fracCells[0].nodeIndices.size() == 4) {
            internalFracType = "FEQUAD";
            break;
        }
    }
    TecplotElemInfo fracTecInfo = GetTecplotInfo(internalFracType);

    for (const auto& frac : fractures) {
        totalFracNodes += frac.fracNodes.size();
        totalFracElems += frac.fracCells.size();
    }

    if (totalFracElems > 0) {
        out << "ZONE T=\"Fracture\"\n";
        out << " STRANDID=2, SOLUTIONTIME=" << time << "\n";
        out << " NODES=" << totalFracNodes << ", ELEMENTS=" << totalFracElems
            << ", DATAPACKING=BLOCK, ZONETYPE=" << fracTecInfo.typeName << "\n";

        if (!allVarNames.empty()) out << " VARLOCATION=([4-" << (3 + allVarNames.size()) << "]=CELLCENTERED)\n";

        // =================================================
        // Fracture Geometry (修复：每10个换行)
        // =================================================
        int col = 0;
        // X
        col = 0; for (const auto& frac : fractures) for (const auto& n : frac.fracNodes) { out << n.coord.m_x << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";
        // Y
        col = 0; for (const auto& frac : fractures) for (const auto& n : frac.fracNodes) { out << n.coord.m_y << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";
        // Z
        col = 0; for (const auto& frac : fractures) for (const auto& n : frac.fracNodes) { out << n.coord.m_z << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";

        // =================================================
        // Fracture Fields (修复：每10个换行)
        // =================================================
        for (const auto& name : allVarNames) {
            auto fieldPtr = fieldMgr_.getFractureScalar(name);
            col = 0;
            if (fieldPtr) {
                for (double val : fieldPtr->data) {
                    out << val << " ";
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
                int lastIdx = nodeOffset + ids.back() + 1;

                for (int k = 0; k < fracTecInfo.targetNodeCount; ++k) {
                    if (k < n) out << (nodeOffset + ids[k] + 1) << " ";
                    else out << lastIdx << " ";
                }
                out << "\n";
            }
            nodeOffset += frac.fracNodes.size();
        }
    }

    out.close();
    std::cout << "[PostProcess] Export completed.\n";
}
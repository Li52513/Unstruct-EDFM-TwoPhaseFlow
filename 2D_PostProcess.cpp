/**
 * @file 2D_PostProcess.cpp
 * @brief 2D-EDFM �������ӻ�ģ�� (Tecplot Exporter) ʵ��
 * @details ���� 2D ������� 1D �ѷ���ļ������˼���������������Ϊ Tecplot BLOCK ��ʽ��
 * ֧�ֻ������������/�ı��Σ�����Ӧ��䣬������������ѷ��������������⡣
 */

#include "2D_PostProcess.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <set>
#include <unordered_map>
#include <vector>
#include <string>

 // ==============================================================================
 // ���캯��
 // ==============================================================================
PostProcess_2D::PostProcess_2D(const MeshManager& meshMgr, const FieldManager_2D& fieldMgr)
    : meshMgr_(meshMgr), fieldMgr_(fieldMgr)
{
}

// ==============================================================================
// ���ĵ�������
// ==============================================================================
void PostProcess_2D::ExportTecplot(const std::string& filename, double time) const
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("[PostProcess_2D Error] Failed to open file for Tecplot export: " + filename);
    }

    std::cout << ">>> Exporting 2D Tecplot data to: " << filename << " (Time = " << time << ") ..." << std::endl;

    // 1. ��ȡȫ��Ψһ������
    std::vector<std::string> varNames = GetAllUniqueFieldNames();

    // 2. д��ȫ�� Header
    file << "TITLE = \"2D EDFM Multiphase Simulation Result\"\n";
    file << "VARIABLES = \"X\", \"Y\"";
    for (const auto& name : varNames) {
        file << ", \"" << name << "\"";
    }
    file << "\n";

    // ���ÿ�ѧ���������߾������
    file << std::scientific << std::setprecision(8);

    // =========================================================
    // Zone 1: Matrix Domain (������ - 2D Cell)
    // =========================================================
    const auto& matrixNodes = meshMgr_.mesh().getNodes();
    const auto& matrixCells = meshMgr_.mesh().getCells();
    size_t numMatrixNodes = matrixNodes.size();
    size_t numMatrixCells = matrixCells.size();

    if (numMatrixCells > 0) {
        // �ж���������� Tecplot ��Ԫ���� (Tri / Quad)
        int maxNodesPerCell = 0;
        for (const auto& cell : matrixCells) {
            maxNodesPerCell = std::max(maxNodesPerCell, static_cast<int>(cell.CellNodeIDs.size()));
        }
        std::string matrixZoneType = GetTecplotElementType(maxNodesPerCell);

        // д�� Matrix Zone Header (BLOCK ��ʽ)
        file << "ZONE T=\"Matrix\", SOLUTIONTIME=" << time << ", DATAPACKING=BLOCK, "
            << "NODES=" << numMatrixNodes << ", ELEMENTS=" << numMatrixCells << ", "
            << "ZONETYPE=" << matrixZoneType;

        // ������Щ�������ڵ�Ԫ���� (X,Y Ϊ Node������Ϊ Cell-Centered)
        if (!varNames.empty()) {
            file << ", VARLOCATION=([3-" << 2 + varNames.size() << "]=CELLCENTERED)";
        }
        file << "\n";

        // (1.1) ������ҽڵ� X ����
        for (const auto& node : matrixNodes) {
            file << node.coord.m_x << "\n";
        }
        // (1.2) ������ҽڵ� Y ����
        for (const auto& node : matrixNodes) {
            file << node.coord.m_y << "\n";
        }

        // (1.3) ������ұ����� (CELLCENTERED)
        for (const auto& name : varNames) {
            auto field = fieldMgr_.getMatrixScalar(name);
            if (field) {
                // �������������������������ʵ����
                for (size_t i = 0; i < numMatrixCells; ++i) {
                    file << field->data[i] << "\n";
                }
            }
            else {
                // ��������ڣ�����ֻ���ѷ����е����ԣ����������
                for (size_t i = 0; i < numMatrixCells; ++i) {
                    file << 0.0 << "\n";
                }
            }
        }

        // (1.4) �������������ͨ�� (Connectivity)
        // ���� Gmsh Global Node ID �� �ֲ���������� 1-based Index ӳ��
        std::unordered_map<int, int> gmshIdToLocalIdx;
        for (size_t i = 0; i < numMatrixNodes; ++i) {
            gmshIdToLocalIdx[matrixNodes[i].id] = static_cast<int>(i) + 1; // Tecplot �� 1-based
        }

        for (const auto& cell : matrixCells) {
            size_t nNodes = cell.CellNodeIDs.size();
            for (int i = 0; i < maxNodesPerCell; ++i) {
                // �����ǰ��Ԫ�������Σ�3���ڵ㣩�������� Zone ���ı��Σ�max=4�������ظ����һ���ڵ�
                int nodeGmshId = (i < nNodes) ? cell.CellNodeIDs[i] : cell.CellNodeIDs.back();

                // ��ȫ��ȡ����� 1-based ����
                auto it = gmshIdToLocalIdx.find(nodeGmshId);
                if (it != gmshIdToLocalIdx.end()) {
                    file << it->second << " ";
                }
                else {
                    file << "1 "; // ���˷����󱸷���
                }
            }
            file << "\n";
        }
    }

    // =========================================================
    // Zone 2: Fracture Domain (�ѷ��� - 1D Segments)
    // =========================================================
    const auto& fractures = meshMgr_.fracture_network().fractures;

    // ͳ���ѷ�΢Ԫ��
    size_t totalFracCells = 0;
    for (const auto& frac : fractures) {
        totalFracCells += frac.elements.size();
    }

    // ���ö����ڵ���ԣ�ÿ��΢Ԫ��̬���� 2 ���˵㣬�ܽڵ���Ϊ΢Ԫ�� * 2
    size_t totalFracNodes = totalFracCells * 2;

    if (totalFracCells > 0) {
        // д�� Fracture Zone Header
        file << "ZONE T=\"Fractures\", SOLUTIONTIME=" << time << ", DATAPACKING=BLOCK, "
            << "NODES=" << totalFracNodes << ", ELEMENTS=" << totalFracCells << ", "
            << "ZONETYPE=" << GetTecplotElementType(2);

        if (!varNames.empty()) {
            file << ", VARLOCATION=([3-" << 2 + varNames.size() << "]=CELLCENTERED)";
        }
        file << "\n";

        // �ռ����е������ڵ����� (������ѭ������ I/O ��λ)
        std::vector<double> fracX;
        std::vector<double> fracY;
        fracX.reserve(totalFracNodes);
        fracY.reserve(totalFracNodes);

        for (const auto& frac : fractures) {
            // �������ѷ�� X, Y ���
            double dx = frac.end.m_x - frac.start.m_x;
            double dy = frac.end.m_y - frac.start.m_y;

            for (const auto& cell : frac.elements) {
                // ���ڹ�һ������ [param0, param1] ���Բ�ֵ����΢Ԫ�˵�����
                fracX.push_back(frac.start.m_x + dx * cell.param0); // �ڵ�1 X
                fracX.push_back(frac.start.m_x + dx * cell.param1); // �ڵ�2 X

                fracY.push_back(frac.start.m_y + dy * cell.param0); // �ڵ�1 Y
                fracY.push_back(frac.start.m_y + dy * cell.param1); // �ڵ�2 Y
            }
        }

        // (2.1) ����ѷ�ڵ� X ����
        for (double x : fracX) file << x << "\n";

        // (2.2) ����ѷ�ڵ� Y ����
        for (double y : fracY) file << y << "\n";

        // (2.3) ����ѷ������ (CELLCENTERED)
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

        // (2.4) ����ѷ���ͨ�� (Connectivity)
        // �����������ˣ��� i ��΢Ԫ�ϸ����ӽڵ� (2*i + 1) �� (2*i + 2) [TecplotҪ��1-based]
        for (size_t i = 0; i < totalFracCells; ++i) {
            file << (2 * i + 1) << " " << (2 * i + 2) << "\n";
        }
    }

    file.close();
    std::cout << ">>> 2D Tecplot export completed successfully." << std::endl;
}

// ==============================================================================
// ������ ParaView (VTK)
// ==============================================================================
void PostProcess_2D::ExportVTK(const std::string& filename, double time) const
{
    std::cout << "[PostProcess_2D] Exporting ParaView VTK file: " << filename << " (t=" << time << ")" << std::endl;

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "[Error] Cannot open file: " << filename << std::endl;
        return;
    }

    out << std::scientific << std::setprecision(6);
    std::vector<std::string> varNames = GetAllUniqueFieldNames();

    const auto& matrixNodes = meshMgr_.mesh().getNodes();
    const auto& matrixCells = meshMgr_.mesh().getCells();
    size_t numMatrixNodes = matrixNodes.size();
    size_t numMatrixCells = matrixCells.size();

    // �������ҽڵ㵽 0-based ������ӳ�� (VTK ����� 0 ��ʼ)
    std::unordered_map<int, int> nodeID2Index;
    for (size_t i = 0; i < numMatrixNodes; ++i) {
        nodeID2Index[matrixNodes[i].id] = static_cast<int>(i);
    }

    const auto& fractures = meshMgr_.fracture_network().fractures;
    size_t totalFracCells = 0;
    for (const auto& frac : fractures) {
        totalFracCells += frac.elements.size();
    }
    // 1D �ѷ��� VTK ����Ҫ�����������ڵ�����Ⱦ�߶�
    size_t totalFracNodes = totalFracCells * 2;

    size_t totalNodes = numMatrixNodes + totalFracNodes;
    size_t totalCells = numMatrixCells + totalFracCells;

    // ���� CELLS �б����ܳ���
    size_t cellListSize = 0;
    for (const auto& cell : matrixCells) {
        cellListSize += (1 + cell.CellNodeIDs.size());
    }
    cellListSize += totalFracCells * 3; // �ѷ��߶�: 1(�ڵ���2) + 2(�ڵ�����) = 3

    // 1. VTK Header
    out << "# vtk DataFile Version 3.0\n";
    out << "2D-EDFM Simulation Result (Time: " << time << ")\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // 2. POINTS
    out << "POINTS " << totalNodes << " double\n";
    int col = 0;
    for (const auto& node : matrixNodes) {
        out << node.coord.m_x << " " << node.coord.m_y << " 0.0 "; // 2D ���� Z=0
        if (++col % 3 == 0) out << "\n";
    }
    for (const auto& frac : fractures) {
        double dx = frac.end.m_x - frac.start.m_x;
        double dy = frac.end.m_y - frac.start.m_y;
        for (const auto& cell : frac.elements) {
            out << (frac.start.m_x + dx * cell.param0) << " " << (frac.start.m_y + dy * cell.param0) << " 0.0 ";
            if (++col % 3 == 0) out << "\n";
            out << (frac.start.m_x + dx * cell.param1) << " " << (frac.start.m_y + dy * cell.param1) << " 0.0 ";
            if (++col % 3 == 0) out << "\n";
        }
    }
    if (col % 3 != 0) out << "\n";

    // 3. CELLS
    out << "CELLS " << totalCells << " " << cellListSize << "\n";
    for (const auto& cell : matrixCells) {
        out << cell.CellNodeIDs.size() << " ";
        for (int id : cell.CellNodeIDs) {
            auto it = nodeID2Index.find(id);
            if (it != nodeID2Index.end()) {
                out << it->second << " ";
            }
            else {
                // [�ϸ�ģʽ] ����ȱʧ�ڵ�ֱ�ӱ����۶ϣ������ڸ��������˴���
                out.close();
                throw std::runtime_error("[PostProcess_2D Error] Node ID " + std::to_string(id) + " not found in Matrix Mesh!");
            }
        }
        out << "\n";
    }

    size_t fracNodeOffset = numMatrixNodes;
    for (size_t i = 0; i < totalFracCells; ++i) {
        out << "2 " << fracNodeOffset << " " << (fracNodeOffset + 1) << "\n";
        fracNodeOffset += 2;
    }

    // 4. CELL_TYPES
    out << "CELL_TYPES " << totalCells << "\n";
    for (const auto& cell : matrixCells) {
        size_t n = cell.CellNodeIDs.size();
        if (n == 3) {
            out << "5\n";      // VTK_TRIANGLE
        }
        else if (n == 4) {
            out << "9\n";      // VTK_QUAD
        }
        else {
            // [�ϸ�ģʽ] �ܾ� 2D Matrix �г��ַ� 3/4 �ڵ���˻������ε�Ԫ����ֹ��������ģ��
            out.close();
            throw std::runtime_error("[PostProcess_2D Error] Invalid matrix cell detected with " + std::to_string(n) + " nodes. Only Triangle (3) and Quad (4) are strictly supported.");
        }
    }
    for (size_t i = 0; i < totalFracCells; ++i) {
        out << "3\n";                  // VTK_LINE
    }

    // 5. CELL_DATA (cell-centered fields for accurate well/fracture visualization)
    out << "CELL_DATA " << totalCells << "\n";

    // Domain tag (0: matrix, 1: fracture)
    out << "SCALARS DomainID int 1\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < numMatrixCells; ++i) out << "0\n";
    for (size_t i = 0; i < totalFracCells; ++i) out << "1\n";

    // Cell-centered physical fields (no interpolation — exact solver values)
    for (const auto& name : varNames) {
        out << "SCALARS cell_" << name << " double 1\nLOOKUP_TABLE default\n";
        auto matField = fieldMgr_.getMatrixScalar(name);
        for (size_t c = 0; c < numMatrixCells; ++c) {
            double val = (matField && c < matField->data.size()) ? matField->data[c] : 0.0;
            if (!std::isfinite(val)) val = 0.0;
            out << val << "\n";
        }
        auto fracField = fieldMgr_.getFractureScalar(name);
        for (size_t c = 0; c < totalFracCells; ++c) {
            double val = (fracField && c < fracField->data.size()) ? fracField->data[c] : 0.0;
            if (!std::isfinite(val)) val = 0.0;
            out << val << "\n";
        }
    }

    // 6. POINT_DATA (node-averaged interpolation)
    out << "POINT_DATA " << totalNodes << "\n";

    for (const auto& name : varNames) {
        out << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";

        // --- 1. ���ҽڵ㳡 (Cell-to-Node ƽ��) ---
        std::vector<double> matVals(numMatrixNodes, 0.0);
        std::vector<int> matCounts(numMatrixNodes, 0);
        auto matField = fieldMgr_.getMatrixScalar(name);

        if (matField) {
            for (size_t c = 0; c < numMatrixCells; ++c) {
                double val = (c < matField->data.size()) ? matField->data[c] : 0.0;
                if (!std::isfinite(val)) val = 0.0;
                // ����Ԫ���ĵ�ֵ�ۼӵ�����������нڵ���
                for (int id : matrixCells[c].CellNodeIDs) {
                    auto it = nodeID2Index.find(id);
                    if (it != nodeID2Index.end()) {
                        int pIdx = it->second;
                        matVals[pIdx] += val;
                        matCounts[pIdx]++;
                    }
                }
            }
            // ��ÿ���ڵ��ֵȡƽ��
            for (size_t p = 0; p < numMatrixNodes; ++p) {
                if (matCounts[p] > 0) matVals[p] /= matCounts[p];
            }
        }
        for (size_t p = 0; p < numMatrixNodes; ++p) {
            out << matVals[p] << "\n";
        }

        // --- 2. �ѷ�ڵ㳡 ---
        // 2D�������ѷ�΢Ԫ����ȫ�����Ķ˵�(totalFracNodes = totalFracCells * 2)
        // ÿ��΢Ԫ�����˵�ֱ�Ӽ̳и�΢Ԫ��ֵ����
        auto fracField = fieldMgr_.getFractureScalar(name);
        if (fracField) {
            for (size_t c = 0; c < totalFracCells; ++c) {
                double val = (c < fracField->data.size()) ? fracField->data[c] : 0.0;
                if (!std::isfinite(val)) val = 0.0;
                // ÿ���ѷ�΢Ԫ��Ӧ��� 2 ���ڵ�����
                out << val << "\n" << val << "\n";
            }
        }
        else {
            for (size_t c = 0; c < totalFracCells; ++c) {
                out << "0.0\n0.0\n";
            }
        }
    }

    out.close();
    std::cout << "[PostProcess_2D] VTK Export completed successfully with Point Interpolation.\n"; }

// ==============================================================================
// ExportVTU  (Step 4: ParaView XML Unstructured Grid for time-series animation)
// ==============================================================================
void PostProcess_2D::ExportVTU(const std::string& filename, double time) const
{
    std::cout << "[PostProcess_2D] Exporting VTU: " << filename << "  (t=" << time << ")" << std::endl;

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "[Error] Cannot open file: " << filename << std::endl;
        return;
    }
    out << std::scientific << std::setprecision(8);

    const std::vector<std::string> varNames = GetAllUniqueFieldNames();
    const auto& matrixNodes = meshMgr_.mesh().getNodes();
    const auto& matrixCells = meshMgr_.mesh().getCells();
    const size_t numMatrixNodes = matrixNodes.size();
    const size_t numMatrixCells = matrixCells.size();

    // 0-based node ID → point-array index (VTK requires 0-based)
    std::unordered_map<int, int> nodeID2Index;
    nodeID2Index.reserve(numMatrixNodes * 2);
    for (size_t i = 0; i < numMatrixNodes; ++i)
        nodeID2Index[matrixNodes[i].id] = static_cast<int>(i);

    const auto& fractures = meshMgr_.fracture_network().fractures;
    size_t totalFracCells = 0;
    for (const auto& frac : fractures) totalFracCells += frac.elements.size();
    const size_t totalFracNodes = totalFracCells * 2;  // 2 endpoints per segment

    const size_t totalPoints = numMatrixNodes + totalFracNodes;
    const size_t totalCells  = numMatrixCells + totalFracCells;

    // ── XML header ───────────────────────────────────────────────────────────
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <!-- Time = " << time << " s -->\n";
    out << "    <Piece NumberOfPoints=\"" << totalPoints
        << "\" NumberOfCells=\"" << totalCells << "\">\n";

    // ── Points ───────────────────────────────────────────────────────────────
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& nd : matrixNodes)
        out << "          " << nd.coord.m_x << " " << nd.coord.m_y << " 0.0\n";
    for (const auto& frac : fractures) {
        const double dx = frac.end.m_x - frac.start.m_x;
        const double dy = frac.end.m_y - frac.start.m_y;
        for (const auto& elem : frac.elements) {
            out << "          "
                << (frac.start.m_x + dx * elem.param0) << " "
                << (frac.start.m_y + dy * elem.param0) << " 0.0\n";
            out << "          "
                << (frac.start.m_x + dx * elem.param1) << " "
                << (frac.start.m_y + dy * elem.param1) << " 0.0\n";
        }
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    // ── Cells ────────────────────────────────────────────────────────────────
    out << "      <Cells>\n";

    // connectivity
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& cell : matrixCells) {
        out << "          ";
        for (int id : cell.CellNodeIDs) {
            auto it = nodeID2Index.find(id);
            out << (it != nodeID2Index.end() ? it->second : 0) << " ";
        }
        out << "\n";
    }
    size_t fracPtBase = numMatrixNodes;
    for (size_t i = 0; i < totalFracCells; ++i) {
        out << "          " << (fracPtBase + 2*i) << " " << (fracPtBase + 2*i + 1) << "\n";
    }
    out << "        </DataArray>\n";

    // offsets (cumulative node count per cell)
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    out << "          ";
    int offset = 0;
    for (const auto& cell : matrixCells) {
        offset += static_cast<int>(cell.CellNodeIDs.size());
        out << offset << " ";
    }
    for (size_t i = 0; i < totalFracCells; ++i) {
        offset += 2;
        out << offset << " ";
    }
    out << "\n        </DataArray>\n";

    // types
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    out << "          ";
    for (const auto& cell : matrixCells) {
        const size_t n = cell.CellNodeIDs.size();
        if      (n == 3) out << "5 ";   // VTK_TRIANGLE
        else if (n == 4) out << "9 ";   // VTK_QUAD
        else             out << "7 ";   // VTK_POLYGON (fallback)
    }
    for (size_t i = 0; i < totalFracCells; ++i) out << "3 "; // VTK_LINE
    out << "\n        </DataArray>\n";

    out << "      </Cells>\n";

    // ── CellData ─────────────────────────────────────────────────────────────
    out << "      <CellData>\n";

    // DomainID (0=matrix, 1=fracture)
    out << "        <DataArray type=\"Int32\" Name=\"DomainID\" format=\"ascii\">\n";
    out << "          ";
    for (size_t i = 0; i < numMatrixCells; ++i) out << "0 ";
    for (size_t i = 0; i < totalFracCells;  ++i) out << "1 ";
    out << "\n        </DataArray>\n";

    // Physical fields
    for (const auto& name : varNames) {
        out << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
        out << "          ";

        auto matField  = fieldMgr_.getMatrixScalar(name);
        auto fracField = fieldMgr_.getFractureScalar(name);

        for (size_t i = 0; i < numMatrixCells; ++i) {
            double v = (matField && i < matField->data.size()) ? matField->data[i] : 0.0;
            if (!std::isfinite(v)) v = 0.0;
            out << v << " ";
        }
        size_t fi = 0;
        for (const auto& frac : fractures) {
            for (size_t ei = 0; ei < frac.elements.size(); ++ei, ++fi) {
                double v = (fracField && fi < fracField->data.size()) ? fracField->data[fi] : 0.0;
                if (!std::isfinite(v)) v = 0.0;
                out << v << " ";
            }
        }
        out << "\n        </DataArray>\n";
    }
    out << "      </CellData>\n";

    // ── Close ────────────────────────────────────────────────────────────────
    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";

    out.close();
    std::cout << "[PostProcess_2D] VTU export completed: " << filename << std::endl;
}

// ==============================================================================
// ExportPVD  (Step 4: ParaView time-series collection file)
// ==============================================================================
/*static*/
void PostProcess_2D::ExportPVD(const std::string& pvd_filename,
                                const std::vector<std::string>& vtu_filenames,
                                const std::vector<double>& times)
{
    if (vtu_filenames.size() != times.size()) {
        std::cerr << "[PostProcess_2D] ExportPVD: vtu_filenames and times size mismatch!\n";
        return;
    }
    std::ofstream out(pvd_filename);
    if (!out.is_open()) {
        std::cerr << "[PostProcess_2D] ExportPVD: cannot open " << pvd_filename << "\n";
        return;
    }
    out << std::scientific << std::setprecision(8);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <Collection>\n";
    for (size_t i = 0; i < vtu_filenames.size(); ++i) {
        out << "    <DataSet timestep=\"" << times[i]
            << "\" group=\"\" part=\"0\" file=\"" << vtu_filenames[i] << "\"/>\n";
    }
    out << "  </Collection>\n";
    out << "</VTKFile>\n";
    out.close();
    std::cout << "[PostProcess_2D] PVD collection written: " << pvd_filename
              << " (" << times.size() << " snapshots)\n";
}

// ==============================================================================
// �ڲ���������
// ==============================================================================

std::vector<std::string> PostProcess_2D::GetAllUniqueFieldNames() const
{
    std::set<std::string> uniqueNames;

    // ��ȡ����������˫���ȱ�����
    for (const auto& pair : fieldMgr_.matrixFields.fields) {
        // ���� dynamic_pointer_cast ��ȫ��� volScalarField (�ų� Vector/ADVar ��)
        if (std::dynamic_pointer_cast<volScalarField>(pair.second)) {
            uniqueNames.insert(pair.first);
        }
    }

    // ��ȡ�ѷ�������˫���ȱ�����
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
        return "FELINESEG";         // 1D �߶� (�ѷ�)
    }
    else if (numNodes == 3) {
        return "FETRIANGLE";        // �������λ���
    }
    else if (numNodes == 4) {
        return "FEQUADRILATERAL";   // �ı��λ��� �� ���(������+�ı���)
    }

    // ��������
    return "FEQUADRILATERAL";
}
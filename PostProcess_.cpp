#include "PostProcess_.h"

#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>

// 2) 便捷封装定义
bool getFaceValueFromCellValue_T(
    MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
    const TemperatureBCAdapter& Tbc,
    const std::string& T_cell, const std::string& T_face_out,
    const std::vector<Vector>* gradT
) {
    return getFaceValueFromCellValue_BC<TemperatureBCAdapter>(
        mgr, reg, freg, Tbc, T_cell, T_face_out, gradT);
}

bool getFaceValueFromCellValue_P(
    MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
    const PressureBCAdapter& Pbc,
    const std::string& p_cell, const std::string& p_face_out,
    const std::vector<Vector>* gradp
) {
    return getFaceValueFromCellValue_BC<PressureBCAdapter>(
        mgr, reg, freg, Pbc, p_cell, p_face_out, gradp);
}

// 3) no-BC cell-to-face averaging
bool getFaceValueFromCellValue_plain(
    MeshManager& mgr,
    const FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const std::string& cellFieldName,
    const std::string& faceFieldName)
{
    Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
    const auto& faces = mesh.getFaces();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();

    auto phiC = reg.get<volScalarField>(cellFieldName);
    if (!phiC)
    {
        std::cerr << "[getFaceValueFromCellValue_plain] missing cell field: " << cellFieldName << "\n";
        return false;
    }
    if (phiC->size != cells.size())
    {
        std::cerr << "[getFaceValueFromCellValue_plain] size mismatch: " << cellFieldName << " vs cells\n";
        return false;
    }

    auto phiF = freg.getOrCreate<faceScalarField>(faceFieldName, faces.size(), 0.0);
    for (const auto& F : faces)
    {
        const int iF = F.id - 1;
        auto itOwner = id2idx.find(F.ownerCell);
        if (itOwner == id2idx.end())
        {
            std::cerr << "[getFaceValueFromCellValue_plain] missing owner cell id " << F.ownerCell << "\n";
            return false;
        }

        const double phiOwner = (*phiC)[itOwner->second];
        if (!F.isBoundary() && F.neighborCell >= 0)
        {
            auto itNei = id2idx.find(F.neighborCell);
            if (itNei == id2idx.end())
            {
                std::cerr << "[getFaceValueFromCellValue_plain] missing neighbor cell id " << F.neighborCell << "\n";
                return false;
            }
            double gamma = F.f_linearInterpolationCoef;
            if (!(gamma > 0.0 && gamma < 1.0)) gamma = 0.5;
            const double phiNei = (*phiC)[itNei->second];
            (*phiF)[iF] = (1.0 - gamma) * phiOwner + gamma * phiNei;
        }
        else
        {
            (*phiF)[iF] = phiOwner;
        }
    }

    return true;
}

// 4) face → node 定义（优先边界面平均）
bool getNodeValueFromFaceValue(
    MeshManager& mgr,
    FaceFieldRegistry& freg,
    const std::string& faceFieldName,
    std::vector<double>& nodeValsOut
)
{
    Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
    const auto& nodes = mesh.getNodes();
    const auto& faces = mesh.getFaces();

    auto phiF = freg.get<faceScalarField>(faceFieldName);
    if (!phiF) { std::cerr << "[getNodeValueFromFaceValue] missing face field: " << faceFieldName << "\n"; return false; }
    if (phiF->size != faces.size()) { std::cerr << "[getNodeValueFromFaceValue] size mismatch in " << faceFieldName << "\n"; return false; }

    nodeValsOut.assign(nodes.size(), 0.0);

    // nodeId → 索引
    std::unordered_map<int, size_t> nodeId2idx;
    nodeId2idx.reserve(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) nodeId2idx[nodes[i].id] = i;

    // node → faces 邻接
    std::vector<std::vector<int>> node2faces(nodes.size());
    for (const auto& F : faces) {
        const int iF = F.id - 1;
        for (int nid : F.FaceNodeIDs) {
            auto it = nodeId2idx.find(nid);
            if (it != nodeId2idx.end()) node2faces[it->second].push_back(iF);
        }
    }

    // 取值
    for (size_t in = 0; in < nodes.size(); ++in)
    {
        const auto& adj = node2faces[in];
        if (adj.empty()) { nodeValsOut[in] = 0.0; continue; }

        std::vector<int> bndAdj;
        bndAdj.reserve(adj.size());
        for (int iF : adj) if (faces[iF].isBoundary()) bndAdj.push_back(iF);

        const std::vector<int>& pick = bndAdj.empty() ? adj : bndAdj;

        double sum = 0.0; int cnt = 0;
        for (int iF : pick) { sum += (*phiF)[iF]; ++cnt; }
        nodeValsOut[in] = (cnt > 0 ? sum / cnt : 0.0);
    }
    return true;
}

// 4) 一条龙：cell→face(含ABC+梯度缓冲)→node→Tecplot(三角)
bool outputTecplot_cellToFaceToNode_BC(
    MeshManager& mgr,
    const FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const TemperatureBCAdapter* Tbc,
    const PressureBCAdapter* Pbc,
    const std::string& cellFieldName,
    const std::string& faceFieldName,
    const std::vector<Vector>* gradBuf,
    const std::string& outFilename
) {
    // 1) cell -> face
    bool ok = false;
    if (Tbc) ok = getFaceValueFromCellValue_T(mgr, reg, freg, *Tbc, cellFieldName, faceFieldName, gradBuf);
    else if (Pbc) ok = getFaceValueFromCellValue_P(mgr, reg, freg, *Pbc, cellFieldName, faceFieldName, gradBuf);
    else ok = getFaceValueFromCellValue_plain(mgr, reg, freg, cellFieldName, faceFieldName);
    if (!ok) return false;

    // 2) face -> node
    std::vector<double> nodeVals;
    if (!getNodeValueFromFaceValue(mgr, freg, faceFieldName, nodeVals)) return false;

    // 3) 写 Tecplot（关键修补：将 CellNodeIDs 的“节点ID”映射为 1-based 顺序下标）
    Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
    const auto& nodes = mesh.getNodes();
    const auto& cells = mesh.getCells();

    if (nodeVals.size() != nodes.size()) {
        std::cerr << "[Tecplot] node value size mismatch.\n";
        return false;
    }

    // —— 建立 nodeId -> 顺序下标（0-based），然后写连通时 +1 —— //
    //    nodes[i].id  可能是任意 ID； Tecplot 需要的是 “第 i 行” 的 1-based 下标
    std::unordered_map<int, int> nodeId2Seq;
    nodeId2Seq.reserve(nodes.size());
    for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
        nodeId2Seq[nodes[i].id] = i;  // 0-based
    }

    std::ofstream ofs(outFilename, std::ios::out);
    if (!ofs) { std::cerr << "[Tecplot] cannot open file: " << outFilename << "\n"; return false; }

    ofs << "Variables = x, y, T\n";
    ofs << "Zone  n=" << nodes.size()
        << "  e=" << cells.size()
        << "  f=fepoint  et=triangle\n";

    // 写节点表：坐标 + 节点值（与 nodes 的顺序一致）
    for (size_t in = 0; in < nodes.size(); ++in) {
        // 你的 Vector 成员可能是 x/y 或 m_x/m_y；按你工程里一致的那个来：
        // 如果是 .x/.y 就改成 nodes[in].coord.x / .y
        ofs << nodes[in].coord.m_x << "  " << nodes[in].coord.m_y << "  " << nodeVals[in] << "\n";
    }

    // 写单元连通：将 “节点ID” 转换为 “顺序下标 + 1（Tecplot 1-based）”
    for (const auto& c : cells) {
        for (size_t k = 0; k < c.CellNodeIDs.size(); ++k) {
            const int nodeId = c.CellNodeIDs[k];
            auto it = nodeId2Seq.find(nodeId);
            if (it == nodeId2Seq.end()) {
                std::cerr << "[Tecplot] cell " << c.id << " references unknown nodeId " << nodeId << "\n";
                ofs.close();
                return false;
            }
            const int tecIndex = it->second + 1; // 1-based
            ofs << tecIndex << (k + 1 == c.CellNodeIDs.size() ? '\n' : ' ');
        }
    }

    ofs.close();
    return true;
}



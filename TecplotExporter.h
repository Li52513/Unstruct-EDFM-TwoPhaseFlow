#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "BCAdapter.h"
#include "FaceFieldRegistry.h"
#include "FieldRegistry.h"
#include "MeshManager.h"
#include "TemperatureBCAdapter.h"

#if __cplusplus >= 201703L
#include <filesystem>
#else
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif
#endif

namespace TecplotExport 
{

    // ---------------- Options ----------------
    struct TecplotOptions {
        std::string folder;
        std::string filePrefix = "field";
        int precision = 16;
    };

    // ---------------- FS helpers ----------------
#if __cplusplus >= 201703L
    inline void ensureParentDirExists(const std::string& filepath) {
        try { std::filesystem::create_directories(std::filesystem::path(filepath).parent_path()); }
        catch (...) {}
    }
    inline void ensureDirExists(const std::string& dirpath) {
        try { std::filesystem::create_directories(dirpath); }
        catch (...) {}
    }
#else
    inline void ensureParentDirExists(const std::string& filepath) {
        std::string path = filepath;
        size_t pos = path.find_last_of("/\\");
        if (pos == std::string::npos) return;
        path = path.substr(0, pos);
#ifdef _WIN32
        std::string acc;
        for (char ch : path) { acc.push_back(ch); if (ch == '\\' || ch == '/') _mkdir(acc.c_str()); }
        _mkdir(path.c_str());
#else
        std::string acc;
        for (char ch : path) { acc.push_back(ch); if (ch == '/') mkdir(acc.c_str(), 0755); }
        mkdir(path.c_str(), 0755);
#endif
    }
    inline void ensureDirExists(const std::string& dirpath) {
#ifdef _WIN32
        std::string acc;
        for (char ch : dirpath) { acc.push_back(ch); if (ch == '\\' || ch == '/') (void)_mkdir(acc.c_str()); }
        _mkdir(dirpath.c_str());
#else
        std::string acc;
        for (char ch : dirpath) { acc.push_back(ch); if (ch == '/') mkdir(acc.c_str(), 0755); }
        mkdir(dirpath.c_str(), 0755);
#endif
    }
#endif

    // ---------------- id maps ----------------
    inline std::unordered_map<int, size_t> buildNodeIdToIndex(const Mesh& mesh) {
        std::unordered_map<int, size_t> id2idx;
        const auto& nodes = mesh.getNodes();
        id2idx.reserve(nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i) id2idx[nodes[i].id] = i;
        return id2idx;
    }

    // ===================================================================
    // 1) cell → face：内部线性插值；边界按 ABC（与离散一致：d⊥=|ownerToNeighbor|）
    // ===================================================================
    inline std::shared_ptr<faceScalarField> computeFaceScalarFromCells(
        MeshManager& mgr,
        const FieldRegistry& reg,
        const std::string& cellField,
        FaceFieldRegistry& freg,
        const std::string& faceField,
        const PressureBCAdapter* Pbc,
        const TemperatureBCAdapter* Tbc)
    {
        const double eps = 1e-12;
        auto cellValues = reg.get<volScalarField>(cellField.c_str());
        if (!cellValues) {
            std::cerr << "[TecplotExport] missing cell field " << cellField << "\n";
            return nullptr;
        }

        Mesh& mesh = mgr.mesh();
        const auto& faces = mesh.getFaces();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        auto facesField = freg.getOrCreate<faceScalarField>(faceField, faces.size(), 0.0);

        for (const auto& F : faces) {
            const size_t fid = static_cast<size_t>(F.id - 1);
            auto itP = id2idx.find(F.ownerCell);
            if (itP == id2idx.end()) { (*facesField)[fid] = 0.0; continue; }
            const size_t iP = static_cast<size_t>(itP->second);
            const double phiP = (*cellValues)[iP];

            if (!F.isBoundary()) {
                auto itN = id2idx.find(F.neighborCell);
                if (itN == id2idx.end()) { (*facesField)[fid] = phiP; continue; }
                const size_t iN = static_cast<size_t>(itN->second);

                double gamma = F.f_linearInterpolationCoef;
                if (!(gamma >= 0.0 && gamma <= 1.0)) {
                    const Vector& CP = cells[iP].center;
                    const Vector& CN = cells[iN].center;
                    const Vector dON = CN - CP;
                    const double D = dON.Mag();
                    if (D > eps) {
                        const Vector eON = dON / D;
                        const double s = (F.midpoint - CP) * eON;
                        gamma = std::min(1.0, std::max(0.0, s / D));
                    }
                    else gamma = 0.5;
                }
                const double phiN = (*cellValues)[iN];
                (*facesField)[fid] = (1.0 - gamma) * phiP + gamma * phiN;
            }
            else {
                double phiF = phiP;
                double a = 0.0, b = 0.0, c = 0.0;
                bool hasBC = false;
                if (Pbc) hasBC = Pbc->getABC(F.id, a, b, c);
                else if (Tbc) hasBC = Tbc->getABC(F.id, a, b, c);

                if (hasBC) {
                    const double dperp = std::max(F.ownerToNeighbor.Mag(), eps);
                    if (std::abs(a) < eps && std::abs(b) < eps) {
                        phiF = phiP;
                    }
                    else if (std::abs(a) < eps) {
                        phiF = phiP + c * dperp;                       // Neumann
                    }
                    else {                                           // Robin
                        const double denom = a + b / dperp;
                        const double safe = (std::abs(denom) < eps) ? (denom >= 0.0 ? eps : -eps) : denom;
                        phiF = (c + (b / dperp) * phiP) / safe;
                    }
                }
                (*facesField)[fid] = phiF;
            }
        }
        return facesField;
    }

    // ===================================================================
    // 2) face → node：节点若连着边界面，只平均边界面；否则平均所有相邻面
    // ===================================================================
    inline std::vector<double> computeNodeScalarFromFaces_boundaryAware(
        const Mesh& mesh,
        const faceScalarField& faceValues)
    {
        const auto& nodes = mesh.getNodes();
        const auto& faces = mesh.getFaces();

        std::unordered_map<int, size_t> nodeId2idx = buildNodeIdToIndex(mesh);
        std::vector<std::vector<int>> nodeFaces(nodes.size());
        std::vector<std::vector<int>> nodeBdryFaces(nodes.size());

        for (const auto& F : faces) {
            const int fid0 = F.id - 1;
            for (int nid : F.FaceNodeIDs) {
                auto it = nodeId2idx.find(nid);
                if (it == nodeId2idx.end()) continue;
                const size_t nidx = it->second;
                nodeFaces[nidx].push_back(fid0);
                if (F.isBoundary()) nodeBdryFaces[nidx].push_back(fid0);
            }
        }

        std::vector<double> nodal(nodes.size(), 0.0);
        for (size_t i = 0; i < nodes.size(); ++i) {
            const auto& pool = (!nodeBdryFaces[i].empty() ? nodeBdryFaces[i] : nodeFaces[i]);
            if (pool.empty()) { nodal[i] = 0.0; continue; }
            double s = 0.0;
            for (int fid0 : pool) s += faceValues.data[static_cast<size_t>(fid0)];
            nodal[i] = s / double(pool.size());
        }
        return nodal;
    }

    // ===================================================================
    // 3) 写 Tecplot FEPOINT（带单元拓扑）
    // ===================================================================
    inline bool writeTecplotSingleVariable(
        const Mesh& mesh,
        const std::vector<double>& nodal,
        const std::string& varName,
        const TecplotOptions& opt,
        int step,
        double time)
    {
        if (nodal.size() != mesh.getNodes().size()) {
            std::cerr << "[TecplotExport] nodal size mismatch for " << varName << "\n";
            return false;
        }

        ensureDirExists(opt.folder);
        std::ostringstream fn;
        fn << opt.folder << "/" << opt.filePrefix << "_" << varName
            << "_step_" << std::setw(6) << std::setfill('0') << step << ".plt";
        const std::string filename = fn.str();
        ensureParentDirExists(filename);

        std::ofstream ofs(filename);
        if (!ofs) { std::cerr << "[TecplotExport] cannot open " << filename << "\n"; return false; }

        const auto& nodes = mesh.getNodes();
        auto nodeId2order = buildNodeIdToIndex(mesh);

        struct Quad { int n1, n2, n3, n4; };
        std::vector<Quad> elements;
        elements.reserve(mesh.getCells().size());

        bool hasDegenerate = false;
        for (const auto& cell : mesh.getCells()) {
            if (cell.id < 0) continue;
            if (cell.CellNodeIDs.size() < 3) continue;

            std::array<int, 4> ids{ 0,0,0,0 };
            const size_t cnt = cell.CellNodeIDs.size();
            for (size_t k = 0; k < std::min<size_t>(cnt, 4); ++k) {
                auto it = nodeId2order.find(cell.CellNodeIDs[k]);
                ids[k] = (it == nodeId2order.end() ? 1 : int(it->second) + 1);
            }
            if (cnt == 3) { ids[3] = ids[2]; hasDegenerate = true; }
            else if (cnt > 4) { hasDegenerate = true; }

            elements.push_back({ ids[0],ids[1],ids[2],ids[3] });
        }

        ofs << std::scientific << std::setprecision(opt.precision);
        ofs << "TITLE = \"Tecplot export for " << varName << "\"\n";
        ofs << "VARIABLES = \"x\" \"y\" \"" << varName << "\"\n";
        ofs << "ZONE T=\"" << varName << "_step_" << step << "\" N=" << nodes.size()
            << " E=" << elements.size() << " F=FEPOINT ET=QUADRILATERAL\n";
        ofs << "# time = " << time << "\n";

        for (size_t i = 0; i < nodes.size(); ++i)
            ofs << nodes[i].coord.m_x << " " << nodes[i].coord.m_y << " " << nodal[i] << "\n";
        for (const auto& e : elements)
            ofs << e.n1 << " " << e.n2 << " " << e.n3 << " " << e.n4 << "\n";

        if (hasDegenerate) ofs << "# warning: degenerate quads present\n";
        return true;
    }

    // ===================================================================
    // 一条龙：cell→face→node→写 FEPOINT
    // ===================================================================
    inline bool export_cellField_to_tecplot_after_step_cell2face2node(
        MeshManager& mgr,
        const FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const TecplotOptions& opt,
        const std::string& cellField,
        const std::string& tmpFaceField,
        const std::string& varName,
        int step,
        double time,
        const PressureBCAdapter* Pbc,
        const TemperatureBCAdapter* Tbc)
    {
        auto faceField = computeFaceScalarFromCells(mgr, reg, cellField, freg, tmpFaceField, Pbc, Tbc);
        if (!faceField) return false;
        auto nodal = computeNodeScalarFromFaces_boundaryAware(mgr.mesh(), *faceField);
        return writeTecplotSingleVariable(mgr.mesh(), nodal, varName, opt, step, time);
    }

    // ===================================================================
    // 兼容旧名：保持老调用不改（内部转到 cell2face2node）
    // ===================================================================
    inline bool export_cellField_to_tecplot_after_step
    (
        MeshManager& mgr,
        const FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const TecplotOptions& opt,
        const std::string& cellField,
        const std::string& tmpFaceField,
        const std::string& varName,
        int step,
        double time,
        const PressureBCAdapter* Pbc,
        const TemperatureBCAdapter* Tbc)
    {
        return export_cellField_to_tecplot_after_step_cell2face2node(
            mgr, reg, freg, opt, cellField, tmpFaceField, varName, step, time, Pbc, Tbc);
    }

    // ===================================================================
    // 同时导出 P/T（可继续被你现有的 runTransient 调用）
    // ===================================================================
    inline bool export_pT_tecplot_after_step(
        MeshManager& mgr,
        const FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& Pbc,
        const TemperatureBCAdapter& Tbc,
        const std::string& pField,
        const std::string& TField,
        int step,
        double time,
        TecplotOptions opt)
    {
        if (opt.folder.empty()) opt.folder = "out_tec";
        if (opt.filePrefix.empty()) opt.filePrefix = "field";

        const bool okP = export_cellField_to_tecplot_after_step(
            mgr, reg, freg, opt, pField, pField + "_face_tmp", "P", step, time, &Pbc, nullptr);

        const bool okT = export_cellField_to_tecplot_after_step(
            mgr, reg, freg, opt, TField, TField + "_face_tmp", "T", step, time, nullptr, &Tbc);

        return okP && okT;
    }

} // namespace TecplotExport

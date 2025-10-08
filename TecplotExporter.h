#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
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

// TecplotExporter.h
// 提供从单元标量场构造面/节点值并导出 Tecplot FEPOINT 文件的工具。

namespace TecplotExport
{
    struct TecplotOptions
    {
        std::string folder;
        std::string filePrefix = "field";
        int precision = 16;
    };

#if __cplusplus >= 201703L
    inline void ensureParentDirExists(const std::string& filepath)
    {
        try
        {
            std::filesystem::create_directories(std::filesystem::path(filepath).parent_path());
        }
        catch (...)
        {
        }
    }

    inline void ensureDirExists(const std::string& dirpath)
    {
        try
        {
            std::filesystem::create_directories(dirpath);
        }
        catch (...)
        {
        }
    }
#else
#ifdef _WIN32
    inline void ensureParentDirExists(const std::string& filepath)
    {
        std::string path = filepath;
        size_t pos = path.find_last_of("/\\");
        if (pos == std::string::npos) return;
        path = path.substr(0, pos);

        std::string acc;
        for (char ch : path)
        {
            acc.push_back(ch);
            if (ch == '\\' || ch == '/') _mkdir(acc.c_str());
        }
        _mkdir(path.c_str());
    }

    inline void ensureDirExists(const std::string& dirpath)
    {
        std::string acc;
        for (char ch : dirpath)
        {
            acc.push_back(ch);
            if (ch == '\\' || ch == '/') (void)_mkdir(acc.c_str());
        }
        _mkdir(dirpath.c_str());
    }
#else
    inline void ensureParentDirExists(const std::string& filepath)
    {
        std::string path = filepath;
        size_t pos = path.find_last_of("/\\");
        if (pos == std::string::npos) return;
        path = path.substr(0, pos);
        std::string acc;
        for (char ch : path)
        {
            acc.push_back(ch);
            if (ch == '/') mkdir(acc.c_str(), 0755);
        }
        mkdir(path.c_str(), 0755);
    }

    inline void ensureDirExists(const std::string& dirpath)
    {
        std::string acc;
        for (char ch : dirpath)
        {
            acc.push_back(ch);
            if (ch == '/') mkdir(acc.c_str(), 0755);
        }
        mkdir(dirpath.c_str(), 0755);
    }
#endif
#endif

    inline std::unordered_map<int, size_t> buildNodeIdToIndex(const Mesh& mesh)
    {
        std::unordered_map<int, size_t> id2idx;
        const auto& nodes = mesh.getNodes();
        id2idx.reserve(nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            id2idx[nodes[i].id] = i;
        }
        return id2idx;
    }

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
        if (!cellValues)
        {
            std::cerr << "[TecplotExport] missing cell field " << cellField << "\n";
            return nullptr;
        }

        Mesh& mesh = mgr.mesh();
        const auto& faces = mesh.getFaces();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        auto facesField = freg.getOrCreate<faceScalarField>(faceField, faces.size(), 0.0);
        for (const auto& F : faces)
        {
            const size_t fid = static_cast<size_t>(F.id - 1);
            auto itP = id2idx.find(F.ownerCell);
            if (itP == id2idx.end())
            {
                (*facesField)[fid] = 0.0;
                continue;
            }
            const size_t idxP = static_cast<size_t>(itP->second);
            const double phiP = (*cellValues)[idxP];

            if (!F.isBoundary())
            {
                auto itN = id2idx.find(F.neighborCell);
                if (itN == id2idx.end())
                {
                    (*facesField)[fid] = phiP;
                    continue;
                }
                const size_t idxN = static_cast<size_t>(itN->second);
                double gamma = F.f_linearInterpolationCoef;

                if (!(gamma >= 0.0 && gamma <= 1.0))
                {
                    const Vector& CP = cells[idxP].center;
                    const Vector& CN = cells[idxN].center;
                    const Vector dON = CN - CP;
                    const double D = dON.Mag();
                    if (D > eps)
                    {
                        const Vector eON = dON / D;
                        const double s = (F.midpoint - CP) * eON;
                        gamma = std::min(1.0, std::max(0.0, s / D));
                    }
                    else
                    {
                        gamma = 0.5;
                    }
                }

                const double phiN = (*cellValues)[idxN];
                (*facesField)[fid] = (1.0 - gamma) * phiP + gamma * phiN;
            }
            else
            {
                double phiF = phiP;
                double a = 0.0, b = 0.0, c = 0.0;
                bool hasBC = false;
                if (Pbc)
                {
                    hasBC = Pbc->getABC(F.id, a, b, c);
                }
                else if (Tbc)
                {
                    hasBC = Tbc->getABC(F.id, a, b, c);
                }

                if (hasBC)
                {
                    const Vector& CP = cells[idxP].center;
                    const double dperp = std::max((F.midpoint - CP).Mag(), eps);
                    if (std::abs(a) < eps && std::abs(b) < eps)
                    {
                        phiF = phiP;
                    }
                    else if (std::abs(a) < eps)
                    {
                        phiF = phiP + c * dperp;
                    }
                    else
                    {
                        const double denom = a + b / dperp;
                        const double safeDenom = (std::abs(denom) < eps) ? (denom >= 0.0 ? eps : -eps) : denom;
                        phiF = (c + (b / dperp) * phiP) / safeDenom;
                    }
                }
                (*facesField)[fid] = phiF;
            }
        }

        return facesField;
    }

    inline std::vector<double> computeNodeScalarFromCells_IDW(
        const Mesh& mesh,
        const FieldRegistry& reg,
        const std::string& cellField,
        double power = 2.0,
        double eps = 1e-12)
    {
        std::vector<double> nodal(mesh.getNodes().size(), 0.0);
        auto cellValues = reg.get<volScalarField>(cellField.c_str());
        if (!cellValues)
        {
            return nodal;
        }

        const auto& nodes = mesh.getNodes();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        auto nodeId2idx = buildNodeIdToIndex(mesh);

        std::vector<std::vector<int>> nodeToCells(nodes.size());
        for (const auto& cell : cells)
        {
            if (cell.id < 0) continue;
            for (int nid : cell.CellNodeIDs)
            {
                auto it = nodeId2idx.find(nid);
                if (it != nodeId2idx.end())
                {
                    nodeToCells[it->second].push_back(cell.id);
                }
            }
        }

        for (size_t i = 0; i < nodes.size(); ++i)
        {
            const auto& adjacent = nodeToCells[i];
            if (adjacent.empty())
            {
                nodal[i] = 0.0;
                continue;
            }

            double numerator = 0.0;
            double denom = 0.0;
            const Vector& nodeCoord = nodes[i].coord;
            for (int cid : adjacent)
            {
                auto it = id2idx.find(cid);
                if (it == id2idx.end()) continue;
                const size_t cidx = static_cast<size_t>(it->second);
                const double dist = (cells[cidx].center - nodeCoord).Mag();
                const double w = 1.0 / std::pow(std::max(dist, eps), power);
                numerator += w * (*cellValues)[cidx];
                denom += w;
            }

            if (denom > 0.0)
            {
                nodal[i] = numerator / denom;
            }
            else
            {
                nodal[i] = 0.0;
            }
        }
        return nodal;
    }

    inline std::vector<double> computeNodeScalarFromFaces(
        const Mesh& mesh,
        const faceScalarField& faceValues,
        const std::vector<double>& fallback)
    {
        const auto& nodes = mesh.getNodes();
        std::vector<double> accum(nodes.size(), 0.0);
        std::vector<double> weights(nodes.size(), 0.0);

        auto nodeId2idx = buildNodeIdToIndex(mesh);
        for (const auto& face : mesh.getFaces())
        {
            const double phiF = faceValues.data[static_cast<size_t>(face.id - 1)];
            const double w = std::max(face.length, 1e-12);
            for (int nid : face.FaceNodeIDs)
            {
                auto it = nodeId2idx.find(nid);
                if (it == nodeId2idx.end()) continue;
                const size_t nidx = it->second;
                accum[nidx] += phiF * w;
                weights[nidx] += w;
            }
        }

        std::vector<double> nodal(nodes.size(), 0.0);
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            if (weights[i] > 0.0)
            {
                nodal[i] = accum[i] / weights[i];
            }
            else if (!fallback.empty())
            {
                nodal[i] = fallback[i];
            }
            else
            {
                nodal[i] = 0.0;
            }
        }
        return nodal;
    }

    inline std::vector<double> reconstructNodalScalar(
        MeshManager& mgr,
        const FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const std::string& cellField,
        const std::string& tmpFaceField,
        const PressureBCAdapter* Pbc,
        const TemperatureBCAdapter* Tbc)
    {
        auto faceField = computeFaceScalarFromCells(mgr, reg, cellField, freg, tmpFaceField, Pbc, Tbc);
        if (!faceField)
        {
            const Mesh& mesh = mgr.mesh();
            return computeNodeScalarFromCells_IDW(mesh, reg, cellField);
        }

        const Mesh& mesh = mgr.mesh();
        auto fallback = computeNodeScalarFromCells_IDW(mesh, reg, cellField);
        return computeNodeScalarFromFaces(mesh, *faceField, fallback);
    }

    inline bool writeTecplotSingleVariable(
        const Mesh& mesh,
        const std::vector<double>& nodal,
        const std::string& varName,
        const TecplotOptions& opt,
        int step,
        double time)
    {
        if (nodal.size() != mesh.getNodes().size())
        {
            std::cerr << "[TecplotExport] nodal values size mismatch for " << varName << "\n";
            return false;
        }

        ensureDirExists(opt.folder);
        std::ostringstream fn;
        fn << opt.folder << "/" << opt.filePrefix << "_" << varName << "_step_"
           << std::setw(6) << std::setfill('0') << step << ".plt";
        const std::string filename = fn.str();
        ensureParentDirExists(filename);

        std::ofstream ofs(filename);
        if (!ofs)
        {
            std::cerr << "[TecplotExport] cannot open " << filename << "\n";
            return false;
        }

        const auto& nodes = mesh.getNodes();
        auto nodeId2order = buildNodeIdToIndex(mesh);

        struct Quad
        {
            int n1;
            int n2;
            int n3;
            int n4;
        };

        std::vector<Quad> elements;
        elements.reserve(mesh.getCells().size());

        bool hasDegenerate = false;
        for (const auto& cell : mesh.getCells())
        {
            if (cell.id < 0) continue;
            if (cell.CellNodeIDs.size() < 3) continue;

            std::array<int, 4> ids{0, 0, 0, 0};
            const size_t count = cell.CellNodeIDs.size();
            for (size_t k = 0; k < std::min<size_t>(count, 4); ++k)
            {
                auto it = nodeId2order.find(cell.CellNodeIDs[k]);
                if (it == nodeId2order.end())
                {
                    ids[k] = 1;
                    continue;
                }
                ids[k] = static_cast<int>(it->second + 1);
            }

            if (count == 3)
            {
                ids[3] = ids[2];
                hasDegenerate = true;
            }
            else if (count >= 4)
            {
                if (count > 4) hasDegenerate = true;
            }
            else
            {
                continue;
            }

            elements.push_back({ids[0], ids[1], ids[2], ids[3]});
        }

        ofs << std::scientific << std::setprecision(opt.precision);
        ofs << "TITLE = \"Tecplot export for " << varName << "\"\n";
        ofs << "VARIABLES = \"x\" \"y\" \"" << varName << "\"\n";
        ofs << "ZONE T=\"" << varName << "_step_" << step << "\" N=" << nodes.size()
            << " E=" << elements.size() << " F=FEPOINT ET=QUADRILATERAL\n";
        ofs << "# time = " << time << "\n";

        for (size_t i = 0; i < nodes.size(); ++i)
        {
            ofs << nodes[i].coord.m_x << " " << nodes[i].coord.m_y << " " << nodal[i] << "\n";
        }

        for (const auto& e : elements)
        {
            ofs << e.n1 << " " << e.n2 << " " << e.n3 << " " << e.n4 << "\n";
        }

        if (hasDegenerate)
        {
            ofs << "# warning: some elements exported as degenerate quads (triangle or >4 nodes)\n";
        }

        ofs.close();
        return true;
    }

    inline bool export_cellField_to_tecplot_after_step(
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
        auto nodal = reconstructNodalScalar(mgr, reg, freg, cellField, tmpFaceField, Pbc, Tbc);
        return writeTecplotSingleVariable(mgr.mesh(), nodal, varName, opt, step, time);
    }

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
        opt.folder = opt.folder.empty() ? std::string("out_tec") : opt.folder;
        if (opt.filePrefix.empty()) opt.filePrefix = "field";

        TecplotOptions optP = opt;
        optP.filePrefix = opt.filePrefix;
        bool okP = export_cellField_to_tecplot_after_step(
            mgr, reg, freg, optP, pField, pField + "_face_tmp", "P", step, time, &Pbc, nullptr);

        TecplotOptions optT = opt;
        optT.filePrefix = opt.filePrefix;
        bool okT = export_cellField_to_tecplot_after_step(
            mgr, reg, freg, optT, TField, TField + "_face_tmp", "T", step, time, nullptr, &Tbc);

        return okP && okT;
    }
}

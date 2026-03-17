#pragma once

#include <stdexcept>
#include <algorithm>
#include "FIM_ConnectionManager.h"
#include "MeshManager.h"
#include "2D_FieldManager.h"

class FIM_TopologyBuilder2D {
public:
    static void LoadAllConnections(FIM_ConnectionManager& connMgr, const MeshManager& meshMgr, const FieldManager_2D& fieldMgr) {
        size_t expected = 0;

        const auto& faces = meshMgr.mesh().getFaces();
        for (const auto& face : faces) {
            if (!face.isBoundary()) ++expected;
        }

        const auto& fracs = meshMgr.fracture_network().fractures;
        for (const auto& frac : fracs) {
            if (frac.elements.size() > 1) expected += (frac.elements.size() - 1);
        }

        const auto& nncMap = meshMgr.getNNCTopologyMap();
        for (const auto& kv : nncMap) expected += kv.second.size();

        expected += fieldMgr.ff_topology.size();
        connMgr.ReserveRawConnections(expected);

        _loadMatrix(connMgr, meshMgr, fieldMgr);
        _loadFI(connMgr, meshMgr, fieldMgr);
        _loadNNC(connMgr, meshMgr, fieldMgr);
        _loadFF(connMgr, meshMgr, fieldMgr);
    }

private:
    static void _loadMatrix(FIM_ConnectionManager& connMgr, const MeshManager& meshMgr, const FieldManager_2D& fieldMgr) {
        auto pFlow = fieldMgr.getMatrixFaceScalar("T_Matrix_Flow");
        auto pHeat = fieldMgr.getMatrixFaceScalar("T_Matrix_Heat");
        if (!pFlow || !pHeat) return;

        const auto& faces = meshMgr.mesh().getFaces();
        const auto& cells = meshMgr.mesh().getCells();

        if (pFlow->data.size() != faces.size() || pHeat->data.size() != faces.size()) {
            throw std::runtime_error("[Build 2D] Matrix face field size mismatch with mesh faces.");
        }

        for (size_t i = 0; i < faces.size(); ++i) {
            if (faces[i].isBoundary()) continue;

            int nI = faces[i].ownerCell_index;
            int nJ = faces[i].neighborCell_index;
            if (nI < 0 || nJ < 0) continue;
            if (nI >= static_cast<int>(cells.size()) || nJ >= static_cast<int>(cells.size())) continue;

            double dist = (faces[i].midpoint - cells[nI].center).Mag()
                + (cells[nJ].center - faces[i].midpoint).Mag();

            connMgr.PushConnection(
                nI, nJ,
                pFlow->data[i], pHeat->data[i],
                faces[i].vectorE.Mag(), std::max(dist, 1e-12),
                ConnectionType::Matrix_Matrix
            );
        }
    }

    static void _loadFI(FIM_ConnectionManager& connMgr, const MeshManager& meshMgr, const FieldManager_2D& fieldMgr) {
        auto pFlow = fieldMgr.getFractureFaceScalar("T_FI_Flow");
        auto pHeat = fieldMgr.getFractureFaceScalar("T_FI_Heat");
        if (!pFlow || !pHeat) return;

        const auto& fracs = meshMgr.fracture_network().fractures;

        size_t expected = 0;
        for (const auto& frac : fracs) {
            if (frac.elements.size() > 1) expected += (frac.elements.size() - 1);
        }

        if (expected != pFlow->data.size() || expected != pHeat->data.size()) {
            throw std::runtime_error("[Build 2D] FI arrays mismatch with fracture topology.");
        }

        size_t fiIdx = 0;
        for (const auto& frac : fracs) {
            if (frac.elements.size() < 2) continue;

            for (size_t i = 0; i + 1 < frac.elements.size(); ++i) {
                int nodeI = frac.elements[i].solverIndex;
                int nodeJ = frac.elements[i + 1].solverIndex;
                double auxDist = 0.5 * (frac.elements[i].length + frac.elements[i + 1].length);

                connMgr.PushConnection(
                    nodeI, nodeJ,
                    pFlow->data[fiIdx], pHeat->data[fiIdx],
                    1.0, std::max(auxDist, 1e-12),
                    ConnectionType::Fracture_Internal
                );
                ++fiIdx;
            }
        }
    }

    static void _loadNNC(FIM_ConnectionManager& connMgr, const MeshManager& meshMgr, const FieldManager_2D& fieldMgr) {
        auto pFlow = fieldMgr.getNNCScalar("T_NNC_Flow");
        auto pHeat = fieldMgr.getNNCScalar("T_NNC_Heat");
        if (!pFlow || !pHeat) return;

        const auto& pairsMap = meshMgr.getNNCTopologyMap();

        size_t totalNNCCount = 0;
        for (const auto& kv : pairsMap) totalNNCCount += kv.second.size();

        if (totalNNCCount != pFlow->data.size() || totalNNCCount != pHeat->data.size()) {
            throw std::runtime_error("[Build 2D] NNC arrays mismatch with topology map.");
        }

        constexpr double thickness = 1.0;
        constexpr double eps = 1e-12;

        size_t idx = 0;
        for (const auto& kv : pairsMap) {
            int nodeI = kv.first;
            for (int nodeJ : kv.second) {
                const FractureElement* pElem = meshMgr.getFractureElementBySolverIndex(nodeJ);
                if (!pElem) {
                    throw std::runtime_error("[Build 2D] NNC solverIndex cannot map to fracture element.");
                }

                const double auxArea = std::max(pElem->length * thickness, eps);
                const double auxDist = std::max(pElem->avgDistance, eps);

                connMgr.PushConnection(
                    nodeI, nodeJ,
                    pFlow->data[idx], pHeat->data[idx],
                    auxArea, auxDist,
                    ConnectionType::Matrix_Fracture
                );
                ++idx;
            }
        }
    }

    static void _loadFF(FIM_ConnectionManager& connMgr, const MeshManager& meshMgr, const FieldManager_2D& fieldMgr) {
        auto pFlow = fieldMgr.getFFScalar("T_FF_Flow");
        auto pHeat = fieldMgr.getFFScalar("T_FF_Heat");
        if (!pFlow || !pHeat) return;

        const auto& ffPairs = fieldMgr.ff_topology;
        if (ffPairs.size() != pFlow->data.size() || ffPairs.size() != pHeat->data.size()) {
            throw std::runtime_error("[Build 2D] FF topology mismatch with FF field arrays.");
        }

        for (size_t i = 0; i < ffPairs.size(); ++i) {
            const int sI = ffPairs[i].first;
            const int sJ = ffPairs[i].second;

            const FractureElement* eI = meshMgr.getFractureElementBySolverIndex(sI);
            const FractureElement* eJ = meshMgr.getFractureElementBySolverIndex(sJ);
            if (!eI || !eJ) {
                throw std::runtime_error("[Build 2D] FF solverIndex cannot map to fracture element.");
            }

            const double mI = std::max(eI->length * std::max(eI->aperture, 1e-12), 1e-12);
            const double mJ = std::max(eJ->length * std::max(eJ->aperture, 1e-12), 1e-12);

            double auxArea = 0.5 * (mI + mJ);
            double auxDist = std::max(0.5 * (eI->length + eJ->length), 1e-12);

            connMgr.PushConnection(
                sI, sJ,
                pFlow->data[i], pHeat->data[i],
                auxArea, auxDist,
                ConnectionType::Fracture_Fracture
            );
        }
    }
};
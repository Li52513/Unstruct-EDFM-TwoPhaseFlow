#pragma once

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include "FIM_ConnectionManager.h"
#include "MeshManager.h"
#include "2D_FieldManager.h"
#include "SolverContrlStrName_op.h"

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
        const auto& faces = meshMgr.mesh().getFaces();
        const auto& cells = meshMgr.mesh().getCells();
        const PhysicalProperties_string_op::TransmissibilityFields tTags;
        auto pFlow = fieldMgr.getMatrixFaceScalar(tTags.matrix_flow);
        auto pHeat = fieldMgr.getMatrixFaceScalar(tTags.matrix_heat);
        if (!pFlow || !pHeat) {
            for (const auto& face : faces) {
                if (!face.isBoundary()) {
                    throw std::runtime_error("[Build 2D] Missing matrix transmissibility fields.");
                }
            }
            return;
        }

        if (pFlow->data.size() != faces.size() || pHeat->data.size() != faces.size()) {
            throw std::runtime_error("[Build 2D] Matrix face field size mismatch with mesh faces.");
        }

        for (size_t i = 0; i < faces.size(); ++i) {
            if (faces[i].isBoundary()) continue;

            int nI = faces[i].ownerCell_index;
            int nJ = faces[i].neighborCell_index;
            if (nI < 0 || nJ < 0) continue;
            if (nI >= static_cast<int>(cells.size()) || nJ >= static_cast<int>(cells.size())) continue;

            const Vector dO = faces[i].midpoint - cells[nI].center;
            const Vector dN = cells[nJ].center - faces[i].midpoint;
            const double nx = faces[i].normal.m_x;
            const double ny = faces[i].normal.m_y;
            const double dist =
                std::max(std::abs(dO.m_x * nx + dO.m_y * ny) +
                         std::abs(dN.m_x * nx + dN.m_y * ny), 1e-12);

            connMgr.PushConnection(
                nI, nJ,
                pFlow->data[i], pHeat->data[i],
                faces[i].vectorE.Mag(), dist,
                ConnectionType::Matrix_Matrix,
                faces[i].vectorT  // non-orthogonal correction vector (Step 2)
            );
        }
    }

    static void _loadFI(FIM_ConnectionManager& connMgr, const MeshManager& meshMgr, const FieldManager_2D& fieldMgr) {
        const auto& fracs = meshMgr.fracture_network().fractures;
        const PhysicalProperties_string_op::TransmissibilityFields tTags;
        auto pFlow = fieldMgr.getFractureFaceScalar(tTags.fi_flow);
        auto pHeat = fieldMgr.getFractureFaceScalar(tTags.fi_heat);

        size_t expected = 0;
        for (const auto& frac : fracs) {
            if (frac.elements.size() > 1) expected += (frac.elements.size() - 1);
        }
        if (!pFlow || !pHeat) {
            if (expected > 0) {
                throw std::runtime_error("[Build 2D] Missing FI transmissibility fields.");
            }
            return;
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
        const auto& pairsMap = meshMgr.getNNCTopologyMap();
        const PhysicalProperties_string_op::TransmissibilityFields tTags;
        auto pFlow = fieldMgr.getNNCScalar(tTags.nnc_flow);
        auto pHeat = fieldMgr.getNNCScalar(tTags.nnc_heat);

        size_t totalNNCCount = 0;
        for (const auto& kv : pairsMap) totalNNCCount += kv.second.size();
        if (!pFlow || !pHeat) {
            if (totalNNCCount > 0) {
                throw std::runtime_error("[Build 2D] Missing NNC transmissibility fields.");
            }
            return;
        }

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
        const auto& ffPairs = fieldMgr.ff_topology;
        const PhysicalProperties_string_op::TransmissibilityFields tTags;
        auto pFlow = fieldMgr.getFFScalar(tTags.ff_flow);
        auto pHeat = fieldMgr.getFFScalar(tTags.ff_heat);
        if (!pFlow || !pHeat) {
            if (!ffPairs.empty()) {
                throw std::runtime_error("[Build 2D] Missing FF transmissibility fields.");
            }
            return;
        }
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

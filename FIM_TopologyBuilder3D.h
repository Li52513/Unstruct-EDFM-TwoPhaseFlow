/**
 * @file FIM_TopologyBuilder3D.h
 * @brief 3D topology to algebraic connections.
 */

#pragma once

#include "FIM_ConnectionManager.h"
#include "3D_MeshManager.h"
#include "3D_FieldManager.h"

class FIM_TopologyBuilder3D {
public:
    static void LoadAllConnections(FIM_ConnectionManager& connMgr, const MeshManager_3D& meshMgr, const FieldManager_3D& fieldMgr) {
        size_t expected = 0;

        const auto& faces = meshMgr.mesh().getFaces();
        for (const auto& face : faces) {
            if (!face.isBoundary()) ++expected;
        }

        const auto& edges = meshMgr.fracture_network().getGlobalEdges();
        for (const auto& edge : edges) {
            if (edge.ownerCell_solverIndex >= 0 && edge.neighborCell_solverIndex >= 0) ++expected;
        }

        expected += meshMgr.getInteractionPairs().size();
        expected += fieldMgr.ff_topology.size();
        connMgr.ReserveRawConnections(expected);

        _loadMatrix(connMgr, meshMgr, fieldMgr);
        _loadFI(connMgr, meshMgr, fieldMgr);
        _loadNNC(connMgr, meshMgr, fieldMgr);
        _loadFF(connMgr, meshMgr, fieldMgr);
    }

private:
    static void _loadMatrix(FIM_ConnectionManager& connMgr, const MeshManager_3D& meshMgr, const FieldManager_3D& fieldMgr) {
        auto pFlow = fieldMgr.getMatrixFaceScalar("T_Matrix_Flow");
        auto pHeat = fieldMgr.getMatrixFaceScalar("T_Matrix_Heat");
        if (!pFlow || !pHeat) return;

        const auto& faces = meshMgr.mesh().getFaces();
        const auto& cells = meshMgr.mesh().getCells();
        const auto& cellId2Idx = meshMgr.mesh().getCellId2Index();

        if (pFlow->data.size() != faces.size() || pHeat->data.size() != faces.size()) {
            throw std::runtime_error("[Build 3D] Matrix face field size mismatch with mesh faces.");
        }

        for (size_t i = 0; i < faces.size(); ++i) {
            if (faces[i].isBoundary()) continue;
            int nodeI = static_cast<int>(cellId2Idx.at(faces[i].ownerCell));
            int nodeJ = static_cast<int>(cellId2Idx.at(faces[i].neighborCell));

            double dist = (faces[i].midpoint - cells[nodeI].center).Mag() +
                (cells[nodeJ].center - faces[i].midpoint).Mag();

            connMgr.PushConnection(
                nodeI, nodeJ,
                pFlow->data[i], pHeat->data[i],
                faces[i].vectorE.Mag(), std::max(dist, 1e-12),
                ConnectionType::Matrix_Matrix,
                faces[i].vectorT
            );
        }
    }

    static void _loadFI(FIM_ConnectionManager& connMgr, const MeshManager_3D& meshMgr, const FieldManager_3D& fieldMgr) {
        auto pFlow = fieldMgr.getFractureEdgeScalar("T_FI_Flow");
        auto pHeat = fieldMgr.getFractureEdgeScalar("T_FI_Heat");
        if (!pFlow || !pHeat) return;

        const auto& edges = meshMgr.fracture_network().getGlobalEdges();
        if (edges.size() != pFlow->data.size() || edges.size() != pHeat->data.size())
            throw std::runtime_error("[Build 3D] FI Transmissibility array size mismatches Geometry Edges.");

        for (size_t i = 0; i < edges.size(); ++i) {
            int nodeI = edges[i].ownerCell_solverIndex;
            int nodeJ = edges[i].neighborCell_solverIndex;
            if (nodeI < 0 || nodeJ < 0) continue;
            connMgr.PushConnection(nodeI, nodeJ, pFlow->data[i], pHeat->data[i], edges[i].length, edges[i].length, ConnectionType::Fracture_Internal);
        }
    }

    static void _loadNNC(FIM_ConnectionManager& connMgr, const MeshManager_3D& meshMgr, const FieldManager_3D& fieldMgr) {
        auto pFlow = fieldMgr.getNNCScalar("T_NNC_Flow");
        auto pHeat = fieldMgr.getNNCScalar("T_NNC_Heat");
        if (!pFlow || !pHeat) return;

        const auto& pairs = meshMgr.getInteractionPairs();
        if (pairs.size() != pFlow->data.size() || pairs.size() != pHeat->data.size()) {
            throw std::runtime_error("[Build 3D] NNC Geometry pairs size mismatches Transmissibility array lengths!");
        }

        for (size_t i = 0; i < pairs.size(); ++i) {
            connMgr.PushConnection(pairs[i].matrixSolverIndex, pairs[i].fracCellSolverIndex,
                pFlow->data[i], pHeat->data[i], pairs[i].intersectionArea,
                pairs[i].distMatrixToFracPlane, ConnectionType::Matrix_Fracture);
        }
    }

    static void _loadFF(FIM_ConnectionManager& connMgr, const MeshManager_3D& meshMgr, const FieldManager_3D& fieldMgr) {
        auto pFlow = fieldMgr.getFFScalar("T_FF_Flow");
        auto pHeat = fieldMgr.getFFScalar("T_FF_Heat");
        if (!pFlow || !pHeat) return;

        const auto& ffPairs = fieldMgr.ff_topology;
        if (ffPairs.size() != pFlow->data.size() || ffPairs.size() != pHeat->data.size()) {
            throw std::runtime_error("[Build 3D] FF Topology pairs mismatch with array lengths!");
        }

        const auto& frNet = meshMgr.fracture_network();
        const auto& elems = frNet.getOrderedFractureElements();
        const int offset = frNet.getSolverIndexOffset();

        for (size_t i = 0; i < ffPairs.size(); ++i) {
            const int sI = ffPairs[i].first;
            const int sJ = ffPairs[i].second;

            const int lI = sI - offset;
            const int lJ = sJ - offset;
            if (lI < 0 || lJ < 0 || lI >= static_cast<int>(elems.size()) || lJ >= static_cast<int>(elems.size()) ||
                elems[lI] == nullptr || elems[lJ] == nullptr) {
                throw std::runtime_error("[Build 3D] FF solverIndex out of ordered fracture cache range.");
            }

            const auto* eI = elems[lI];
            const auto* eJ = elems[lJ];

            double auxArea = 0.5 * (std::max(eI->area, 0.0) + std::max(eJ->area, 0.0));
            if (auxArea <= 1e-14) auxArea = 1.0;

            double auxDist = (eI->centroid - eJ->centroid).Mag();
            if (auxDist <= 1e-14) auxDist = 1e-12;

            connMgr.PushConnection(
                sI, sJ,
                pFlow->data[i], pHeat->data[i],
                auxArea, auxDist,
                ConnectionType::Fracture_Fracture
            );
        }
    }
};

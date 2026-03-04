/**
 * @file BoundaryAssembler.cpp
 * @brief Boundary and leakoff assembly facade implementation
 */
#include "BoundaryAssembler.h"
#include "FVM_Ops_AD.h"
#include "ADVar.hpp"
#include <cmath>
#include <iostream>

namespace {
constexpr double kBoundaryCoeffEps = 1.0e-14;
}

BoundaryAssemblyStats BoundaryAssembler::Assemble_2D(
    MeshManager& mgr,
    const BoundarySetting::BoundaryConditionManager& bcMgr,
    int dofOffset,
    FieldManager_2D& fm,
    const std::string& fieldName,
    std::vector<double>& residual,
    std::vector<double>& jacobianDiag)
{
    BoundaryAssemblyStats stats;
    if (dofOffset < 0 || dofOffset >= 3) {
        std::cerr << "[BoundaryAssembler::Assemble_2D] Invalid dofOffset = " << dofOffset << " (expected 0..2)." << std::endl;
        return stats;
    }
    if (residual.size() != jacobianDiag.size()) {
        std::cerr << "[BoundaryAssembler::Assemble_2D] residual/jacobianDiag size mismatch: "
                  << residual.size() << " vs " << jacobianDiag.size() << std::endl;
        return stats;
    }

    auto pField = fm.getMatrixScalar(fieldName);
    if (!pField) return stats;

    const auto& mesh = mgr.mesh();
    for (size_t i = 0; i < mesh.getCells().size(); ++i) {
        const auto& cell = mesh.getCells()[i];

        int solverIdx = static_cast<int>(i);
        int eqIdx = mgr.getEquationIndex(solverIdx, dofOffset);
        if (eqIdx < 0 || eqIdx >= static_cast<int>(residual.size()) || eqIdx >= static_cast<int>(jacobianDiag.size())) {
            continue;
        }

        ADVar<3> var_cell;
        var_cell.val = (*pField)[i];
        for (int k = 0; k < 3; ++k) var_cell.grad(k) = 0.0;
        var_cell.grad(dofOffset) = 1.0;

        for (int globalFaceID : cell.CellFaceIDs) {
            int faceIdx = mesh.getFaceIndex(globalFaceID);
            if (faceIdx < 0) continue;
            const auto& face = mesh.getFaces()[faceIdx];

            if (face.isBoundary() && bcMgr.HasBC(face.physicalGroupId)) {
                auto bc = bcMgr.GetBCCoefficients(face.physicalGroupId, face.midpoint);
                ADVar<3> q_flux;
                q_flux.val = 0.0;
                for (int k = 0; k < 3; ++k) q_flux.grad(k) = 0.0;

                if (bc.type == BoundarySetting::BoundaryType::Dirichlet) {
                    if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                        std::cerr << "[BoundaryAssembler::Assemble_2D] Skip Dirichlet face due to near-zero BC coefficient a on tag "
                                  << face.physicalGroupId << std::endl;
                        continue;
                    }
                    const double dist = (cell.center - face.midpoint).Mag();
                    const double T_geom = (dist > 1e-12) ? (face.length / dist) : 0.0;
                    q_flux = FVM_Ops::Op_Boundary_Dirichlet_AD<3, ADVar<3>>(T_geom, var_cell, bc.c / bc.a);
                }
                else if (bc.type == BoundarySetting::BoundaryType::Neumann) {
                    q_flux = FVM_Ops::Op_Boundary_Neumann_AD<3, ADVar<3>>(face.length, bc.c);
                }
                else if (bc.type == BoundarySetting::BoundaryType::Robin) {
                    if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                        std::cerr << "[BoundaryAssembler::Assemble_2D] Skip Robin face due to near-zero BC coefficient a on tag "
                                  << face.physicalGroupId << std::endl;
                        continue;
                    }
                    const double C_L = -bc.a;
                    const double far_val = bc.c / bc.a;

                    // Kernel gives leakoff per unit area; multiply by geometric boundary measure.
                    q_flux = FVM_Ops::Op_Leakoff_Source_AD<3, ADVar<3>>(true, C_L, var_cell, far_val);
                    const double A_face = face.length; // 2D: edge length with unit thickness
                    q_flux = A_face * q_flux;
                }

                residual[eqIdx] += q_flux.val;
                jacobianDiag[eqIdx] += q_flux.grad(dofOffset);

                stats.sumResidual += q_flux.val;
                stats.sumJacobianDiag += q_flux.grad(dofOffset);
                stats.matrixBCCount++;
            }
        }
    }

    for (const auto* elem : mgr.fracture_network().getOrderedFractureElements()) {
        (void)elem;
        stats.fractureBCCount += 0;
    }

    return stats;
}

BoundaryAssemblyStats BoundaryAssembler::Assemble_3D(
    MeshManager_3D& mgr,
    const BoundarySetting::BoundaryConditionManager& bcMgr,
    int dofOffset,
    FieldManager_3D& fm,
    const std::string& fieldName,
    std::vector<double>& residual,
    std::vector<double>& jacobianDiag)
{
    BoundaryAssemblyStats stats;
    if (dofOffset < 0 || dofOffset >= 3) {
        std::cerr << "[BoundaryAssembler::Assemble_3D] Invalid dofOffset = " << dofOffset << " (expected 0..2)." << std::endl;
        return stats;
    }
    if (residual.size() != jacobianDiag.size()) {
        std::cerr << "[BoundaryAssembler::Assemble_3D] residual/jacobianDiag size mismatch: "
                  << residual.size() << " vs " << jacobianDiag.size() << std::endl;
        return stats;
    }

    auto pField = fm.getMatrixScalar(fieldName);
    if (!pField) return stats;

    const auto& mesh = mgr.mesh();
    for (size_t i = 0; i < mesh.getCells().size(); ++i) {
        const auto& cell = mesh.getCells()[i];

        int solverIdx = static_cast<int>(i);
        int eqIdx = mgr.getEquationIndex(solverIdx, dofOffset);
        if (eqIdx < 0 || eqIdx >= static_cast<int>(residual.size()) || eqIdx >= static_cast<int>(jacobianDiag.size())) {
            continue;
        }

        ADVar<3> var_cell;
        var_cell.val = (*pField)[i];
        for (int k = 0; k < 3; ++k) var_cell.grad(k) = 0.0;
        var_cell.grad(dofOffset) = 1.0;

        for (int globalFaceID : cell.CellFaceIDs) {
            int faceIdx = mesh.getFaceIndex(globalFaceID);
            if (faceIdx < 0) continue;
            const auto& face = mesh.getFaces()[faceIdx];

            if (face.isBoundary() && bcMgr.HasBC(face.physicalGroupId)) {
                auto bc = bcMgr.GetBCCoefficients(face.physicalGroupId, face.midpoint);
                ADVar<3> q_flux;
                q_flux.val = 0.0;
                for (int k = 0; k < 3; ++k) q_flux.grad(k) = 0.0;

                if (bc.type == BoundarySetting::BoundaryType::Dirichlet) {
                    if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                        std::cerr << "[BoundaryAssembler::Assemble_3D] Skip Dirichlet face due to near-zero BC coefficient a on tag "
                                  << face.physicalGroupId << std::endl;
                        continue;
                    }
                    const double dist = (cell.center - face.midpoint).Mag();
                    const double T_geom = (dist > 1e-12) ? (face.length / dist) : 0.0;
                    q_flux = FVM_Ops::Op_Boundary_Dirichlet_AD<3, ADVar<3>>(T_geom, var_cell, bc.c / bc.a);
                }
                else if (bc.type == BoundarySetting::BoundaryType::Neumann) {
                    q_flux = FVM_Ops::Op_Boundary_Neumann_AD<3, ADVar<3>>(face.length, bc.c);
                }
                else if (bc.type == BoundarySetting::BoundaryType::Robin) {
                    if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                        std::cerr << "[BoundaryAssembler::Assemble_3D] Skip Robin face due to near-zero BC coefficient a on tag "
                                  << face.physicalGroupId << std::endl;
                        continue;
                    }
                    const double C_L = -bc.a;
                    const double far_val = bc.c / bc.a;

                    // Kernel gives leakoff per unit area; multiply by geometric boundary measure.
                    q_flux = FVM_Ops::Op_Leakoff_Source_AD<3, ADVar<3>>(true, C_L, var_cell, far_val);
                    const double A_face = face.length; // 3D: project stores boundary face area in length
                    q_flux = A_face * q_flux;
                }

                residual[eqIdx] += q_flux.val;
                jacobianDiag[eqIdx] += q_flux.grad(dofOffset);

                stats.sumResidual += q_flux.val;
                stats.sumJacobianDiag += q_flux.grad(dofOffset);
                stats.matrixBCCount++;
            }
        }
    }

    for (const auto* elem : mgr.fracture_network().getOrderedFractureElements()) {
        (void)elem;
        stats.fractureBCCount += 0;
    }

    return stats;
}

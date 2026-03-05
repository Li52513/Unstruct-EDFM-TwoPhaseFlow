/**
 * @file BoundaryAssembler.cpp
 * @brief Boundary and leakoff assembly facade implementation
 */
#include "BoundaryAssembler.h"
#include "FVM_Ops_AD.h"
#include "Well_WellControlTypes.h"
#include "SolverContrlStrName_op.h"
#include "ADVar.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

namespace {
    constexpr double kBoundaryCoeffEps = 1.0e-14;
    constexpr double kEps = 1.0e-12; // [Patch 2] 分母保护极小值

    struct WellFieldTags {
        std::string pressure;
        std::string k_xx;
        std::string k_yy;
        std::string k_zz;
        std::string rho_w;
        std::string rho_g;
        std::string lambda_w_mob;
        std::string lambda_g_mob;
        std::string h_w;
        std::string h_g;
    };

    inline WellFieldTags BuildWellFieldTags() {
        WellFieldTags tags;
        const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const PhysicalProperties_string_op::Rock rock;
        const PhysicalProperties_string_op::Water water;
        const PhysicalProperties_string_op::CO2 gas;
        tags.pressure = pCfg.pressure_field;
        tags.k_xx = rock.k_xx_tag;
        tags.k_yy = rock.k_yy_tag;
        tags.k_zz = rock.k_zz_tag;
        tags.rho_w = water.rho_tag;
        tags.rho_g = gas.rho_tag;
        tags.lambda_w_mob = water.lambda_w_tag;
        tags.lambda_g_mob = gas.lambda_g_tag;
        tags.h_w = water.h_tag;
        tags.h_g = gas.h_tag;
        return tags;
    }

    // [Patch 1] 新增 Helper 区
    inline bool isValidEqIdx(int eqIdx, const std::vector<double>& residual, const std::vector<double>& jacobianDiag) {
        return eqIdx >= 0 && eqIdx < static_cast<int>(residual.size()) && eqIdx < static_cast<int>(jacobianDiag.size());
    }

    inline bool isValidLocalIdx(int localIdx, const std::shared_ptr<volScalarField>& field) {
        return field && localIdx >= 0 && localIdx < static_cast<int>(field->data.size());
    }

    inline double safeGetFieldValue(const std::shared_ptr<volScalarField>& field, int localIdx, double default_val = 0.0) {
        if (isValidLocalIdx(localIdx, field)) return field->data[localIdx];
        return default_val;
    }

    inline int get3DMatrixOffset(const MeshManager_3D& mgr) {
        return mgr.fracture_network().getSolverIndexOffset();
    }

    inline void ResolveIndex(int input_id, WellTargetDomain domain, int n_matrix_dofs, int& solverIdx, int& localIdx) {
        if (domain == WellTargetDomain::Matrix) {
            localIdx = input_id;
            solverIdx = input_id;
        }
        else {
            if (input_id >= n_matrix_dofs) {
                solverIdx = input_id;
                localIdx = input_id - n_matrix_dofs;
            }
            else {
                localIdx = input_id;
                solverIdx = input_id + n_matrix_dofs;
            }
        }
    }
} // namespace end

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

BoundaryAssemblyStats BoundaryAssembler::Assemble_Wells_2D(
    MeshManager& mgr, FieldManager_2D& fm, const std::vector<WellScheduleStep>& active_steps,
    int dofOffset_P, int dofOffset_W, int dofOffset_G, int dofOffset_E,
    std::vector<double>& residual, std::vector<double>& jacobianDiag)
{
    BoundaryAssemblyStats stats;
    if (residual.size() != jacobianDiag.size() || dofOffset_P < 0 || dofOffset_P >= 3) {
        std::cerr << "[Assemble_Wells_2D] Error: Invalid Vector Size or P-Offset.\n";
        return stats;
    }
    if ((dofOffset_W >= 3) || (dofOffset_G >= 3) || (dofOffset_E >= 3)) {
        std::cerr << "[Assemble_Wells_2D] Error: One or more equation offsets exceed ADVar<3> bounds.\n";
        return stats;
    }

    const auto& mesh = mgr.mesh();
    int nMat = mgr.getMatrixDOFCount();
    const WellFieldTags tags = BuildWellFieldTags();

    for (const auto& step : active_steps) {
        int solverIdx = -1, localIdx = -1;
        ResolveIndex(step.completion_id, step.domain, nMat, solverIdx, localIdx);

        std::shared_ptr<volScalarField> pField_P, pField_Kxx, pField_Kyy;
        std::shared_ptr<volScalarField> pField_MobDenW, pField_MobDenG, pField_RhoW, pField_RhoG, pField_LamW, pField_LamG;
        std::shared_ptr<volScalarField> pField_Hw, pField_Hg;

        if (step.domain == WellTargetDomain::Matrix) {
            if (localIdx < 0 || localIdx >= static_cast<int>(mesh.getCells().size())) {
                std::cerr << "[WellSkip][2D][" << step.well_name << "] Invalid Matrix completion_id=" << localIdx << "\n";
                continue;
            }
            pField_P = fm.getMatrixScalar(tags.pressure);
            pField_Kxx = fm.getMatrixScalar(tags.k_xx);
            pField_Kyy = fm.getMatrixScalar(tags.k_yy);
            pField_MobDenW = fm.getMatrixScalar("mob_density_w");
            pField_MobDenG = fm.getMatrixScalar("mob_density_g");
            pField_RhoW = fm.getMatrixScalar(tags.rho_w);
            pField_RhoG = fm.getMatrixScalar(tags.rho_g);
            pField_LamW = fm.getMatrixScalar(tags.lambda_w_mob);
            pField_LamG = fm.getMatrixScalar(tags.lambda_g_mob);
            pField_Hw = fm.getMatrixScalar(tags.h_w);
            pField_Hg = fm.getMatrixScalar(tags.h_g);
        }
        else {
            pField_P = fm.getFractureScalar(tags.pressure);
            pField_MobDenW = fm.getFractureScalar("mob_density_w");
            pField_MobDenG = fm.getFractureScalar("mob_density_g");
            pField_RhoW = fm.getFractureScalar(tags.rho_w);
            pField_RhoG = fm.getFractureScalar(tags.rho_g);
            pField_LamW = fm.getFractureScalar(tags.lambda_w_mob);
            pField_LamG = fm.getFractureScalar(tags.lambda_g_mob);
            pField_Hw = fm.getFractureScalar(tags.h_w);
            pField_Hg = fm.getFractureScalar(tags.h_g);
        }

        if (!isValidLocalIdx(localIdx, pField_P)) {
            std::cerr << "[WellSkip][2D][" << step.well_name << "] Missing pressure field.\n";
            continue;
        }

        double mob_den_w = 0.0, mob_den_g = 0.0;
        if (pField_MobDenW) {
            mob_den_w = safeGetFieldValue(pField_MobDenW, localIdx, 0.0);
        }
        else if (pField_RhoW && pField_LamW) {
            mob_den_w = safeGetFieldValue(pField_RhoW, localIdx) * safeGetFieldValue(pField_LamW, localIdx);
        }
        else {
            std::cerr << "[WellSkip][2D][" << step.well_name << "] Missing water mobility data.\n";
            continue;
        }

        if (pField_MobDenG) {
            mob_den_g = safeGetFieldValue(pField_MobDenG, localIdx, 0.0);
        }
        else if (pField_RhoG && pField_LamG) {
            mob_den_g = safeGetFieldValue(pField_RhoG, localIdx) * safeGetFieldValue(pField_LamG, localIdx);
        }
        else {
            mob_den_g = 0.0;
        }

        double WI = step.wi_override;
        if (WI < 0.0) {
            if (step.domain == WellTargetDomain::Fracture) {
                std::cerr << "[WellSkip][2D][" << step.well_name << "] Fracture requires explicit wi_override.\n";
                continue;
            }
            double V = mesh.getCells()[localIdx].volume;
            double kxx = safeGetFieldValue(pField_Kxx, localIdx, 1e-13);
            double kyy = safeGetFieldValue(pField_Kyy, localIdx, 1e-13);
            double k_eff = std::sqrt(kxx * kyy);
            double r_eq = std::sqrt(V / 3.141592653589793);

            if (V <= kEps || step.rw <= 0.0 || r_eq <= step.rw || k_eff <= 0.0) {
                std::cerr << "[WellSkip][2D][" << step.well_name << "] Invalid geom for auto-WI (V=" << V << ", r_eq=" << r_eq << ").\n";
                continue;
            }
            double denom = std::log(r_eq / step.rw) + step.skin;
            if (denom <= kEps) {
                std::cerr << "[WellSkip][2D][" << step.well_name << "] WI denominator <= eps.\n";
                continue;
            }
            WI = (2.0 * 3.141592653589793 * k_eff * 1.0) / denom;
        }

        ADVar<3> P_cell; P_cell.val = pField_P->data[localIdx];
        for (int k = 0; k < 3; ++k) P_cell.grad(k) = 0.0;
        P_cell.grad(dofOffset_P) = 1.0;

        ADVar<3> q_mass_w = { 0 }, q_mass_g = { 0 };

        if (step.control_mode == WellControlMode::BHP) {
            if (step.component_mode == WellComponentMode::Water || step.component_mode == WellComponentMode::Total) {
                double frac = (step.component_mode == WellComponentMode::Total) ? step.frac_w : 1.0;
                q_mass_w = FVM_Ops::Op_Well_BHP_Source_AD<3, ADVar<3>>(WI * frac * mob_den_w, P_cell, step.target_value);
            }
            if (step.component_mode == WellComponentMode::Gas || step.component_mode == WellComponentMode::Total) {
                double frac = (step.component_mode == WellComponentMode::Total) ? step.frac_g : 1.0;
                q_mass_g = FVM_Ops::Op_Well_BHP_Source_AD<3, ADVar<3>>(WI * frac * mob_den_g, P_cell, step.target_value);
            }
        }
        else if (step.control_mode == WellControlMode::Rate) {
            if (step.component_mode == WellComponentMode::Water || step.component_mode == WellComponentMode::Total) {
                double frac = (step.component_mode == WellComponentMode::Total) ? step.frac_w : 1.0;
                q_mass_w = FVM_Ops::Op_Well_Rate_Source_AD<3, ADVar<3>>(step.target_value * frac);
            }
            if (step.component_mode == WellComponentMode::Gas || step.component_mode == WellComponentMode::Total) {
                double frac = (step.component_mode == WellComponentMode::Total) ? step.frac_g : 1.0;
                q_mass_g = FVM_Ops::Op_Well_Rate_Source_AD<3, ADVar<3>>(step.target_value * frac);
            }
        }

        if (dofOffset_W >= 0) {
            int eqIdx_W = mgr.getEquationIndex(solverIdx, dofOffset_W);
            if (isValidEqIdx(eqIdx_W, residual, jacobianDiag)) {
                residual[eqIdx_W] += q_mass_w.val;
                jacobianDiag[eqIdx_W] += q_mass_w.grad(dofOffset_P);
                stats.sumResidual += q_mass_w.val;
                stats.sumJacobianDiag += q_mass_w.grad(dofOffset_P); // [Patch 2] 完整统计
            }
        }
        if (dofOffset_G >= 0) {
            int eqIdx_G = mgr.getEquationIndex(solverIdx, dofOffset_G);
            if (isValidEqIdx(eqIdx_G, residual, jacobianDiag)) {
                residual[eqIdx_G] += q_mass_g.val;
                jacobianDiag[eqIdx_G] += q_mass_g.grad(dofOffset_P);
                stats.sumResidual += q_mass_g.val;
                stats.sumJacobianDiag += q_mass_g.grad(dofOffset_P);
            }
        }
        if (dofOffset_E >= 0 && isValidLocalIdx(localIdx, pField_Hw)) {
            double hw = pField_Hw->data[localIdx];
            double hg = safeGetFieldValue(pField_Hg, localIdx, 0.0);
            ADVar<3> q_energy = FVM_Ops::Op_Well_Energy_Source_AD<3, ADVar<3>>(q_mass_w, hw, q_mass_g, hg);
            int eqIdx_E = mgr.getEquationIndex(solverIdx, dofOffset_E);
            if (isValidEqIdx(eqIdx_E, residual, jacobianDiag)) {
                residual[eqIdx_E] += q_energy.val;
                jacobianDiag[eqIdx_E] += q_energy.grad(dofOffset_P);
                stats.sumResidual += q_energy.val;
                stats.sumJacobianDiag += q_energy.grad(dofOffset_P);
            }
        }

        if (step.domain == WellTargetDomain::Matrix) stats.matrixBCCount++;
        else stats.fractureBCCount++;
    }
    return stats;
}

// [Patch 3 & 9] 3D 井装配
BoundaryAssemblyStats BoundaryAssembler::Assemble_Wells_3D(
    MeshManager_3D& mgr, FieldManager_3D& fm, const std::vector<WellScheduleStep>& active_steps,
    int dofOffset_P, int dofOffset_W, int dofOffset_G, int dofOffset_E,
    std::vector<double>& residual, std::vector<double>& jacobianDiag)
{
    BoundaryAssemblyStats stats;
    if (residual.size() != jacobianDiag.size() || dofOffset_P < 0 || dofOffset_P >= 3) {
        std::cerr << "[Assemble_Wells_3D] Error: Invalid Vector Size or P-Offset.\n";
        return stats;
    }
    if ((dofOffset_W >= 3) || (dofOffset_G >= 3) || (dofOffset_E >= 3)) {
        std::cerr << "[Assemble_Wells_3D] Error: One or more equation offsets exceed ADVar<3> bounds.\n";
        return stats;
    }

    const auto& mesh = mgr.mesh();
    int nMat = get3DMatrixOffset(mgr); // [Patch 3] 修复 3D Matrix Offset
    const WellFieldTags tags = BuildWellFieldTags();

    for (const auto& step : active_steps) {
        int solverIdx = -1, localIdx = -1;
        ResolveIndex(step.completion_id, step.domain, nMat, solverIdx, localIdx);

        std::shared_ptr<volScalarField> pField_P, pField_Kxx, pField_Kyy, pField_Kzz;
        std::shared_ptr<volScalarField> pField_MobDenW, pField_MobDenG, pField_RhoW, pField_RhoG, pField_LamW, pField_LamG;
        std::shared_ptr<volScalarField> pField_Hw, pField_Hg;

        if (step.domain == WellTargetDomain::Matrix) {
            if (localIdx < 0 || localIdx >= static_cast<int>(mesh.getCells().size())) {
                std::cerr << "[WellSkip][3D][" << step.well_name << "] Invalid Matrix completion_id.\n";
                continue;
            }
            pField_P = fm.getMatrixScalar(tags.pressure);
            pField_Kxx = fm.getMatrixScalar(tags.k_xx);
            pField_Kyy = fm.getMatrixScalar(tags.k_yy);
            pField_Kzz = fm.getMatrixScalar(tags.k_zz);
            pField_MobDenW = fm.getMatrixScalar("mob_density_w");
            pField_MobDenG = fm.getMatrixScalar("mob_density_g");
            pField_RhoW = fm.getMatrixScalar(tags.rho_w);
            pField_RhoG = fm.getMatrixScalar(tags.rho_g);
            pField_LamW = fm.getMatrixScalar(tags.lambda_w_mob);
            pField_LamG = fm.getMatrixScalar(tags.lambda_g_mob);
            pField_Hw = fm.getMatrixScalar(tags.h_w);
            pField_Hg = fm.getMatrixScalar(tags.h_g);
        }
        else {
            pField_P = fm.getFractureScalar(tags.pressure);
            pField_MobDenW = fm.getFractureScalar("mob_density_w");
            pField_MobDenG = fm.getFractureScalar("mob_density_g");
            pField_RhoW = fm.getFractureScalar(tags.rho_w);
            pField_RhoG = fm.getFractureScalar(tags.rho_g);
            pField_LamW = fm.getFractureScalar(tags.lambda_w_mob);
            pField_LamG = fm.getFractureScalar(tags.lambda_g_mob);
            pField_Hw = fm.getFractureScalar(tags.h_w);
            pField_Hg = fm.getFractureScalar(tags.h_g);
        }

        if (!isValidLocalIdx(localIdx, pField_P)) {
            std::cerr << "[WellSkip][3D][" << step.well_name << "] Missing pressure field.\n";
            continue;
        }

        double mob_den_w = 0.0, mob_den_g = 0.0;
        if (pField_MobDenW) {
            mob_den_w = safeGetFieldValue(pField_MobDenW, localIdx, 0.0);
        }
        else if (pField_RhoW && pField_LamW) {
            mob_den_w = safeGetFieldValue(pField_RhoW, localIdx) * safeGetFieldValue(pField_LamW, localIdx);
        }
        else {
            std::cerr << "[WellSkip][3D][" << step.well_name << "] Missing mobility data.\n";
            continue;
        }

        if (pField_MobDenG) {
            mob_den_g = safeGetFieldValue(pField_MobDenG, localIdx, 0.0);
        }
        else if (pField_RhoG && pField_LamG) {
            mob_den_g = safeGetFieldValue(pField_RhoG, localIdx) * safeGetFieldValue(pField_LamG, localIdx);
        }

        double WI = step.wi_override;
        if (WI < 0.0) {
            if (step.domain == WellTargetDomain::Fracture) {
                std::cerr << "[WellSkip][3D][" << step.well_name << "] Fracture requires wi_override.\n";
                continue;
            }
            const auto& cell = mesh.getCells()[localIdx];
            double V = cell.volume;

            double L = step.L_override;
            WellAxis axis = step.well_axis;

            // [Patch 3] None 轴向保护
            if (L < 0.0 && axis == WellAxis::None) {
                std::cerr << "[Assemble_Wells_3D] Warning: well_axis=None with L_override<0. Defaulting to Z-axis.\n";
                axis = WellAxis::Z;
            }

            if (L < 0.0) {
                double min_ax = 1e20, max_ax = -1e20;
                for (int nodeID : cell.CellNodeIDs) {
                    const auto& pt = mesh.getNodesMap().at(nodeID).coord;
                    double val = (axis == WellAxis::X) ? pt.m_x : (axis == WellAxis::Y) ? pt.m_y : pt.m_z;
                    if (val < min_ax) min_ax = val;
                    if (val > max_ax) max_ax = val;
                }
                L = max_ax - min_ax;
            }

            if (V <= kEps || step.rw <= 0.0 || L <= kEps) {
                std::cerr << "[WellSkip][3D][" << step.well_name << "] Invalid geometry (V/L<=eps).\n";
                continue;
            }

            double r_eq = std::sqrt(V / (3.141592653589793 * L));
            double kxx = safeGetFieldValue(pField_Kxx, localIdx, 1e-13);
            double kyy = safeGetFieldValue(pField_Kyy, localIdx, 1e-13);
            double kzz = safeGetFieldValue(pField_Kzz, localIdx, 1e-13);

            double k_plane = 1e-13;
            if (axis == WellAxis::X) k_plane = std::sqrt(kyy * kzz);
            else if (axis == WellAxis::Y) k_plane = std::sqrt(kxx * kzz);
            else k_plane = std::sqrt(kxx * kyy);

            if (r_eq <= step.rw || k_plane <= 0.0) {
                std::cerr << "[WellSkip][3D][" << step.well_name << "] Auto-WI degenerated (r_eq <= rw or k<=0).\n";
                continue;
            }

            double denom = std::log(r_eq / step.rw) + step.skin;
            if (denom <= kEps) {
                std::cerr << "[WellSkip][3D][" << step.well_name << "] WI denominator <= eps.\n";
                continue;
            }
            WI = (2.0 * 3.141592653589793 * k_plane * L) / denom;
        }

        ADVar<3> P_cell; P_cell.val = pField_P->data[localIdx];
        for (int k = 0; k < 3; ++k) P_cell.grad(k) = 0.0;
        P_cell.grad(dofOffset_P) = 1.0;

        ADVar<3> q_mass_w = { 0 }, q_mass_g = { 0 };

        if (step.control_mode == WellControlMode::BHP) {
            if (step.component_mode == WellComponentMode::Water || step.component_mode == WellComponentMode::Total) {
                double frac = (step.component_mode == WellComponentMode::Total) ? step.frac_w : 1.0;
                q_mass_w = FVM_Ops::Op_Well_BHP_Source_AD<3, ADVar<3>>(WI * frac * mob_den_w, P_cell, step.target_value);
            }
            if (step.component_mode == WellComponentMode::Gas || step.component_mode == WellComponentMode::Total) {
                double frac = (step.component_mode == WellComponentMode::Total) ? step.frac_g : 1.0;
                q_mass_g = FVM_Ops::Op_Well_BHP_Source_AD<3, ADVar<3>>(WI * frac * mob_den_g, P_cell, step.target_value);
            }
        }
        else if (step.control_mode == WellControlMode::Rate) {
            if (step.component_mode == WellComponentMode::Water || step.component_mode == WellComponentMode::Total) {
                double frac = (step.component_mode == WellComponentMode::Total) ? step.frac_w : 1.0;
                q_mass_w = FVM_Ops::Op_Well_Rate_Source_AD<3, ADVar<3>>(step.target_value * frac);
            }
            if (step.component_mode == WellComponentMode::Gas || step.component_mode == WellComponentMode::Total) {
                double frac = (step.component_mode == WellComponentMode::Total) ? step.frac_g : 1.0;
                q_mass_g = FVM_Ops::Op_Well_Rate_Source_AD<3, ADVar<3>>(step.target_value * frac);
            }
        }

        if (dofOffset_W >= 0) {
            int eqIdx_W = mgr.getEquationIndex(solverIdx, dofOffset_W);
            if (isValidEqIdx(eqIdx_W, residual, jacobianDiag)) {
                residual[eqIdx_W] += q_mass_w.val;
                jacobianDiag[eqIdx_W] += q_mass_w.grad(dofOffset_P);
                stats.sumResidual += q_mass_w.val;
                stats.sumJacobianDiag += q_mass_w.grad(dofOffset_P);
            }
        }
        if (dofOffset_G >= 0) {
            int eqIdx_G = mgr.getEquationIndex(solverIdx, dofOffset_G);
            if (isValidEqIdx(eqIdx_G, residual, jacobianDiag)) {
                residual[eqIdx_G] += q_mass_g.val;
                jacobianDiag[eqIdx_G] += q_mass_g.grad(dofOffset_P);
                stats.sumResidual += q_mass_g.val;
                stats.sumJacobianDiag += q_mass_g.grad(dofOffset_P);
            }
        }
        if (dofOffset_E >= 0 && isValidLocalIdx(localIdx, pField_Hw)) {
            double hw = pField_Hw->data[localIdx];
            double hg = safeGetFieldValue(pField_Hg, localIdx, 0.0);
            ADVar<3> q_energy = FVM_Ops::Op_Well_Energy_Source_AD<3, ADVar<3>>(q_mass_w, hw, q_mass_g, hg);
            int eqIdx_E = mgr.getEquationIndex(solverIdx, dofOffset_E);
            if (isValidEqIdx(eqIdx_E, residual, jacobianDiag)) {
                residual[eqIdx_E] += q_energy.val;
                jacobianDiag[eqIdx_E] += q_energy.grad(dofOffset_P);
                stats.sumResidual += q_energy.val;
                stats.sumJacobianDiag += q_energy.grad(dofOffset_P);
            }
        }

        if (step.domain == WellTargetDomain::Matrix) stats.matrixBCCount++;
        else stats.fractureBCCount++;
    }
    return stats;
}

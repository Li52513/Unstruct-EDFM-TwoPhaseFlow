#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "ADVar.hpp"
#include "AD_FluidEvaluator.h"
#include "BoundaryAssembler.h"
#include "BoundaryConditionManager.h"
#include "MeshManager.h"
#include "SolverContrlStrName_op.h"
#include "2D_FieldManager.h"

namespace Test_Day3_BoundaryFixes {

inline double ComputeBoundaryMeasure2D(const Mesh& mesh) {
    double sum_area = 0.0;
    for (const auto& face : mesh.getFaces()) {
        if (!face.isBoundary()) continue;
        double area = std::max(face.vectorE.Mag(), 0.0);
        if (area <= 1.0e-20) area = std::max(face.length, 0.0);
        sum_area += area;
    }
    return sum_area;
}

inline void Test_N2_TemperatureJacobianSlot_2D() {
    const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
    const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();

    MeshManager mgr(10.0, 10.0, 0.0, 2, 2, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OrthogonalCorrection);
    mgr.BuildGlobalSystemIndexing();
    mgr.setNumDOFs(2);

    const int totalEq = mgr.getTotalEquationDOFs();
    std::vector<double> res(totalEq, 0.0);
    std::vector<double> diag(totalEq, 0.0);

    BoundarySetting::BoundaryConditionManager bcMgrT;
    for (const auto& face : mgr.mesh().getFaces()) {
        if (face.isBoundary()) {
            bcMgrT.SetDirichletBC(face.physicalGroupId, 280.0);
        }
    }

    FieldManager_2D fm;
    fm.InitSizes(mgr.mesh().getGridCount(), 0, 0, 0, mgr.mesh().getFaces().size(), 0);
    auto pField = fm.matrixFields.create<volScalarField>(pCfg.pressure_field, mgr.mesh().getGridCount());
    auto tField = fm.matrixFields.create<volScalarField>(tCfg.temperatue_field, mgr.mesh().getGridCount());
    for (int i = 0; i < pField->size; ++i) {
        (*pField)[i] = 2.0e5;
        (*tField)[i] = 320.0;
    }

    const auto stats = BoundaryAssembler::Assemble_2D(
        mgr, bcMgrT, 1, fm, tCfg.temperatue_field, res, diag);
    assert(stats.matrixBCCount > 0 && "N=2 temperature BC should touch boundary rows");

    bool has_nonzero_t_diag = false;
    for (int bi = 0; bi < static_cast<int>(mgr.mesh().getGridCount()); ++bi) {
        const int eqT = mgr.getEquationIndex(bi, 1);
        if (eqT < 0 || eqT >= totalEq) continue;
        if (std::abs(diag[eqT]) > 1.0e-20) {
            has_nonzero_t_diag = true;
            assert(diag[eqT] > 0.0 && "N=2 temperature Dirichlet diagonal should be positive");
        }
    }
    assert(has_nonzero_t_diag && "N=2 wrapper must map temperature diagonal to AD slot T");
    std::cout << "  [PASS] Day3 fix: N=2 temperature Jacobian slot mapping\n";
}

inline void Test_EnergyNeumannPhaseSplit_2D() {
    const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
    const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
    const auto sCfg = PhysicalProperties_string_op::SaturationEquation_String::FIM();

    MeshManager mgr(10.0, 10.0, 0.0, 1, 1, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OrthogonalCorrection);
    mgr.BuildGlobalSystemIndexing();
    mgr.setNumDOFs(3);

    const int totalEq = mgr.getTotalEquationDOFs();
    std::vector<double> res(totalEq, 0.0);
    std::vector<std::array<double, 3>> jac(totalEq, std::array<double, 3>{ 0.0, 0.0, 0.0 });

    constexpr double p_neumann_mass_flux = 1.0e-5;
    BoundarySetting::BoundaryConditionManager bcMgrP;
    BoundarySetting::BoundaryConditionManager bcMgrS;
    BoundarySetting::BoundaryConditionManager bcMgrT;
    for (const auto& face : mgr.mesh().getFaces()) {
        if (!face.isBoundary()) continue;
        bcMgrP.SetNeumannBC(face.physicalGroupId, p_neumann_mass_flux);
        bcMgrS.SetNeumannBC(face.physicalGroupId, 0.0);
        bcMgrT.SetNeumannBC(face.physicalGroupId, 0.0);
    }

    FieldManager_2D fm;
    fm.InitSizes(mgr.mesh().getGridCount(), 0, 0, 0, mgr.mesh().getFaces().size(), 0);
    auto pField = fm.matrixFields.create<volScalarField>(pCfg.pressure_field, mgr.mesh().getGridCount());
    auto tField = fm.matrixFields.create<volScalarField>(tCfg.temperatue_field, mgr.mesh().getGridCount());
    auto swField = fm.matrixFields.create<volScalarField>(sCfg.saturation, mgr.mesh().getGridCount());

    constexpr double p_val = 2.0e5;
    constexpr double t_val = 330.0;
    for (int i = 0; i < pField->size; ++i) {
        (*pField)[i] = p_val;
        (*tField)[i] = t_val;
        (*swField)[i] = 0.4;
    }

    const auto stats = BoundaryAssembler::Assemble_2D_FullJac(
        mgr, bcMgrT, 2, fm, tCfg.temperatue_field, res, jac,
        &bcMgrP, &bcMgrS, FluidPropertyEvalConfig(), CapRelPerm::VGParams(), CapRelPerm::RelPermParams());
    assert(stats.matrixBCCount > 0 && "Energy boundary assembly should visit matrix boundary faces");

    ADVar<1> p_ad(p_val);
    ADVar<1> t_ad(t_val);
    const auto props_w = AD_Fluid::Evaluator::evaluateWater<1>(p_ad, t_ad);
    const double expected = ComputeBoundaryMeasure2D(mgr.mesh()) * p_neumann_mass_flux * props_w.h.val;

    double actual = 0.0;
    for (int bi = 0; bi < static_cast<int>(mgr.mesh().getGridCount()); ++bi) {
        const int eqT = mgr.getEquationIndex(bi, 2);
        if (eqT >= 0 && eqT < totalEq) actual += res[eqT];
    }

    const double tol = 1.0e-10 * (std::abs(expected) + 1.0);
    assert(std::abs(actual - expected) < tol &&
        "Energy Neumann convection must match qE=(area*c_w)*h_w+(area*c_g)*h_g with c_g=0");
    std::cout << "  [PASS] Day3 fix: Neumann energy phase split consistency\n";
}

inline void Test_FractureSwClampGradient_2D() {
    const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
    const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
    const auto sCfg = PhysicalProperties_string_op::SaturationEquation_String::FIM();

    MeshManager mgr(10.0, 10.0, 0.0, 8, 8, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OrthogonalCorrection);
    mgr.addFracture(Vector(1.0e-9, 5.0, 0.0), Vector(10.0 - 1.0e-9, 5.0, 0.0));
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    mgr.BuildGlobalSystemIndexing();
    mgr.setNumDOFs(3);

    const int totalEq = mgr.getTotalEquationDOFs();
    std::vector<double> res(totalEq, 0.0);
    std::vector<std::array<double, 3>> jac(totalEq, std::array<double, 3>{ 0.0, 0.0, 0.0 });

    BoundarySetting::BoundaryConditionManager bcMgrSw;
    for (const auto& face : mgr.mesh().getFaces()) {
        if (face.isBoundary()) {
            bcMgrSw.SetDirichletBC(face.physicalGroupId, 0.2);
        }
    }

    const int nFrac = static_cast<int>(mgr.fracture_network().getOrderedFractureElements().size());
    FieldManager_2D fm;
    fm.InitSizes(mgr.mesh().getGridCount(), nFrac, 0, 0, mgr.mesh().getFaces().size(), 0);

    auto pMat = fm.matrixFields.create<volScalarField>(pCfg.pressure_field, mgr.mesh().getGridCount());
    auto tMat = fm.matrixFields.create<volScalarField>(tCfg.temperatue_field, mgr.mesh().getGridCount());
    auto sMat = fm.matrixFields.create<volScalarField>(sCfg.saturation, mgr.mesh().getGridCount());
    for (int i = 0; i < pMat->size; ++i) {
        (*pMat)[i] = 2.0e5;
        (*tMat)[i] = 320.0;
        (*sMat)[i] = 0.5;
    }

    auto pFrac = fm.fractureFields.create<volScalarField>(pCfg.pressure_field, nFrac);
    auto tFrac = fm.fractureFields.create<volScalarField>(tCfg.temperatue_field, nFrac);
    auto sFrac = fm.fractureFields.create<volScalarField>(sCfg.saturation, nFrac);
    for (int i = 0; i < nFrac; ++i) {
        (*pFrac)[i] = 2.0e5;
        (*tFrac)[i] = 320.0;
        (*sFrac)[i] = 1.0 + 1.0e-4;
    }

    const auto stats = BoundaryAssembler::Assemble_2D_FullJac(
        mgr, bcMgrSw, 1, fm, sCfg.saturation, res, jac,
        nullptr, nullptr, FluidPropertyEvalConfig(), CapRelPerm::VGParams(), CapRelPerm::RelPermParams());
    assert(stats.fractureBCCount > 0 && "Fracture boundary saturation test must touch fracture BC rows");

    bool found_nonzero_drdsw = false;
    for (const auto* elemPtr : mgr.fracture_network().getOrderedFractureElements()) {
        if (!elemPtr) continue;
        const int eqSw = mgr.getEquationIndex(elemPtr->solverIndex, 1);
        if (eqSw < 0 || eqSw >= totalEq) continue;
        if (std::abs(jac[eqSw][1]) > 1.0e-20) {
            found_nonzero_drdsw = true;
            break;
        }
    }

    assert(found_nonzero_drdsw &&
        "Fracture Sw clamp must preserve a non-zero dR/dSw when Sw is slightly out of [0,1]");
    std::cout << "  [PASS] Day3 fix: fracture Sw clamp gradient safety\n";
}

inline void RunPatchChecks() {
    Test_N2_TemperatureJacobianSlot_2D();
    Test_EnergyNeumannPhaseSplit_2D();
}

inline void RunLeakoffChecks() {
    Test_FractureSwClampGradient_2D();
}

} // namespace Test_Day3_BoundaryFixes

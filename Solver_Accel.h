#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <vector>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "PhysicalPropertiesManager.h"
#include "BCAdapter.h"
#include "TemperatureBCAdapter.h"
#include "WellConfig.h"
#include "FVM_WellDOF.h"

#include "DiffusionCentral.h"
#include "ConvectionUpwind_Flux.h"
#include "ConvectionUpwind.h"
#include "Timeterm_BDF.h"

#include "Solver_AssemblerCOO.h"
#include "Solver_TimeLoopUtils.h"
#include "Solver_PostChecks.h"
#include "LinearSolver_Eigen.h"

namespace Solver {
namespace Accel {

struct OuterIterRuntime
{
    struct FieldHistory
    {
        std::vector<double> prev;      // x^{k}
        std::vector<double> prevPrev;  // x^{k-1}
        double omega = 1.0;            // 上一次的 Aitken 系数/占位符
    };

    FieldHistory pHist;
    FieldHistory THist;

    std::vector<int8_t> faceFluxSign;  // 记录面质量通量符号
    std::vector<int>    flippedFaces;  // 最近一次检测出的翻转面

    std::vector<double> scratchP;      // 通用工作缓冲
    std::vector<double> scratchT;

    int pressurePatternStamp = -1;
    int temperaturePatternStamp = -1;
    int precondCountdownP = 0;
    int precondCountdownT = 0;

    double lastDpInner = std::numeric_limits<double>::infinity();
    double lastCFL_T = 0.0;
    double lastDt = -1.0;
};

inline void ensureOuterRuntimeCapacity(const Mesh& mesh, OuterIterRuntime& rt)
{
    const size_t nCells = mesh.getCells().size();
    const size_t nFaces = mesh.getFaces().size();

    auto ensureHistory = [nCells](OuterIterRuntime::FieldHistory& h)
    {
        if (h.prev.size() != nCells)      h.prev.assign(nCells, 0.0);
        if (h.prevPrev.size() != nCells)  h.prevPrev.assign(nCells, 0.0);
    };

    ensureHistory(rt.pHist);
    ensureHistory(rt.THist);

    if (rt.scratchP.size() != nCells) rt.scratchP.assign(nCells, 0.0);
    if (rt.scratchT.size() != nCells) rt.scratchT.assign(nCells, 0.0);

    if (rt.faceFluxSign.size() != nFaces) rt.faceFluxSign.assign(nFaces, 0);
    if (!rt.flippedFaces.empty()) rt.flippedFaces.clear();
}

//------------------------------------------------------------------------------
// Pressure assemble+solve: 统一处理有井/无井并返回 dp_inf
//------------------------------------------------------------------------------
inline bool assembleSolvePressureInner(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const PressureBCAdapter& Pbc,
    const Vector& g,
    double dt,
    const SolverControls& ctrl,
    const OperatorFieldNames& nmP,
    const std::vector<WellConfig>& wellsCfg_in,
    std::vector<WellConfig>& wellsCfg_work,
    std::vector<WellDOF>& wellsWork,
    double& dp_inf_out,
    SparseSystemCOO* lastSys = nullptr)
{
    Mesh& mesh = mgr.mesh();
    const bool hasWells = !wellsCfg_in.empty();

    if (!FVM::Diffusion::build_FaceCoeffs_Central(
        mgr, reg, freg,
        nmP.a_f_diff, nmP.s_f_diff,
        "p_g",
        { "kxx:kxx","kyy:kyy","kzz:kzz","/mu_g","rho:rho_g" },
        "rho_g",
        FVM::Diffusion::RhoFaceMethod::Linear,
        g, Pbc,
        /*massForm=*/false, /*alpha_anisotropy=*/0))
        return false;

    if (!FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Flow(
        mgr, reg, dt,
        "c_phi", "phi",
        "p_g_old", "rho_g",
        "p_g_prev", "rho_g",
        "Drho_Dp_g",
        nmP.a_time, nmP.b_time))
        return false;

    SparseSystemCOO sysp;
    if (!assemble_COO(mgr, reg, freg, "ddt+diffusion", nmP, &sysp)) return false;
    if (lastSys) *lastSys = sysp;

    int Np = 0; auto lid_p = buildUnknownMap(mesh, Np);
    auto pvec = gatherFieldToVec(reg, mesh, "p_g", lid_p, Np);

    if (hasWells) {
        wellsCfg_work = wellsCfg_in;
        build_masks_and_PI_for_all(mgr, reg, wellsCfg_work);

        wellsWork.clear();
        const int Nc = (int)mesh.getCells().size();
        const int Ntot = register_well_dofs_for_all(Nc, wellsCfg_work, wellsWork);
        extend_linear_system_size(sysp, Ntot);

        std::vector<int> lid_cell(Nc);
        for (int i = 0; i < Nc; ++i) lid_cell[i] = i;

        for (const auto& w : wellsWork) {
            add_peaceman_coupling_cell_rows(
                sysp, mesh, reg, w.PI_field, w.mask_field,
                lid_cell, w.lid, /*scaleByCellMeasure=*/false);
            add_well_row(
                sysp, mesh, reg, w.PI_field, w.mask_field,
                lid_cell, w.lid, w.mode, w.target,
                /*scaleByCellMeasure=*/false);
        }

        pvec.resize(Ntot, 0.0);
        for (const auto& w : wellsWork) {
            pvec[w.lid] = (w.mode == WellDOF::Mode::Pressure) ? w.target : pvec[0];
        }
    }

    auto optP = ctrl.lin_p;
    if (optP.tol <= 0.0) optP.tol = ctrl.tol_p_abs;

    double resP = 0.0; int itP = 0;
    if (!solveCOO_Eigen(sysp, pvec, optP, &itP, &resP)) return false;

    if (hasWells) {
        std::vector<double> p_cells(Np);
        std::copy(pvec.begin(), pvec.begin() + Np, p_cells.begin());
        scatterVecToField(reg, mesh, "p_g", lid_p, p_cells);
        writeback_pw_fields_for_all(reg, wellsWork, pvec);
    }
    else {
        scatterVecToField(reg, mesh, "p_g", lid_p, pvec);
    }

    underRelaxInPlace(reg, "p_g", "p_g_prev", ctrl.urf_p);
    dp_inf_out = maxAbsDiff(reg, "p_g", "p_g_prev");
    updatePrevIterates(reg, { {"p_g","p_g_prev"} });
    return true;
}

struct PressureStageResult
{
    double dpInf = std::numeric_limits<double>::infinity();
    int    sweepsPerformed = 0;
    bool   metInnerTol = false;
};

inline bool runPressureSweeps(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const PressureBCAdapter& Pbc,
    const Vector& g,
    double dt,
    const SolverControls& ctrl,
    const OperatorFieldNames& nmP,
    const std::vector<WellConfig>& wellsCfg_in,
    std::vector<WellConfig>& wellsCfg_work,
    std::vector<WellDOF>& wellsWork,
    PressureStageResult& outRes,
    SparseSystemCOO* lastSys = nullptr)
{
    outRes = {};
    const int maxSweeps = std::max(1, ctrl.NsweepP_max);
    for (int sweep = 0; sweep < maxSweeps; ++sweep) {
        double dpInner = 0.0;
        if (!assembleSolvePressureInner(
                mgr, reg, freg, Pbc, g, dt,
                ctrl, nmP, wellsCfg_in,
                wellsCfg_work, wellsWork,
                dpInner, lastSys))
            return false;

        outRes.dpInf = dpInner;
        outRes.sweepsPerformed = sweep + 1;
        if (dpInner < ctrl.tol_p_inner) {
            outRes.metInnerTol = true;
            break;
        }
    }
    return true;
}

inline bool rebuildDarcyFluxes(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const OperatorFieldNames& nmP,
    const PressureBCAdapter& Pbc)
{
    if (!FVM::Convection::buildFlux_Darcy_Mass(
            mgr, reg, freg,
            nmP.a_f_diff, nmP.s_f_diff,
            "p_g", "rho_g",
            "mf_g", "Qf_g", "ufn_g",
            &Pbc, true))
        return false;

    debugCheckMassFlux(mgr, freg, "mf_g", 1e-20);
    return true;
}

inline bool solveTemperatureStage(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const TemperatureBCAdapter& Tbc,
    const Vector& g,
    double dt,
    const SolverControls& ctrl,
    const OperatorFieldNames& nmT,
    double& dT_inf_out,
    SparseSystemCOO* lastSys = nullptr)
{
    if (!FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Temperature(
            mgr, reg, dt,
            "C_eff", "T_old",
            nmT.a_time, nmT.b_time))
        return false;

    if (!FVM::Diffusion::build_FaceCoeffs_Central(
            mgr, reg, freg,
            nmT.a_f_diff, nmT.s_f_diff,
            "T",
            { "iso:lambda_eff" },
            "",
            FVM::Diffusion::RhoFaceMethod::Linear,
            g, Tbc,
            /*massForm=*/false, /*alpha_anisotropy=*/0))
        return false;

    if (!FVM::Convection::build_FaceCoeffs_Upwind(
            mgr, reg, freg,
            "T", "mf_g",
            { "cp_g" },
            nmT, Tbc))
        return false;

    SparseSystemCOO sysT;
    if (!assemble_COO(mgr, reg, freg, "ddt+diffusion+convection", nmT, &sysT))
        return false;
    if (lastSys) *lastSys = sysT;

    Mesh& mesh = mgr.mesh();
    int Nt = 0; auto lid_t = buildUnknownMap(mesh, Nt);
    auto Tvec = gatherFieldToVec(reg, mesh, "T", lid_t, Nt);

    auto optT = ctrl.lin_T;
    if (optT.tol <= 0.0) optT.tol = ctrl.tol_T_abs;

    double resT = 0.0; int itT = 0;
    if (!solveCOO_Eigen(sysT, Tvec, optT, &itT, &resT)) return false;

    scatterVecToField(reg, mesh, "T", lid_t, Tvec);
    underRelaxInPlace(reg, "T", "T_prev", ctrl.urf_T);
    dT_inf_out = maxAbsDiff(reg, "T", "T_prev");
    updatePrevIterates(reg, { {"T","T_prev"} });
    return true;
}

inline bool outerIter_constProperties_singlePhase_CO2_T_H_withWell_accel(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    PhysicalPropertiesManager& ppm,
    const TemperatureBCAdapter& Tbc,
    const PressureBCAdapter& Pbc,
    const Vector& g,
    const std::vector<WellConfig>& wellsCfg_in,
    double dt,
    const SolverControls& ctrl,
    OuterIterRuntime& runtime)
{
    Mesh& mesh = mgr.mesh();
    ensureOuterRuntimeCapacity(mesh, runtime);

    double prev_dp_g = std::numeric_limits<double>::infinity();
    double prev_dT = std::numeric_limits<double>::infinity();

    std::vector<WellConfig> wellsCfg_work;
    std::vector<WellDOF> wellsWork;

    for (int it = 0; it < ctrl.maxOuter; ++it) {
        if (!startOuterIteration_scatter(reg, "p_g", "p_g_prev")) return false;
        if (!startOuterIteration_T(reg, "T", "T_prev")) return false;

        ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
        ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);
        ppm.ComputeEffectiveThermalProperties_constProperties_singlePhase_CO2_T_H(mgr, reg);

        const OperatorFieldNames nmP = makeNames("p_g");
        const OperatorFieldNames nmT = makeNames("T");

        SparseSystemCOO lastP;
        SparseSystemCOO lastT;

        PressureStageResult presRes;
        if (!runPressureSweeps(
                mgr, reg, freg, Pbc, g, dt, ctrl,
                nmP, wellsCfg_in, wellsCfg_work, wellsWork,
                presRes,
                ctrl.dumpMMOnLastIter ? &lastP : nullptr))
            return false;
        runtime.lastDpInner = presRes.dpInf;

        if (!rebuildDarcyFluxes(mgr, reg, freg, nmP, Pbc)) return false;

        double dT = 0.0;
        if (!solveTemperatureStage(
                mgr, reg, freg, Tbc, g, dt,
                ctrl, nmT, dT,
                ctrl.dumpMMOnLastIter ? &lastT : nullptr))
            return false;

        runtime.lastCFL_T = 0.0;

        std::cout << "[Outer " << it << "]  sweeps=" << presRes.sweepsPerformed
                  << "  |Δp|_inf=" << presRes.dpInf
                  << "  |ΔT|_inf=" << dT << "\n";

        if (it > 0) {
            double rT = dT / std::max(prev_dT, 1e-30);
            double uT = (rT < 0.7) ? +0.05 : (rT > 0.95 ? -0.05 : 0.0);
            auto& ctrl_mut = const_cast<SolverControls&>(ctrl);
            ctrl_mut.urf_T = std::min(0.7, std::max(0.15, ctrl_mut.urf_T + uT));

            double rp = presRes.dpInf / std::max(prev_dp_g, 1e-30);
            double up = (rp < 0.7) ? +0.05 : (rp > 0.95 ? -0.05 : 0.0);
            ctrl_mut.urf_p = std::min(0.7, std::max(0.15, ctrl_mut.urf_p + up));
        }

        prev_dp_g = presRes.dpInf;
        prev_dT = dT;

        auto maxAbsField = [&](const std::string& fld) -> double {
            double m = 0.0;
            auto f = reg.get<volScalarField>(fld);
            if (!f) return 1.0;
            const auto& cells = mesh.getCells();
            const auto& id2idx = mesh.getCellId2Index();
            for (const auto& c : cells) {
                if (c.id < 0) continue;
                size_t i = id2idx.at(c.id);
                m = std::max(m, std::abs((*f)[i]));
            }
            return std::max(1.0, m);
        };

        const double PScale = maxAbsField("p_g");
        const double TScale = maxAbsField("T");

        const bool convP = (presRes.dpInf < std::max(ctrl.tol_p_abs, ctrl.tol_p_rel * PScale));
        const bool convT = (dT < std::max(ctrl.tol_T_abs, ctrl.tol_T_rel * TScale));

        if (convP && convT) {
            std::cout << "Converged at outer iter " << it << "\n";

            if (ctrl.dumpMMOnLastIter) {
    #if __cplusplus >= 201703L
                try { std::filesystem::create_directories("mm"); }
                catch (...) {}
    #endif
                PostChecks::dumpCOO_to_matrix_market(
                    lastP,
                    "mm/A_P_CO2_pTH_accel.mtx",
                    "mm/b_P_CO2_pTH_accel.txt",
                    /*sym=*/false);
                PostChecks::dumpCOO_to_matrix_market(
                    lastT,
                    "mm/A_T_CO2_pTH_accel.mtx",
                    "mm/b_T_CO2_pTH_accel.txt",
                    /*sym=*/false);
            }
            return true;
        }

        if (it == ctrl.maxOuter - 1) {
            std::cout << "Reached maxOuter without meeting P/T tolerances.\n";
        }
    }

    return true;
}

} // namespace Accel
} // namespace Solver

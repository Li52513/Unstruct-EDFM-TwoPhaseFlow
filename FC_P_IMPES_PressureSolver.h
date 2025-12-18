#pragma once

#include <iostream>
#include <vector>

#include "FC_P_IMPES_PressureEqAssembler.h"
#include "Solver_TimeLoopUtils.h" // underRelaxInPlace / maxAbsDiff / updatePrevIterates
#include "PressureEqSolver.h" // reuse LinearSolverOptions, PressureSolveControls, gather/scatter, under-relax utilities

namespace FC_P_IMPES_I
{
    inline bool solver_FC_IMPES_I_PressureEq(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& Pbc,
        std::vector<WellDOF_TwoPhase>& wells,
        double dt,
        const IMPES_Iteration::PressureSolveControls& ctrl,
        double& dp_inf,
        double& lin_residual,
        int& lin_iterations,
        PhaseOperatorConfig fccfg = PhaseOperatorConfig{})
    {
        if (dt <= 0.0)
        {
            std::cerr << "[FC-IMPES-I][Pressure] invalid dt.\n";
            return false;
        }

        IMPES_Iteration::PressureAssemblyResult asmb;
        if (!assemblePressureTwoPhase_FC(mgr, reg, freg, Pbc, wells, dt, ctrl.assembly, asmb, fccfg))
        {
            std::cerr << "[FC-IMPES-I][Pressure] assemblePressureTwoPhase_FC failed.\n";
            return false;
        }

        const int nUnknowns = asmb.system.n;

        std::vector<double> pvec;
        IMPES_Iteration::gatherFieldToVector(
            reg,
            mgr.mesh(),
            ctrl.assembly.pressure_field,
            asmb.cell_lid,
            nUnknowns,
            pvec);

        double linRes = 0.0;
        int    linIters = 0;
        if (!solveCOO_Eigen(asmb.system, pvec, ctrl.linear, &linIters, &linRes))
        {
            std::cerr << "[FC-IMPES-I][Pressure] linear solver failed.\n";
            return false;
        }

        IMPES_Iteration::scatterVectorToField(
            reg,
            mgr.mesh(),
            ctrl.assembly.pressure_field,
            asmb.cell_lid,
            pvec);

        for (auto& w : wells)
        {
            if (w.lid >= 0 && w.lid < nUnknowns)
            {
                w.p_bh = pvec[w.lid];
            }
        }

        underRelaxInPlace(reg, ctrl.assembly.pressure_field, ctrl.assembly.pressure_prev_field, ctrl.under_relax);
        const double dpInf = maxAbsDiff(reg, ctrl.assembly.pressure_field, ctrl.assembly.pressure_prev_field);
        updatePrevIterates(reg, { { ctrl.assembly.pressure_field, ctrl.assembly.pressure_prev_field } });

        dp_inf = dpInf;
        lin_residual = linRes;
        lin_iterations = linIters;

        // Build phase mass fluxes directly (no split).
        fccfg.pressure_field = ctrl.assembly.pressure_field;
        fccfg.Pc_field = ctrl.assembly.Pc_field;
        if (!buildPhaseMassFluxes(mgr, reg, freg, Pbc, fccfg))
        {
            std::cerr << "[FC-IMPES-I][Flux] buildPhaseMassFluxes failed.\n";
            return false;
        }

        return true;
    }
}

#pragma once

#include <algorithm>
#include <cmath>
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

        // Save previous well BHPs for consistent under-relaxation (Rate wells only).
        std::vector<double> pbh_prev;
        pbh_prev.reserve(wells.size());
        for (const auto& w : wells) pbh_prev.push_back(w.p_bh);

        auto compute_row_residual = [&](int row) -> double
            {
                if (row < 0 || row >= asmb.system.n) return 0.0;
                double r = 0.0;
                if (row < static_cast<int>(asmb.system.b.size()))
                {
                    r -= asmb.system.b[static_cast<size_t>(row)];
                }
                for (const auto& t : asmb.system.A)
                {
                    if (t.r != row) continue;
                    if (t.c < 0 || static_cast<size_t>(t.c) >= pvec.size()) continue;
                    r += t.v * pvec[static_cast<size_t>(t.c)];
                }
                return r;
            };

        auto compute_max_rate_well_rel_residual = [&](double target_floor) -> double
            {
                double max_rel = 0.0;
                for (const auto& w : wells)
                {
                    if (w.mode != WellDOF_TwoPhase::Mode::Rate) continue;
                    const double denom = std::max(std::abs(w.target), target_floor);
                    const double r = compute_row_residual(w.lid); // units: kg/s (strict-rate equation)
                    const double rel = std::abs(r) / denom;
                    if (rel > max_rel) max_rel = rel;
                }
                return max_rel;
            };

        double linRes = 0.0;
        int    linIters = 0;
        LinearSolverOptions opt = ctrl.linear;
        if (!solveCOO_Eigen(asmb.system, pvec, opt, &linIters, &linRes))
        {
            std::cerr << "[FC-IMPES-I][Pressure] linear solver failed.\n";
            return false;
        }

        // If strict-rate wells are present, the global convergence criterion (relative residual) can
        // still leave a noticeable absolute residual on the well rate equation (kg/s).
        // Fall back to a direct solve when the rate-well equation is not accurate enough.
        {
            constexpr double rate_well_rel_tol = 1e-6;  // relative to target (or target_floor)
            constexpr double rate_well_target_floor = 1e-6; // kg/s, avoids blow-up when target ~ 0
            const double max_rel = compute_max_rate_well_rel_residual(rate_well_target_floor);
            if (max_rel > rate_well_rel_tol)
            {
                if (ctrl.verbose)
                {
                    std::cout << "[FC-IMPES-I][Pressure] rate-well residual max_rel=" << max_rel
                        << " (tol=" << rate_well_rel_tol << "), fallback to SparseLU.\n";
                }

                LinearSolverOptions opt2 = opt;
                opt2.type = LinearSolverOptions::Type::SparseLU;

                double linRes2 = 0.0;
                int    linIters2 = 0;
                if (solveCOO_Eigen(asmb.system, pvec, opt2, &linIters2, &linRes2))
                {
                    linRes = linRes2;
                    linIters = linIters2;

                    if (ctrl.verbose)
                    {
                        const double max_rel2 = compute_max_rate_well_rel_residual(rate_well_target_floor);
                        std::cout << "[FC-IMPES-I][Pressure] rate-well residual after SparseLU max_rel=" << max_rel2 << "\n";
                    }
                }
                else
                {
                    std::cerr << "[FC-IMPES-I][Pressure] SparseLU fallback failed (keeping iterative solution).\n";
                }
            }
        }

        IMPES_Iteration::scatterVectorToField(
            reg,
            mgr.mesh(),
            ctrl.assembly.pressure_field,
            asmb.cell_lid,
            pvec);

        for (size_t iw = 0; iw < wells.size(); ++iw)
        {
            auto& w = wells[iw];
            if (w.lid >= 0 && w.lid < nUnknowns)
            {
                w.p_bh = pvec[w.lid];
            }
            if (w.mode == WellDOF_TwoPhase::Mode::Pressure)
            {
                // Pressure-controlled wells: keep BHP exactly at target for downstream modules.
                w.p_bh = w.target;
            }
            // Apply the same under-relaxation to Rate wells' BHP to keep consistency with relaxed cell pressures.
            if (w.mode == WellDOF_TwoPhase::Mode::Rate)
            {
                const double urf = std::max(0.0, std::min(1.0, ctrl.under_relax));
                w.p_bh = pbh_prev[iw] + urf * (w.p_bh - pbh_prev[iw]);
            }
        }

        GeneralTools::underRelaxInPlace(reg, ctrl.assembly.pressure_field, ctrl.assembly.pressure_prev_field, ctrl.under_relax);
        const double dpInf = GeneralTools::maxAbsDiff(reg, ctrl.assembly.pressure_field, ctrl.assembly.pressure_prev_field);
        GeneralTools::updatePrevIterates(reg, { { ctrl.assembly.pressure_field, ctrl.assembly.pressure_prev_field } });

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

```cpp`r`n#pragma once

#include <vector>
#include <iostream>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "LinearSolver_Eigen.h"
#include "Solver_TimeLoopUtils.h"
#include "IMPES_PressureAssembler.h"
#include "IMPES_FluxSplitter.h"
#include "IMPES_SaturationTransport.h"
#include "IMPES_TemperatureAssembler.h"
#include "BCAdapter.h"
#include "TemperatureBCAdapter.h"
#include "FVM_WellDOF_TwoPhase.h"
#include "FaceSignMask.hpp"

namespace Solver
{
namespace IMPES
{
    using ::IMPES::PressureAssemblyConfig;
    using ::IMPES::PressureAssemblyResult;
    using ::IMPES::FluxSplitConfig;
    using ::IMPES::FluxSplitResult;
    using ::IMPES::SaturationTransportConfig;
    using ::IMPES::SaturationTransportReport;
    using ::IMPES::TemperatureSolveControls;
    using ::IMPES::TemperatureSolveReport;

    struct PressureSolveControls
    {
        PressureAssemblyConfig assembly;
        LinearSolverOptions linear;
        double under_relax = 1.0;
    };

    struct PressureSolveReport
    {
        double lin_residual = 0.0;
        int lin_iterations = 0;
        double dp_inf = 0.0;
        std::shared_ptr<faceScalarField> total_mass_flux;
        std::shared_ptr<faceScalarField> total_vol_flux;
        std::shared_ptr<faceScalarField> total_face_velocity;
    };

    inline void gatherFieldToVector(
        const FieldRegistry& reg,
        Mesh& mesh,
        const std::string& name,
        const std::vector<int>& lid_of_cell,
        int nUnknowns,
        std::vector<double>& out)
    {
        out.assign(nUnknowns, 0.0);
        auto fld = reg.get<volScalarField>(name);
        if (!fld) return;

        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const int r = lid_of_cell[i];
            if (r >= 0 && r < nUnknowns)
            {
                out[r] = (*fld)[i];
            }
        }
    }

    inline void scatterVectorToField(
        FieldRegistry& reg,
        Mesh& mesh,
        const std::string& name,
        const std::vector<int>& lid_of_cell,
        const std::vector<double>& vec)
    {
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        auto fld = reg.getOrCreate<volScalarField>(name, cells.size(), 0.0);

        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const int r = lid_of_cell[i];
            if (r >= 0 && static_cast<size_t>(r) < vec.size())
            {
                (*fld)[i] = vec[r];
            }
        }
    }

    inline bool assemblePressure_IMPES(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& Pbc,
        std::vector<WellDOF_TwoPhase>& wells,
        double dt,
        const PressureSolveControls& ctrl,
        PressureSolveReport* report = nullptr)
    {
        PressureAssemblyResult asmb;
        if (!::IMPES::assemblePressureSystem(
                mgr, reg, freg, Pbc,
                wells, dt,
                ctrl.assembly,
                asmb))
        {
            return false;
        }

        const int nUnknowns = asmb.system.n;
        std::vector<double> pvec;
        gatherFieldToVector(reg, mgr.mesh(), ctrl.assembly.pressure_field, asmb.cell_lid, nUnknowns, pvec);

        double linRes = 0.0;
        int linIters = 0;
        if (!solveCOO_Eigen(asmb.system, pvec, ctrl.linear, &linIters, &linRes))
        {
            return false;
        }

        scatterVectorToField(reg, mgr.mesh(), ctrl.assembly.pressure_field, asmb.cell_lid, pvec);

        if (ctrl.under_relax > 0.0 && ctrl.under_relax < 1.0 && !ctrl.assembly.pressure_prev_field.empty())
        {
            underRelaxInPlace(reg, ctrl.assembly.pressure_field, ctrl.assembly.pressure_prev_field, ctrl.under_relax);
        }

        double dpInf = 0.0;
        if (!ctrl.assembly.pressure_prev_field.empty())
        {
            dpInf = maxAbsDiff(reg, ctrl.assembly.pressure_field, ctrl.assembly.pressure_prev_field);
            updatePrevIterates(reg, { { ctrl.assembly.pressure_field, ctrl.assembly.pressure_prev_field } });
        }

        if (report)
        {
            report->lin_residual = linRes;
            report->lin_iterations = linIters;
            report->dp_inf = dpInf;
            report->total_mass_flux = asmb.total_mass_flux;
            report->total_vol_flux = asmb.total_vol_flux;
            report->total_face_velocity = asmb.total_face_velocity;
        }
        return true;
    }

    inline bool runTransient_IMPES(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& Pbc,
        const TemperatureBCAdapter& Tbc,
        std::vector<WellDOF_TwoPhase>& wells,
        int nSteps,
        double dt,
        const PressureSolveControls& pressureCtrl,
        const FluxSplitConfig& fluxCfg,
        const SaturationTransportConfig& satCfg,
        const TemperatureSolveControls& tempCtrl)
    {
        if (!(dt > 0.0))
        {
            std::cerr << "[IMPES] invalid dt.\n";
            return false;
        }

        FaceSignMask fluxMask;

        for (int step = 0; step < nSteps; ++step)
        {
            const int stepId = step + 1;

            if (!startTimeStep_scalar(mgr.mesh(), reg,
                    pressureCtrl.assembly.pressure_field,
                    pressureCtrl.assembly.pressure_old_field,
                    pressureCtrl.assembly.pressure_prev_field))
            {
                return false;
            }

            PressureSolveReport pRpt;
            if (!assemblePressure_IMPES(mgr, reg, freg, Pbc, wells, dt, pressureCtrl, &pRpt))
            {
                std::cerr << "[IMPES] pressure solve failed at step " << stepId << ".\n";
                return false;
            }

            FluxSplitResult fluxRpt;
            if (!::IMPES::splitTwoPhaseMassFlux(mgr, reg, freg, fluxCfg, &fluxMask, &fluxRpt))
            {
                std::cerr << "[IMPES] flux splitting failed at step " << stepId << ".\n";
                return false;
            }

            if (!startTimeStep_scalar(mgr.mesh(), reg,
                    satCfg.saturation,
                    satCfg.saturation_old,
                    satCfg.saturation_prev))
            {
                return false;
            }

            SaturationTransportReport satRpt;
            if (!::IMPES::advanceSaturationExplicit(mgr, reg, freg, wells, dt, satCfg, &satRpt))
            {
                std::cerr << "[IMPES] saturation update failed at step " << stepId << ".\n";
                return false;
            }

            if (!startTimeStep_scalar(mgr.mesh(), reg,
                    tempCtrl.assembly.temperature_field,
                    tempCtrl.assembly.temperature_old_field,
                    tempCtrl.assembly.temperature_prev_field))
            {
                return false;
            }

            TemperatureSolveReport tRpt;
            if (!::IMPES::assembleAndSolveTemperature_IMPES(
                    mgr, reg, freg, Tbc, wells, dt, tempCtrl, &tRpt))
            {
                std::cerr << "[IMPES] temperature solve failed at step " << stepId << ".\n";
                return false;
            }

            std::cout << "[IMPES] Step " << stepId
                      << " | Pressure iters=" << pRpt.lin_iterations
                      << ", |锟斤拷p|_inf=" << pRpt.dp_inf
                      << " | Sat CFL=" << satRpt.max_CFL
                      << " | Temp iters=" << tRpt.lin_iterations
                      << ", |锟斤拷T|_inf=" << tRpt.dT_inf
                      << std::endl;

            if (satRpt.max_CFL > satCfg.CFL_limit && satRpt.suggested_dt > 0.0)
            {
                std::cout << "[IMPES]   CFL warning: max=" << satRpt.max_CFL
                          << ", suggested 锟斤拷t=" << satRpt.suggested_dt << "\n";
            }
        }
        return true;
    }
} 
} 
`r`n```

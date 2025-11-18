#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <functional>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <limits>
#if __cplusplus >= 201703L
#include <filesystem>
#endif

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "LinearSolver_Eigen.h"
#include "Solver_TimeLoopUtils.h"
#include "IMPES_PressureEqAssemblerandSolver.h"
#include "IMPES_FluxSplitterandSolver.h"
#include "IMPES_SaturationTransportEqAssemblerandSolver.h"
#include "IMPES_TimeTermAssemblerandSolver.h"
#include  "IMPES_TemperatureAssemblerandSolver.h"
#include "IMPES_PostProcessIO.h"
#include "BCAdapter.h"
#include "TemperatureBCAdapter.h"
#include "FVM_WellDOF_TwoPhase.h"
#include "FaceSignMask.hpp"
#include "PostProcess_.h"
#include "Diff_TPFA_GradientsOperation.h"

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
            double tol_abs = 1e-6;
            double tol_rel = 1e-6;
            double urf_min = 0.15;
            double urf_max = 0.7;
            double urf_step = 0.05;
            int max_outer_iterations = 10;
            bool report_outer_iterations = false;
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
            const TemperatureSolveControls& tempCtrl,
            int writeEveryP = 0,
            int writeEveryT = 0,
            int writeEverySw = 0,
            const std::string& outPrefixP = "",
            const std::string& outPrefixT = "",
            const std::string& outPrefixSw = "",
            const std::string& timeSeriesFile = "",
            int snapshotEveryCsv = 0,
            const std::string& snapshotPrefix = "",
            std::vector<std::string> snapshotFields = {})
        {
            if (!(dt > 0.0))
            {
                std::cerr << "[IMPES] invalid dt.\n";
                return false;
            }
#if __cplusplus >= 201703L
            try
            {
                if (!outPrefixP.empty())
                {
                    auto dirP = std::filesystem::path(outPrefixP).parent_path();
                    if (!dirP.empty()) std::filesystem::create_directories(dirP);
                }
                if (!outPrefixT.empty())
                {
                    auto dirT = std::filesystem::path(outPrefixT).parent_path();
                    if (!dirT.empty()) std::filesystem::create_directories(dirT);
                }
                if (!outPrefixSw.empty())
                {
                    auto dirSw = std::filesystem::path(outPrefixSw).parent_path();
                    if (!dirSw.empty()) std::filesystem::create_directories(dirSw);
                }
                if (!timeSeriesFile.empty())
                {
                    auto dirTS = std::filesystem::path(timeSeriesFile).parent_path();
                    if (!dirTS.empty()) std::filesystem::create_directories(dirTS);
                }
                if (!snapshotPrefix.empty())
                {
                    auto dirSnap = std::filesystem::path(snapshotPrefix).parent_path();
                    if (!dirSnap.empty()) std::filesystem::create_directories(dirSnap);
                }
            }
            catch (...)
            {
                std::cerr << "[IMPES] cannot create directories for Tecplot output prefixes.\n";
                return false;
            }
#endif

            PressureSolveControls pressureState = pressureCtrl;
            TemperatureSolveControls tempState = tempCtrl;
            FaceSignMask fluxMask;
            FluxSplitConfig fluxState = fluxCfg;
            if (!fluxState.pressure_bc)
            {
                fluxState.pressure_bc = &Pbc;
            }
            if (snapshotFields.empty())
            {
                snapshotFields = {
                    pressureState.assembly.pressure_field,
                    tempState.assembly.temperature_field,
                    satCfg.saturation
                };
            }

            const auto maxAbsField = [&](const std::string& fld) -> double
            {
                auto field = reg.get<volScalarField>(fld);
                if (!field) return 1.0;
                const auto& cells = mgr.mesh().getCells();
                const auto& id2idx = mgr.mesh().getCellId2Index();
                double m = 0.0;
                for (const auto& c : cells)
                {
                    if (c.id < 0) continue;
                    const size_t i = id2idx.at(c.id);
                    m = std::max(m, std::abs((*field)[i]));
                }
                return std::max(1.0, m);
            };

            for (int step = 0; step < nSteps; ++step)
            {
                const int stepId = step + 1;

                if (!startTimeStep_scalar(mgr.mesh(), reg,
                    pressureState.assembly.pressure_field,
                    pressureState.assembly.pressure_old_field,
                    pressureState.assembly.pressure_prev_field))
                {
                    return false;
                }
                if (!startTimeStep_scalar(mgr.mesh(), reg,
                    satCfg.saturation,
                    satCfg.saturation_old,
                    satCfg.saturation_prev))
                {
                    return false;
                }
                if (!startTimeStep_scalar(mgr.mesh(), reg,
                    tempState.assembly.temperature_field,
                    tempState.assembly.temperature_old_field,
                    tempState.assembly.temperature_prev_field))
                {
                    return false;
                }

                PressureSolveReport pRpt;
                TemperatureSolveReport tRpt;
                SaturationTransportReport satRpt;

                const int maxOuter = std::max(1, pressureState.max_outer_iterations);
                double prev_dp = std::numeric_limits<double>::infinity();
                double prev_dT = std::numeric_limits<double>::infinity();
                bool convergedOuter = false;
                int lastOuterIters = 0;

                for (int outer = 0; outer < maxOuter; ++outer)
                {
                    lastOuterIters = outer + 1;
                    if (!pressureState.assembly.pressure_prev_field.empty())
                    {
                        if (!startOuterIteration_scatter(
                            reg,
                            pressureState.assembly.pressure_field,
                            pressureState.assembly.pressure_prev_field))
                        {
                            return false;
                        }
                    }
                    if (!tempState.assembly.temperature_prev_field.empty())
                    {
                        if (!startOuterIteration_scatter(
                            reg,
                            tempState.assembly.temperature_field,
                            tempState.assembly.temperature_prev_field))
                        {
                            return false;
                        }
                    }

                    if (!assemblePressure_IMPES(mgr, reg, freg, Pbc, wells, dt, pressureState, &pRpt))
                    {
                        std::cerr << "[IMPES] pressure solve failed at step " << stepId
                            << ", outer iter " << (outer + 1) << ".\n";
                        return false;
                    }

                    FluxSplitResult fluxRpt;
                    if (!::IMPES::splitTwoPhaseMassFlux(mgr, reg, freg, fluxState, &fluxMask, &fluxRpt))
                    {
                        std::cerr << "[IMPES] flux splitting failed at step " << stepId
                            << ", outer iter " << (outer + 1) << ".\n";
                        return false;
                    }

                    debugCheckMassFlux(mgr, freg, "mf_g", 1e-20);

                    if (!::IMPES::advanceSaturationExplicit(mgr, reg, freg, wells, dt, satCfg, &satRpt))
                    {
                        std::cerr << "[IMPES] saturation update failed at step " << stepId
                            << ", outer iter " << (outer + 1) << ".\n";
                        return false;
                    }

                    if (!::IMPES::assembleAndSolveTemperature_IMPES(
                        mgr, reg, freg, Tbc, wells, dt, tempState, &tRpt))
                    {
                        std::cerr << "[IMPES] temperature solve failed at step " << stepId
                            << ", outer iter " << (outer + 1) << ".\n";
                        return false;
                    }

                    if (pressureState.report_outer_iterations)
                    {
                        std::cout << "[IMPES] Step " << stepId
                            << " Outer " << (outer + 1)
                            << " |dp|_inf=" << pRpt.dp_inf
                            << " |dT|_inf=" << tRpt.dT_inf
                            << " |Sat CFL=" << satRpt.max_CFL << "\n";
                    }

                    if (outer > 0)
                    {
                        if (pressureState.urf_step > 0.0 && pressureState.urf_max > pressureState.urf_min)
                        {
                            const double ratioP = (std::isfinite(prev_dp) && prev_dp > 0.0)
                                ? pRpt.dp_inf / std::max(prev_dp, 1e-30)
                                : 1.0;
                            double deltaP = 0.0;
                            if (ratioP < 0.7) deltaP = pressureState.urf_step;
                            else if (ratioP > 0.95) deltaP = -pressureState.urf_step;
                            if (deltaP != 0.0)
                            {
                                double newUrfP = pressureState.under_relax + deltaP;
                                newUrfP = std::max(pressureState.urf_min, std::min(pressureState.urf_max, newUrfP));
                                pressureState.under_relax = newUrfP;
                            }
                        }

                        if (tempState.urf_step > 0.0 && tempState.urf_max > tempState.urf_min)
                        {
                            const double ratioT = (std::isfinite(prev_dT) && prev_dT > 0.0)
                                ? tRpt.dT_inf / std::max(prev_dT, 1e-30)
                                : 1.0;
                            double deltaT = 0.0;
                            if (ratioT < 0.7) deltaT = tempState.urf_step;
                            else if (ratioT > 0.95) deltaT = -tempState.urf_step;
                            if (deltaT != 0.0)
                            {
                                double newUrfT = tempState.under_relax + deltaT;
                                newUrfT = std::max(tempState.urf_min, std::min(tempState.urf_max, newUrfT));
                                tempState.under_relax = newUrfT;
                            }
                        }
                    }

                    prev_dp = pRpt.dp_inf;
                    prev_dT = tRpt.dT_inf;

                    const bool checkP = (pressureState.tol_abs > 0.0 || pressureState.tol_rel > 0.0);
                    const bool checkT = (tempState.tol_abs > 0.0 || tempState.tol_rel > 0.0);
                    const double pScale = checkP ? maxAbsField(pressureState.assembly.pressure_field) : 1.0;
                    const double tScale = checkT ? maxAbsField(tempState.assembly.temperature_field) : 1.0;
                    const double dpTol = checkP ? std::max(pressureState.tol_abs, pressureState.tol_rel * pScale) : 0.0;
                    const double dTTol = checkT ? std::max(tempState.tol_abs, tempState.tol_rel * tScale) : 0.0;
                    const bool convP = !checkP || (pRpt.dp_inf < dpTol);
                    const bool convT = !checkT || (tRpt.dT_inf < dTTol);

                    if (convP && convT)
                    {
                        convergedOuter = true;
						cout << "[IMPES] Step " << stepId
							<< " converged in " << (outer + 1) << " outer iterations.\n";
                        break;
                    }
                }

                if (!convergedOuter && maxOuter > 1)
                {
                    std::cout << "[IMPES] Step " << stepId
                        << " reached max outer iterations (" << maxOuter
                        << ") without satisfying pressure/temperature tolerances.\n";
                }
                if (lastOuterIters == 0) lastOuterIters = maxOuter;

                std::cout << "[IMPES] Step " << stepId
                    << " | Pressure iters=" << pRpt.lin_iterations
                    << ", |dp|_inf=" << pRpt.dp_inf
                    << " | Sat CFL=" << satRpt.max_CFL
                    << " | Temp iters=" << tRpt.lin_iterations
                    << ", |dT|_inf=" << tRpt.dT_inf
                    << std::endl;

                if (satRpt.max_CFL > satCfg.CFL_limit && satRpt.suggested_dt > 0.0)
                {
                    std::cout << "[IMPES]   CFL warning: max=" << satRpt.max_CFL
                        << ", suggested dt=" << satRpt.suggested_dt << "\n";
                }

                ::IMPES::Output::MassFluxSummary fluxSummary;
                ::IMPES::Output::computeMassFluxSummary(mgr, freg, fluxState.gas_mass_flux, fluxSummary);

                const double simTime = stepId * dt;
                if (!timeSeriesFile.empty())
                {
                    ::IMPES::Output::TimeSeriesRecord rec;
                    rec.step = stepId;
                    rec.simTime = simTime;
                    rec.outerIterations = lastOuterIters;
                    rec.dp_inf = pRpt.dp_inf;
                    rec.dT_inf = tRpt.dT_inf;
                    rec.satCFL = satRpt.max_CFL;
                    rec.satSuggestedDt = satRpt.suggested_dt;
                    rec.pressureLinIters = pRpt.lin_iterations;
                    rec.temperatureLinIters = tRpt.lin_iterations;
                    rec.flux = fluxSummary;
                    if (!::IMPES::Output::appendTimeSeriesRecord(timeSeriesFile, rec))
                    {
                        return false;
                    }
                }

                if (snapshotEveryCsv > 0 && !snapshotPrefix.empty() && (stepId % snapshotEveryCsv == 0))
                {
                    if (!::IMPES::Output::writeFieldSnapshotCSV(snapshotPrefix, stepId, simTime, mgr, reg, snapshotFields))
                    {
                        return false;
                    }
                }

                const bool dumpP = !outPrefixP.empty() && ((writeEveryP <= 0) || (stepId % writeEveryP == 0));
                const bool dumpT = !outPrefixT.empty() && ((writeEveryT <= 0) || (stepId % writeEveryT == 0));
                const bool dumpSw = !outPrefixSw.empty() && ((writeEverySw <= 0) || (stepId % writeEverySw == 0));

                if (dumpP)
                {
                    const std::vector<Vector> gradP =
                        computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, pressureState.assembly.pressure_field.c_str(), 0);

                    std::ostringstream fnP;
                    fnP << outPrefixP << "_step_" << std::setw(5) << std::setfill('0') << stepId << ".plt";

                    const std::string faceFieldName = pressureState.assembly.pressure_field + "_face_tmp";

                    const bool okPltP = outputTecplot_cellToFaceToNode_BC(
                        mgr, reg, freg,
                        nullptr, &Pbc,
                        pressureState.assembly.pressure_field,
                        faceFieldName,
                        &gradP,
                        fnP.str());

                    if (!okPltP)
                    {
                        std::cerr << "[IMPES] Tecplot export (pressure) failed at step " << stepId << ".\n";
                        return false;
                    }
                }

                if (dumpT)
                {
                    const std::vector<Vector> gradT =
                        computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, tempState.assembly.temperature_field.c_str(), 0);

                    std::ostringstream fnT;
                    fnT << outPrefixT << "_step_" << std::setw(5) << std::setfill('0') << stepId << ".plt";

                    const std::string faceFieldName = tempState.assembly.temperature_field + "_face_tmp";

                    const bool okPltT = outputTecplot_cellToFaceToNode_BC(
                        mgr, reg, freg,
                        &Tbc, nullptr,
                        tempState.assembly.temperature_field,
                        faceFieldName,
                        &gradT,
                        fnT.str());

                    if (!okPltT)
                    {
                        std::cerr << "[IMPES] Tecplot export (temperature) failed at step " << stepId << ".\n";
                        return false;
                    }
                }

                if (dumpSw)
                {
                    const std::vector<Vector> gradSw =
                        computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, satCfg.saturation.c_str(), 0);

                    std::ostringstream fnSw;
                    fnSw << outPrefixSw << "_step_" << std::setw(5) << std::setfill('0') << stepId << ".plt";

                    const std::string faceFieldName = satCfg.saturation + "_face_tmp";

                    const bool okPltSw = outputTecplot_cellToFaceToNode_BC(
                        mgr, reg, freg,
                        nullptr, nullptr,
                        satCfg.saturation,
                        faceFieldName,
                        &gradSw,
                        fnSw.str());

                    if (!okPltSw)
                    {
                        std::cerr << "[IMPES] Tecplot export (saturation) failed at step " << stepId << ".\n";
                        return false;
                    }
                }
        }
            return true;
        }

        inline bool runTransient_IMPES_P_s(
            MeshManager& mgr,
            FieldRegistry& reg,
            FaceFieldRegistry& freg,
            const PressureBCAdapter& Pbc,
            std::vector<WellDOF_TwoPhase>& wells,
            int nSteps,
            double dt,
            const PressureSolveControls& pressureCtrl,
            const FluxSplitConfig& fluxCfg,
            const SaturationTransportConfig& satCfg,
            std::function<bool()> updateMobility = nullptr,
            int writeEveryP = 0,
            int writeEverySw = 0,
            const std::string& outPrefixP = "",
            const std::string& outPrefixSw = "",
            const std::string& timeSeriesFile = "",
            int snapshotEveryCsv = 0,
            const std::string& snapshotPrefix = "",
            std::vector<std::string> snapshotFields = {})
        {
            if (!(dt > 0.0))
            {
                std::cerr << "[IMPES][P-s] invalid dt.\n";
                return false;
            }
#if __cplusplus >= 201703L
            try
            {
                if (!outPrefixP.empty())
                {
                    auto dirP = std::filesystem::path(outPrefixP).parent_path();
                    if (!dirP.empty()) std::filesystem::create_directories(dirP);
                }
                if (!outPrefixSw.empty())
                {
                    auto dirSw = std::filesystem::path(outPrefixSw).parent_path();
                    if (!dirSw.empty()) std::filesystem::create_directories(dirSw);
                }
                if (!timeSeriesFile.empty())
                {
                    auto dirTS = std::filesystem::path(timeSeriesFile).parent_path();
                    if (!dirTS.empty()) std::filesystem::create_directories(dirTS);
                }
                if (!snapshotPrefix.empty())
                {
                    auto dirSnap = std::filesystem::path(snapshotPrefix).parent_path();
                    if (!dirSnap.empty()) std::filesystem::create_directories(dirSnap);
                }
            }
            catch (...)
            {
                std::cerr << "[IMPES][P-s] cannot create directories for Tecplot/CSV output prefixes.\n";
                return false;
            }
#endif

            PressureSolveControls pressureState = pressureCtrl;
            FaceSignMask fluxMask;
            FluxSplitConfig fluxState = fluxCfg;
            if (!fluxState.pressure_bc)
            {
                fluxState.pressure_bc = &Pbc;
            }
            if (snapshotFields.empty())
            {
                snapshotFields = {
                    pressureState.assembly.pressure_field,
                    satCfg.saturation
                };
            }

            const auto maxAbsField = [&](const std::string& fld) -> double
            {
                auto field = reg.get<volScalarField>(fld);
                if (!field) return 1.0;
                const auto& cells = mgr.mesh().getCells();
                const auto& id2idx = mgr.mesh().getCellId2Index();
                double m = 0.0;
                for (const auto& c : cells)
                {
                    if (c.id < 0) continue;
                    const size_t i = id2idx.at(c.id);
                    m = std::max(m, std::abs((*field)[i]));
                }
                return std::max(1.0, m);
            };

            for (int step = 0; step < nSteps; ++step)
            {
                const int stepId = step + 1;

                if (!startTimeStep_scalar(mgr.mesh(), reg,
                    pressureState.assembly.pressure_field,
                    pressureState.assembly.pressure_old_field,
                    pressureState.assembly.pressure_prev_field))
                {
                    return false;
                }
                if (!startTimeStep_scalar(mgr.mesh(), reg,
                    satCfg.saturation,
                    satCfg.saturation_old,
                    satCfg.saturation_prev))
                {
                    return false;
                }

                PressureSolveReport pRpt;
                SaturationTransportReport satRpt;

                const int maxOuter = std::max(1, pressureState.max_outer_iterations);
                double prev_dp = std::numeric_limits<double>::infinity();
                bool convergedOuter = false;
                int lastOuterIters = 0;
               
                if (updateMobility)
                {
                    if (!updateMobility())
                    {
                        std::cerr << "[IMPES][P-s] mobility refresh failed at step " << stepId
                            << ".\n";
                        return false;
                    }
                }

                for (int outer = 0; outer < maxOuter; ++outer)
                {
                    lastOuterIters = outer + 1;
                    if (!pressureState.assembly.pressure_prev_field.empty())
                    {
                        if (!startOuterIteration_scatter(
                            reg,
                            pressureState.assembly.pressure_field,
                            pressureState.assembly.pressure_prev_field))
                        {
                            return false;
                        }
                    }

                    if (!assemblePressure_IMPES(mgr, reg, freg, Pbc, wells, dt, pressureState, &pRpt))
                    {
                        std::cerr << "[IMPES][P-s] pressure solve failed at step " << stepId
                            << ", outer iter " << (outer + 1) << ".\n";
                        return false;
                    }

                    FluxSplitResult fluxRpt;
                    if (!::IMPES::splitTwoPhaseMassFlux(mgr, reg, freg, fluxState, &fluxMask, &fluxRpt))
                    {
                        std::cerr << "[IMPES][P-s] flux splitting failed at step " << stepId
                            << ", outer iter " << (outer + 1) << ".\n";
                        return false;
                    }

                    debugCheckMassFlux(mgr, freg, fluxState.gas_mass_flux.empty() ? "mf_g" : fluxState.gas_mass_flux, 1e-20);

                    if (!::IMPES::advanceSaturationExplicit(mgr, reg, freg, wells, dt, satCfg, &satRpt))
                    {
                        std::cerr << "[IMPES][P-s] saturation update failed at step " << stepId
                            << ", outer iter " << (outer + 1) << ".\n";
                        return false;
                    }

                    

                    if (pressureState.report_outer_iterations)
                    {
                        std::cout << "[IMPES][P-s] Step " << stepId
                            << " Outer " << (outer + 1)
                            << " |dp|_inf=" << pRpt.dp_inf
                            << " |Sat CFL=" << satRpt.max_CFL << "\n";
                    }

                    if (outer > 0)
                    {
                        if (pressureState.urf_step > 0.0 && pressureState.urf_max > pressureState.urf_min)
                        {
                            const double ratioP = (std::isfinite(prev_dp) && prev_dp > 0.0)
                                ? pRpt.dp_inf / std::max(prev_dp, 1e-30)
                                : 1.0;
                            double deltaP = 0.0;
                            if (ratioP < 0.7) deltaP = pressureState.urf_step;
                            else if (ratioP > 0.95) deltaP = -pressureState.urf_step;
                            if (deltaP != 0.0)
                            {
                                double newUrfP = pressureState.under_relax + deltaP;
                                newUrfP = std::max(pressureState.urf_min, std::min(pressureState.urf_max, newUrfP));
                                pressureState.under_relax = newUrfP;
                            }
                        }
                    }

                    prev_dp = pRpt.dp_inf;

                    const bool checkP = (pressureState.tol_abs > 0.0 || pressureState.tol_rel > 0.0);
                    const double pScale = checkP ? maxAbsField(pressureState.assembly.pressure_field) : 1.0;
                    const double dpTol = checkP ? std::max(pressureState.tol_abs, pressureState.tol_rel * pScale) : 0.0;
                    const bool convP = !checkP || (pRpt.dp_inf < dpTol);

                    if (convP)
                    {
                        convergedOuter = true;
                        break;
                    }
                }

                if (!convergedOuter && maxOuter > 1)
                {
                    std::cout << "[IMPES][P-s] Step " << stepId
                        << " reached max outer iterations (" << maxOuter
                        << ") without satisfying pressure tolerance.\n";
                }
                if (lastOuterIters == 0) lastOuterIters = maxOuter;

                std::cout << "[IMPES][P-s] Step " << stepId
                    << " | Pressure iters=" << pRpt.lin_iterations
                    << ", |dp|_inf=" << pRpt.dp_inf
                    << " | Sat CFL=" << satRpt.max_CFL
                    << std::endl;

                ::IMPES::Output::MassFluxSummary fluxSummary;
                ::IMPES::Output::computeMassFluxSummary(mgr, freg, fluxState.gas_mass_flux.empty() ? "mf_g" : fluxState.gas_mass_flux, fluxSummary);

                const double simTime = stepId * dt;
                if (!timeSeriesFile.empty())
                {
                    ::IMPES::Output::TimeSeriesRecord rec;
                    rec.step = stepId;
                    rec.simTime = simTime;
                    rec.outerIterations = lastOuterIters;
                    rec.dp_inf = pRpt.dp_inf;
                    rec.dT_inf = 0.0;
                    rec.satCFL = satRpt.max_CFL;
                    rec.satSuggestedDt = satRpt.suggested_dt;
                    rec.pressureLinIters = pRpt.lin_iterations;
                    rec.temperatureLinIters = 0;
                    rec.flux = fluxSummary;
                    if (!::IMPES::Output::appendTimeSeriesRecord(timeSeriesFile, rec))
                    {
                        return false;
                    }
                }

                if (snapshotEveryCsv > 0 && !snapshotPrefix.empty() && (stepId % snapshotEveryCsv == 0))
                {
                    if (!::IMPES::Output::writeFieldSnapshotCSV(snapshotPrefix, stepId, simTime, mgr, reg, snapshotFields))
                    {
                        return false;
                    }
                }

                const bool dumpP = !outPrefixP.empty() && ((writeEveryP <= 0) || (stepId % writeEveryP == 0));
                const bool dumpSw = !outPrefixSw.empty() && ((writeEverySw <= 0) || (stepId % writeEverySw == 0));

                if (dumpP)
                {
                    const std::vector<Vector> gradP =
                        computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, pressureState.assembly.pressure_field.c_str(), 0);

                    std::ostringstream fnP;
                    fnP << outPrefixP << "_step_" << std::setw(5) << std::setfill('0') << stepId << ".plt";

                    const std::string faceFieldName = pressureState.assembly.pressure_field + "_face_tmp";

                    const bool okPltP = outputTecplot_cellToFaceToNode_BC(
                        mgr, reg, freg,
                        nullptr, &Pbc,
                        pressureState.assembly.pressure_field,
                        faceFieldName,
                        &gradP,
                        fnP.str());

                    if (!okPltP)
                    {
                        std::cerr << "[IMPES][P-s] Tecplot export (pressure) failed at step " << stepId << ".\n";
                        return false;
                    }
                }

                if (dumpSw)
                {
                    const std::vector<Vector> gradSw =
                        computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, satCfg.saturation.c_str(), 0);

                    std::ostringstream fnSw;
                    fnSw << outPrefixSw << "_step_" << std::setw(5) << std::setfill('0') << stepId << ".plt";

                    const std::string faceFieldName = satCfg.saturation + "_face_tmp";

                    const bool okPltSw = outputTecplot_cellToFaceToNode_BC(
                        mgr, reg, freg,
                        nullptr, nullptr,
                        satCfg.saturation,
                        faceFieldName,
                        &gradSw,
                        fnSw.str());

                    if (!okPltSw)
                    {
                        std::cerr << "[IMPES][P-s] Tecplot export (saturation) failed at step " << stepId << ".\n";
                        return false;
                    }
                }
            }
            return true;
        }
    } // namespace IMPES
} // namespace Solver

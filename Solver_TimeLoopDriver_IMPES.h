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
                snapshotFields = 
                {
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


namespace IMPES_revised
{
    ///压力求解器控制参数
    struct PressureSolveControls
    {
        PressureAssemblyConfig assembly;   ///< 里面指定 p 字段名、c_t 字段名等
        LinearSolverOptions    linear;     ///< 线性求解器设置（tol, maxIter, 等）

        double under_relax = 1.0;   ///< 欠松弛系数 (1.0 = 不欠松弛)
        double tol_abs = 1e-6;  ///< 外迭代绝对收敛阈值
        double tol_rel = 1e-6;  ///< 外迭代相对收敛阈值（对 dp_inf ）
        double urf_min = 0.15;  ///< URF 自适应下界（如需）
        double urf_max = 0.7;   ///< URF 自适应上界（如需）
        double urf_step = 0.05;  ///< 每步调整幅度（如需）
        int    max_outer_iterations = 10;  ///< 时间步内最大外迭代次数
        bool   report_outer_iterations = false;
    };
    
    /// 单次压力组装+线性求解的报告
    struct PressureSolveReport
    {
        double lin_residual = 0.0; ///< 线性求解器残差
        int    lin_iterations = 0;   ///< 线性求解器迭代步数
        double dp_inf = 0.0; ///< 本次 outer 中 p 的无穷范数变化

        std::shared_ptr<faceScalarField> total_mass_flux;
        std::shared_ptr<faceScalarField> total_vol_flux;
        std::shared_ptr<faceScalarField> total_face_velocity;
    };

    /// 把 cell 场收集到未知向量
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

    /// 把未知向量散到 cell 场
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

    /**
     * @brief  单次“组装两相流压力方程 + 线性求解器 + 欠松弛 + 井底压力回写”。
     *
     * 假设:
     * - 物性 (mobility, density, c_t) 已在时间步开始被冻结，不在此函数内更新；
     * - assemblePressureTwoPhase 已包含井 DOF，并构造好 sys.n >= max(cell, well lid)；
     */
    inline bool solver_IMPES_Pressure(
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
        if (!assemblePressureTwoPhase(mgr, reg, freg, Pbc, wells, dt, ctrl.assembly, asmb))
        {
            std::cerr << "[IMPES_revised][Pressure] assemblePressureTwoPhase failed.\n";
            return false;
        }

        const int nUnknowns = asmb.system.n;

        // 1. 用当前压力场作为初始猜测
        std::vector<double> pvec;
        gatherFieldToVector(reg, mgr.mesh(), ctrl.assembly.pressure_field,
            asmb.cell_lid, nUnknowns, pvec);

        double linRes = 0.0;
        int    linIters = 0;
        if (!solveCOO_Eigen(asmb.system, pvec, ctrl.linear, &linIters, &linRes))
        {
            std::cerr << "[IMPES_revised][Pressure] linear solver failed.\n";
            return false;
        }

        // 2. 把 cell 部分解写回压力场
        scatterVectorToField(reg, mgr.mesh(), ctrl.assembly.pressure_field,
            asmb.cell_lid, pvec);

        // 3. 更新井底压力: p_bh = 解向量中对应的 DOF
        for (auto& w : wells)
        {
            if (w.lid >= 0 && w.lid < nUnknowns)
            {
                w.p_bh = pvec[w.lid];
            }
        }

        // 4. 对 cell 压力做欠松弛（注意：目前没有对井 DOF 欠松弛，如果要严格 Rate，可设 under_relax=1）
        if (ctrl.under_relax > 0.0 && ctrl.under_relax < 1.0 &&
            !ctrl.assembly.pressure_prev_field.empty())
        {
            underRelaxInPlace(reg,
                ctrl.assembly.pressure_field,
                ctrl.assembly.pressure_prev_field,
                ctrl.under_relax);
        }

        // 5. 计算本次 outer 的 dp_inf，并更新 prev 场
        double dpInf = 0.0;
        if (!ctrl.assembly.pressure_prev_field.empty())
        {
            dpInf = maxAbsDiff(reg,
                ctrl.assembly.pressure_field,
                ctrl.assembly.pressure_prev_field);
            updatePrevIterates(reg,
                { { ctrl.assembly.pressure_field,
                    ctrl.assembly.pressure_prev_field } });
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



    /**
     * @brief  主时间步推进函数：两相流 IMPES (P-S)。
     *
     * 时间步循环结构：
     *  - 对每个时间步 n→n+1:
     *    1. startTimeStep_scalar 更新 p_old, p_prev, Sw_old, Sw_prev；
     *    2. 复制 rho_t → rho_t_old，为压力时间项提供 ρ_t^n；
     *    3. 调用 updateMobility 冻结物性（两相岩石 + CO2 + 水，含 rho_t, drho_t_dp 等）；
     *    4. 仅对压力做外迭代：调用 solver_IMPES_Pressure，直到 dp_inf 收敛或 outer 达到上限；
     *    5. 压力收敛后，拆总质量通量 → 相质量通量，并显式更新饱和度；
     *    6. 基于 CFL / 饱和度抖动调整下一时间步 dt，输出/快照。
     */
    inline bool runTransient_IMPES(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& Pbc,
        std::vector<WellDOF_TwoPhase>& wells,
        int    nSteps,
        double dt0,
        const PressureSolveControls& pressureCtrl,
        const FluxSplitConfig& fluxCfg,              ///< 目前不再在本函数内使用
        const SaturationTransportConfig& satCfg,
        std::function<bool()> updateMobility = nullptr,
        int writeEveryP = 0,
        int writeEverySw = 0,
        const std::string& outPrefixP = "",
        const std::string& outPrefixSw = "",
        const std::string& timeSeriesFile = "",
        int  snapshotEveryCsv = 0,
        const std::string& snapshotPrefix = "",
        std::vector<std::string> snapshotFields = {})
    {
        if (!(dt0 > 0.0))
        {
            std::cerr << "[IMPES_revised] invalid initial dt.\n";
            return false;
        }

        // 当前版本不再在这里直接使用 fluxCfg（flux splitting 交给 Sat 模块内部）。
        (void)fluxCfg;

#if __cplusplus >= 201703L
        // 创建输出目录（如果需要）
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

        Mesh& mesh = mgr.mesh();

        double time = 0.0;
        double dt_cur = dt0;

        if (snapshotFields.empty())
        {
            snapshotFields =
            {
                pressureCtrl.assembly.pressure_field,
                satCfg.saturation
            };
        }

        // ===== 预处理：确保 rho_total_old_field 存在，并在 t=0 做一次初始化 =====
        {
            const auto& assm = pressureCtrl.assembly;
            if (!assm.rho_total_field.empty() && !assm.rho_total_old_field.empty())
            {
                auto rho_t = reg.get<volScalarField>(assm.rho_total_field);
                if (rho_t)
                {
                    auto rho_t_old = reg.get<volScalarField>(assm.rho_total_old_field);
                    if (!rho_t_old)
                    {
                        // 创建 rho_t_old，并初始化为 rho_t（首时间层）
                        reg.getOrCreate<volScalarField>(
                            assm.rho_total_old_field,
                            rho_t->data.size(),
                            0.0);

                        if (!copyField(reg,
                            assm.rho_total_field,
                            assm.rho_total_old_field))
                        {
                            std::cerr << "[IMPES_revised] failed to init rho_total_old_field.\n";
                            return false;
                        }
                    }
                }
                else
                {
                    std::cerr << "[IMPES_revised] warning: rho_total_field '"
                        << assm.rho_total_field
                        << "' not found when initializing rho_total_old_field.\n";
                }
            }
        }

        // 小工具：量级估计，用于构造相对收敛判据
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

        // ==================== ?????? ====================
        const double dt_min = 1e-6;
        const double dt_max = 1e9;
        const double grow_max = 1.5;
        const double shrink_min = 0.05;
        const double retry_shrink = 0.5;
        const double safety = 0.9;

        const auto rollbackPrimaryFields = [&]() -> bool
            {
                bool ok = true;
                if (!pressureCtrl.assembly.pressure_field.empty() &&
                    !pressureCtrl.assembly.pressure_old_field.empty())
                {
                    ok = ok && copyField(reg,
                        pressureCtrl.assembly.pressure_old_field,
                        pressureCtrl.assembly.pressure_field);
                }
                if (!pressureCtrl.assembly.pressure_prev_field.empty() &&
                    !pressureCtrl.assembly.pressure_old_field.empty())
                {
                    ok = ok && copyField(reg,
                        pressureCtrl.assembly.pressure_old_field,
                        pressureCtrl.assembly.pressure_prev_field);
                }
                if (!satCfg.saturation.empty() &&
                    !satCfg.saturation_old.empty())
                {
                    ok = ok && copyField(reg,
                        satCfg.saturation_old,
                        satCfg.saturation);
                }
                if (!satCfg.saturation_prev.empty() &&
                    !satCfg.saturation_old.empty())
                {
                    ok = ok && copyField(reg,
                        satCfg.saturation_old,
                        satCfg.saturation_prev);
                }
                return ok;
            };

        const auto shrinkTimeStepForRetry = [&](const std::string& reason) -> bool
            {
                const double newDt = std::max(dt_min, dt_cur * retry_shrink);
                if (newDt >= dt_cur && dt_cur <= dt_min)
                {
                    std::cerr << "[IMPES_revised] " << reason
                        << ", but dt already reached minimum (" << dt_cur << ").\n";
                    return false;
                }
                dt_cur = newDt;
                std::cout << "  [TimeStep] " << reason
                    << ", retry with dt = " << dt_cur << "\n";
                return true;
            };

        int step = 0;
        while (step < nSteps)
        {
            bool stepAccepted = false;
            int attempt = 0;

            while (!stepAccepted)
            {
                ++attempt;
                const int stepId = step + 1;
                const double trialTime = time + dt_cur;

                std::cout << "\n[IMPES_revised] ===== Time step " << stepId
                    << " (t = " << trialTime << " s, dt = " << dt_cur << " s) =====\n";
                if (attempt > 1)
                {
                    std::cout << "  [TimeStep] retry attempt " << attempt << "\n";
                }

                // 1) ??????????? p_old/p_prev, Sw_old/Sw_prev
                if (!pressureCtrl.assembly.pressure_field.empty() &&
                    !pressureCtrl.assembly.pressure_old_field.empty() &&
                    !pressureCtrl.assembly.pressure_prev_field.empty())
                {
                    startTimeStep_scalar(
                        mesh,
                        reg,
                        pressureCtrl.assembly.pressure_field,
                        pressureCtrl.assembly.pressure_old_field,
                        pressureCtrl.assembly.pressure_prev_field);
                }

                if (!satCfg.saturation.empty() &&
                    !satCfg.saturation_old.empty() &&
                    !satCfg.saturation_prev.empty())
                {
                    startTimeStep_scalar(
                        mesh,
                        reg,
                        satCfg.saturation,
                        satCfg.saturation_old,
                        satCfg.saturation_prev);
                }

                // 2) ?????????????????????? (mobility, density, c_t, drho_dp, ...)
                if (updateMobility)
                {
                    if (!updateMobility())
                    {
                        std::cerr << "[IMPES_revised] updateMobility() failed at step "
                            << stepId << ".\n";
                        return false;
                    }
                }

                // 3) ??????????????
                IMPES_revised::PressureSolveControls ctrl = pressureCtrl;
                IMPES_revised::PressureSolveReport   pRpt;
                bool   converged = false;

                const std::string pFieldName = ctrl.assembly.pressure_field;
                const double PScale = maxAbsField(pFieldName);

                for (int outer = 0; outer < ctrl.max_outer_iterations; ++outer)
                {
                    const int outerId = outer + 1;

                    if (!solver_IMPES_Pressure(
                        mgr, reg, freg, Pbc, wells, dt_cur, ctrl, &pRpt))
                    {
                        std::cerr << "[IMPES_revised] pressure solver failed at step "
                            << stepId << ", outer iter " << outerId << ".\n";
                        return false;
                    }

                    std::cout << "  [Pressure] outer " << outerId
                        << ": dp_inf = " << pRpt.dp_inf
                        << ", lin res = " << pRpt.lin_residual
                        << ", iters = " << pRpt.lin_iterations << "\n";

                    const double tolEff = std::max(ctrl.tol_abs, ctrl.tol_rel * PScale);
                    const bool   convP = (pRpt.dp_inf < tolEff);

                    if (convP)
                    {
                        std::cout << "  [Pressure] converged at outer " << outerId
                            << ", dp_inf = " << pRpt.dp_inf
                            << ", PScale = " << PScale
                            << ", tolEff = " << tolEff << "\n";
                        converged = true;
                        break;
                    }
                }

                if (!converged)
                {
                    std::cerr << "[IMPES_revised] WARNING: pressure not fully converged at step "
                        << stepId << ", last dp_inf = " << pRpt.dp_inf << ".\n";

                    if (!rollbackPrimaryFields()) return false;
                    if (!shrinkTimeStepForRetry("pressure not converged")) return false;
                    continue;
                }

                std::cout << "  [IMPES_revised] Advancing saturation explicitly.\n";
                SaturationTransportReport satRpt;

                if (!IMPES_revised::advanceSaturationExplicit(
                    mgr, reg, freg, wells, dt_cur, satCfg, &satRpt))
                {
                    std::cerr << "[IMPES_revised][P-s] saturation update failed at step "
                        << stepId << ".\n";
                    return false;
                }

                std::cout << "  [Sat] max_CFL = " << satRpt.max_CFL
                    << ", max_dS = " << satRpt.max_dS
                    << ", suggested_dt = " << satRpt.suggested_dt << "\n";

                bool needRetry = false;
                double retryDt = dt_cur;
                double dt_next = dt_cur;

                if (satRpt.suggested_dt > 0.0 &&
                    (satCfg.track_cfl || satCfg.track_saturation_jitter))
                {
                    double target = safety * satRpt.suggested_dt;
                    target = std::max(dt_min, std::min(dt_max, target));

                    if (target < dt_cur)
                    {
                        needRetry = true;
                        retryDt = std::max(dt_min, std::max(target, dt_cur * shrink_min));
                    }

                    double ratio = target / dt_cur;
                    if (ratio > grow_max)   ratio = grow_max;
                    if (ratio < shrink_min) ratio = shrink_min;

                    dt_next = dt_cur * ratio;
                }

                if (needRetry)
                {
                    std::cout << "  [TimeStep] CFL limit violated, retry with dt = "
                        << retryDt << "\n";
                    if (!rollbackPrimaryFields()) return false;
                    dt_cur = retryDt;
                    continue;
                }

                // 5) ???????
                const double dt_used = dt_cur;
                time = trialTime;
                stepAccepted = true;

                {
                    const auto& assm = pressureCtrl.assembly;
                    if (!assm.rho_total_field.empty() && !assm.rho_total_old_field.empty())
                    {
                        if (!copyField(reg,
                            assm.rho_total_field,
                            assm.rho_total_old_field))
                        {
                            std::cerr << "[IMPES_revised] copy rho_total -> rho_total_old failed at step "
                                << stepId << ".\n";
                            return false;
                        }
                    }
                }

                std::cout << "  [TimeStep] dt_old = " << dt_used
                    << ", dt_new = " << dt_next << "\n";

                // 6) ??? / ????
                if (snapshotEveryCsv > 0 && !snapshotPrefix.empty() &&
                    (stepId % snapshotEveryCsv == 0))
                {
                    if (!::IMPES::Output::writeFieldSnapshotCSV(
                        snapshotPrefix, stepId, time, mgr, reg, snapshotFields))
                    {
                        return false;
                    }
                }

                const bool dumpP = !outPrefixP.empty() &&
                    ((writeEveryP <= 0) || (stepId % writeEveryP == 0));
                const bool dumpSw = !outPrefixSw.empty() &&
                    ((writeEverySw <= 0) || (stepId % writeEverySw == 0));

                if (dumpP)
                {
                    const std::vector<Vector> gradP =
                        computeCellGradients_LSQ_with_GG(
                            mgr.mesh(), reg,
                            pressureCtrl.assembly.pressure_field.c_str(), 0);

                    std::ostringstream fnP;
                    fnP << outPrefixP << "_step_"
                        << std::setw(5) << std::setfill('0') << stepId << ".plt";

                    const std::string faceFieldName =
                        pressureCtrl.assembly.pressure_field + "_face_tmp";

                    const bool okPltP = outputTecplot_cellToFaceToNode_BC(
                        mgr, reg, freg,
                        nullptr, nullptr,
                        pressureCtrl.assembly.pressure_field,
                        faceFieldName,
                        &gradP,
                        fnP.str());

                    if (!okPltP)
                    {
                        std::cerr << "[IMPES][P-s] Tecplot export (pressure) failed at step "
                            << stepId << ".\n";
                        return false;
                    }
                }

                if (dumpSw)
                {
                    const std::vector<Vector> gradSw =
                        computeCellGradients_LSQ_with_GG(
                            mgr.mesh(), reg, satCfg.saturation.c_str(), 0);

                    std::ostringstream fnSw;
                    fnSw << outPrefixSw << "_step_"
                        << std::setw(5) << std::setfill('0') << stepId << ".plt";

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
                        std::cerr << "[IMPES][P-s] Tecplot export (saturation) failed at step "
                            << stepId << ".\n";
                        return false;
                    }
                }

                dt_cur = dt_next;
                ++step;
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
            snapshotFields =
            {
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

                if (!solver_IMPES_Pressure(
                    mgr, reg, freg, Pbc, wells, dt, pressureState, &pRpt))
                {
                    std::cerr << "[IMPES][P-s] pressure solve failed at step " << stepId
                        << ", outer iter " << (outer + 1) << ".\n";
                    return false;
                }

                if (!IMPES_revised::advanceSaturationExplicit(
                    mgr, reg, freg, wells, dt, satCfg, &satRpt))
                {
                    std::cerr << "[IMPES_revised][P-s] saturation update failed at step "
                        << stepId << ".\n";
                    return false;
                }

                debugCheckMassFlux(mgr, freg, fluxState.gas_mass_flux.empty() ? "mf_g" : fluxState.gas_mass_flux, 1e-20);


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
}


namespace IMPES_Iteration
{
    struct PressureSolveControls
    {
        PressureAssemblyConfig assembly;
		
        
		//Todo add more controls if needed

    };

	//=======小工具函数=======//
    /**
    *\brief 
    * 计算两个 cell 标量场的最大无穷范数差值
    */ 
    inline double maxAbsDifference_volField
    (
        MeshManager& mgr,
        FieldRegistry& reg,
        const std::string& a_name,
        const std::string& b_name
    )
    {
        auto a = reg.get<volScalarField>(a_name);
        auto b = reg.get<volScalarField>(b_name);
        if (!a || !b) return 0.0;

        const auto& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        double m = 0.0;
        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const double diff = std::abs((*a)[i] - (*b)[i]);
            if (diff > m) m = diff;
        }
        return m;
    }


    /// 锐利阶跃前沿：
    ///   x <= x_front : S_w = Swr          （CO2 区，只剩残余水）
    ///   x >  x_front : S_w = s_w_inj      （未被扫到的水区，通常取 1.0）
    /// x_front = xmin + min(Lx, front_speed * t)
    inline bool applySaturationFrontPattern
    (
        MeshManager& mgr,
        FieldRegistry& reg,
        const SaturationTransportConfig& satCfg,
        double xmin,
        double Lx,
        int    stepId,
        double dt,
        double front_speed,
        double s_w_inj
    )
    {
        if (dt <= 0.0) {
            std::cerr << "[Analytic] applySaturationFrontPattern: invalid dt.\n";
            return false;
        }
        if (Lx <= 0.0) {
            std::cerr << "[Analytic] applySaturationFrontPattern: invalid Lx.\n";
            return false;
        }

        auto s_w = reg.get<volScalarField>(satCfg.saturation);
        if (!s_w) {
            std::cerr << "[Analytic] applySaturationFrontPattern: saturation field '"
                << satCfg.saturation << "' not found.\n";
            return false;
        }

        const auto& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        const double t = stepId * dt;
        const double x_front = xmin + std::min(Lx, front_speed * t);
        const double Sw_res = satCfg.vg_params.Swr;  // 残余水饱和度

        for (const auto& c : cells) {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const double x = c.center.m_x;

            // 饱和度：前沿左侧 Sw=Swr（残余饱和度），右侧 Sw=s_w_inj
            (*s_w)[i] = (x <= x_front) ? Sw_res : s_w_inj;
        }

        return true;
    }


    /// 光滑 tanh 型前沿：
    ///   x << x_front : S_w ≈ Swr
    ///   x >> x_front : S_w ≈ s_w_inj
    ///   过渡区厚度 w(t) = max(w_min, w0 + k * sqrt(t))
    inline bool applySaturationFront_SmoothTanh(
        MeshManager& mgr,
        FieldRegistry& reg,
        const SaturationTransportConfig& satCfg,
        double xmin,
        double Lx,
        int    stepId,
        double dt,
        double front_speed,
        double s_w_inj,
        double w0 = 1.0,   ///< 初始前沿半厚度 [m]
        double w_sqrtcoef = 0.0    ///< 随 sqrt(t) 增长的系数 [m / sqrt(s)]
    )
    {
        if (dt <= 0.0) {
            std::cerr << "[Analytic] applySaturationFront_SmoothTanh: invalid dt.\n";
            return false;
        }
        if (Lx <= 0.0) {
            std::cerr << "[Analytic] applySaturationFront_SmoothTanh: invalid Lx.\n";
            return false;
        }

        auto s_w = reg.get<volScalarField>(satCfg.saturation);
        if (!s_w) {
            std::cerr << "[Analytic] applySaturationFront_SmoothTanh: saturation field '"
                << satCfg.saturation << "' not found.\n";
            return false;
        }

        const auto& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        const double t = stepId * dt;
        const double x_front = xmin + std::min(Lx, front_speed * t);
        const double Sw_res = satCfg.vg_params.Swr;

        double w = w0 + w_sqrtcoef * std::sqrt(std::max(0.0, t));
        w = std::max(w, 1e-6); // 数值兜底，避免除零

        for (const auto& c : cells) {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const double x = c.center.m_x;

            const double xi = (x - x_front) / w;
            const double weight = 0.5 * (1.0 + std::tanh(xi)); // ∈(0,1)
            (*s_w)[i] = Sw_res + (s_w_inj - Sw_res) * weight;
        }

        return true;
    }

    /// 有限长度 CO2 slug：
    ///   x < x_tail      : S_w = s_w_inj
    ///   x_tail ≤ x ≤ x_head : S_w = Swr (slug 区, CO2 主导)
    ///   x > x_head      : S_w = s_w_inj
    /// slug 以 front_speed 向 +x 方向运动，长度 = slug_length
    inline bool applySaturationFront_Slug(
        MeshManager& mgr,
        FieldRegistry& reg,
        const SaturationTransportConfig& satCfg,
        double xmin,
        double Lx,
        int    stepId,
        double dt,
        double front_speed,
        double s_w_inj,
        double slug_length   ///< slug 的物理长度 [m]
    )
    {
        if (dt <= 0.0) {
            std::cerr << "[Analytic] applySaturationFront_Slug: invalid dt.\n";
            return false;
        }
        if (Lx <= 0.0) {
            std::cerr << "[Analytic] applySaturationFront_Slug: invalid Lx.\n";
            return false;
        }
        if (slug_length <= 0.0) {
            std::cerr << "[Analytic] applySaturationFront_Slug: invalid slug_length.\n";
            return false;
        }

        auto s_w = reg.get<volScalarField>(satCfg.saturation);
        if (!s_w) {
            std::cerr << "[Analytic] applySaturationFront_Slug: saturation field '"
                << satCfg.saturation << "' not found.\n";
            return false;
        }

        const auto& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        const double t = stepId * dt;
        const double x_head = xmin + std::min(Lx, front_speed * t);
        const double x_tail = std::max(xmin, x_head - slug_length);
        const double Sw_res = satCfg.vg_params.Swr;

        for (const auto& c : cells) {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const double x = c.center.m_x;

            if (x >= x_tail && x <= x_head) {
                // slug 区域：CO2 主导
                (*s_w)[i] = Sw_res;
            }
            else {
                // 其余区域：仍为水（或背景饱和度）
                (*s_w)[i] = s_w_inj;
            }
        }

        return true;
    }
  
    /*
    * /brief: 解析模式下的 IMPES 时间步推进测试：
    * - 不求解压力/饱和度方程；
    * - 每个时间步直接赋值  s_w 的解析模式；
    * - 调用两相物性更新与时间项 TimeTerm_IMPES_Pressure；
    * - 保留 Tecplot / CSV 等后处理输出。
    */
    inline bool runTransient_IMPES_AnalyticTest(
        MeshManager& mgr,
        FieldRegistry& reg,
		FaceFieldRegistry& freg,
        const PressureBCAdapter& Pbc,
        int nSteps,
        double dt,
        const PressureSolveControls& pressureCtrl,
		const SaturationTransportConfig& satCfg,
        int writeEveryP = 0,
        int writeEverySw = 0,
        const std::string& outPrefixP = "",
        const std::string& outPrefixSw = "",
        int snapshotEveryCsv = 0,
        const std::string& snapshotPrefix = "",
        std::vector<std::string> snapshotFields = {}
    )
    {
        if (!(dt > 0.0))
        {
            std::cerr << "[IMPES][AnalyticTest] invalid dt.\n";
            return false;
        }

        // ===== 1. Tecplot / CSV 输出目录准备（原样保留） =====
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
            std::cerr << "[IMPES][AnalyticTest] cannot create directories for Tecplot output prefixes.\n";
            return false;
        }
#endif



        // ===== 2. 控制参数与默认场名 =====
        const auto& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

		auto p_w = reg.get<volScalarField>(pressureCtrl.assembly.pressure_field);
        auto s_w = reg.get<volScalarField>(satCfg.saturation);

		if (!p_w || !s_w)
		{
			std::cerr << "[IMPES][AnalyticTest] pressure or saturation field not found.\n";
			return false;
		}

		auto p_w_old = reg.get<volScalarField>(pressureCtrl.assembly.pressure_old_field);
		auto p_w_prev = reg.get<volScalarField>(pressureCtrl.assembly.pressure_prev_field);
		auto s_w_old = reg.get<volScalarField>(satCfg.saturation_old);
		auto s_w_prev = reg.get<volScalarField>(satCfg.saturation_prev);
		
        if (!p_w_old || !p_w_prev || !s_w_old || !s_w_prev)
        {
            std::cerr << "[IMPES][AnalyticTest] missing old/prev fields for p_w or s_w.\n";
            return false;
        }

        const OperatorFieldNames nm = makeNames(pressureCtrl.assembly.operator_tag);

        if (snapshotFields.empty())
        {
            snapshotFields = {
                pressureCtrl.assembly.pressure_field,
                satCfg.saturation,
                nm.a_time,                             // 时间项对角系数
                nm.b_time                              // 时间项源项
            };
        }

        //===== 3. 解析模式几何信息（x_min, x_max, Lx） =====
        double xmin = std::numeric_limits<double>::max(); 
        double xmax = -std::numeric_limits<double>::max();
        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            xmin = std::min(xmin, c.center.m_x);
            xmax = std::max(xmax, c.center.m_x);
        }
        const double Lx = std::max(xmax - xmin, 1e-12);

        // ====4. 初始化物性：先用 t=0, p_w(x), Sw=1.0 更新一次，并把 basic props 写入 *_old
        TwoPhase::updateTwoPhasePropertiesAtTimeStep(mgr, reg, satCfg.saturation, satCfg.vg_params, satCfg.rp_params);  //计算各相的有效饱和度、各相相对渗透率、毛细压力
        TwoPhase::updateWaterBasicPropertiesAtStep(mgr, reg, pressureCtrl.assembly.pressure_field, "T"); //基于相对渗透率更新流度以及其他基本物性参数
        
        ///额外计算CO2相的压力
        auto Pc = reg.get<volScalarField>(TwoPhase::Auxiliaryparameters().Pc_tag);
		auto p_g = reg.getOrCreate<volScalarField>(pressureCtrl.assembly.pressure_g,p_w->data.size(),0.0);
        for (const auto& c : cells)
        {
			if (c.id < 0) continue;
			const size_t i = id2idx.at(c.id);
            
			(*p_g)[i] = (*p_w)[i] + (*Pc)[i];
        }

		TwoPhase::updateCO2BasicPropertiesAtStep(mgr, reg, pressureCtrl.assembly.pressure_g, "T"); //基于相对渗透率更新流度以及其他基本物性参数

		TwoPhase::copyBasicPropertiesToOldLayer(reg);  //把基本物性参数拷贝到旧时间层

       
        // lambda: 每个时间步，根据 stepId 赋值 t^{n+1} 时刻的解析 p_w / s_w
        ///解析模型参数设置及模型

        //////解析解case1 ：锐利前沿
        // //      const double front_speed = 0.1;     ///< 饱和度前沿速度 [m/s]
        // //      const double s_w_inj = 1.0;

        ////解析解case2:平滑 tanh 前沿案例（applySaturationFront_SmoothTanh） 
        //double front_speed = 1.0;      // [m/s]
        //double s_w_inj = 1.0;      // 右侧纯水
        //double w0 = 3.0;      // 起始半厚度 ~ 一个单元
        //double w_sqrtcoef = 0.5;      // 随时间略微变宽

        //解析解case3 :  Slug 案例（applySaturationFront_Slug）
        double front_speed = 0.5;      // [m/s]
        double s_w_inj = 1.0;      // slug 外部为纯水
        double slug_length = 20.0;     // slug 长度 20 m


		// ===== 5. 时间步循环 =====
        for (int step = 0; step < nSteps; ++step)
        {
            const int    stepId = step + 1;
            const double simTime = stepId * dt;  //记录模拟时间

            std::cout << "\n[IMPES][AnalyticTest] Time step " << stepId
                << " (t = " << simTime << " s, dt = " << dt << " s) =====\n";

            //1） 时间步启动，将 t^n 的 p_w / s_w 赋值到 old / prev  这里压力是常数其实都一样
            if (!startTimeStep_scalar(mgr.mesh(), reg,
                pressureCtrl.assembly.pressure_field,
                pressureCtrl.assembly.pressure_old_field,
                pressureCtrl.assembly.pressure_prev_field))
            {
                return false;
            }
            if (!startTimeStep_scalar(mgr.mesh(), reg, satCfg.saturation, satCfg.saturation_old, satCfg.saturation_prev))
            {
                std::cerr << "[IMPES][AnalyticTest] startTimeStep for saturation failed.\n";
                return false;
            }
            
            // 2) 解析赋值 t^{n+1} 的 S_w
            if (!applySaturationFront_Slug(
                mgr, reg, satCfg,
                xmin, Lx,
                stepId, dt,
                front_speed,
                s_w_inj,
                slug_length))
            {
                std::cerr << "[IMPES][AnalyticTest] applySaturationFrontPattern failed at step "
                    << stepId << ".\n";
                return false;
            }

			// 3) 得到新的两相物性参数
			TwoPhase::updateTwoPhasePropertiesAtTimeStep(mgr, reg, satCfg.saturation, satCfg.vg_params, satCfg.rp_params);  //计算各相的有效饱和度、各相相对渗透率、毛细压力
			TwoPhase::updateWaterBasicPropertiesAtStep(mgr, reg, pressureCtrl.assembly.pressure_field, "T"); //基于相对渗透率更新流度以及其他基本物性参数
            ///额外计算CO2相的压力
            auto Pc = reg.get<volScalarField>(TwoPhase::Auxiliaryparameters().Pc_tag);
            auto p_g = reg.getOrCreate<volScalarField>(pressureCtrl.assembly.pressure_g, p_w->data.size(), 0.0);
            for (const auto& c : cells)
            {
                if (c.id < 0) continue;
                const size_t i = id2idx.at(c.id);

                (*p_g)[i] = (*p_w)[i] + (*Pc)[i];
            }

            TwoPhase::updateCO2BasicPropertiesAtStep(mgr, reg, pressureCtrl.assembly.pressure_g, "T"); //基于相对渗透率更新流度以及其他基本物性参数


			// 4) 计算时间项系数

            //组装时间项

            if (!IMPES_Iteration::TimeTerm_IMPES_Pressure(mgr, reg, dt, pressureCtrl.assembly.pressure_old_field, pressureCtrl.assembly.pressure_field, nm.a_time, nm.b_time))
            {
                std::cerr << "[IMPES][AnalyticTest] TimeTerm_IMPES_Pressure failed at step "
                    << stepId << ".\n";
                return false;
            }

            // ---- 6.6：CSV snapshot 输出（若启用） ----
            if (snapshotEveryCsv > 0 && !snapshotPrefix.empty() && (stepId % snapshotEveryCsv == 0))
            {
                if (!::IMPES::Output::writeFieldSnapshotCSV(
                    snapshotPrefix, stepId, simTime, mgr, reg, snapshotFields))
                {
                    std::cerr << "[IMPES][AnalyticTest] CSV snapshot failed at step "
                        << stepId << ".\n";
                    return false;
                }
            }
            // ---- 6.7：Tecplot 输出 (P  / Sw) ----
            const bool dumpP = !outPrefixP.empty() &&
                ((writeEveryP <= 0) || (stepId % writeEveryP == 0));
            const bool dumpSw = !outPrefixSw.empty() &&
                ((writeEverySw <= 0) || (stepId % writeEverySw == 0));
            if (dumpP)
            {
				const std::vector<Vector>gradP =computeCellGradients_LSQ_with_GG(mgr.mesh(), reg,pressureCtrl.assembly.pressure_field.c_str(), 0);
                std::ostringstream fnP;
                fnP << outPrefixP << "_step_" << std::setw(5) << std::setfill('0')
                    << stepId << ".plt";
				const std::string faceFieldName = pressureCtrl.assembly.pressure_field + "_face_tmp";
                const bool okPltP = outputTecplot_cellToFaceToNode_BC(
                    mgr, reg, freg,
                    nullptr, &Pbc,
                    pressureCtrl.assembly.pressure_field,
                    faceFieldName,
                    &gradP,
                    fnP.str());

                if (!okPltP)
                {
                    std::cerr << "[IMPES][AnalyticTest] Tecplot export (pressure) failed at step "
                        << stepId << ".\n";
                    return false;
                }
            }
            if (dumpSw)
            {
                const std::vector<Vector> gradSw =
                    computeCellGradients_LSQ_with_GG(
                        mgr.mesh(), reg,
                        satCfg.saturation.c_str(), 0);

                std::ostringstream fnSw;
                fnSw << outPrefixSw << "_step_" << std::setw(5) << std::setfill('0')
                    << stepId << ".plt";

                const std::string faceFieldName =
                    satCfg.saturation + "_face_tmp";

                const bool okPltSw = outputTecplot_cellToFaceToNode_BC(
                    mgr, reg, freg,
                    nullptr, nullptr,
                    satCfg.saturation,
                    faceFieldName,
                    &gradSw,
                    fnSw.str());

                if (!okPltSw)
                {
                    std::cerr << "[IMPES][AnalyticTest] Tecplot export (saturation) failed at step "
                        << stepId << ".\n";
                    return false;
                }
            }
			TwoPhase::copyBasicPropertiesToOldLayer(reg);  //把基本物性参数拷贝到旧时间层
        }
        std::cout << "===== [IMPES][AnalyticTest] finished " << nSteps
            << " time steps. =====\n";
		return true;
    };
}

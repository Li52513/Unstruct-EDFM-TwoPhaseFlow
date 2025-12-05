#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#if __cplusplus >= 201703L
#include <filesystem>
#endif
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "PressureEqSolver.h"
#include "SaturationTransportEqAssemblerandSolver.h"
#include "0_PhysicalParametesCalculateandUpdata.h"
#include "PostProcess_.h"
#include "FluxSplitterandSolver.h"
#include "TwoPhaseWells_StrictRate.h"
#include "Solver_TimeLoopUtils.h"
#include "IMPES_PostProcessIO.h"

namespace IMPES_Iteration
{
    /// 时间步控制参数
    struct TimeStepControl
    {
        double dt_min = 1e-9;   ///< 允许的最小时间步 [s]
        double dt_max = 1000;   ///< 允许的最大时间步 [s]
        double grow_factor = 5;    ///< 接受时间步后，最大放大倍数
        double shrink_factor = 0.7;    ///< 拒绝时间步时收缩倍数
        double safety_factor = 0.9;    ///< 相对于 Redondo 建议步长的安全系数
        int    max_retries = 8;      ///< 同一物理步最多重试次数
    };
    inline bool runTransient_IMPES_Iteration(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& Pbc,
        std::vector<WellDOF_TwoPhase>& wells,
        int                        nSteps,
        double                     dt_initial,
        const PressureSolveControls& pressureCtrl,
        const SaturationTransportConfig& satCfg,
        const FluxSplitConfig& fluxCfg,
        int                        writeEveryP = 0,
        int                        writeEverySw = 0,
        const std::string& outPrefixP = "",
        const std::string& outPrefixSw = "",
        int                        snapshotEveryCsv = 0,
        const std::string& snapshotPrefix = "",
        std::vector<std::string>   snapshotFields = {}
    ) 
    {
        //================Step0 初始化设置，确认字符串以及导出文件是否合理===================//
        if (nSteps <= 0) {
            std::cerr << "[IMPES][Iteration] invalid nSteps.\n";
            return false;
        }
        if (!(dt_initial > 0.0)) {
            std::cerr << "[IMPES][Iteration] invalid initial dt.\n";
            return false;
        }
        // ==== 1. Tecplot / CSV output directory prep ====
#if __cplusplus >= 201703L
        try {
            if (!outPrefixP.empty()) {
                auto dirP = std::filesystem::path(outPrefixP).parent_path();
                if (!dirP.empty()) std::filesystem::create_directories(dirP);
            }
            if (!outPrefixSw.empty()) {
                auto dirSw = std::filesystem::path(outPrefixSw).parent_path();
                if (!dirSw.empty()) std::filesystem::create_directories(dirSw);
            }
            if (!snapshotPrefix.empty()) {
                auto dirSnap = std::filesystem::path(snapshotPrefix).parent_path();
                if (!dirSnap.empty()) std::filesystem::create_directories(dirSnap);
            }
        }
        catch (...) {
            std::cerr << "[IMPES][Iteration] cannot create directories for Tecplot / CSV output.\n";
            return false;
        }
#else
        // std::filesystem not available: skip directory creation
#endif
        TimeStepControl Tsc;
        Mesh& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        auto p_w = reg.get<volScalarField>(pressureCtrl.assembly .pressure_field);
        auto s_w = reg.get<volScalarField>(satCfg.saturation);
        auto p_w_old = reg.get<volScalarField>(pressureCtrl.assembly.pressure_old_field);
        auto p_w_prev = reg.get<volScalarField>(pressureCtrl.assembly.pressure_prev_field);
        auto s_w_old = reg.get<volScalarField>(satCfg.saturation_old);

        if (!p_w || !s_w || !p_w_old || !p_w_prev || !s_w_old ) 
        {
            std::cerr << "[IMPES][Iteration] missing pressure/saturation or old/prev fields.\n";
            return false;
        }
        TwoPhase::updateTwoPhasePropertiesAtTimeStep(mgr, reg, "s_w", satCfg.VG_Parameter.vg_params, satCfg.VG_Parameter.relperm_params);
        auto p_g = reg.get<volScalarField>("p_g");
        auto p_g_old = reg.get<volScalarField>("p_g_old");
        auto p_g_prev = reg.get<volScalarField>("p_g_prev");
        auto Pc = reg.get<volScalarField>(TwoPhase::Auxiliaryparameters().Pc_tag);
        if (!p_g || !p_g_old || !p_g_prev || !Pc)
        {
            std::cerr << "[IMPES][Iteration] missing p_g/Pc fields.\n";
            return false;
        }
        if (!satCfg.water_source_field.empty())
        {
            const size_t nCells = cells.size();
            reg.getOrCreate<volScalarField>(satCfg.water_source_field, nCells, 0.0);
        }
        //===========================================================================================//


       //====================================Step1 进入时间步循环========================================//
        double time = 0.0;
        double dt = dt_initial;
        std::cout << "[IMPES][Iteration] Starting transient run: " << "nSteps = " << nSteps<< ", dt_initial = " << dt_initial << " s\n";  
        for (int step = 0; step < nSteps; /* 手动控制 step++ */)
        {
            const int    stepId = step + 1;
            const double t_next = time + dt;
            std::cout << "\n[IMPES][Iteration] Time step " << stepId
                << " (t^n = " << time
                << " s, dt = " << dt
                << " s, t^{n+1} = " << t_next << " s) =====\n";

            auto rollback_after_failure = [&](const char* stage)->bool
            {
                copyField(reg,
                    pressureCtrl.assembly.pressure_old_field,
                    pressureCtrl.assembly.pressure_field);
                copyField(reg,
                    satCfg.saturation_old,
                    satCfg.saturation);
                TwoPhase::updateTwoPhasePropertiesAtTimeStep(
                    mgr, reg,
                    satCfg.saturation,
                    satCfg.VG_Parameter.vg_params,
                    satCfg.VG_Parameter.relperm_params
                );
                for (const auto& c : cells)
                {
                    if (c.id < 0) continue;
                    const size_t i = id2idx.at(c.id);
                    (*p_g)[i] = (*p_w)[i] + (*Pc)[i];
                }
                TwoPhase::updateWaterBasicPropertiesAtStep(
                    mgr, reg,
                    pressureCtrl.assembly.pressure_field,
                    "T"
                );
                TwoPhase::updateCO2BasicPropertiesAtStep(
                    mgr, reg,
                    pressureCtrl.assembly.pressure_g,
                    "T"
                );
                TwoPhase::copyBasicPropertiesToOldLayer(reg);

                dt *= Tsc.shrink_factor;
                if (dt < Tsc.dt_min) {
                    std::cerr << "[IMPES][Iteration] dt fell below dt_min after "
                        << stage << ".\n";
                    return false;
                }
                return true;
            };

            /// 1.1 备份当前时间层的 p_w / S_w（用 *_old 做回滚）
            if (!startTimeStep_scalar(
                mesh, reg,
                pressureCtrl.assembly.pressure_field,
                pressureCtrl.assembly.pressure_old_field,
                pressureCtrl.assembly.pressure_prev_field))
            {
                std::cerr << "[IMPES][Iteration] startTimeStep for pressure failed.\n";
                return false;
            }
            if (!startTimeStep_scalar(
                mesh, reg,
                satCfg.saturation,
                satCfg.saturation_old,
                satCfg.saturation_prev))
            {
                std::cerr << "[IMPES][Iteration] startTimeStep for saturation failed.\n";
                return false;
            }
           
            bool accept_step = true;                // 标记本时间步是否接受

            /// 1.2 总压力方程外迭代求解 --------
            double dp_prev = 0.0;
            double linRes_last = 0.0;
            int    linIters_last = 0;
            int    used_outer = 0;
            bool   p_converged = false;
            for (int it = 0; it < pressureCtrl.max_outer; ++it) //外迭代开始
            {
                //// 1.2.1 更新 p_g = p_w + Pc
                for (const auto& c : cells) {
                    if (c.id < 0) continue;
                    const size_t i = id2idx.at(c.id);
                    (*p_g)[i] = (*p_w)[i] + (*Pc)[i];
                }
                //// 1.2.2 更新基本物性（水 + CO2），温度字段统一用 "T"
                TwoPhase::updateWaterBasicPropertiesAtStep(
                    mgr, reg,
                    pressureCtrl.assembly.pressure_field,
                    "T"
                );
                TwoPhase::updateCO2BasicPropertiesAtStep(
                    mgr, reg,
                    pressureCtrl.assembly.pressure_g,
                    "T"
                );
                //// 1.2.3 单次外迭代：组装 + 线性求解 + 欠松弛
                double dp_inf = 0.0;
                double linRes = 0.0;
                int    linIts = 0;

                const bool okP = solver_IMPES_Iteration_PressureEq
                (
                    mgr, reg, freg,
                    Pbc,
                    wells,
                    dt,
                    pressureCtrl,
                    dp_inf,
                    linRes,
                    linIts
                );
                if (!okP) 
                {
                    std::cerr << "[IMPES][PressureStep] solver_IMPES_Iteration_PressureEq failed at outer it="
                        << it << ".\n";
                    accept_step = false;
                    break;
                }
                used_outer = it + 1;
                dp_prev = dp_inf;
                linRes_last = linRes;
                linIters_last = linIts;
                if (pressureCtrl.verbose) 
                {
                    std::cout << "[IMPES][Pressure] outer " << used_outer
                        << " dp_inf=" << dp_prev
                        << " linRes=" << linRes_last
                        << " linIters=" << linIters_last
                        << " (URF=" << pressureCtrl.under_relax << ")\n";
                }
                const bool conv_abs = (dp_prev <= pressureCtrl.tol_abs);
                if (conv_abs) 
                {
                    std::cout << "[IMPES][Pressure] converged at outer iter " << it
                        << " with dp_inf=" << dp_prev << "\n";
                    p_converged = true;
                    break;
                }
                if (it == pressureCtrl.max_outer - 1) 
                {
                    std::cout << "[IMPES][Pressure] reached max_outer_iterations without convergence.\n";
                }
            }//外迭代结束
            
            /// 1.3 外迭代结束后的收敛处理
            if (!p_converged || !accept_step)
            {
                std::cerr << "[IMPES][Iteration] Pressure not converged at time step "
                    << stepId << ", rollback and shrink dt.\n";

                // 使用统一的回滚逻辑：恢复 *_old，重算 p_g / 物性，缩小 dt
                if (!rollback_after_failure("pressure convergence failure"))
                {
                    return false;   // dt 已经小于 dt_min
                }
                continue;           // 重新计算该时间段（step 不自增）
            }

            /// 1.4 得到收敛的压力后，计算主扩散通量、毛细通量以及重力通量

            // -------- 4.2 两相总通量 → 水/气相通量拆分 --------
            {
                FluxSplitResult fluxRep;
                
                //if (!splitTwoPhaseMassFlux(mgr, reg, freg, fluxCfg,pressureCtrl.assembly, nullptr, &fluxRep))
                //{
                //    std::cerr << "[IMPES][Iteration] splitTwoPhaseMassFlux failed at time step "
                //        << stepId << ".\n";
                //    if (!rollback_after_failure("flux split failure"))
                //    {
                //        return false;
                //    }
                //    continue;
                //}
            }
            if (!satCfg.water_source_field.empty())
            {
                auto qField = reg.getOrCreate<volScalarField>(
                    satCfg.water_source_field,
                    cells.size(), 0.0);
                if (!wells.empty())
                {
                    std::vector<double> wellSources;
                    if (!FVM::TwoPhaseWellsStrict::build_saturation_well_sources_strict(
                        mgr,
                        reg,
                        wells,
                        satCfg.VG_Parameter.vg_params,
                        satCfg.VG_Parameter.relperm_params,
                        pressureCtrl.assembly.pressure_field,
                        wellSources))
                    {
                        std::cerr << "[IMPES][Iteration] failed to evaluate saturation well sources at time step "
                            << stepId << ".\n";
                        if (!rollback_after_failure("well source update"))
                        {
                            return false;
                        }
                        continue;
                    }
                    if (wellSources.size() != qField->data.size())
                    {
                        qField->data.assign(cells.size(), 0.0);
                    }
                    const size_t nWrite = std::min(qField->data.size(), wellSources.size());
                    for (size_t ic = 0; ic < nWrite; ++ic)
                    {
                        (*qField)[ic] = wellSources[ic];
                    }
                    for (size_t ic = nWrite; ic < qField->data.size(); ++ic)
                    {
                        (*qField)[ic] = 0.0;
                    }
                }
                else
                {
                    qField->data.assign(qField->data.size(), 0.0);
                }
            }
            // -------- 4.3 显式饱和度步（Euler） + Redondo / SimpleCFL 时间控制 --------
            SaturationStepStats satStats;
            if (!advanceSaturationExplicit_Euler(mgr, reg, freg, satCfg, dt, satStats))
            {
                std::cerr << "[IMPES][Iteration] splitTwoPhaseMassFlux failed at time step "
                    << stepId << ".\n";
                if (!rollback_after_failure("saturation update failure"))
                {
                    return false;
                }
                continue;
            }
            std::cout << "[IMPES][Saturation] max_CFL=" << satStats.max_CFL
                << ", max_dS=" << satStats.max_dS
                << ", dt_suggest=" << satStats.suggested_dt << "\n";
            
           
            
            // ---- 6.7：Tecplot 输出 (P  / Sw) ----
            const bool dumpP = !outPrefixP.empty() &&
                ((writeEveryP <= 0) || (stepId % writeEveryP == 0));
            const bool dumpSw = !outPrefixSw.empty() &&
                ((writeEverySw <= 0) || (stepId % writeEverySw == 0));
            if (dumpP)
            {
                const std::vector<Vector>gradP = computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, pressureCtrl.assembly.pressure_field.c_str(), 0);
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
            if (snapshotFields.empty())
            {
                snapshotFields =
                {
                    pressureCtrl.assembly.pressure_field,
                    satCfg.saturation
                };
            }
            // CSV snapshot（cell 标量场）
            if (snapshotEveryCsv > 0 &&
                (stepId % snapshotEveryCsv) == 0 &&
                !snapshotPrefix.empty())
            {
                const double simTime = t_next; // 本步的 t^{n+1}
                if (!IMPES::Output::writeFieldSnapshotCSV(
                    snapshotPrefix, stepId, simTime,
                    mgr, reg, snapshotFields))
                {
                    std::cerr << "[IMPES][Iteration] CSV snapshot failed at step "
                        << stepId << ".\n";
                    return false;
                }
            }


            // -------- 4.4 更新两相物性到新时间层，并写入 *_old --------
            TwoPhase::updateTwoPhasePropertiesAtTimeStep(
                mgr, reg,
                satCfg.saturation,
                satCfg.VG_Parameter.vg_params,
                satCfg.VG_Parameter.relperm_params
            );
            for (const auto& c : cells) {
                if (c.id < 0) continue;
                const size_t i = id2idx.at(c.id);
                (*p_g)[i] = (*p_w)[i] + (*Pc)[i];
            }
            TwoPhase::updateWaterBasicPropertiesAtStep(
                mgr, reg,
                pressureCtrl.assembly.pressure_field,
                "T"
            );
            TwoPhase::updateCO2BasicPropertiesAtStep(
                mgr, reg,
                pressureCtrl.assembly.pressure_g,
                "T"
            );
            TwoPhase::copyBasicPropertiesToOldLayer(reg);


            // -------- 4.6 接受时间步，更新时间和 dt --------
            time += dt;
            ++step;

            // 使用 Redondo / SimpleCFL 建议的 dt 作为主要约束，
            // 并允许适度增长（dt_growth_factor）
            double dt_target = satStats.suggested_dt;
            if (!(dt_target > 0.0)) {
                dt_target = dt; // 容错：如果算法返回了非正值，就保持原 dt
            }
            // 不能比当前 dt 增长太多，但也不超过建议值
            double dt_candidate = std::min({dt_target, Tsc.grow_factor * dt, Tsc.dt_max});
            dt = dt_candidate;

            std::cout << "[IMPES][Iteration] step " << stepId
                << " accepted. Next dt = " << dt << " s, simTime = "
                << time << " s\n";
        }

        std::cout << "\n[IMPES][Iteration] Transient run finished. Final time = "
            << time << " s\n";
        return true;
    }

} // namespace IMPES_Iteration

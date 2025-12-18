```C++
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
#include "PhysicalPropertiesManager_TwoPhase.h"

#include "Solver_TimeLoopUtils.h"
#include "PressureEqSolver.h"
#include "FaceMassRateCalculate.h"
#include "FluxSplitterandSolver.h"
#include "SaturationTransportEqAssemblerandSolver.h"

#include "PostProcess_.h"
#include "IMPES_PostProcessIO.h"

#include "TwoPhaseWells_StrictRate.h"
#include "IMPES_Iteration_WellHelpers.h"
#include "PostProcess_OutputReport.h"

namespace IMPES_Iteration
{
    /// 时间步控制参数
    struct TimeStepControl
    {
        double dt_min = 1e-5;   ///< 允许的最小时间步 [s]
        double dt_max = 100;   ///< 允许的最大时间步 [s]
        double grow_factor = 100;    ///< 接受时间步后，最大放大倍数
        double shrink_factor = 0.7;    ///< 拒绝时间步时收缩倍数
        double safety_factor = 0.9;    ///< 步长的安全系数
        int    max_retries = 8;      ///< 同一物理步最多重试次数
    };
    inline bool runTransient_IMPES_Iteration
    (
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& Pbc,
        std::vector<WellDOF_TwoPhase>& wells,
        int                        nSteps,
        double                     dt_initial,
         PressureSolveControls& pressureCtrl,
         SaturationTransportConfig& satCfg,
        const FluxSplitConfig& fluxCfg,
        const FaceMassRateConfig& FaceMassRateCfg,   //新增
        int                        writeEveryP = 0,
        int                        writeEverySw = 0,
        const std::string& outPrefixP = "",
        const std::string& outPrefixSw = "",
        int                        snapshotEveryCsv = 0,
        const std::string& snapshotPrefix = ""
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
        TwoPhase::updateTwoPhasePropertiesAtTimeStep(mgr, reg, satCfg.saturation, satCfg.VG_Parameter.vg_params, satCfg.VG_Parameter.relperm_params);
        auto p_g = reg.get<volScalarField>("p_g");
        auto p_g_old = reg.get<volScalarField>("p_g_old");
        auto p_g_prev = reg.get<volScalarField>("p_g_prev");
        auto Pc = reg.get<volScalarField>(PhysicalProperties_string::TwoPhase_case().Pc_tag);
        if (!p_g || !p_g_old || !p_g_prev || !Pc)
        {
            std::cerr << "[IMPES][Iteration] missing p_g/Pc fields.\n";
            return false;
        }
        // 初始化 p_g = p_w + Pc，并更新一次物性 old 层
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

        if (!satCfg.water_source_field.empty())
        {
            const size_t nCells = cells.size();
            reg.getOrCreate<volScalarField>(satCfg.water_source_field, nCells, 0.0);
        }
        // 定义一个水相井源场名
        constexpr const char* Qw_well_name = "Qw_well";

        //===========================================================================================//


       //====================================Step1 进入时间步循环========================================//
        double time = 0.0;
        double dt = dt_initial;
        std::cout << "[IMPES][Iteration] Starting transient run: " << "nSteps = " << nSteps<< ", dt_initial = " << dt_initial << " s\n";  
        std::vector<std::string> snapshotFields;
        const std::string fluxTimeSeriesFile =
            snapshotPrefix.empty()
            ? "boundary_mass_flux_timeseries.csv"
            : (snapshotPrefix + "_boundary_mass_flux_timeseries.csv");

        const std::string WaterfluxTimeSeriesFile =
            snapshotPrefix.empty()
            ? "boundary_mass_flux_timeseries.csv"
            : (snapshotPrefix + "_boundary_water_mass_flux_timeseries.csv");

        const std::string CO2TimeSeriesFile =
            snapshotPrefix.empty()
            ? "boundary_mass_flux_timeseries.csv"
            : (snapshotPrefix + "_boundary_CO2_mass_flux_timeseries.csv");
        const std::string wellBdryCsvFile =
            snapshotPrefix.empty()
            ? "well_boundary_timeseries.csv"
            : (snapshotPrefix + "__well_boundary_timeseries.csv");

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
           
            bool accept_step = true; // 标记本时间步是否接受

            /// 1.2 总压力方程外迭代求解 --------
            double dp_prev = 0.0;
            double linRes_last = 0.0;
            int    linIters_last = 0;
            int    used_outer = 0;
            bool   p_converged = false;

            for (int it = 0; it < pressureCtrl.max_outer; ++it) //外迭代开始
            {
                //// 1.2.1 更新 p_g = p_w + Pc
                for (const auto& c : cells) 
                {
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
                //if (!rollback_after_failure("pressure convergence failure"))
                //{
                //    return false;   // dt 已经小于 dt_min
                //}
                continue;           // 重新计算该时间段（step 不自增）
            }

            if (accept_step)
            {
                // 遍历 wells，找出所有 Pressure 模式的生产井，依次打印
                for (const auto& w : wells)
                {
                    if (w.role == WellDOF_TwoPhase::Role::Producer &&
                        w.mode == WellDOF_TwoPhase::Mode::Pressure)
                    {
                        IMPES_Iteration::debugPrintPeacemanPIAndRate(
                            mgr,
                            reg,
                            w,
                            pressureCtrl.assembly.pressure_field  // 一般是 "p_w"
                        );
                    }
                    if (w.role == WellDOF_TwoPhase::Role::Producer &&
                        w.mode == WellDOF_TwoPhase::Mode::Rate)
                    {
                        IMPES_Iteration::debugPrintPeacemanPIAndRate(
                            mgr,
                            reg,
                            w,
                            pressureCtrl.assembly.pressure_field  // 一般是 "p_w"
                        );
                    }
                    if (w.role == WellDOF_TwoPhase::Role::Injector &&
                        w.mode == WellDOF_TwoPhase::Mode::Pressure)
                    {
                        IMPES_Iteration::debugPrintPeacemanPIAndRate(
                            mgr,
                            reg,
                            w,
                            pressureCtrl.assembly.pressure_field  // 一般是 "p_w"
                        );
                    }
                    if (w.role == WellDOF_TwoPhase::Role::Injector &&
                        w.mode == WellDOF_TwoPhase::Mode::Rate)
                    {
                        IMPES_Iteration::debugPrintPeacemanPIAndRate(
                            mgr,
                            reg,
                            w,
                            pressureCtrl.assembly.pressure_field  // 一般是 "p_w"
                        );
                    }

                }
            }
            /// 1.4 压力收敛后：构建面质量通量
            if (accept_step)
            {
                if (!buildFaceMassRates(mgr, reg, freg, pressureCtrl.assembly, Pbc, FaceMassRateCfg))
                {
                    std::cerr << "[IMPES_Iteration] buildFaceMassRates failed after pressure convergence.\n";
                    accept_step = false;
                }
            }
            /// 1.5 两相通量分配，mf_total -> mf_w / mf_g
            FluxSplitResult fluxRep; // 可选：用于收集诊断信息（max|fw-?, ...）
            FaceSignMask faceMask;   // 可选：记录面方向 + 上游信息
            if (accept_step)
            {
                if (!splitTwoPhaseMassFlux(mgr, reg, freg, fluxCfg, FaceMassRateCfg, &faceMask, &fluxRep))
                {
                    std::cerr << "[IMPES_Iteration] splitTwoPhaseMassFlux failed.\n";
                    if (!rollback_after_failure("flux split")) return false;
                    accept_step = false;
                }
                //// 3.5) 调试：检查边界 inflow/outflow 通量与 fw 因为边界上忽略了毛细通量和重力通量
                //IMPES_Iteration::diagnoseBoundaryFluxWithFw(
                //    mgr,
                //    freg,
                //    FaceMassRateCfg.total_mass_flux,        // 对应 mf_total
                //    fluxCfg.fractional_flow_face,   // 对应 fw_face，一般是 "fw_face"
                //    fluxCfg.flux_sign_epsilon       // 与 splitTwoPhaseMassFlux 使用同一阈值
                //);
            }

            if (!IMPES_Iteration::buildWaterSourceFieldFromWells(mgr, reg, wells, satCfg.VG_Parameter.vg_params, satCfg.VG_Parameter.relperm_params, pressureCtrl.assembly.pressure_field, Qw_well_name))
            {
                std::cerr << "[IMPES][Well] failed to build water well sources, abort step.\n";
                return false;
            }
            satCfg.water_source_field = Qw_well_name;
            /// 1.5 显式推进水相饱和度 S_w
            SaturationStepStats satStats;
            if (accept_step)
            {
                bool sat_ok = false;

                switch (satCfg.time_integration_scheme)
                {
                case SatTimeIntegrationScheme::ExplicitEuler:
                    sat_ok = advanceSaturationExplicit_Euler(
                        mgr, reg, freg, satCfg, fluxCfg, dt, satStats);
                    break;

                case SatTimeIntegrationScheme::HeunRK2:
                    sat_ok = advanceSaturationHeun_RK2(
                        mgr, reg, freg, satCfg, fluxCfg, dt, satStats);
                    break;

                default:
                    std::cerr << "[IMPES][Saturation] Unknown time integration scheme.\n";
                    sat_ok = false;
                    break;
                }

                if (!sat_ok)
                {
                    std::cerr << "[IMPES][Iteration] advanceSaturation("
                        << (satCfg.time_integration_scheme ==
                            SatTimeIntegrationScheme::HeunRK2
                            ? "Heun_RK2" : "Explicit_Euler")
                        << ") failed at step " << stepId << ".\n";

                    if (!rollback_after_failure("saturation"))
                        return false;

                    accept_step = false;
                }
                else
                {
                    std::cout << "[IMPES][Saturation] step " << stepId
                        << " (" << (satCfg.time_integration_scheme ==
                            SatTimeIntegrationScheme::HeunRK2
                            ? "Heun/RK2" : "Euler")
                        << "): max_CFL=" << satStats.max_CFL
                        << ", max_dS=" << satStats.max_dS
                        << ", dt_suggest=" << satStats.suggested_dt << "\n";
                }
            }

            //====================================Step2 当前时间步计算结束，时步更新并结果输出========================================//           
            if (accept_step)
            {
                // ---- 2.0 统计边界总质量 inflow / outflow，并写入时间序列 CSV ----
                {
                    IMPES::Output::MassFluxSummary mfSummary;
                    if (IMPES::Output::computeMassFluxSummary(
                        mgr,
                        freg,
                        FaceMassRateCfg.total_mass_flux,  // 总质量通量面场名
                        mfSummary))
                    {
                        IMPES::Output::TimeSeriesRecord rec;
                        rec.step = stepId;
                        rec.simTime = t_next;            // 当前步的 t^{n+1}
                        rec.outerIterations = used_outer;        // 压力外迭代次数
                        rec.dp_inf = dp_prev;           // 最后一轮 dp_inf
                        rec.dT_inf = 0.0;               // 当前没有温度外迭代，用 0 占位
                        rec.satCFL = satStats.max_CFL;  // 饱和度 CFL 诊断
                        rec.satSuggestedDt = satStats.suggested_dt;
                        rec.pressureLinIters = linIters_last;     // 线性求解器迭代数
                        rec.temperatureLinIters = 0;                 // 没有温度方程，写 0
                        rec.flux = mfSummary;         // 这里面包含 boundaryInflow/outflow

                        IMPES::Output::appendTimeSeriesRecord(
                            fluxTimeSeriesFile, rec);
                    }
                    else
                    {
                        std::cerr << "[IMPES][Iteration] computeMassFluxSummary failed at step "
                            << stepId << ".\n";
                    }

                    IMPES::Output::MassFluxSummary mfSummary_water;
                    if (IMPES::Output::computeMassFluxSummary(
                        mgr,
                        freg,
                        fluxCfg.water_mass_flux,  // 总质量通量面场名
                        mfSummary_water))
                    {
                        IMPES::Output::TimeSeriesRecord rec;
                        rec.step = stepId;
                        rec.simTime = t_next;            // 当前步的 t^{n+1}
                        rec.outerIterations = used_outer;        // 压力外迭代次数
                        rec.dp_inf = dp_prev;           // 最后一轮 dp_inf
                        rec.dT_inf = 0.0;               // 当前没有温度外迭代，用 0 占位
                        rec.satCFL = satStats.max_CFL;  // 饱和度 CFL 诊断
                        rec.satSuggestedDt = satStats.suggested_dt;
                        rec.pressureLinIters = linIters_last;     // 线性求解器迭代数
                        rec.temperatureLinIters = 0;                 // 没有温度方程，写 0
                        rec.flux = mfSummary_water;         // 这里面包含 boundaryInflow/outflow

                        IMPES::Output::appendTimeSeriesRecord(
                            WaterfluxTimeSeriesFile, rec);
                    }
                    else
                    {
                        std::cerr << "[IMPES][Iteration] compute_Water_MassFluxSummary failed at step "
                            << stepId << ".\n";
                    }

                    IMPES::Output::MassFluxSummary mfSummary_CO2;
                    if (IMPES::Output::computeMassFluxSummary(
                        mgr,
                        freg,
                        fluxCfg.gas_mass_flux,  // 总质量通量面场名
                        mfSummary_CO2))
                    {
                        IMPES::Output::TimeSeriesRecord rec;
                        rec.step = stepId;
                        rec.simTime = t_next;            // 当前步的 t^{n+1}
                        rec.outerIterations = used_outer;        // 压力外迭代次数
                        rec.dp_inf = dp_prev;           // 最后一轮 dp_inf
                        rec.dT_inf = 0.0;               // 当前没有温度外迭代，用 0 占位
                        rec.satCFL = satStats.max_CFL;  // 饱和度 CFL 诊断
                        rec.satSuggestedDt = satStats.suggested_dt;
                        rec.pressureLinIters = linIters_last;     // 线性求解器迭代数
                        rec.temperatureLinIters = 0;                 // 没有温度方程，写 0
                        rec.flux = mfSummary_CO2;         // 这里面包含 boundaryInflow/outflow

                        IMPES::Output::appendTimeSeriesRecord(
                            CO2TimeSeriesFile, rec);
                    }
                    else
                    {
                        std::cerr << "[IMPES][Iteration] compute_Water_MassFluxSummary failed at step "
                            << stepId << ".\n";
                    }
                }
                // ---- 2.1：Tecplot 输出 (P  / Sw) ----
                const bool dumpP = !outPrefixP.empty() &&((writeEveryP <= 0) || (stepId % writeEveryP == 0));
                const bool dumpSw = !outPrefixSw.empty() &&((writeEverySw <= 0) || (stepId % writeEverySw == 0));
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

                // =======================================================================
                // [BEGIN] 井 + 边界 时间序列 CSV 输出
                // =======================================================================
                using IMPES_Iteration::WellMassReport;
                using IMPES_Iteration::BoundaryMassReport;
                using IMPES_Iteration::WellBoundaryCSVRow;

                // 3.1 收集每口井的质量流量报告
                std::vector<WellMassReport> wellReports;
                {
                    Mesh& mesh = mgr.mesh();
                    const auto& cells = mesh.getCells();
                    const auto& id2idx = mesh.getCellId2Index();
                    const int Nc = static_cast<int>(cells.size());

                    std::vector<double> Mw_cell, Mg_cell; // 临时缓存，不同井复用

                    for (const auto& w : wells)
                    {

                        WellMassReport rep;
                        rep.name = w.name;
                        rep.role = w.role;
                        rep.mode = w.mode;
                        rep.Tin = w.Tin;

                        // Pressure 模式：BHP = target；Rate 模式：BHP = 解出的 p_bh
                        if (w.mode == WellDOF_TwoPhase::Mode::Pressure) {
                            rep.p_bh = w.target; // Pa
                        }
                        else {
                            rep.p_bh = w.p_bh;   // 由线性系统给出
                        }

                        // 调用你已经实现好的 严格 Peaceman 两相质量率计算函数
                        if (!FVM::TwoPhaseWellsStrict::compute_well_phase_mass_rates_strict(
                            mgr,
                            reg,
                            w,
                            pressureCtrl.assembly.VG_Parameter.vg_params,
                            pressureCtrl.assembly.VG_Parameter.relperm_params,
                            pressureCtrl.assembly.pressure_field, // 一般是 "p_w"
                            Mw_cell,
                            Mg_cell))
                        {
                            std::cerr << "[WellReport] compute_well_phase_mass_rates_strict failed for well "
                                << w.name << ".\n";
                            continue;
                        }

                        // 汇总这口井的总 Mw / Mg
                        double Mw_tot = 0.0;
                        double Mg_tot = 0.0;
                        for (int ic = 0; ic < Nc; ++ic) {
                            Mw_tot += Mw_cell[ic];
                            Mg_tot += Mg_cell[ic];
                        }
                        rep.Mw_tot = Mw_tot;
                        rep.Mg_tot = Mg_tot;

                        const double denom = std::abs(Mw_tot) + std::abs(Mg_tot);
                        rep.fw_mass = (denom > 1e-30) ? (Mw_tot / denom) : 0.0; // 水相质量分率 fw

                        // 找一个代表性的 hostCellId：第一个完井单元 ID
                        rep.hostCellId = -1;
                        auto mask_field = reg.get<volScalarField>(w.mask_field.c_str());
                        if (mask_field)
                        {
                            for (int ic = 0; ic < Nc; ++ic) {
                                const auto& c = cells[ic];
                                if (c.id < 0) continue;
                                auto it = id2idx.find(c.id);
                                if (it == id2idx.end()) continue;
                                size_t idx = it->second;
                                if (idx < mask_field->data.size()
                                    && (*mask_field)[idx] > 0.5)  // 视为完井单元
                                {
                                    rep.hostCellId = c.id;
                                    break;
                                }
                            }
                        }

                        wellReports.push_back(rep);
                    }
                }

                // 3.2 诊断边界外流，并同时汇总边界水/气质量通量
                BoundaryMassReport bdryRep;
                IMPES_Iteration::diagnoseBoundaryFluxWithFw(
                    mgr,
                    freg,
                    FaceMassRateCfg.total_mass_flux,    // 总质量通量面场名
                    fluxCfg.fractional_flow_face,      // fw_face 面场名，比如 "fw_face"
                    fluxCfg.flux_sign_epsilon,         // 判定流向阈值
                    &bdryRep                            // 额外：把结果写到这个结构里
                );

                // 3.3 基于 wellReports + bdryRep 构造一行 CSV 的数据
                WellBoundaryCSVRow row;
                {
                    // 这里的 time 变量应该是当前物理时间 t^{n+1}
                    // 如果你的代码用的是别的名字（比如 t 或 tn1），替换一下即可
                    double t_phys = time;

                    IMPES_Iteration::buildWellBoundaryRowFromReports(
                        mgr,
                        reg,
                        wellReports,
                        pressureCtrl.assembly.pressure_field, // 用来取井所在基岩单元压力
                        bdryRep,
                        t_phys,
                        row
                    );
                }
                ensureWellBoundaryCSVHeader(wellBdryCsvFile);
                // 3.4 追加写入 CSV 文件
                IMPES_Iteration::appendWellBoundaryCSVRow(
                    wellBdryCsvFile,
                    row
                );
            }
            if (accept_step)
            {
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

            }
            //====================================Step3 更新时间步========================================//  
            if (accept_step) 
            {
                double dt_new = dt;
                if (satStats.suggested_dt > 0.0 &&
                    std::isfinite(satStats.suggested_dt))
                {
                    double dt_suggest = satStats.suggested_dt;
                    dt_suggest = std::max(Tsc.dt_min,
                                   std::min(Tsc.dt_max, dt_suggest));

                    const double dt_grow = dt * Tsc.grow_factor;
                    dt_new = std::min(dt_grow,
                                      Tsc.safety_factor * dt_suggest);
                }
                else
                {
                    // 如果饱和度没给出靠谱的建议步长，就保守缩小一点
                    dt_new = dt * Tsc.shrink_factor;
                }
                dt = std::max(Tsc.dt_min, std::min(Tsc.dt_max, dt_new));

                // ---- 真正推进物理时间 & 步号 ----
                time = t_next;  // 本步物理时间走到 t^{n+1}
                ++step;
            }
        }

        std::cout << "\n[IMPES][Iteration] Transient run finished. Final time = "
            << time << " s\n";
        return true;
    }

} // namespace IMPES_Iteration
```


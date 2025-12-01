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
#include "Solver_TimeLoopUtils.h"

namespace IMPES_Iteration
{
    /// 时间步控制参数
    struct TimeStepControl
    {
        double dt_min = 1e-5;   ///< 允许的最小时间步 [s]
        double dt_max = 1e10;   ///< 允许的最大时间步 [s]
        double grow_factor = 1.3;    ///< 接受时间步后，最大放大倍数
        double shrink_factor = 0.5;    ///< 拒绝时间步时收缩倍数
        double safety_factor = 0.9;    ///< 相对于 Redondo 建议步长的安全系数
        int    max_retries = 8;      ///< 同一物理步最多重试次数
    };

    /**
 * \brief 基于 IMPES 的两相 CO2–H2O 主时间步推进（迭代版）。
 *
 * 时间推进顺序：
 *  - t = 0：根据初始 p_w, S_w, T 计算两相 VG / kr / Pc 与水/CO2 物性，并写入 *_old；
 *  - 每个时间步：
 *      1) 启动时间步：p_w, S_w → *_old, *_prev；
 *      2) 基于当前 S_w^n 更新 VG/kr/Pc；
 *      3) 压力步：外迭代 + 欠松弛，调用 solver_IMPES_Iteration_PressureEq，得到收敛的 p_w^{n+1}
 *         和总质量通量（face field，由压力装配模块写入 FaceFieldRegistry）；
 *      4) 调用两相通量拆分函数，将总质量通量拆成水相/气相质量通量（face field）；
 *      5) 显式 Euler 饱和度推进，调用 advanceSaturationExplicit_Euler，统计 CFL / ΔS / dt_suggest；
 *      6) 基于 p_w^{n+1}, S_w^{n+1}, T 更新两相物性，并拷贝到 *_old 供下一时间步使用；
 *      7) 按需输出 Tecplot / CSV / time-series。
 *
 * \note
 *  - 时间步控制：
 *      - 在 advanceSaturationExplicit_Euler 内部已经根据 cfg.time_control_scheme
 *        计算了 Redondo 风格或 SimpleCFL 风格的 dt_suggest；
 *      - 本函数仅用 stats.suggested_dt 和简单增长因子来更新 dt。
 *  - 所有温度字段均统一用 "T" 字符串。
 *  - 本函数依赖：
 *      - solver_IMPES_Iteration_PressureEq(...)
 *      - splitTwoPhaseMassFlux(...)
 *      - advanceSaturationExplicit_Euler(...)
 *      - TwoPhase::updateTwoPhasePropertiesAtTimeStep(...)
 *      - TwoPhase::updateWaterBasicPropertiesAtStep(...)
 *      - TwoPhase::updateCO2BasicPropertiesAtStep(...)
 *      - TwoPhase::copyBasicPropertiesToOldLayer(...)
 */
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
        if (nSteps <= 0) {
            std::cerr << "[IMPES][Iteration] invalid nSteps.\n";
            return false;
        }
        if (!(dt_initial > 0.0)) {
            std::cerr << "[IMPES][Iteration] invalid initial dt.\n";
            return false;
        }
        // ==== 1. Tecplot / CSV 输出目录准备 ====
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
#endif

        // ===== 2. 基本 field / mesh 句柄 =====
        TimeStepControl Tsc;
        Mesh& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        auto p_w = reg.get<volScalarField>(pressureCtrl.assembly .pressure_field);
        auto s_w = reg.get<volScalarField>(satCfg.saturation);
        auto p_w_old = reg.get<volScalarField>(pressureCtrl.assembly.pressure_old_field);
        auto p_w_prev = reg.get<volScalarField>(pressureCtrl.assembly.pressure_prev_field);
        auto s_w_old = reg.get<volScalarField>(satCfg.saturation_old);
        auto s_w_prev = reg.get<volScalarField>(satCfg.saturation_prev);

        if (!p_w || !s_w || !p_w_old || !p_w_prev || !s_w_old || !s_w_prev) {
            std::cerr << "[IMPES][Iteration] missing pressure/saturation or old/prev fields.\n";
            return false;
        }

        // 操作符字段名（时间项系数等），用于 snapshot 时导出
        const OperatorFieldNames opNames = makeNames(pressureCtrl.assembly.operator_tag);

        if (snapshotFields.empty()) 
        {
            snapshotFields = {
                pressureCtrl.assembly.pressure_field,
                satCfg.saturation
            };
        }

        // VG 中的 Pc 字段：若不存在则按照 cell 数量创建
        auto Pc = reg.get<volScalarField>(TwoPhase::Auxiliaryparameters().Pc_tag);
        if (!Pc)
        {
            Pc = reg.getOrCreate<volScalarField>(
                TwoPhase::Auxiliaryparameters().Pc_tag,
                p_w->data.size(),  // 与压力场同长度
                0.0                // 初值 0，后面 Step 3 会被 VG 更新覆盖
            );
            if (!Pc)
            {
                std::cerr << "[IMPES][Iteration] failed to create Pc field '"
                    << TwoPhase::Auxiliaryparameters().Pc_tag << "'.\n";
                return false;
            }
        }

        // 气相压力场 p_g（可能在别处创建，这里确保存在）
        auto p_g = reg.getOrCreate<volScalarField>(
            pressureCtrl.assembly.pressure_g,
            p_w->data.size(), 0.0
        );

        // ===== 3. t = 0：初始化两相物性，并写入 *_old =====
        {
            // 3.1: 基于初始 S_w^0 更新 VG / kr / Pc
            TwoPhase::updateTwoPhasePropertiesAtTimeStep
            (
                mgr, reg,
                satCfg.saturation,
                satCfg.VG_Parameter.vg_params,
                satCfg.VG_Parameter.relperm_params
            );

            // 3.2: p_g^0 = p_w^0 + Pc(S_w^0)
            for (const auto& c : cells) 
            {
                if (c.id < 0) continue;
                const size_t i = id2idx.at(c.id);
                (*p_g)[i] = (*p_w)[i] + (*Pc)[i];
            }

            // 3.3: 初始基本物性（水 + CO2），温度字段统一用 "T"
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
            // 3.4: 把当前基本物性写入 *_old
            TwoPhase::copyBasicPropertiesToOldLayer(reg);
        }
        // ===== 4. 时间循环 =====
        double time = 0.0;
        double dt = dt_initial;
        std::cout << "[IMPES][Iteration] Starting transient run: "
            << "nSteps = " << nSteps
            << ", dt_initial = " << dt_initial << " s\n";
        for (int step = 0; step < nSteps; /* 手动控制 step++ */)
        {
            const int    stepId = step + 1;
            const double t_next = time + dt;
            std::cout << "\n[IMPES][Iteration] Time step " << stepId
                << " (t^n = " << time
                << " s, dt = " << dt
                << " s, t^{n+1} = " << t_next << " s) =====\n";

            // 4.0 备份当前时间层的 p_w / S_w（用 *_old 做回滚）
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
            // 标记本时间步是否接受
            bool accept_step = true;

            // -------- 4.1 压力步：外迭代 --------
            double dp_prev = std::numeric_limits<double>::infinity();
            double dp_last = dp_prev;
            double linRes_last = 0.0;
            int    linIters_last = 0;
            int    used_outer = 0;
            bool   p_converged = false;

            for (int it = 0; it < pressureCtrl.max_outer; ++it)
            {
                // 4.1.1: 更新 p_g = p_w + Pc
                for (const auto& c : cells) {
                    if (c.id < 0) continue;
                    const size_t i = id2idx.at(c.id);
                    (*p_g)[i] = (*p_w)[i] + (*Pc)[i];
                }
                // 4.1.2: 更新基本物性（水 + CO2），温度字段统一用 "T"
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
                // 4.1.3: 单次外迭代：组装 + 线性求解 + 欠松弛
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

                if (!okP) {
                    std::cerr << "[IMPES][PressureStep] solver_IMPES_Iteration_PressureEq failed at outer it="
                        << it << ".\n";
                    accept_step = false;
                    break;
                }

                used_outer = it + 1;
                dp_prev = dp_last;
                dp_last = dp_inf;
                linRes_last = linRes;
                linIters_last = linIts;

                if (pressureCtrl.verbose) 
                {
                    std::cout << "[IMPES][Pressure] outer " << used_outer
                        << " dp_inf=" << dp_last
                        << " linRes=" << linRes_last
                        << " linIters=" << linIters_last
                        << " (URF=" << pressureCtrl.under_relax << ")\n";
                }

                const bool conv_abs = (dp_last <= pressureCtrl.tol_abs);
                const bool conv_rel = (it > 0 && dp_last <= pressureCtrl.tol_rel * dp_prev);

                if (conv_abs || conv_rel) {
                    std::cout << "[IMPES][Pressure] converged at outer iter " << it
                        << " with dp_inf=" << dp_last << "\n";
                    p_converged = true;
                    break;
                }
                if (it == pressureCtrl.max_outer - 1) 
                {
                    std::cout << "[IMPES][Pressure] reached max_outer_iterations without convergence.\n";
                }
            }

            if (!p_converged || !accept_step)
            {
                // 压力步没收敛：回滚 + 缩小 dt 重算该步
                std::cerr << "[IMPES][Iteration] Pressure not converged at time step "
                    << stepId << ", rollback and shrink dt.\n";

                // 回滚 p_w, S_w 到 old 层
                copyField(reg,
                    pressureCtrl.assembly.pressure_old_field,
                    pressureCtrl.assembly.pressure_field);
                copyField(reg,
                    satCfg.saturation_old,
                    satCfg.saturation);

                // 回滚物性：简单做法是重新基于回滚后的 p_w, S_w 计算
                TwoPhase::updateTwoPhasePropertiesAtTimeStep
                (
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
                    std::cerr << "[IMPES][Iteration] dt fell below dt_min after pressure rollback.\n";
                    return false;
                }
                // 不推进时间，不加 step，直接重算本时间步
                continue;
            }

            // -------- 4.2 两相总通量 → 水/气相通量拆分 --------
            {
                FluxSplitResult fluxRep;
                
                if (!splitTwoPhaseMassFlux(mgr, reg, freg, fluxCfg, nullptr, &fluxRep))
                {
                    std::cerr << "[IMPES][Iteration] splitTwoPhaseMassFlux failed at time step "
                        << stepId << ".\n";
                    // 保守一点：回滚 + 减小 dt
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

                    dt *= Tsc.shrink_factor;
                    if (dt < Tsc.dt_min) {
                        std::cerr << "[IMPES][Iteration] dt fell below dt_min after flux split failure.\n";
                        return false;
                    }
                    continue;
                }
            }
            // -------- 4.3 显式饱和度步（Euler） + Redondo / SimpleCFL 时间控制 --------
            SaturationStepStats satStats;
            if (!advanceSaturationExplicit_Euler(mgr, reg, freg, satCfg, dt, satStats))
            {
                std::cerr << "[IMPES][Iteration] splitTwoPhaseMassFlux failed at time step "
                    << stepId << ".\n";
                // 保守一点：回滚 + 减小 dt
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

                dt *= Tsc.shrink_factor;
                if (dt < Tsc.dt_min) {
                    std::cerr << "[IMPES][Iteration] dt fell below dt_min after flux split failure.\n";
                    return false;
                }
                continue;
            }
            std::cout << "[IMPES][Saturation] max_CFL=" << satStats.max_CFL
                << ", max_dS=" << satStats.max_dS
                << ", dt_suggest=" << satStats.suggested_dt << "\n";

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
            double dt_candidate = std::min(dt_target, Tsc.grow_factor * dt);
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

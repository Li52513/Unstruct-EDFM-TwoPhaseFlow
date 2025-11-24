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


#include"PressureEqAssemblerandSolver.h"
#include "IMPES_SaturationTransportEqAssemblerandSolver.h"
#include "IMPES_PostProcessIO.h"
#include "TimeTermAssemblerandSolver.h"
#include "Diff_TPFA_GradientsOperation.h"
#include "PostProcess_.h"



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
    };
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
    };
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
    };
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
    };

//===========================================================================================//
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
        auto p_g = reg.getOrCreate<volScalarField>(pressureCtrl.assembly.pressure_g, p_w->data.size(), 0.0);
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

		////解析解case1 ：锐利前沿
  //      const double front_speed = 0.1;     ///< 饱和度前沿速度 [m/s]
  //      const double s_w_inj = 1.0;

        //解析解case2:平滑 tanh 前沿案例（applySaturationFront_SmoothTanh） 
        double front_speed = 1.0;      // [m/s]
        double s_w_inj = 1.0;      // 右侧纯水
        double w0 = 3.0;      // 起始半厚度 ~ 一个单元
        double w_sqrtcoef = 0.5;      // 随时间略微变宽

        ////解析解case3 :  Slug 案例（applySaturationFront_Slug）
        //double front_speed = 0.5;      // [m/s]
        //double s_w_inj = 1.0;      // slug 外部为纯水
        //double slug_length = 20.0;     // slug 长度 20 m


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

            ////2.1) case 1- 尖端前沿推进
            //if (!applySaturationFrontPattern(
            //    mgr, reg, satCfg,
            //    xmin, Lx,
            //    stepId, dt,
            //    front_speed,
            //    s_w_inj))
            //{
            //    std::cerr << "[IMPES][AnalyticTest] applySaturationFrontPattern failed at step "
            //        << stepId << ".\n";
            //    return false;
            //}


            // 2-2) case2- 光滑 tanh 型前沿
            if (!applySaturationFront_SmoothTanh(
                mgr, reg, satCfg,
                xmin, Lx,
                stepId, dt,
                front_speed,
                s_w_inj,
                w0,
                w_sqrtcoef))
            {
                std::cerr << "[IMPES][AnalyticTest] applySaturationFrontPattern failed at step "
                    << stepId << ".\n";
                return false;
            }

			//// 2-3) case3- Slug 前沿推进
   //         if (!applySaturationFront_SmoothTanh(
   //             mgr,
   //             reg,
   //             satCfg,
   //             xmin,
   //             Lx,
   //             stepId,
   //             dt,
   //             front_speed,  // [m/s]
   //             s_w_inj,      // 右侧水饱和度，一般 1.0
   //             w0,           // 初始前沿半厚度 [m]
   //             w_sqrtcoef    // 随 sqrt(t) 增长的厚度系数 [m/sqrt(s)]
   //         ))
   //         {
   //             std::cerr << "[IMPES][AnalyticTest] applySaturationFrontPattern failed at step "
   //                 << stepId << ".\n";
   //             return false;
   //         }


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
            TwoPhase::copyBasicPropertiesToOldLayer(reg);  //把基本物性参数拷贝到旧时间层
        }
        std::cout << "===== [IMPES][AnalyticTest] finished " << nSteps
            << " time steps. =====\n";
        return true;
    };
}

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
#include "SaturationTransportEqAssemblerandSolver.h"
#include "IMPES_PostProcessIO.h"
#include "TimeTermAssemblerandSolver.h"
#include "Diff_TPFA_GradientsOperation.h"
#include "PostProcess_.h"



namespace IMPES_Iteration
{
    struct PressureSolveControls_analytic
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
        const double Sw_res = satCfg.VG_Parameter.vg_params.Swr;  // 残余水饱和度

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
        const double Sw_res = satCfg.VG_Parameter.vg_params.Swr;

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
        const double Sw_res = satCfg.VG_Parameter.vg_params.Swr;

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
        const PressureSolveControls_analytic& pressureCtrl,
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
        TwoPhase::updateTwoPhasePropertiesAtTimeStep(mgr, reg, satCfg.saturation, satCfg.VG_Parameter.vg_params, satCfg.VG_Parameter.relperm_params);  //计算各相的有效饱和度、各相相对渗透率、毛细压力
        TwoPhase::updateWaterBasicPropertiesAtStep(mgr, reg, pressureCtrl.assembly.pressure_field, "T"); //基于相对渗透率更新流度以及其他基本物性参数

        ///额外计算CO2相的压力
        auto Pc = reg.get<volScalarField>(PhysicalProperties_string::TwoPhase_case().Pc_tag);
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
            TwoPhase::updateTwoPhasePropertiesAtTimeStep(mgr, reg, satCfg.saturation, satCfg.VG_Parameter.vg_params, satCfg.VG_Parameter.relperm_params);  //计算各相的有效饱和度、各相相对渗透率、毛细压力
            TwoPhase::updateWaterBasicPropertiesAtStep(mgr, reg, pressureCtrl.assembly.pressure_field, "T"); //基于相对渗透率更新流度以及其他基本物性参数
            ///额外计算CO2相的压力
            auto Pc = reg.get<volScalarField>(PhysicalProperties_string::TwoPhase_case().Pc_tag);
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


	//===================================测试重力离散系数的准确性========================================//
    /**
     * @brief 在给定网格和物性场上做“静水压 + 重力”算子测试。
     *
     * 测试思路：
     * - 先调用 computeDerivedMobilityFields 构造 rho_coeff_mass(=λ_mass) 与 rho_gravity_mass(=ρ_buoy)；
     * - 构造解析静水压场 p(x) = p_ref + ρ_buoy(x) * (g · x)，使得 ∇p = ρ_buoy g；
     * - 用 buildAndApplyOperator + assemble_COO 得到离散算子 L；
     * - 计算 y = L[p]，检查 max|y| / max|p| 是否足够小。
     *
     * @param mgr  网格管理器
     * @param reg  单元场仓库
     * @param freg 面场仓库
     * @param pbc  压力边界条件适配器（建议已按静水压解析解设置好 Dirichlet/Neumann）
     * @param cfg_in 压力装配配置（重力方向等可在外部设定，这里会复制一份）
     * @return true  测试成功执行（不代表残差为零，只代表流程正确跑完）
     * @return false 流程中出现错误（缺场或装配失败等）
     */
    inline bool runHydrostaticGravityOperatorTest(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& pbc,
        const PressureAssemblyConfig& cfg_in)
    {
        Mesh& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        // ---- 0) 准备配置：开启重力，若外部未给重力方向则默认竖直向下 ----
        PressureAssemblyConfig cfg = cfg_in;
        cfg.enable_buoyancy = true;
        if (cfg.gravity.Mag() <= 0.0) {
            // 这里假定 y 方向为竖直向上轴，重力向下
            cfg.gravity = Vector{ 0.0, -9.81, 0.0 };
        }

        // ---- 1) 计算 λ_mass / ρ_buoy 等派生场 ----
        if (!detailed::computeDerivedMobilityFields(mgr, reg, cfg)) {
            std::cerr << "[HydrostaticTest] computeDerivedMobilityFields failed.\n";
            return false;
        }

        auto rhoBuoyF = reg.get<volScalarField>(cfg.rho_gravity_field); // 这里存的是 ρ_buoy
        if (!rhoBuoyF) {
            std::cerr << "[HydrostaticTest] missing buoyancy field '"
                << cfg.rho_gravity_field << "'.\n";
            return false;
        }

        // ---- 2) 构造解析静水压场 p(x) = p_ref + ρ_buoy(x) * (g · x) ----
        //
        //     梯度：∇p = ρ_buoy g （g 为常向量） ⇒ (∇p - ρ_buoy g) = 0
        //     达西质量通量：q_m = -k λ_mass (∇p - ρ_buoy g) = 0
        //
        auto pF = reg.getOrCreate<volScalarField>(
            cfg.pressure_field, cells.size(), 0.0);

        const double p_ref = 1.0e7; // 10 MPa，任意常数，只影响绝对水平不影响通量
        const Vector& g = cfg.gravity;

        for (const auto& c : cells) {
            if (c.id < 0) continue; // 跳过 ghost 等
            const size_t i = id2idx.at(c.id);
            const double rhoB = (*rhoBuoyF)[i];  // ρ_buoy(x)
            const double gx = g * c.center;    // g · x  （假定 Vector::operator* 为点积）
            (*pF)[i] = p_ref + rhoB * gx;
        }

        // ---- 3) 构造 unknown 映射，并把单元场 p_w 收集到未知向量 p_vec 中 ----
        int nUnknowns = 0;
        std::vector<int> lid_of_cell = buildUnknownMap(mesh, nUnknowns);
        if (nUnknowns <= 0) {
            std::cerr << "[HydrostaticTest] no unknowns.\n";
            return false;
        }

        std::vector<double> p_vec;
        if (!detailed::gatherFieldToVector(
            mesh, reg,
            cfg.pressure_field,
            lid_of_cell, nUnknowns,
            p_vec))
        {
            std::cerr << "[HydrostaticTest] gatherFieldToVector failed.\n";
            return false;
        }

        // ---- 4) 用“扩散 + 重力”模板构建算子，并对 p_vec 应用 ----
        //
        // mobility_tokens 只给渗透率张量：kxx/kyy/kzz
        // λ_mass 完全通过 rho_coeff_field (= cfg.rho_coeff_field) 进入，
        // 重力项的 ρ_buoy 通过 rho_buoy_field (= cfg.rho_gravity_field) 进入，
        // 与 computeDerivedMobilityFields 的定义严格对应。
        //
        OperatorFieldNames nm = makeNames(cfg.operator_tag);
        std::vector<std::string> mobility_tokens = {
            "kxx:kxx",
            "kyy:kyy",
            "kzz:kzz"
        };

        std::vector<double> Lp; // 存 y = L[p]
        if (!detailed::buildAndApplyOperator(
            mgr, reg, freg,
            pbc,
            nm,
            mobility_tokens,
            cfg.rho_coeff_field,    // λ_mass 场名
            cfg.rho_gravity_field,  // ρ_buoy 场名
            cfg.pressure_field,     // 这里对 p_w 做 L[p]
            cfg,
            /*enable_buoyancy =*/ true,
            p_vec,
            Lp))
        {
            std::cerr << "[HydrostaticTest] buildAndApplyOperator failed.\n";
            return false;
        }

        // ---- 5) 计算残差范数，输出诊断信息 ----
        double maxRes = 0.0;
        double l2Res = 0.0;
        for (int i = 0; i < nUnknowns; ++i) {
            const double r = std::abs(Lp[i]);
            maxRes = std::max(maxRes, r);
            l2Res += Lp[i] * Lp[i];
        }
        l2Res = std::sqrt(l2Res / std::max(1, nUnknowns));

        double maxP = 0.0;
        for (double v : p_vec) {
            maxP = std::max(maxP, std::abs(v));
        }
        const double relMax = (maxP > 0.0) ? (maxRes / maxP) : 0.0;

        std::cout << "[HydrostaticTest] N = " << nUnknowns
            << ", max|L[p]| = " << maxRes
            << ", L2|L[p]| = " << l2Res
            << ", relMax = " << relMax << std::endl;

        // 可选：打印前几个未知的 p 和 L[p]，帮助 debug
        int printed = 0;
        for (size_t ic = 0; ic < cells.size() && printed < 5; ++ic) {
            const int lid = lid_of_cell[ic];
            if (lid < 0) continue;
            const auto& c = cells[ic];
            std::cout << "  cell id = " << c.id
                << ", p = " << p_vec[lid]
                << ", L[p] = " << Lp[lid]
                << std::endl;
            ++printed;
        }

        return true;
    }


    /// \brief 测试毛细扩散算子：解 ∇·(k λ_g ρ_g ∇p_c) = 0，解析解 p_c = p0 + α y
    /// 测试毛细扩散算子：L[Pc] = ∇·(k λ_g ρ_g ∇Pc)
/// 这里取 k=1, λ_g=1, ρ_g=1 => L[Pc] = ∇² Pc
    inline bool runCapillaryOperatorTest
    (
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& PbcA,
        const IMPES_Iteration::PressureAssemblyConfig& cfg
    )
    {
        using namespace IMPES_Iteration::detailed;

        Mesh& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        const int Nc = static_cast<int>(cells.size());

        // 1) 准备 mobility / ρ_coeff 场
        //    这里假定你已经提前把 lambda_w, lambda_g, rho_w, rho_g 填好了
        //    对于测试，可以直接填：lambda_w=0, lambda_g=1, rho_w=0, rho_g=1, kxx=kyy=1
        if (!computeDerivedMobilityFields(mgr, reg, cfg)) {
            std::cerr << "[CapillaryTest] computeDerivedMobilityFields failed.\n";
            return false;
        }

        // 2) 构造一个解析解：Pc(x,y) = P0 + alpha * y
        const double P0 = 2.0e5;   // 基础 capillary 压力 [Pa]
        const double alpha = 1.0e3;   // 梯度 [Pa/m]

        auto PcF = reg.getOrCreate<volScalarField>(cfg.Pc_field.c_str(), cells.size(), 0.0);
        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const Vector& ctr = c.center;
            const double y = ctr.m_y;
            (*PcF)[i] = P0 + alpha * y;
        }

        // 3) 设置与解析解一致的边界条件：
        //    - 左/右：Neumann (no flux) =>  a=0, b=1, c=0
        //    - 下：Dirichlet Pc = P0   =>  a=1, b=0, c=P0
        //    - 上：Dirichlet Pc = P0 + alpha * Ly  (Ly=domain size in y)
        //       这里假设你已经用同样的 BoxBCs 给 p_w 设过一次；
        //       我们只需保证对于 Pc 的 BC 也是同样形式即可。
        //    实现上：你可以单独建一个 PcBC::Registry；为了简单演示，这里直接假设
        //    当前 PbcA 里已经按上述方式配置好；否则你需要和 HydrostaticTest 一样写一遍 BC 注册。

        // 4) 构建 lid 映射（与压力/重力测试相同）
        std::vector<int> lid_of_cell(cells.size(), -1);
        int lid = 0;
        for (size_t ic = 0; ic < cells.size(); ++ic)
        {
            const auto& c = cells[ic];
            if (c.id < 0) continue;
            lid_of_cell[ic] = lid++;
        }
        const int N = lid;

        // 5) 从 Pc 场收集成向量 eval_vec
        std::vector<double> eval_vec, result_vec;
        if (!gatherFieldToVector(mesh, reg, cfg.Pc_field, lid_of_cell, N, eval_vec)) {
            std::cerr << "[CapillaryTest] gatherFieldToVector(Pc) failed.\n";
            return false;
        }

        // 6) 调用 build + assemble + applyCOO，构建 L[Pc]
        OperatorFieldNames nm = makeNames("Pc_capillary_test");

        // mobility_tokens: 这里只需要 kxx,kyy,kzz => 几何部分 k
        std::vector<std::string> mobility_tokens = {
            "kxx:kxx", "kyy:kyy", "kzz:kzz"
        };

        bool ok = buildAndApplyOperator(
            mgr, reg, freg,
            PbcA,
            nm,
            mobility_tokens,
            cfg.rho_capillary_field,   // ρ_coeff_field = λ_g * ρ_g
            "",                        // rho_buoy_field (不用重力，留空)
            cfg.Pc_field,              // x_field = Pc
            cfg,
            /*enable_buoyancy*/ false,
            eval_vec,
            result_vec
        );

        if (!ok) {
            std::cerr << "[CapillaryTest] buildAndApplyOperator failed.\n";
            return false;
        }

        // 7) 统计残差范数
        double maxAbs = 0.0;
        double l2 = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const double v = std::abs(result_vec[i]);
            maxAbs = std::max(maxAbs, v);
            l2 += v * v;
        }
        l2 = std::sqrt(l2 / std::max(N, 1));

        std::cout << "[CapillaryTest] N = " << N
            << ", max|L[Pc]| = " << maxAbs
            << ", L2|L[Pc]| = " << l2 << "\n";

        // 取一个样本单元输出
        int sample_cell_id = -1;
        double Pc_sample = 0.0, L_sample = 0.0;
        for (size_t ic = 0; ic < cells.size(); ++ic)
        {
            const auto& c = cells[ic];
            if (c.id < 0) continue;
            const int li = lid_of_cell[ic];
            if (li < 0) continue;

            sample_cell_id = c.id;
            Pc_sample = (*PcF)[ic];
            L_sample = result_vec[li];
            break;
        }

        std::cout << "  sample cell id = " << sample_cell_id
            << ", Pc = " << Pc_sample
            << ", L[Pc] = " << L_sample << "\n";

        // 期望：max|L[Pc]| 在 1e-10~1e-12 量级（取决于网格、非正交等）
        std::cout << "[CapillaryTest] finished. Expect residuals near machine precision if capillary operator is assembled correctly.\n";
        return true;
    }







}

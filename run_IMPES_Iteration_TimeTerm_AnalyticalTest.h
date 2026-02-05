#pragma once
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include   "PhysicalPropertiesManager_TwoPhase.h"
#include    "TimeLoopDriver.h" 

int run_IMPES_Iteration_TimeTerm_AnalyticalTest()
{
    // ---------------- 0. 网格与场注册 ----------------
    const double lengthX = 100.0;
    const double lengthY = 100.0;
    const double lengthZ = 0.0;
    const int sectionNumX = 50;
    const int sectionNumY = 50;
    const int sectionNumZ = 0;
    const bool usePrism = true;
    const bool useQuadBase = false;

    std::cout << "\n===== [TEST] IMPES two-phase time-term assembly =====\n";

    std::cout << "--- IMPES: building mesh (CO2 displacement) ---\n";
    MeshManager mgr(lengthX, lengthY, lengthZ,
        sectionNumX, sectionNumY, sectionNumZ,
        usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    FieldRegistry     reg;
    FaceFieldRegistry freg;

    // ---------- 1. 初始主变量场 ----------
    std::cout << "--- IMPES: initializing primary fields ---\n";

    // 填充初始场
        ///水相压力场，p_w
    InitFields ic_pw; //创建一个新的 InitFields 结构体用于水相压力场初始化
    ic_pw.x0 = 1e7;     // 初始水相压力 10 MPa
    ic_pw.x_dx = 0.0;   // x方向梯度为 0
    ic_pw.x_dy = 0.0;   // y方向梯度为 0
    ic_pw.x_dz = 0.0;   // z方向梯度为 0
    Initializer::fillBaseDistributions1(mgr.mesh(), reg, ic_pw, "p_w"); //填充水相压力场

    ///水相饱和度场，s_w
    InitFields ic_sw; //创建一个新的 InitFields 结构体用于水相饱和度场初始化
    ic_sw.x0 = 0.95;    // 初始水相饱和度 0.95
    ic_sw.x_dx = 0.0;   // x方向梯度为 0
    ic_sw.x_dy = 0.0;   // y方向梯度为 0
    ic_sw.x_dz = 0.0;   // z方向梯度为 0
    Initializer::fillBaseDistributions1(mgr.mesh(), reg, ic_sw, "s_w"); //填充水相饱和度场

    ///温度场，T
    InitFields ic_T; //创建一个新的 InitFields 结构体用于温度场初始化
    ic_T.x0 = 300.0;    // 初始温度 300 K
    ic_T.x_dx = 0.0;   // x方向梯度为 0
    ic_T.x_dy = 0.0;   // y方向梯度为 0
    ic_T.x_dz = 0.0;   // z方向梯度为 0
    Initializer::fillBaseDistributions1(mgr.mesh(), reg, ic_T, "T"); //填充温度场;

    // 时间相关辅助场 p_old/prev, s_old/prev, T_old/prev
    GeneralTools::ensureTransientFields_scalar1(mgr.mesh(), reg, "p_w", "p_w_old", "p_w_prev");
    GeneralTools::ensureTransientFields_scalar1(mgr.mesh(), reg, "s_w", "s_w_old", "s_w_prev");
    GeneralTools::ensureTransientFields_scalar1(mgr.mesh(), reg, "T", "T_old", "T_prev");

    // 区域分类（这里只用一个 RegionType::Medium）
    PhysicalPropertiesManager ppm;
    ppm.classifyRockRegionsByGeometry(mgr, {}, Cell::RegionType::Medium);
    ppm.UpdateRockProperties(mgr, reg, "p_w", "T");

    // VG / 相对渗透率参数
    VGParams vg_params;
    vg_params.Swr = 0.20;
    vg_params.Sgr = 0.00;
    vg_params.alpha = 2.0e-6;
    vg_params.n = 2.5;

    RelPermParams rp_params;
    rp_params.L = 0.5;

    const auto& bfaces = mgr.boundaryFaces();
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 1.0, 0.0, 6.5e6 };
    PressureBC::BoundaryCoefficient P_Right{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Down{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Up{ 0.0, 1.0, 0.0 };
    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);
    PressureBCAdapter PbcA{ pbc_pw };

    // ---------- 6. 时间项装配测试 ----------

    const int    nSteps = 1000;
    const double totalTime = 1000.0;
    const double dt = totalTime / nSteps;

    IMPES_Iteration::PressureSolveControls_analytic pCtrl;
    pCtrl.assembly.pressure_field = "p_w";
    pCtrl.assembly.pressure_old_field = "p_w_old";
    pCtrl.assembly.pressure_prev_field = "p_w_prev";


    IMPES_Iteration::SaturationTransportConfig satCfg;
    satCfg.saturation = "s_w";
    satCfg.saturation_old = "s_w_old";
    satCfg.saturation_prev = "s_w_prev";
    satCfg.VG_Parameter.vg_params = vg_params;
    satCfg.VG_Parameter.relperm_params = rp_params;

    const std::string outPrefixP = "./Postprocess_Data/IMPES_Iteration_Test/Case3/p_impes_ps_revised/p_ps";
    const std::string outPrefixSw = "./Postprocess_Data/IMPES_Iteration_Test/Case3/s_impes_ps_revised/s_ps";
    const std::string timeSeriesFn = "./Postprocess_Data/IMPES_Iteration_Test/Case3/time_series/impes_ps_revised_stats.csv";  // 如需写出时间序列可填路径

    if (!IMPES_Iteration::runTransient_IMPES_AnalyticTest(
        mgr, reg, freg, PbcA,
        nSteps, dt,
        pCtrl, satCfg,
        10, 10,
        outPrefixP, outPrefixSw,
        10, timeSeriesFn))
    {
        std::cerr << "[TEST] runTransient_IMPES_AnalyticTest failed.\n";
        return 1;
    }

    std::cout << "===== [TEST] IMPES two-phase time-term assembly done. =====\n";
    return 0;
}
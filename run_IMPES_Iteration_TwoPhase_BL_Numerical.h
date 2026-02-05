#pragma once
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>

#include    "PhysicalPropertiesManager_TwoPhase.h"
#include    "IMPES_Iteration_Loop.h"
/**
 * @brief 1D Buckley–Leverett 风格的两相 IMPES 数值测试（仅数值解，不构造解析解）。
 *
 * 设定：
 *  - 1D 区域：x ∈ [0, Lx]，y 方向只有 1 层单元；
 *  - 左端高压、右端低压，形成近似常总通量的两相驱替；
 *  - 初场为典型 Riemann 型：x≈0 处 S_w ≈ S_w_inj，其余区域 S_w ≈ S_w_i；
 *  - VG 模型设置为 Swr = Sgr = 0，且 alpha 较大 → Pc 数值很小，对应经典 Buckley–Leverett 假设。
 *
 * 流程：
 *  1) 构建 1D 网格 + 主变量场 (p_w, s_w, T)；
 *  2) 基于 p_w, T 更新基岩物性；
 *  3) 施加 BL 风格饱和度初场；
 *  4) 设置 VG/Mualem 参数、压力边界条件与 IMPES 控制参数；
 *  5) 调用 IMPES_Iteration::runTransient_IMPES_Iteration 做主时间步推进；
 *  6) 输出 Tecplot 结果用于后处理。
 */
int run_IMPES_Iteration_TwoPhase_BL_Numerical()
{
    // ---------------- 0. 网格与场注册 ----------------
    const double lengthX = 1.0;   // [m]
    const double lengthY = 1.0;    // 几乎 1D
    const double lengthZ = 0.0;

    const int sectionNumX = 50;    // 细一点，方便看前沿
    const int sectionNumY = 50;
    const int sectionNumZ = 0;

    const bool usePrism = true;
    const bool useQuadBase = false;

    std::cout << "\n===== [TEST] IMPES two-phase (1D Buckley–Leverett, numerical only) =====\n";
    std::cout << "--- IMPES: building mesh (1D BL test) ---\n";

    MeshManager mgr(lengthX, lengthY, lengthZ,
        sectionNumX, sectionNumY, sectionNumZ,
        usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OverRelaxed);

    FieldRegistry     reg;
    FaceFieldRegistry freg;

    // ---------- 1. 初始主变量场 ----------
    std::cout << "--- IMPES: initializing primary fields ---\n";

    // 4. VG / 相对渗透率参数：进一步减小 Pc
    VGParams vg_params;
    vg_params.Swr = 0.1;  // 少量束缚水
    vg_params.Sgr = 0.15;  // 适当残余气，避免气相端点过头
    vg_params.alpha = 1e-4;   // 保持 Pc 小
    vg_params.n = 2.0;   // 曲线更平缓，前沿不至于过陡

    RelPermParams rp_params;
    rp_params.L = 0.5;       // Mualem 默认，减缓相对渗透率陡峭度

    double P_left = 9.0e6;
    double P_right = 8e6;
    
    InitFields ic;  // 默认 p0 / T0 / Sw0，如有需要可以修改 ic.p0 / ic.T0 等
    ic.p_w0 = 8e6;                        // x=0 处的压力
    ic.dp_wdx = 0.0;  // 与 lengthX 匹配的梯度
    ic.dp_wdy = 0.0;
    ic.dp_wdz = 0.0;
    
    ic.s_w = 0.2;   // Sw_i
    // ------------ 2.变量字符设置 ------------
    IMPES_Iteration::PressureEquation_String P_Eq;
    IMPES_Iteration::SaturationEquation_String S_Eq;
    IMPES_Iteration::FaceMassRate_String Mf_Eq;
    IMPES_Iteration::FluxSplitConfig_String Fsp_Eq;

    // 填充初始场
        ///水相压力场，p_w
    InitFields ic_pw; //创建一个新的 InitFields 结构体用于水相压力场初始化
    ic_pw.x0 = 1e7;     // 初始水相压力 10 MPa
    ic_pw.x_dx = 0.0;   // x方向梯度为 0
    ic_pw.x_dy = 0.0;   // y方向梯度为 0
    ic_pw.x_dz = 0.0;   // z方向梯度为 0
    Initializer::fillBaseDistributions1(mgr.mesh(), reg, ic_pw, P_Eq.pressure_field); //填充水相压力场

    ///水相饱和度场，s_w
    InitFields ic_sw; //创建一个新的 InitFields 结构体用于水相饱和度场初始化
    ic_sw.x0 = 0.95;    // 初始水相饱和度 0.95
    ic_sw.x_dx = 0.0;   // x方向梯度为 0
    ic_sw.x_dy = 0.0;   // y方向梯度为 0
    ic_sw.x_dz = 0.0;   // z方向梯度为 0
    Initializer::fillBaseDistributions1(mgr.mesh(), reg, ic_sw, S_Eq.saturation); //填充水相饱和度场

    ///温度场，T
    InitFields ic_T; //创建一个新的 InitFields 结构体用于温度场初始化
    ic_T.x0 = 300.0;    // 初始温度 300 K
    ic_T.x_dx = 0.0;   // x方向梯度为 0
    ic_T.x_dy = 0.0;   // y方向梯度为 0
    ic_T.x_dz = 0.0;   // z方向梯度为 0
    Initializer::fillBaseDistributions1(mgr.mesh(), reg, ic_T, "T"); //填充温度场

    // 时间层：old / prev（p_w、s_w、p_g 都准备好）
    GeneralTools::ensureTransientFields_scalar1(mgr.mesh(), reg, P_Eq.pressure_field, P_Eq.pressure_old_field, P_Eq.pressure_prev_field);
    GeneralTools::ensureTransientFields_scalar1(mgr.mesh(), reg, S_Eq.saturation, S_Eq.saturation_old, S_Eq.saturation_prev);
    GeneralTools::ensureTransientFields_scalar1(mgr.mesh(), reg, "p_g", "p_g_old", "p_g_prev");

    // ---------- 2. 基岩区域分类与固相物性 ----------
    PhysicalPropertiesManager ppm;
    ppm.classifyRockRegionsByGeometry(mgr, {}, Cell::RegionType::Medium);
    ppm.UpdateRockProperties(mgr, reg, P_Eq.pressure_field, "T");

    // ---------- 5. 压力边界条件（左高右低，其他 no-flow） ----------
    const auto& bfaces = mgr.boundaryFaces();
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 1.0, 0.0, P_left };
    PressureBC::BoundaryCoefficient P_Right{ 1.0, 0.0, P_right };
    PressureBC::BoundaryCoefficient P_Down{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Up{ 0.0, 1.0, 0.0 };

    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);
    PressureBCAdapter PbcA{ pbc_pw };

    // ---------- 6. 井配置（BL 测试：无井源） ----------
    std::vector<WellConfig_TwoPhase> wells_cfg;
    build_masks_and_WI_for_all(mgr, reg, wells_cfg);
    const int Ncells = static_cast<int>(mgr.mesh().getCells().size());
    std::vector<WellDOF_TwoPhase> wells_dof;
    register_well_dofs_for_all_TwoPhase(Ncells, wells_cfg, wells_dof);

    // ---------- 7. 饱和度方程配置 ----------
    IMPES_Iteration::SaturationTransportConfig satCfg;
    satCfg.VG_Parameter.vg_params = vg_params;
    satCfg.VG_Parameter.relperm_params = rp_params;
    satCfg.dS_max = 0.1;
    satCfg.CFL_safety = 0.8;

    // ---------- 8. 压力方程组装与求解控制参数 ----------   
    IMPES_Iteration::PressureSolveControls pCtrl;
    //组装
    pCtrl.assembly.VG_Parameter.vg_params= vg_params;
    pCtrl.assembly.VG_Parameter.relperm_params = rp_params;
    pCtrl.assembly.enable_buoyancy = false;
    pCtrl.assembly.gradient_smoothing = 6;
    pCtrl.assembly.gravity = Vector{ 0.0, 0.0, 0.0 };

    pCtrl.max_outer = 1000;        // 常物性可保持很小
    pCtrl.tol_abs = 1e8;
    pCtrl.tol_rel = 1e-4;
    pCtrl.under_relax = 1;
    pCtrl.verbose = true;

    // ----------- 9. 面通量计算配置---------------
    IMPES_Iteration::FaceMassRateConfig m_FCtrl;
    m_FCtrl.clamp_dirichlet_backflow = true;
    m_FCtrl.dirichlet_zero_flux_tol = 1e-10;
    m_FCtrl.enable_conservation_check = true;
    m_FCtrl.flux_check_tol = 1e-10;
    
    // ---------- 9. 两相通量拆分配置 ----------
    IMPES_Iteration::FluxSplitConfig fluxCfg;
    fluxCfg.boundary_inflow_fw = 1;
    fluxCfg.enforce_boundary_inflow_fw = true;
    fluxCfg.flux_sign_epsilon = 1e-10;
    fluxCfg.pressure_bc = &PbcA;

    // ---------- 11. IMPES 主时间推进 ----------
    const int    nSteps = 500;
    double       dt_initial = 1e-7;   // 更保守的初始时间步

    const int writeEveryP = 1;
    const int writeEverySw = 1;

    const std::string outPrefixP = "./Postprocess_Data/IMPES_Iteration_Test/Case7/p_impes_ps_revised/p_ps";
    const std::string outPrefixSw = "./Postprocess_Data/IMPES_Iteration_Test/Case7/s_impes_ps_revised/s_ps";
    const int snapshotEveryCsv = 50;
    const std::string snapshotPrefix = "./Postprocess_Data/csv_snapshots/Case7/ps_state_reviesed";

    std::cout << "--- IMPES: start transient run (BL numerical test) ---\n";
    const bool ok = IMPES_Iteration::runTransient_IMPES_Iteration(mgr, reg, freg, PbcA, wells_dof, nSteps, dt_initial, pCtrl, satCfg, fluxCfg, m_FCtrl, writeEveryP, writeEverySw, outPrefixP, outPrefixSw, snapshotEveryCsv, snapshotPrefix);
    if (!ok)
    {
        std::cerr << "[BL-TEST] IMPES transient run failed.\n";
        return EXIT_FAILURE;
    }

    std::cout << "[BL-TEST] IMPES two-phase Buckley–Leverett numerical test finished successfully.\n";
    return EXIT_SUCCESS;
}

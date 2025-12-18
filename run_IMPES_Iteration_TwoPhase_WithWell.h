#pragma once
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "PhysicalPropertiesManager_TwoPhase.h"
#include "IMPES_Iteration_Loop.h"
#include "WellConfig_TwoPhase.h"  
/**
 * @brief Two-phase IMPES test with Peaceman wells (1 injector + 1 producer).
 *
 * 工况：
 *  - 域：1 m × 1 m，厚度 H = 1 m
 *  - 网格：50 × 50
 *  - 边界：上下左右 全部无流（Neumann，∂p/∂n = 0）
 *  - 初始：p_w = 9 MPa，s_w = 0.2
 *          VG: Swr = 0.1, Sgr = 0.15
 *  - 井：
 *      注入井 INJ：pos = (0.25, 0.25)，rw=0.001, H=1.0, maxHitCells=1
 *                  Mode = Rate, target = 5 kg/s，总质量注入为纯水（s_w_bh=1）
 *      采出井 PROD：pos = (0.75, 0.75)，rw=0.001, H=1.0, maxHitCells=1
 *                  Mode = Pressure, target = 8 MPa
 */

int run_IMPES_Iteration_TwoPhase_WellCase()
{
    // ---------------- 0. 网格设置与场注册 ----------------
    const double lengthX = 1.0;   // [m]
    const double lengthY = 1.0;    // 几乎 1D
    const double lengthZ = 0.0;

    const int sectionNumX = 50;    // 细一点，方便看前沿
    const int sectionNumY = 50;
    const int sectionNumZ = 0;

    const bool usePrism = true;
    const bool useQuadBase = false;
    std::cout << "\n ==========Mesh: BUILDING ==========\n";
    MeshManager mgr(lengthX, lengthY, lengthZ,sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OverRelaxed);
	std::cout << " \n ==========Mesh :BUILDING COMPLETED ==========\n";




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

    // 初始场
    InitFields ic;
    ic.p_w0 = 7e6;      // 初始水相压力 7 MPa
    ic.dp_wdx = 0.0;
    ic.dp_wdy = 0.0;
    ic.dp_wdz = 0.0;
    ic.s_w = 0.2;      // 初始水相饱和度

    IMPES_Iteration::PressureEquation_String P_Eq;
    IMPES_Iteration::SaturationEquation_String S_Eq;
    IMPES_Iteration::FaceMassRate_String Mf_Eq;
    IMPES_Iteration::FluxSplitConfig_String Fsp_Eq;

    // 主变量场：水相压力 p_w、水相饱和度 s_w、温度 T
    Initializer::createPrimaryFields_TwoPhase_HT_IMPES(mgr.mesh(), reg, P_Eq.pressure_field, S_Eq.saturation, "T");
    Initializer::fillBaseDistributions_TwoPhase_HT_IMPES(mgr.mesh(), reg, ic);

    // 时间层：old / prev（p_w、s_w、p_g 都准备好）
    ensureTransientFields_scalar(mgr.mesh(), reg, P_Eq.pressure_field, P_Eq.pressure_old_field, P_Eq.pressure_prev_field);
    ensureTransientFields_scalar(mgr.mesh(), reg, S_Eq.saturation, S_Eq.saturation_old, S_Eq.saturation_prev);
    ensureTransientFields_scalar(mgr.mesh(), reg, "p_g", "p_g_old", "p_g_prev");

    // ---------- 2. 基岩区域分类与固相物性 ----------
    PhysicalPropertiesManager ppm;
    ppm.classifyRockRegionsByGeometry(mgr, {}, Cell::RegionType::Medium);
    ppm.UpdateRockProperties(mgr, reg, P_Eq.pressure_field, "T");

    // ---------- 3. 压力边界条件（左高右低，其他 no-flow） ----------
    const auto& bfaces = mgr.boundaryFaces();
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 0.0, 1.0, 0.0 }; // ∂p/∂n = 0
    PressureBC::BoundaryCoefficient P_Right{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Down{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Up{ 0.0, 1.0, 0.0 };

    PressureBC::setBoxBCs2D(pbc_pw, bfaces,
        P_Left, P_Right, P_Down, P_Up);
    PressureBCAdapter PbcA{ pbc_pw };

    // ---------- 4. 井配置（BL 测试：无井源） ----------
    std::vector<WellConfig_TwoPhase> wells_cfg;
    // 4.1 注入井 INJ
    {
        WellConfig_TwoPhase inj;
        inj.name = "INJ";
        inj.role = WellDOF_TwoPhase::Role::Injector;
        inj.mode = WellDOF_TwoPhase::Mode::Rate;    // 质量流量控制
        inj.target = 0.03;                           // 5 kg/s 总质量注入

        inj.Tin = 300.0;                            // 注入温度，随意设置一个值，占位
        inj.s_w_bh = 1.0;                           // 纯水注入
		inj.mu_w_inj = 1e-4;                       // 注入水相粘度
		inj.rho_w_inj = 1000.0;                     // 注入水相密度

        inj.geom.pos = Vector{ 0.25 * lengthX, 0.25 * lengthY, 0.0 };
        inj.geom.rw = 0.001;                        // 井筒半径
        inj.geom.skin = 0.0;
        inj.geom.H = 1.0;                           // 有效厚度
        inj.geom.perfRadius = 0.0;                  // 0 => 用最近 maxHitCells 个单元
        inj.geom.maxHitCells = 1;                   // 只完井最近的 1 个单元
        
        wells_cfg.push_back(inj);
    }

    //4.2 采出井
    {
        WellConfig_TwoPhase prod;
        prod.name = "PROD";
        prod.role = WellDOF_TwoPhase::Role::Producer;
        prod.mode = WellDOF_TwoPhase::Mode::Pressure; // BHP 控制
        prod.target = 5e6;                             // 8 MPa

        prod.Tin = 300.0;  // 对产井无实际意义，仅占位
        prod.s_w_bh = 1.0;    // 占位

        prod.geom.pos = Vector{ 0.43 * lengthX, 0.43 * lengthY, 0.0 };
        prod.geom.rw = 0.01;
        prod.geom.skin = 0.0;
        prod.geom.H = 1.0;
        prod.geom.perfRadius = 0.0;
        prod.geom.maxHitCells = 5;

        wells_cfg.push_back(prod);
    }

    // 4.3 基于几何构建完井掩码 + WI
    build_masks_and_WI_for_all(mgr, reg, wells_cfg);

    const int Ncells = static_cast<int>(mgr.mesh().getCells().size());
    std::vector<WellDOF_TwoPhase> wells_dof;
    register_well_dofs_for_all_TwoPhase(Ncells, wells_cfg, wells_dof);

    // ---------- 7. 饱和度方程配置 ----------
    IMPES_Iteration::SaturationTransportConfig satCfg;
    satCfg.VG_Parameter.vg_params = vg_params;
    satCfg.VG_Parameter.relperm_params = rp_params;
    satCfg.dS_max = 0.15;
    satCfg.CFL_safety = 0.9;
    satCfg.time_control_scheme = IMPES_Iteration::SatTimeControlScheme::SimpleCFL;
    satCfg.time_integration_scheme = IMPES_Iteration::SatTimeIntegrationScheme::ExplicitEuler;

    // ---------- 8. 压力方程组装与求解控制参数 ----------   
    IMPES_Iteration::PressureSolveControls pCtrl;
    //组装
    pCtrl.assembly.VG_Parameter.vg_params = vg_params;
    pCtrl.assembly.VG_Parameter.relperm_params = rp_params;
    pCtrl.assembly.enable_buoyancy = false;
    pCtrl.assembly.gradient_smoothing = 6;
    pCtrl.assembly.gravity = Vector{ 0.0, 0.0, 0.0 };

    pCtrl.max_outer = 1000;        // 常物性可保持很小
    pCtrl.tol_abs = 1e6;
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
    double       dt_initial = 1e-5; 

	IMPES_Iteration::TimeStepControl timeCtrl;
    timeCtrl.dt_min = 1e-5;
	timeCtrl.dt_max = 1.0;
    timeCtrl.grow_factor = 100;

    const int writeEveryP = 10;
    const int writeEverySw = 10;

    const std::string outPrefixP = "./Postprocess_Data/IMPES_Iteration_Test/IMPES_ExplicitEuler/p_impes_ps_withwell/p_ps";
    const std::string outPrefixSw = "./Postprocess_Data/IMPES_Iteration_Test/IMPES_ExplicitEuler/s_impes_ps_withwell/s_ps";
    const int snapshotEveryCsv = 10;
    const std::string snapshotPrefix = "./Postprocess_Data/csv_snapshots/IMPES_ExplicitEuler/ps_state_withwell_case2";

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

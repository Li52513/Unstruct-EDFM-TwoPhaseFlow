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
    const double lengthX = 100.0;   // [m]
    const double lengthY = 100.0;    // 几乎 1D
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
    vg_params.Sgr = 0.05;  // 适当残余气，避免气相端点过头
    vg_params.alpha = 1e-4;   // 保持 Pc 小
    vg_params.n = 2.0;   // 曲线更平缓，前沿不至于过陡
    RelPermParams rp_params;
    rp_params.L = 0.5;       // Mualem 默认，减缓相对渗透率陡峭度

    // 初始场
    InitFields ic;
    ic.p_w0 = 1e7;      // 初始水相压力 7 MPa
    ic.dp_wdx = 0.0;
    ic.dp_wdy = 0.0;
    ic.dp_wdz = 0.0;
    ic.s_w = 0.95;      // 初始水相饱和度

    IMPES_Iteration::PressureEquation_String P_Eq;
    IMPES_Iteration::SaturationEquation_String S_Eq;
    IMPES_Iteration::FaceMassRate_String Mf_Eq;
    IMPES_Iteration::FluxSplitConfig_String Fsp_Eq;

    // 主变量场：水相压力 p_w
    Initializer::createPrimaryFields(mgr.mesh(), reg, P_Eq.pressure_field);
	// 主变量场：水相饱和度 s_w
	Initializer::createPrimaryFields(mgr.mesh(), reg, S_Eq.saturation);
	// 主变量场：温度 T
    Initializer::createPrimaryFields(mgr.mesh(), reg, "T");

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
	GeneralTools::ensureTransientFields_scalar1(mgr.mesh(), reg, P_Eq.pressure_g, P_Eq.pressure_g_old, P_Eq.pressure_g_prev);


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

    // ---------- 4. 井配置 ----------
    std::vector<WellConfig_TwoPhase> wells_cfg;
    // 4.1 注入井 INJ
    {
        WellConfig_TwoPhase inj;
        inj.name = "INJ";
        inj.role = WellDOF_TwoPhase::Role::Injector;
        inj.mode = WellDOF_TwoPhase::Mode::Rate;
        inj.target = 10; // kg/s total mass rate
        inj.Tin = 300.0;
        inj.s_w_bh = 0.05; // pure water injection
        inj.geom.pos = Vector{ 0.0, 0.0, 0.0 };
        inj.geom.rw = 0.1;
        inj.geom.skin = 0.0;
        inj.geom.H = 10;
        inj.geom.perfRadius = 0.0;
        inj.geom.maxHitCells = 1;
        wells_cfg.push_back(inj);
    }

    //4.2 采出井
    {
        WellConfig_TwoPhase prod;
        prod.name = "PROD";
        prod.role = WellDOF_TwoPhase::Role::Producer;
        prod.mode = WellDOF_TwoPhase::Mode::Pressure;
        prod.target = 8e6;
        prod.Tin = 300.0;
        prod.s_w_bh = 1.0;
        prod.geom.pos = Vector{  lengthX, lengthY, 0.0 };
        prod.geom.rw = 0.1;
        prod.geom.skin = 0.0;
        prod.geom.H = 10;
        prod.geom.perfRadius = 0.0;
        prod.geom.maxHitCells = 1;
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
	satCfg.useTraditionMethod = false;
    satCfg.primary_phase = IMPES_Iteration::SaturationPrimaryPhase::Water;

    // ---------- 8. 压力方程组装与求解控制参数 ----------   
    IMPES_Iteration::PressureSolveControls pCtrl;
    //组装
    pCtrl.assembly.VG_Parameter.vg_params = vg_params;
    pCtrl.assembly.VG_Parameter.relperm_params = rp_params;
    pCtrl.assembly.enable_buoyancy = false;
    pCtrl.assembly.gradient_smoothing = 6;
    pCtrl.assembly.gravity = Vector{ 0.0, 0.0, 0.0 };

    pCtrl.max_outer = 1000;        // 常物性可保持很小
    pCtrl.tol_abs = 1e2;
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
    fluxCfg.use_phase_potential_upwind = true;
	fluxCfg.gravity = Vector{ 0.0, 0.0, 0.0 };

    // ---------- 11. IMPES 主时间推进 ----------
    const int    nSteps = 5000;
    double       dt_initial = 1e-5; 

	IMPES_Iteration::TimeStepControl timeCtrl;
    timeCtrl.dt_min = 1e-5;
	timeCtrl.dt_max = 1.0;
    timeCtrl.grow_factor = 100;


    const int writeEveryP = 10;
    const int writeEverySw = 10;

    const std::string outPrefixP = "./Postprocess_Data/IMPES_Iteration_Test/IMPES_ExplicitEuler_withoutGravity_withoutPotentila_Upwind/p_impes_ps_withwell_case2/p_ps";
    const std::string outPrefixSw = "./Postprocess_Data/IMPES_Iteration_Test/IMPES_ExplicitEuler_withoutGravity_withoutPotentila_Upwind/s_impes_ps_withwell_case2/s_ps";
    const int snapshotEveryCsv = 10;
    const std::string snapshotPrefix = "./Postprocess_Data/csv_snapshots/IMPES_ExplicitEuler_withoutGravity_withoutPotentila_Upwind_case2/ps_state_withwell";

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

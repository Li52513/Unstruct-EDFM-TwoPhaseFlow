#pragma once
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>

#include    "0_PhysicalParametesCalculateandUpdata.h"
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
    const double lengthX = 50.0;   // [m]
    const double lengthY = 50.0;    // 几乎 1D
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
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    FieldRegistry     reg;
    FaceFieldRegistry freg;

    // ---------- 1. 初始主变量场 ----------
    std::cout << "--- IMPES: initializing primary fields ---\n";

    // 4. VG / 相对渗透率参数：进一步减小 Pc
    VGParams vg_params;
    vg_params.Swr = 0.2;  // 少量束缚水
    vg_params.Sgr = 0.25;  // 适当残余气，避免气相端点过头
    vg_params.alpha = 1e-4;   // 保持 Pc 小
    vg_params.n = 2.0;   // 曲线更平缓，前沿不至于过陡

    RelPermParams rp_params;
    rp_params.L = 0.5;       // Mualem 默认，减缓相对渗透率陡峭度

    double P_left = 9.0e6;
    double P_right = 8e6;
    
    InitFields ic;  // 默认 p0 / T0 / Sw0，如有需要可以修改 ic.p0 / ic.T0 等
    ic.p_w0 = P_left;                        // x=0 处的压力
    ic.dp_wdx = (P_right - P_left) / lengthX;  // 与 lengthX 匹配的梯度
    ic.dp_wdy = 0.0;
    ic.dp_wdz = 0.0;
    
    ic.s_w = 1-vg_params.Sgr-0.2;   // Sw_i

    // 主变量场：水相压力 p_w、水相饱和度 s_w、温度 T
    Initializer::createPrimaryFields_TwoPhase_HT_IMPES(
        mgr.mesh(), reg, "p_w", "s_w", "T");
    Initializer::fillBaseDistributions_TwoPhase_HT_IMPES(
        mgr.mesh(), reg, ic);

    // 时间层：old / prev（p_w、s_w、p_g 都准备好）
    ensureTransientFields_scalar(mgr.mesh(), reg, "p_w", "p_w_old", "p_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "s_w", "s_w_old", "s_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "p_g", "p_g_old", "p_g_prev");

    // ---------- 2. 基岩区域分类与固相物性 ----------
    PhysicalPropertiesManager ppm;
    ppm.classifyRockRegionsByGeometry(mgr, {}, Cell::RegionType::Medium);
    ppm.UpdateMatrixRockAt(mgr, reg, "p_w", "T");

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

    // ---------- 7. 饱和度输运配置 ----------
    IMPES_Iteration::SaturationTransportConfig satCfg;
    satCfg.saturation = "s_w";
    satCfg.saturation_old = "s_w_old";
    satCfg.saturation_prev = "s_w_prev";
    satCfg.VG_Parameter.vg_params = vg_params;
    satCfg.VG_Parameter.relperm_params = rp_params;

    // ---------- 8. 压力求解控制参数 ----------
    IMPES_Iteration::PressureSolveControls pCtrl;
    pCtrl.max_outer = 1;        // 常物性可保持很小
    pCtrl.tol_abs = 1e15;
    pCtrl.tol_rel = 1e-4;
    pCtrl.under_relax = 1;
    pCtrl.verbose = true;

    pCtrl.assembly.enable_buoyancy = false;
    pCtrl.assembly.gradient_smoothing = 1;
    pCtrl.assembly.gravity = Vector{ 0.0, 0.0, 0.0 };
    pCtrl.assembly.dirichlet_zero_flux_tol = 1.0; // Pa tolerance to treat boundary flux as zero

    // ---------- 9. 两相通量拆分配置 ----------
    IMPES_Iteration::FluxSplitConfig fluxCfg;
    fluxCfg.rho_water = TwoPhase::Water().rho_tag; // "rho_w"
    fluxCfg.rho_gas = TwoPhase::CO2().rho_tag;   // "rho_g"
    fluxCfg.pressure_bc = &PbcA;

    // ---------- 10. 初场后处理：同步 Pc / p_g / 物性到 old 层 ----------
    {
        // 更新 Pc, krw, krg, dPc/dSw
        TwoPhase::updateTwoPhasePropertiesAtTimeStep(mgr, reg, "s_w", vg_params, rp_params);

        // 设置 p_g = p_w + Pc，并同步到 old/prev
        auto p_w = reg.get<volScalarField>("p_w");
        auto p_g = reg.get<volScalarField>("p_g");
        auto p_w_old = reg.get<volScalarField>("p_w_old");
        auto p_w_prev = reg.get<volScalarField>("p_w_prev");
        auto p_g_old = reg.get<volScalarField>("p_g_old");
        auto p_g_prev = reg.get<volScalarField>("p_g_prev");
        auto Pc = reg.get<volScalarField>(TwoPhase::Auxiliaryparameters().Pc_tag);
        const auto& cells = mgr.mesh().getCells();
        const auto& id2idx = mgr.mesh().getCellId2Index();
        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const double pg_val = (*p_w)[i] + (*Pc)[i];
            (*p_g)[i] = pg_val;
            (*p_g_old)[i] = pg_val;
            (*p_g_prev)[i] = pg_val;
            // 确保 p_w_old/p_w_prev 也和 p_w 同步
            (*p_w_old)[i] = (*p_w)[i];
            (*p_w_prev)[i] = (*p_w)[i];
        }

        // 更新水/CO2 物性并拷贝到 *_old 层
        TwoPhase::updateWaterBasicPropertiesAtStep(mgr, reg, "p_w", "T");
        TwoPhase::updateCO2BasicPropertiesAtStep(mgr, reg, "p_g", "T");
        TwoPhase::copyBasicPropertiesToOldLayer(reg);
    }

    // ---------- 11. IMPES 主时间推进 ----------
    const int    nSteps = 200;
    double       dt_initial = 1e-4;   // 更保守的初始时间步

    const int writeEveryP = 1;
    const int writeEverySw = 1;

    const std::string outPrefixP = "./Postprocess_Data/IMPES_Iteration_Test/Case4/p_impes_ps_revised/p_ps";
    const std::string outPrefixSw = "./Postprocess_Data/IMPES_Iteration_Test/Case4/s_impes_ps_revised/s_ps";
    const int snapshotEveryCsv = 1;
    const std::string snapshotPrefix = "./Postprocess_Data/csv_snapshots/Case4/ps_state_reviesed";
    std::vector<std::string> snapshotFields;

    std::cout << "--- IMPES: start transient run (BL numerical test) ---\n";

    const bool ok = IMPES_Iteration::runTransient_IMPES_Iteration
    (
        mgr,
        reg,
        freg,
        PbcA,
        wells_dof,
        nSteps,
        dt_initial,
        pCtrl,
        satCfg,
        fluxCfg,
        writeEveryP,
        writeEverySw,
        outPrefixP,
        outPrefixSw,
        snapshotEveryCsv,
        snapshotPrefix,
        snapshotFields
    );

    if (!ok)
    {
        std::cerr << "[BL-TEST] IMPES transient run failed.\n";
        return EXIT_FAILURE;
    }

    std::cout << "[BL-TEST] IMPES two-phase Buckley–Leverett numerical test finished successfully.\n";
    return EXIT_SUCCESS;
}

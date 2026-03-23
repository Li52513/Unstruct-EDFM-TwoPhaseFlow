/**
 * @file Test_Day6_TransientSolver.h
 * @brief Day6 全隐式(FIM)牛顿瞬态求解器测试入口
 */

#pragma once

namespace Test_Day6 {

    /** @brief 运行 2D 单相 注入/生产 瞬态测试 (>= 50步) */
    void Run_Day6_Transient_2D_SP_InjProd();

    /** @brief 运行 2D 两相 注入/生产 瞬态测试 (>= 50步) */
    void Run_Day6_Transient_2D_TP_InjProd();

    /** @brief 运行 2D 两相 多井干扰 瞬态测试 (>= 50步) */
    void Run_Day6_Transient_2D_TP_Multiwell();

    /** @brief 运行 3D 单相 注入/生产 瞬态测试 (>= 50步) */
    void Run_Day6_Transient_3D_SP_InjProd();

    /** @brief 运行 3D 两相 注入/生产 瞬态测试 (>= 50步) */
    void Run_Day6_Transient_3D_TP_InjProd();

    /** @brief 运行 3D 两相 多井干扰 瞬态测试 (>= 50步) */
    void Run_Day6_Transient_3D_TP_Multiwell();

    /** @brief 运行 2D 矩阵审计测试，检查 NNC/FF 系数装配是否正常 */
    void Run_Day6_MatrixAudit_2D_EDFM();

    /** @brief 运行 3D 矩阵审计测试，检查 NNC/FF 系数装配是否正常 */
    void Run_Day6_MatrixAudit_3D_EDFM();

    /** @brief Day6 campaign (2D): 16 physical scenarios × F0/F1/F2 topology axis with fail-fast gating */
    void Run_Day6_Campaign_2D_All_TxxFy();

    /** @brief Day6 campaign (3D mirror): 16 physical scenarios × F0/F1/F2 topology axis with fail-fast gating */
    void Run_Day6_Campaign_3D_All_TxxFy();

    /** @brief Immediate T1 split gates (run after T1_F0): no-fracture/single-fracture/cross-fracture */
    void Run_Day6_Campaign_2D_T01_F0();
    void Run_Day6_Campaign_2D_T01_F1();
    void Run_Day6_Campaign_2D_T01_F2();

    /** @brief T1 baseline: 2D no-well single-phase pressure diffusion (constant fluid + constant rock) with Fourier analytical check */
    void Run_Day6_T1_2D_SP_NoWell_Analytical();

} // namespace Test_Day6
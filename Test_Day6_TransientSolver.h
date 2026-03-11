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

} // namespace Test_Day6
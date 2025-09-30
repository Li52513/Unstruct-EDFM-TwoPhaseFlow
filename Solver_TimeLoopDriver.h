#pragma once
#include <functional>
#include <iostream>
#include "Solver_TimeLoopSkeleton.h"  

// 可选的导出回调
using WriteCallback = std::function<void(int step, double time)>;

/**
 * @brief CO2 单相 渗流-传热：多步时间推进驱动
 * @param nSteps     时间步数
 * @param dt         时间步长
 * @param ctrl       求解控制参数（外迭代/欠松弛/Jacobi 迭代设置等）
 * @param writeEvery 每多少步调用一次回调（<=0 表示不回调）
 * @param onWrite    回调：用于导出场（可为空）
 */

inline bool runTransient_CO2_singlePhase
(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    PhysicalPropertiesManager& ppm,
    const PressureBCAdapter& Pbc,
    const TemperatureBCAdapter& Tbc,
    const GravUpwind& gu,
    const RockDefaults& rock,
    int nSteps,
    double dt,
    const SolverControls& ctrl,
    int writeEvery = 0,
    WriteCallback onWrite = nullptr
)
{
    if (nSteps <= 0 || dt <= 0.0) {
        std::cerr << "[runTransient] invalid nSteps/dt.\n";
        return false;
    }

    double t = 0.0;
    for (int step = 0; step < nSteps; ++step)
    {
        const int step1 = step + 1;
        t += dt;
        // 单步推进（内部含：复制 *_old，初始化 *_prev = *_old，外迭代，装配与求解，欠松弛，提交 p^{n+1},T^{n+1}）
        bool ok = outerIter_OneStep_singlePhase(mgr, reg, freg, ppm,  Pbc, Tbc, gu, rock, dt, ctrl,"CO2");
        if (!ok) {
            std::cerr << "[runTransient] step " << step1 << " failed.\n";
            return false;
        }

        // 可选写出
        if (writeEvery > 0 && (step1 % writeEvery == 0)) {
            if (onWrite) onWrite(step1, t);
            std::cout << "[runTransient] wrote step " << step1 << " at t=" << t << "\n";
        }

    }
    // 最后一步也写一次（如果没有刚好命中 writeEvery）
    if (writeEvery > 0 && onWrite && (nSteps % writeEvery != 0)) {
        onWrite(nSteps, nSteps * dt);
        std::cout << "[runTransient] wrote final step " << nSteps << " at t=" << (nSteps * dt) << "\n";
    }
    return true;
}
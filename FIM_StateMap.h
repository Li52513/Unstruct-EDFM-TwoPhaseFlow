/**
 * @file FIM_StateMap.h
 * @brief 全隐式求解器 (FIM) 全局状态容器 (轻量级工业版)
 * @details
 * 彻底摒弃全局驻留 ADVar，仅存储当前牛顿步的纯标量状态值 (double)。
 * AD 变量的升维与局部播种 (Seeding) 将限制在 Jacobian 组装器的局部循环中，
 * 以确保内存极简与多线程带宽优化。
 */

#pragma once
#include <cstddef>
#include <vector>

template<int N>
struct FIM_StateMap {
    std::vector<double> P;  ///< 主变量 1: 压力 (单相) 或 水相压力 (两相)
    std::vector<double> T;  ///< 主变量 2: 体系温度
    std::vector<double> Sw; ///< 主变量 3: 水相饱和度 (仅 N=3 时分配内存)

    /**
     * @brief 初始化轻量化全局状态数组
     * @param totalBlocks 系统总自由度块数
     */
    void InitSizes(size_t totalBlocks) {
        P.assign(totalBlocks, 0.0);
        T.assign(totalBlocks, 0.0);
        if constexpr (N == 3) {
            Sw.assign(totalBlocks, 0.0);
        }
    }
};
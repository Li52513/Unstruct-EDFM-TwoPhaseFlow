/**
 * @file FVM_Ops.h
 * @brief 有限体积法 (FVM) 基础通用算子库 - Production Ready (Enhanced)
 * @details
 * 包含插值、梯度计算、传导率调和平均等核心离散算子。
 * 采用“两层算子架构”：
 * 1. 数学物理层 (Math Layer)：剥离拓扑，提供纯粹的公式计算 (供 NNC / FF 直接调用)。
 * 2. 拓扑解析层 (Topology Layer)：基于 Mesh 面结构，自动提取几何并调用数学层 (供 Matrix Face 调用)。
 * * 设计为 Header-only 模板库，以支持标量、向量及张量场。
 * @version 1.1.0
 */

#pragma once

#include <cmath>        // std::abs, std::sqrt
#include <algorithm>    // std::max, std::min
#include <limits>       // std::numeric_limits
#include <iostream>     // std::cerr

#include "mesh.h"           // 包含 Mesh, Face, Cell, Node
#include "VolField.h"       // 包含 VolField<T>
#include "UserDefineVarType.h"   // 包含 Vector, Tensor

namespace FVM_Ops
{
    // =========================================================
    // 常量定义
    // =========================================================
    static constexpr double kEpsilon = 1e-20; ///! 防止除零的严格小量，与传导率求解器的极小值标准对齐

    // =========================================================
    // 第一组：插值类算子 (Interpolation Operators)
    // =========================================================

    /**
     * @brief 算术平均插值 (Arithmetic Interpolation)
     * @details 基于几何权重将单元中心值插值到面中心：
     * phi_f = w * phi_O + (1 - w) * phi_N
     * @tparam T 物理量类型 (double, Vector, etc.)
     * @param mesh 网格对象
     * @param field 待插值的体心场
     * @param faceIdx 面索引，在输入前需要确保索引正确
     * @return T 插值后的面心值
     */
    template<typename T>
    inline T Op_Interp_Arithmetic(const Mesh& mesh, const VolField<T>& field, int faceIdx)
    {
        const Face& face = mesh.getFaces()[faceIdx];
        int owner = face.ownerCell;
        size_t idxO = mesh.getCellId2Index().at(owner);
        const T& phiO = field[idxO];

        if (face.isBoundary())
        {
            return phiO; // 边界面直接取 owner 值
        }
        else
        {
            int neighbor = face.neighborCell;
            size_t idxN = mesh.getCellId2Index().at(neighbor);
            const T& phiN = field[idxN];
            double w = face.f_linearInterpolationCoef;
            return w * phiN + (1.0 - w) * phiO;
        }
    }

    // =========================================================
    // 第二组：纯数学物理传导算子 (Math & Physics Core Operators)
    // =========================================================

    /**
     * @brief 纯数学串联热阻/流阻计算 (Series Resistance)
     * @details
     * 适用于非标准拓扑 (如 EDFM 的 NNC 和 FF)，剥离了网格拓扑，仅基于物理距离和传导率。
     * 公式: R = d1 / K1 + d2 / K2
     * @param d1 侧1的投影物理距离
     * @param K1 侧1的物理传导系数 (如渗透率、导热系数)
     * @param d2 侧2的投影物理距离
     * @param K2 侧2的物理传导系数
     * @return double 总物理阻力
     */
    inline double Op_Math_SeriesResistance(double d1, double K1, double d2, double K2)
    {
        // 确保物性不小于 1e-30，防止出现非物理的负值或极端奇异
        double term1 = d1 / (std::max(K1, 1e-30) + kEpsilon);
        double term2 = d2 / (std::max(K2, 1e-30) + kEpsilon);
        return term1 + term2;
    }

    /**
     * @brief 纯数学两点传导率计算 (Two-Point Transmissibility)
     * @details 直接返回可用于乘以势能差的传导率因子：T = Area / R
     * @param d1 侧1的投影距离
     * @param K1 侧1的物理传导系数
     * @param d2 侧2的投影距离
     * @param K2 侧2的物理传导系数
     * @param area 截面积 (Face Area 或 Intersection Area)
     * @return double 界面静态传导率 T
     */
    inline double Op_Math_Transmissibility(double d1, double K1, double d2, double K2, double area)
    {
        double resistance = Op_Math_SeriesResistance(d1, K1, d2, K2);
        return area / (resistance + kEpsilon);
    }

    /**
     * @brief 纯数学调和平均等效系数计算 (Harmonic Mean)
     * @details 公式: K_eff = (d1 + d2) / (d1/K1 + d2/K2)
     * @return double 界面等效传导系数 K_eff
     */
    inline double Op_Math_HarmonicMean(double d1, double K1, double d2, double K2)
    {
        double resistance = Op_Math_SeriesResistance(d1, K1, d2, K2);
        return (d1 + d2) / (resistance + kEpsilon);
    }

    // =========================================================
    // 第三组：拓扑级面心调和算子 (Topological Face Operators)
    // =========================================================

    /**
     * @brief 面心调和平均插值 (Harmonic Interpolation at Matrix Face)
     * @details
     * 适用于标准基岩网格边界 (Matrix-Matrix Face)。
     * 自动提取体心到面心的法向投影距离，并隐式调用底层数学算子。
     * @tparam T 物理量类型 (当前支持标量 double)
     * @param mesh 网格对象
     * @param field 待插值的体心场 (如基岩导热系数场)
     * @param faceIdx 面索引
     * @return T 插值后的面心等效物理量
     */
    template<typename T>
    inline T Op_Interp_Harmonic_Face(const Mesh& mesh, const VolField<T>& field, int faceIdx)
    {
        const Face& face = mesh.getFaces()[faceIdx];
        int owner = face.ownerCell;
        size_t idxO = mesh.getCellId2Index().at(owner);
        const T& phiO = field[idxO];

        // 边界处理：直接返回内部侧属性
        if (face.isBoundary())
        {
            return phiO;
        }

        int neighbor = face.neighborCell;
        size_t idxN = mesh.getCellId2Index().at(neighbor);
        const T& phiN = field[idxN];

        // 提取拓扑几何中心
        const Cell& cellO = mesh.getCells()[idxO];
        const Cell& cellN = mesh.getCells()[idxN];

        // 统一法向量方向 (Owner -> Neighbor)
        Vector n = face.normal;
        if (face.ownerCell_index != static_cast<int>(idxO)) {
            n = Vector(-n.m_x, -n.m_y, -n.m_z);
        }

        // 计算投影距离 d = |(Center - FaceMid) . n|
        Vector dVecO = face.midpoint - cellO.center;
        Vector dVecN = cellN.center - face.midpoint;

        double dO = std::abs(dVecO.m_x * n.m_x + dVecO.m_y * n.m_y + dVecO.m_z * n.m_z);
        double dN = std::abs(dVecN.m_x * n.m_x + dVecN.m_y * n.m_y + dVecN.m_z * n.m_z);

        // 防御性截断，防止网格高度扭曲导致投影距离为0
        dO = std::max(dO, 1e-6);
        dN = std::max(dN, 1e-6);

        // 调用底层纯数学算子完成计算
        return static_cast<T>(Op_Math_HarmonicMean(dO, static_cast<double>(phiO), dN, static_cast<double>(phiN)));
    }

} // namespace FVM_Ops
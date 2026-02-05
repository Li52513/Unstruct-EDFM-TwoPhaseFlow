#pragma once

#ifndef BOX_8DOP_HEADER
#define BOX_8DOP_HEADER

#include <algorithm>
#include <limits>
#include "UserDefineVarType.h" // 包含 Vector 定义
#include "AABB.h"              // 包含 AABB 定义

/**
 * @class Box8DOP
 * @brief 2.5D 离散定向多面体 (8-Discrete Oriented Polytope) 包围盒
 * @details
 * 该结构体用于 2D/3D EDFM 裂缝求交的几何精筛。
 * 几何形态为：XY 平面上的八边形 (8-DOP) 拉伸至 Z 轴高度 (AABB)。
 * 相比纯 AABB，它能更紧密地包围对角线方向的裂缝线段，显著减少
 * "AABB相交但实际分离" (Type II Error) 的情况。
 * * 包含 5 个投影轴的区间：
 * 1. X 轴: [minX, maxX]
 * 2. Y 轴: [minY, maxY]
 * 3. Z 轴: [minZ, maxZ] (即 AABB 的 Z 范围)
 * 4. D1轴 (x+y): [minD1, maxD1]
 * 5. D2轴 (x-y): [minD2, maxD2]
 */
class Box8DOP
{
public:
    // =========================================================
    // 成员变量 (Member Variables)
    // =========================================================

    // 标准笛卡尔坐标轴范围 (兼容 AABB)
    double minX, maxX;
    double minY, maxY;
    double minZ, maxZ;

    // 对角线轴范围
    // D1 = x + y
    double minD1, maxD1;
    // D2 = x - y
    double minD2, maxD2;

    // =========================================================
    // 构造与初始化 (Construction & Initialization)
    // =========================================================

    /**
     * @brief 默认构造函数
     * @details 初始化为无效的反向区间，确保第一次合并时能正确更新。
     */
    Box8DOP();

    /**
     * @brief 构造函数：由线段的两个端点构建
     * @param p1 线段端点 1
     * @param p2 线段端点 2
     */
    Box8DOP(const Vector& p1, const Vector& p2);

    /**
     * @brief 从线段重新构建 8-DOP
     * @param p1 线段端点 1
     * @param p2 线段端点 2
     */
    void fromSegment(const Vector& p1, const Vector& p2);

    // =========================================================
    // 几何检测 (Intersection Tests)
    // =========================================================

    /**
     * @brief 判断 8-DOP 是否与基岩网格的 AABB 相交
     * @details
     * 使用分离轴定理 (SAT)。基岩网格仅需存储 AABB，无需构建 8-DOP。
     * 算法将 AABB 投影到 8-DOP 的 5 个轴上进行区间重叠测试。
     * * @param otherMatrixBox 基岩网格面的 AABB
     * @return true 如果两者在所有轴上都重叠（即可能相交）
     * @return false 如果至少在一个轴上分离（即一定不相交）
     */
    bool overlaps(const AABB& otherMatrixBox) const;
};

#endif // BOX_8DOP_HEADER
#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

#include "UserDefineVarType.h" // 包含 Vector 定义
#include "AABB.h"              // 包含 AABB 定义 (成员: min, max)

/**
 * @class Box14DOP
 * @brief 14-Discrete Oriented Polytope (14-DOP) 包围盒
 * @details
 * 用于 3D 空间的高效几何剔除。
 * 包含 7 个轴向 (每轴 2 个面，共 14 面):
 * Axis 0-2: (1,0,0), (0,1,0), (0,0,1) [即 AABB]
 * Axis 3-6: (1,1,1), (1,-1,1), (1,1,-1), (1,-1,-1) [体对角线方向]
 */
class Box14DOP
{
public:
    // =========================================================
    // 成员变量
    // =========================================================
    // 存储 7 个轴向上的 [min, max] 区间
    double minValues[7];
    double maxValues[7];

    // =========================================================
    // 构造与初始化
    // =========================================================
    Box14DOP();
    void reset();

    // =========================================================
    // 构建方法
    // =========================================================
    /**
     * @brief 将一个 3D 点加入包围盒
     */
    void addPoint(const Vector& p);

    /**
     * @brief 从一个 3D 三角形构建 14-DOP
     */
    void fromTriangle(const Vector& p1, const Vector& p2, const Vector& p3);

    // =========================================================
    // 核心检测方法
    // =========================================================
    /**
     * @brief 检测 14-DOP 与 AABB 是否重叠 (利用 SAT 分离轴定理)
     * @param box 待检测的 AABB (基岩单元包围盒)
     * @return true 如果重叠
     */
    bool overlaps(const AABB& box) const;

private:
    // =========================================================
    // 内部辅助
    // =========================================================
    double getProjection(const Vector& p, int axisIdx) const;
    void getAABBProjectionInterval(const AABB& box, int axisIdx, double& outMin, double& outMax) const;
};
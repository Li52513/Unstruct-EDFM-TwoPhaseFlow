#include "8_DOP.h"
#include <cmath>

// =========================================================
// 构造与初始化实现
// =========================================================

Box8DOP::Box8DOP()
{
    // 初始化为最大可能的反向区间
    double maxVal = std::numeric_limits<double>::max();
    double minVal = std::numeric_limits<double>::lowest(); // -max

    minX = maxVal; maxX = minVal;
    minY = maxVal; maxY = minVal;
    minZ = maxVal; maxZ = minVal;
    minD1 = maxVal; maxD1 = minVal;
    minD2 = maxVal; maxD2 = minVal;
}

Box8DOP::Box8DOP(const Vector& p1, const Vector& p2)
{
    fromSegment(p1, p2);
}

void Box8DOP::fromSegment(const Vector& p1, const Vector& p2)
{
    // 1. 标准轴 (X, Y, Z)
    if (p1.m_x < p2.m_x) { minX = p1.m_x; maxX = p2.m_x; }
    else { minX = p2.m_x; maxX = p1.m_x; }

    if (p1.m_y < p2.m_y) { minY = p1.m_y; maxY = p2.m_y; }
    else { minY = p2.m_y; maxY = p1.m_y; }

    if (p1.m_z < p2.m_z) { minZ = p1.m_z; maxZ = p2.m_z; }
    else { minZ = p2.m_z; maxZ = p1.m_z; }

    // 2. 对角轴 D1 (x + y)
    double d1_p1 = p1.m_x + p1.m_y;
    double d1_p2 = p2.m_x + p2.m_y;

    if (d1_p1 < d1_p2) { minD1 = d1_p1; maxD1 = d1_p2; }
    else { minD1 = d1_p2; maxD1 = d1_p1; }

    // 3. 对角轴 D2 (x - y)
    double d2_p1 = p1.m_x - p1.m_y;
    double d2_p2 = p2.m_x - p2.m_y;

    if (d2_p1 < d2_p2) { minD2 = d2_p1; maxD2 = d2_p2; }
    else { minD2 = d2_p2; maxD2 = d2_p1; }
}

// =========================================================
// 几何检测实现 (SAT Algorithm)
// =========================================================

bool Box8DOP::overlaps(const AABB& otherMatrixBox) const
{
    // ---------------------------------------------------------
    // 阶段 1: 检查标准笛卡尔轴 (X, Y, Z)
    // 这等同于标准 AABB vs AABB 检测
    // ---------------------------------------------------------

    // X 轴投影不重叠
    if (maxX < otherMatrixBox.min.m_x || minX > otherMatrixBox.max.m_x) return false;

    // Y 轴投影不重叠
    if (maxY < otherMatrixBox.min.m_y || minY > otherMatrixBox.max.m_y) return false;

    // Z 轴投影不重叠 (对于 2D 问题，Z通常相等或接近，此检查兼容 3D)
    if (maxZ < otherMatrixBox.min.m_z || minZ > otherMatrixBox.max.m_z) return false;

    // ---------------------------------------------------------
    // 阶段 2: 检查对角线轴 (D1, D2)
    // 这是 8-DOP 的核心优势，能剔除对角线附近的虚空区域
    // ---------------------------------------------------------

    // --- D1 轴: V = x + y ---
    // 计算 AABB 在 D1 轴上的投影区间 [aabbMinD1, aabbMaxD1]
    // min(x+y) 发生在 (minX, minY)
    // max(x+y) 发生在 (maxX, maxY)
    double aabbMinD1 = otherMatrixBox.min.m_x + otherMatrixBox.min.m_y;
    double aabbMaxD1 = otherMatrixBox.max.m_x + otherMatrixBox.max.m_y;

    // 检查 D1 轴是否分离
    if (maxD1 < aabbMinD1 || minD1 > aabbMaxD1) return false;


    // --- D2 轴: V = x - y ---
    // 计算 AABB 在 D2 轴上的投影区间 [aabbMinD2, aabbMaxD2]
    // min(x-y) 发生在 x 最小且 y 最大处 -> (minX - maxY)
    // max(x-y) 发生在 x 最大且 y 最小处 -> (maxX - minY)
    double aabbMinD2 = otherMatrixBox.min.m_x - otherMatrixBox.max.m_y;
    double aabbMaxD2 = otherMatrixBox.max.m_x - otherMatrixBox.min.m_y;

    // 检查 D2 轴是否分离
    if (maxD2 < aabbMinD2 || minD2 > aabbMaxD2) return false;

    // 所有 5 个轴均重叠，判定为相交
    return true;
}
#include "14_DOP.h"

// =========================================================
// 构造与初始化
// =========================================================

Box14DOP::Box14DOP()
{
    reset();
}

void Box14DOP::reset()
{
    for (int i = 0; i < 7; ++i) {
        minValues[i] = std::numeric_limits<double>::max();
        maxValues[i] = std::numeric_limits<double>::lowest();
    }
}

// =========================================================
// 内部辅助: 投影计算
// =========================================================
double Box14DOP::getProjection(const Vector& p, int axisIdx) const
{
    switch (axisIdx) {
    case 0: return p.m_x;
    case 1: return p.m_y;
    case 2: return p.m_z;
    case 3: return p.m_x + p.m_y + p.m_z;
    case 4: return p.m_x - p.m_y + p.m_z;
    case 5: return p.m_x + p.m_y - p.m_z;
    case 6: return p.m_x - p.m_y - p.m_z;
    default: return 0.0;
    }
}

// =========================================================
// 构建方法
// =========================================================

void Box14DOP::addPoint(const Vector& p)
{
    for (int i = 0; i < 7; ++i) {
        double val = getProjection(p, i);
        if (val < minValues[i]) minValues[i] = val;
        if (val > maxValues[i]) maxValues[i] = val;
    }
}

void Box14DOP::fromTriangle(const Vector& p1, const Vector& p2, const Vector& p3)
{
    reset();
    addPoint(p1);
    addPoint(p2);
    addPoint(p3);
}

// =========================================================
// 核心检测: Overlaps AABB
// =========================================================

void Box14DOP::getAABBProjectionInterval(const AABB& box, int axisIdx, double& outMin, double& outMax) const
{
    // [修正] 使用 AABB.h 中定义的 min 和 max 成员
    const double& x0 = box.min.m_x; const double& x1 = box.max.m_x;
    const double& y0 = box.min.m_y; const double& y1 = box.max.m_y;
    const double& z0 = box.min.m_z; const double& z1 = box.max.m_z;

    // 轴对齐情况优化
    if (axisIdx <= 2) {
        if (axisIdx == 0) { outMin = x0; outMax = x1; }
        else if (axisIdx == 1) { outMin = y0; outMax = y1; }
        else { outMin = z0; outMax = z1; }
        return;
    }

    // 针对 4 个对角轴，计算 AABB 8 个角点的投影极值
    // Corners: (x0,y0,z0) ... (x1,y1,z1)

    // 预计算 8 个角点的投影
    double corners[8];
    corners[0] = getProjection(Vector(x0, y0, z0), axisIdx);
    corners[1] = getProjection(Vector(x0, y0, z1), axisIdx);
    corners[2] = getProjection(Vector(x0, y1, z0), axisIdx);
    corners[3] = getProjection(Vector(x0, y1, z1), axisIdx);
    corners[4] = getProjection(Vector(x1, y0, z0), axisIdx);
    corners[5] = getProjection(Vector(x1, y0, z1), axisIdx);
    corners[6] = getProjection(Vector(x1, y1, z0), axisIdx);
    corners[7] = getProjection(Vector(x1, y1, z1), axisIdx);

    outMin = corners[0];
    outMax = corners[0];
    for (int i = 1; i < 8; ++i) {
        if (corners[i] < outMin) outMin = corners[i];
        if (corners[i] > outMax) outMax = corners[i];
    }
}

bool Box14DOP::overlaps(const AABB& box) const
{
    // SAT 检测：只要有一个轴分离，即为不重叠
    for (int i = 0; i < 7; ++i) {
        double boxMin, boxMax;
        getAABBProjectionInterval(box, i, boxMin, boxMax);

        if (maxValues[i] < boxMin || minValues[i] > boxMax) {
            return false;
        }
    }
    return true;
}
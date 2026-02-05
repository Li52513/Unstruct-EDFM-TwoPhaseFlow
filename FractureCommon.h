#pragma once

/// @brief 交点来源类型枚举
enum class IntersectionOrigin
{
    FracFrac,       // 裂缝–裂缝交点
    FracFace,       // 裂缝–网格面交点
    FracStart,      // 裂缝起点
    FracEnd         // 裂缝终点
};

/// @brief 特征距离计算指标枚举
enum class DistanceMetric
{
    CellCenter,     // 使用单元中心点距离
    NodeAverage,    // 使用单元节点平均位置距离
    AreaWeight,     // 使用面积加权的数值积分平均
    CrossAwareGauss // 考虑穿越情形的高斯积分法
};

/// @brief 裂缝段物理类型枚举
enum class FractureElementType
{
    Conductive,     // 导流
    Blocking         // 阻塞
};

// ==========================================
//  求交策略与性能统计
// ==========================================

/// @brief 裂缝-网格面求交搜索策略
enum class IntersectionSearchStrategy_2D
{
    BruteForce,                     // 方法1：暴力遍历 (无 AABB 预筛)
    GlobalAABB,                     // 方法2：全局遍历 + AABB 预筛
    GridIndexing,                   // 方法3：背景网格索引 (Grid-based Indexing)
    GridIndexing_BasedOn8DOP,       // 方法4：背景网格索引 + 8-DOP 几何精筛 
    GridIndexing_BasedOn8DOP_DDA    // 方法5：DDA 射线追踪索引 + 8-DOP 几何精筛 
};

/// @brief 性能统计指标
struct IntersectionStatistics_2D
{
    long long candidateCount = 0;      // 1. 进入循环的候选面总数
    long long aabbCheckCount = 0;      // 2. 执行 AABB 碰撞检测的次数
    long long geometryCheckCount = 0;  // 3. 执行精确几何求交的次数 (耗时大户)
    int trueIntersectionCount = 0;     // 4. 最终命中的真实交点数

    void reset() {
        candidateCount = 0;
        aabbCheckCount = 0;
        geometryCheckCount = 0;
        trueIntersectionCount = 0;
    }
};
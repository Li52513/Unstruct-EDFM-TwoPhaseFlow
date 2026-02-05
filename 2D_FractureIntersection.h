#pragma once

#include "UserDefineVarType.h" // 包含 Vector 定义

// =========================================================
// 3D-EDFM 裂缝-基岩交点类型定义
// 参考文献: An improved embedded discrete fracture model... (Wang et al., 2022)
// =========================================================

/// @brief 3D 空间中的交点类型
enum class FracIntersectionType
{
    None,
    Type1_FracVertex,   // Type 1: 裂缝网格的顶点 (位于基岩单元内部)
    Type2_FracEdge,     // Type 2: 裂缝网格的边 与 基岩单元的面 相交
    Type3_MatrixEdge    // Type 3: 裂缝网格的面 与 基岩单元的棱(Edge) 相交
};

/// @brief 描述一个 3D 几何交点
class IntersectionPoint3D
{
    Vector coord;               // 交点坐标
    FracIntersectionType type;  // 交点类型

    // --- 拓扑溯源信息 (用于构建交互多边形时的连通性分析) ---
    int fracNodeID = -1;        // 若为 Type 1，对应裂缝节点的 ID
    int matrixFaceID = -1;      // 若为 Type 2，对应基岩 Face ID
    int matrixEdgeID = -1;      // 若为 Type 3，对应基岩 Edge ID (需配合 MatrixEdge 使用)

    // 有参构造函数
    IntersectionPoint3D(const Vector& p, FracIntersectionType t)
        : coord(p), type(t) {
    }
};
#pragma once
#include <vector>
#include "UserDefineVarType.h"  // 定义 Vector 类型
#include "FractureCommon.h"

struct FractureIntersectionPoint
{
    int id; // 交点编号
    Vector point; // 交点坐标
    int edgeID; // 与之相交的边（Face）的编号；若无则为 -1
    double param; // 裂缝起点到交点在裂缝线段上的归一化参数（0～1）
    bool    isFF;         // <<< 是否是裂缝–裂缝交点
    int     globalFFID;   // <<< 全局 FF 交点的 ID（face 交点填 0）
    IntersectionOrigin origin;   ///交点的"来源"

    // 裂缝交点的有参构造函数声明
    FractureIntersectionPoint(
        int _id,
        const Vector& _pt,
        int _edgeID,
        double _param,
        bool _isFF = false,
        int  _globalID = 0,
        IntersectionOrigin _orig = IntersectionOrigin::FracFace);
};
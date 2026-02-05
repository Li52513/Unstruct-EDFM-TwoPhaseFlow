#pragma once

#include <vector>

// 引入现有基础结构
#include "UserDefineVarType.h" // 包含 Vector 定义
#include "AABB.h"              // 包含 AABB 定义

// 前向声明 
class Node;

// =========================================================
// 裂缝网格单元 (Fracture Mesh Element / Sub-grid)
// 对应拓扑层级: Level 2
// 物理含义: 宏观裂缝经过独立离散后的网格
// =========================================================

class FractureElement_2D
{
public:
    // =========================================================
    // 1. 基础属性
    // =========================================================
    int id;                         // 单元全局唯一编号
    std::vector<int> nodeIndices;   // 构成该单元的节点索引  对应Fracture_2D类中的fracNodes 的局部下标

    // 索引管理体系 (Index Management)
    int localIndex;                 /// 在所属 Fracture::fracCells 中的局部下标 (0-based Local Index)
	int solverIndex;				/// 在全局线性方程组 AX=B 中的行号 (0-based Solver Index), 初始值为 -1      
    // =========================================================
    // 2. 几何属性 (Geometry)
    // =========================================================
    Vector centroid;                // 几何中心
    Vector normal;                  // 单元局部平均法向 (对角线叉乘)
    double area;                    // 单元面积 (双三角形面积之和)
    AABB boundingBox;               // 包围盒 (用于八叉树快速求交 AABB Test)

    // =========================================================
    // 3. 构造函数声明
    // =========================================================
    
    // 默认构造
    FractureElement_2D();

    // 有参构造
    FractureElement_2D(int _id, const std::vector<int>& _nodes);

    // =========================================================
    // 4. 成员行为声明
    // =========================================================
    
    /**
     * @brief 计算几何属性 (支持扭曲四边形)
     * @param allFracNodes 裂缝面节点列表
     * @details
     * 1. 质心: 顶点平均
     * 2. 面积: Area(Tri_012) + Area(Tri_023)
     * 3. 法向: (Area1 * n1 + Area2 * n2) / TotalArea
     */
    void computeGeometry(const std::vector<Node>& allFracNodes);

};
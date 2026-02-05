#pragma once

#include <vector>
#include <map>
#include <unordered_map>

#include "UserDefineVarType.h"
#include "Node.h"
#include "AABB.h"


class Cell 
{
public:
    
    // =========================================================
    // 1. 类型定义 (Type Definitions)
    // =========================================================
    // 网格单元位置分类
    enum class LocationType
    {
        Inner,      // 内部单元
        Boundary,   // 边界单元
    };

    // 区域分类
    enum class RegionType
    {
        Low,
        Medium,
        High
    };
    // =========================================================
    // 2. 成员属性 (Member Properties)
    // =========================================================

    // --- 基础几何信息 ---
    int id;                                 // 单元编号
    vector<int> CellFaceIDs;                // 构成单元的面编号
    vector<int> CellNodeIDs;                // 构成单元的节点编号
    Vector center;                          // 单元中心（几何重心）
    double volume;                          // 单元面积（对三角形）/ 体积

   
    AABB boundingBox;                       //包围盒，用于加速几何求交

    // --- 分类状态 ---
    LocationType location = LocationType::Inner; // 默认内部单元
    RegionType region = RegionType::Medium;      // 默认中等

    // =========================================================
    // 3. 成员行为 (Member Behaviors)
    // =========================================================
    
    //构造函数
    Cell(int id, const vector<int>& nodeIDs);  // 构造函数

    //计算2D多边形面积与Z平均值
    void calculatePolygonProps(const std::vector<int>& ids,
        const std::unordered_map<int, Node>& allNodes, double& outArea, Vector& outCentroid);


    //计算单元的几何中心和面积
    ///2D基岩网格
	void computeCenterAndVolume_2D(const unordered_map<int, Node>& allNodes); 

    ///3D基岩网格
    void computeCenterAndVolume_3D(const unordered_map<int, Node>& allNodes);

    //获取单元的局部面信息，在mesh.h的BuildMesh函数中进行了调用，用于构建face与cell和node之间的拓扑关系
    ///2D基岩网格
    vector<vector<int>> getLocalFaces_2D(const std::unordered_map<int, Node>& allNodes) const; // 获取单元的局部面信息 
   
    ///3D基岩网格
    vector<vector<int>> getLocalFaces_3D() const;

    // 计算单元的 AABB 包围盒
    // 需传入全局节点表以查询坐标
    void computeAABB(const std::unordered_map<int, Node>& allNodes);
};
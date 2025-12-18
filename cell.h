#pragma once
#include <vector>
#include <map>
#include "UserDefineVarType.h"
#include "Node.h"
#include "PropertiesSummary.h"
#include <unordered_map>


class Cell 
{
public:
   
    
    //------------------网格单元信息----------------------//
    int id;                          // 单元编号
    vector<int> CellFaceIDs;             // 构成单元的面编号
    vector<int> CellNodeIDs;             // 构成单元的节点编号
    Vector center;                   // 单元中心（几何重心）
    double volume;                   // 单元面积（对三角形）
    

	//-----------------网格单元分类------------------//
    enum class LocationType
    {
		Inner,       // 内部单元
		Boundary,    // 边界单元
    };
	LocationType location = LocationType::Inner; // 默认内部单元

    //-----------------区域分类-------------------------//
    enum class RegionType
    {
        Low,
		Medium,
		High
    };
    RegionType region = RegionType::Medium; // 默认中等
    
    //------------------裂缝交互信息----------------------//
    vector<int> fractureIDs;          // 与该单元交互的裂缝编号
    vector<double> fracturePressure;  // 存储与该单元交互的裂缝的压力
    vector<Vector> fracturePositions; // 存储裂缝交点的位置
    vector<double> fractureFlux;      // 裂缝流量影响

    //------------------离散化系数----------------------//
    double faceDiscreCoef;            // 面的离散化系数（控制体与面之间的流量系数）
    double sourceTerm;                // 源项，可能来自裂缝或外部来源

    //------------------计算辅助信息----------------------//
    vector<double> CI;                // 离散化系数，用于求解矩阵
    vector<int> CI_belongFraction;   // 存储相邻单元ID
    double error;                     // 误差，用于收敛判断

    Cell(int id, const vector<int>& nodeIDs);  // 构造函数

	void computeCenterAndVolume(const unordered_map<int, Node>& allNodes); //计算单元的几何中心和面积

    vector<vector<int>> getLocalFaces(const std::unordered_map<int, Node>& allNodes) const; // 获取单元的局部面信息
   

};

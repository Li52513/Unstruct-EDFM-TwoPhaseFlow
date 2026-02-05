#pragma once
#include <vector>
#include <algorithm>
#include <cmath>

#include "UserDefineVarType.h"  // 定义 Vector 类型（即 Variable3D<double>）

#include "AABB.h"
#include "Mesh.h"
#include "Cell.h"
#include "Node.h"



#include "FractureCommon.h"
#include "FractureElement.h"
#include "FractureIntersectionPoint.h"


using namespace std;

class Fracture  
{
public:
	
    //构造函数
	Fracture(const Vector& s, const Vector& e);   // 裂缝起点、终点坐标初始化构造函数

    //成员属性
    ///基础id与几何信息
    int id  ;               // 裂缝编号
    Vector start;           // 裂缝起点坐标
	Vector end;             // 裂缝终点坐标
    AABB fractureBox;       // 裂缝的包围盒（用于碰撞检测等）

    ///离散化信息
	vector<FractureIntersectionPoint> intersections;    // 裂缝交点信息
	vector<FractureElement> elements;		            // 裂缝段信息
    
    ///辅助统计信息
	int candidateFaceCount_ = 0;        //用于记录精确检测的候选面               

    //静态成员函数
    ///计算裂缝段长度
	static double computeSegmentLength(const Vector& a, const Vector& b); 

	///计算裂缝段中点坐标
	static Vector computeMidpoint(const Vector& a, const Vector& b);
	
    ///精确判断两线段是否有交点且在两线段内
    static bool lineSegmentIntersection(const Vector& p, const Vector& q, 
        const Vector& r, const Vector& s, Vector& ip); 

	///定位裂缝段所在基岩网格单元
    static int findContainingCell(const Vector& point, const std::vector<Cell>& cells, 
        const std::unordered_map<int, Node>& nodes);

	///计算裂缝段到基岩网格单元的距离--方法1利用 Cell 中的各节点到裂缝段的距离计算平均距离
	static double computeAverageDistanceFromNodes(const Cell& cell, const std::unordered_map<int, Node>& meshNodes, 
        const Vector& segStart, const Vector& segEnd); 

	///计算裂缝段到基岩网格单元的距离--方法2利用 Cell 中点计算到裂缝段的距离
	static double computeDistanceFromCenter(const Cell& cell, const Vector& segStart, const Vector& segEnd); 

    ///计算裂缝段到基岩网格单元的距离--方法3 高精度方法。将单元划分为多个三角形，通过面积加权积分计算平均距离
    static double computeAreaWeightedDistance(const Cell& cell, const std::unordered_map<int, Node>& meshNodes, 
        const Vector& segStart, const Vector& segEnd);

    ///计算裂缝段到基岩网格单元的距离--方法4 最高精度方法。考虑裂缝是否完全贯穿单元，若贯穿则分割多边形分别积分（高斯积分+扇形分割）
    static double computeCrossAwareAverageDistance(const Cell& cell, const std::unordered_map<int, Node>& meshNodes, 
        const Vector& segStart, const Vector& segEnd);

    //成员行为
    ///根据该裂缝的起点和终点，计算并更新 fracture的AABBbox 属性
    void computeAABB(); 

	// 检测裂缝与基岩网格面的交点（改进版，支持多种搜索策略）
    void DetectFracturetoMeshFaceIntersections_improved
    (
        const Mesh& mesh,
        const std::vector<Cell>& meshCells,
        const std::unordered_map<int, Node>& meshNodes,
        IntersectionSearchStrategy_2D strategy,      // <--- 策略选择
        IntersectionStatistics_2D& stats             // <--- 统计输出
    );

    ///排序交点并重新编号
    void sortAndRenumberIntersections();
    
    ///裂缝段划分
    void subdivide(const vector<Cell>& meshCells, const unordered_map<int, Node>& meshNodes, DistanceMetric metric);

    ///计算几何耦合系数 geomCI 和 geomAlpha
    void computeGeometryCouplingCoefficientgeomCIandgeomAlpha();
    
    ///给定 param，定位它属于哪一段
    int locateSegment(double param) const; 
};
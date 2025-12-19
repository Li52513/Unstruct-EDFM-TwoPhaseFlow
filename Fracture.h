#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include "UserDefineVarType.h"  // 定义 Vector 类型（即 Variable3D<double>）
#include "Cell.h"
#include "Face.h"
#include "Node.h"
#include "Mesh.h"
#include "FractureSolidProperties.h"

/// 交点来源类型进行枚举
enum class IntersectionOrigin 
{
    FracFrac,       // 裂缝–裂缝交点
    FracFace,       // 裂缝–网格面交点
    FracStart,      // 裂缝起点
    FracEnd         // 裂缝终点
};


enum class DistanceMetric 
{
	CellCenter, // 使用单元中心点距离
	NodeAverage, // 使用单元节点平均位置距离
	AreaWeight,   // 使用面积加权的“数值积分” <d> = sum(d_j * ΔS_j) / S
    CrossAwareGauss 
};

struct FractureIntersectionPointByMatrixMesh  //// 用于描述裂缝与基岩网格面交点的结构体
{
	int id; // 交点编号
	Vector point; // 交点坐标
	int edgeID; // 与之相交的边（Face）的编号；若无则为 -1
	double param; // 裂缝起点到交点在裂缝线段上的归一化参数（0～1）
    bool    isFF;         // <<< 是否是裂缝–裂缝交点
    int     globalFFID;   // <<< 全局 FF 交点的 ID（face 交点填 0）
    IntersectionOrigin origin;   ///交点的“来源”
    
    FractureIntersectionPointByMatrixMesh(int _id,
        const Vector& _pt,
        int _edgeID,
        double _param,
        bool _isFF = false,
        int  _globalID = 0,
        IntersectionOrigin _orig = IntersectionOrigin::FracFace)
        : id(_id)
        , point(_pt)
        , edgeID(_edgeID)
        , param(_param)
        , isFF(_isFF)
        , globalFFID(_globalID)
        , origin(_orig)
    {
	} //有参构造函数
};

struct FractureElement  ///// 描述裂缝单元（裂缝段）的结构体
{
    
	// ————— 编号信息 —————  （在哪里）
    int id;                  ///< 裂缝单元序号
    int cellID;              ///< 裂缝单元所在基岩网格单元编号
    // —————— 几何信息 —————— （有多大）
	double aperture = 1.0;         ///< 裂缝单元开度 w
    double length = 0.0;           ///< 裂缝单元长度 L
    double avgDistance =0.0;      ///< 平均距离 d
    double geomCI = 0.0;   ///< 纯几何耦合系数  CI_geom = (L·w) / d
    double geomAlpha = 0.0;   ///< 纯几何 α_geom = 2·w / L  可以进一步修正
    //——————裂缝段类型———————
    FractureElementType type = FractureElementType::Conductive;
	///< 裂缝段类型（阻塞/导流）

    // —————— 物性耦合 & 离散系数 —————— （储存离散系数）
    double alpha_fr = 0.0;   ///< TI 计算用 α_phys = 2·w·k / (μ·L)
    double k_eff = 0.0;   ///< 并联系统有效渗透率
    double aW_fr = 0.0;   ///< 左邻接段导流系数
    double aE_fr = 0.0;   ///< 右邻接段导流系数
    double b_fr = 0.0;   ///< 裂缝–基岩源项系数 CI_phys = θ·A / d
    double aP_fr = 0.0;   ///< 本段总离散系数 aW+aE+b
    // —————— 归一化参数区间 ——————
    double param0 = 0.0;   ///< 本段起始归一化位置
    double param1 = 0.0;   ///< 本段终止归一化位置

    // —— 裂缝–裂缝 TI 交换信息 ——
    struct FFF_Exchange
    {
		int peerFracID{ -1 }; // 对端裂缝编号
		int peerSegID{ -1 };  // 对端裂缝段编号
        int peerGlobalSeg{ -1 };   // 对端全局段索引（便于快速定位）
		int atGlobalFF{ -1 };	  // 所在全局 FF 交点编号
        double TIw{ 0.0 };         // 水相交汇导纳（Star–Delta 后的“边”）
        double TIg{ 0.0 };         // 气相交汇导纳
    };
    vector<FFF_Exchange> ffEx;

    bool isFFatStart = false, isFFatEnd = false;
    int  gIDstart = 0, gIDend = 0;

    FractureElement(int id_, int cellID_, double len = 0.0, double avgDist = 0.0)
        : id(id_)
        , cellID(cellID_)
        , length(len)
        , avgDistance(avgDist)
    {
        // 其余成员已经在声明处给出默认值
    }
};

class Fracture  
{
public:
	
    /*=====构造函数=====*/
    Fracture(const Vector& s, const Vector& e);
	/*===裂缝编号===*/ 
	int id  ; // 裂缝编号
        
    /*===裂缝起点终点信息===*/
    Vector start; // 裂缝起点坐标
	Vector end;   // 裂缝终点坐标

	/*===存储裂缝交点及裂缝段信息===*/
	vector<FractureIntersectionPointByMatrixMesh> intersections; 
	vector<FractureElement> elements;		  
 
    //构建每条裂缝的AABB
	AABB fractureBox; // 裂缝的包围盒（用于碰撞检测等）
    int candidateFaceCount_ = 0; // 可选加 mutable
    
	// 计算裂缝与网格面交点
    void DetectFracturetoMeshFaceIntersections(const Mesh& mesh, const std::vector<Cell>& meshCells, const std::unordered_map<int, Node>& meshNodes, bool useAABBFilter);

    void computeAABB();

	// 排序交点并重新编号
    void sortAndRenumberIntersections();

    //裂缝段划分
    void subdivide(const vector<Cell>& meshCells, const unordered_map<int, Node>& meshNodes, DistanceMetric metric);

/// -/*=== 计算几何耦合系数 geomCI 和 geomAlpha======*/
    void computeGeometryCouplingCoefficientgeomCIandgeomAlpha();
        
    /*==裂缝段离散====*/
  

    /*== 给定 param，定位它属于哪一段====*/
    int locateSegment(double param) const;

    static bool lineSegmentIntersection(const Vector& p, const Vector& q, const Vector& r, const Vector& s, Vector& ip);
    


private:
    
    double computeSegmentLength(const Vector& a, const Vector& b);
    Vector computeMidpoint(const Vector& a, const Vector& b);
    int findContainingCell(const Vector& point, const std::vector<Cell>& cells, const std::unordered_map<int, Node>& nodes);
    static double computeAverageDistanceFromNodes(const Cell& cell, const std::unordered_map<int, Node>& meshNodes, const Vector& segStart, const Vector& segEnd);
    static double computeDistanceFromCenter(const Cell& cell, const Vector& segStart, const Vector& segEnd);
    static double computeAreaWeightedDistance(const Cell& cell, const std::unordered_map<int, Node>& meshNodes, const Vector& segStart, const Vector& segEnd);
    static double computeCrossAwareAverageDistance(const Cell& cell, const std::unordered_map<int, Node>& meshNodes, const Vector& segStart, const Vector& segEnd);
};


//#pragma once
//using namespace std;
//#include <vector>
//#include <algorithm>
//#include "UserDefineVarType.h"
//
//class Fraction
//{
//public:
//	//the Geometry&Mesh information of the fraction/////////
//	
//	int fracNum;				// number of fractions
//	int fraction_ID;			// fraction ID
//	Vector startPoint_pos;		// start point position
//	Vector endPoint_pos;		// end point position
//	int belongMatrix_ID;		// matrix ID to which the fraction belongs
//	double A_fr;				// area of the fraction
//	double d_fr;				// distance between the fraction to the center of the cell
//	vector<int>inter_fraction_ID;// the ID of the intersection fraction
//	vector<Fraction>fraction_cell;// store discrete fracture segments
//	vector<Vector>node;		// store the node position of the fracture
//
//	//the Physical information of the fraction/////////
//	double k_f;					// permeability of the fraction
//	double k_eff;				// effective permeability of the fraction
//	double phie_f;				// porosity of the fraction
//	double k_width;				// width of the fraction
//	vector <double> interPointofOtherFraction_pressure; // pressure of the intersection point of other fractions
//	vector<Vector>interPointwithFractions_pos; // position of the intersection point with other fractions
//	int interFractionCellwithOthers_ID; // the ID of the intersection fraction cell with other fractions
//
//	//the Numerical information of the fraction/////////
//	double aP_fr;				// aP coefficient of the fraction
//	double aE_fr;				// aE coefficient of the fraction
//	double aW_fr;				// aW coefficient of the fraction
//	double aP0_fr;				// aP0 coefficient of the fraction
//	double b_fr;				// b coefficient of the fraction
//	double b_fr_TI_alpha;		// TI_coefficient of the fraction
//	vector <double> b_fr_TI;	// TI_coefficient of the fraction
//	double Erroe_fr;			// Error of the fraction
//
//};
#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include "UserDefineVarType.h"  // 定义 Vector 类型（即 Variable3D<double>）
#include "Cell.h"
#include "Face.h"
#include "Node.h"
#include "Fluid.h"
#include "Matrix.h"
#include "Mesh.h"
#include "FractureTypes.h"

#pragma message("Fracture.h included!")

/// 交点来源类型进行枚举
enum class IntersectionOrigin 
{
    FracFrac,       // 裂缝–裂缝交点
    FracFace,       // 裂缝–网格面交点
    FracStart,      // 裂缝起点
    FracEnd         // 裂缝终点
};

struct FractureIntersectionPoint  //// 用于描述裂缝与基岩网格边界交点的结构体
{
	int id; // 交点编号
	Vector point; // 交点坐标
	int edgeID; // 与之相交的边（Face）的编号；若无则为 -1
	double param; // 裂缝起点到交点在裂缝线段上的归一化参数（0～1）
    bool    isFF;         // <<< 新增：是否是裂缝–裂缝交点
    int     globalFFID;   // <<< 新增：全局 FF 交点的 ID（face 交点填 0）
    IntersectionOrigin origin;   ///交点的“来源”
    
    FractureIntersectionPoint(int _id,
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
	double aperture = 0.0;         ///< 裂缝单元开度 w
    double length = 0.0;           ///< 裂缝单元长度 L
    double avgDistance =0.0;      ///< 平均距离 d
    double geomCI = 0.0;   ///< 纯几何耦合系数  CI_geom = (L·w) / d
    double geomAlpha = 0.0;   ///< 纯几何 α_geom = 2·w / L  可以进一步修正

    //——————裂缝段类型———————
    FractureElementType type = FractureElementType::Conductive;
	///< 裂缝段类型（阻塞/导流）
	
	// —————— 储存物性参数 —————— （特性是什么）
    SolidProperties_Frac solidProps;  ///< 固相骨架物性
    WaterProperties fluidProps;  ///< 水相物性
    CO2Properties   gasProps;    ///< CO₂ 相物性

    // —————— 物性耦合 & 离散系数 —————— （储存离散系数）
    double alpha_fr = 0.0;   ///< TI 计算用 α_phys = 2·w·k / (μ·L)
    double k_eff = 0.0;   ///< 并联系统有效渗透率

    double aW_fr = 0.0;   ///< 左邻接段导流系数
    double aE_fr = 0.0;   ///< 右邻接段导流系数
    double b_fr = 0.0;   ///< 裂缝–基岩源项系数 CI_phys = θ·A / d
    double aP_fr = 0.0;   ///< 本段总离散系数 aW+aE+b

    // —————— 待求解变量 ——————  （仍待补充 T 和 s）
    double p_fr = 5.5e6;   ///< 裂缝段压力
	double T_fr = 383.15;   ///< 裂缝段温度
	double s_water_fr = 0.8;   ///< 裂缝段水相饱和度
	double s_CO2_fr = 0.2;   ///< 裂缝段气相饱和度

    // —————— 归一化参数区间 ——————
    double param0 = 0.0;   ///< 本段起始归一化位置
    double param1 = 0.0;   ///< 本段终止归一化位置

    // —— 裂缝–裂缝 TI 交换信息 ——
    struct FFF_Exchange
    {
        int    peerFracID;      ///< 对端裂缝索引
        int    peerSegID;       ///< 对端裂缝段索引
        double TI;              ///< 交换系数
        double* peerPressure;   ///< 指向对端段的 p_fr
    };
    vector<FFF_Exchange> ffEx;

    bool isFFatStart = false, isFFatEnd = false;
    int  gIDstart = 0, gIDend = 0;

    FractureElement(int id_, int cellID_, double len = 0.0, double avgDist = 0.0)
        : id(id_)
        , cellID(cellID_)
        , length(len)
        , avgDistance(avgDist)
		, solidProps()
        , fluidProps()
        , gasProps()
    {
        // 其余成员已经在声明处给出默认值
    }
};

class Fracture  
{
public:
	
	/*===裂缝编号===*/ 
	int id  ; // 裂缝编号
        
    /*===裂缝起点终点信息===*/
    Vector start; // 裂缝起点坐标
	Vector end;   // 裂缝终点坐标

	/*===存储裂缝交点及裂缝段信息===*/
	vector<FractureIntersectionPoint> intersections; 
	vector<FractureElement> elements;		  

	///*===裂缝物性参数（考虑到流体物性的变化以及多物理场过程对裂缝物性的影响，将物性赋值在裂缝段内）===*/
   
     /*==裂缝与基岩交点计算====*/
    void DetectFracturetoMeshFaceIntersections(const vector<Face>& meshFaces, const std::vector<Cell>& meshCells, const std::map<int, Node>& meshNodes);
  
 /*=====构造函数=====*/
    Fracture(const Vector& s, const Vector& e);

/// -/*=== 计算几何耦合系数 geomCI 和 geomAlpha======*/
/// 仅基于几何：计算每段的 geomCI 和 geomAlpha
/// geomCI   = (段长 * aperture) / avgDistance
/// geomAlpha= 2 * aperture / 段长
    void computeGeometryCouplingCoefficientgeomCIandgeomAlpha();
        
    /*==裂缝段离散====*/
    void subdivide(const vector<Cell>& meshCells, const map<int, Node>& meshNodes,  bool useCenterDistance = false);    

    /*== 给定 param，定位它属于哪一段====*/
    int locateSegment(double param) const;

    static bool lineSegmentIntersection(const Vector& p, const Vector& q, const Vector& r, const Vector& s, Vector& ip);
    void sortAndRenumberIntersections();

    

private:
    static double pointToSegmentDistance(const Vector& p, const Vector& s, const Vector& e);
    static bool pointInTriangle(const Vector& p, const Vector& a, const Vector& b, const Vector& c);
    double computeSegmentLength(const Vector& a, const Vector& b);
    Vector computeMidpoint(const Vector& a, const Vector& b);
    int findContainingCell(const Vector& point, const std::vector<Cell>& cells, const std::map<int, Node>& nodes);
    static double computeAverageDistanceFromNodes(const Cell& cell, const std::map<int, Node>& meshNodes, const Vector& segStart, const Vector& segEnd);
    static double computeDistanceFromCenter(const Cell& cell, const Vector& segStart, const Vector& segEnd);
};


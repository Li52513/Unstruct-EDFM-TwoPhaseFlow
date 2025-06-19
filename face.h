#pragma once
#include <vector>
#include <cmath>
#include "UserDefineVarType.h"
using namespace std;

// ---------- 非正交校正方法枚举 ----------
enum class NormalVectorCorrectionMethod 
{
    MinimumCorrection,    // 最小修正值法（Minimum Correction Approach）
    OrthogonalCorrection, // 正交修正法（Orthogonal Correction Approach）
    OverRelaxed          // 超松弛修正法（Over-relaxed Approach）
};


class Face 
{
public:
    int id;                              // 面的编号
    vector<int> nodeIDs;            // 构成该面节点的编号（2D 为2个，3D 为3个）
    vector<Vector> nodeCoords;      // 构成面节点的坐标
    double length;                       // 面的长度（2D边长 or 3D面积）
    Vector normal;                       // 面的法向量（单位向量）
    Vector midpoint;                     // 面的几何中心
    int ownerCell;                       // 所属的单元编号
    int neighborCell;                    // 相邻的单元编号（无为 -1）

    Vector ownerToNeighbor;            // 从 ownerCell.center 指向 neighborCell.center 的向量
    Vector vectorE;                    // 分解得到的正交分量（沿 d 方向）
    Vector vectorT;                    // 分解得到的非正交分量
    double faceDiscreCoef = 0.0;       // 面的离散化系数 a_face
    double faceSource = 0.0;        // 面的源项

    // —— 下面这些字段用于压力插值和梯度（待修正） ——  
    Vector f_pressureGrad;             // 面上压力梯度
    double f_linearInterpolationCoef;  // 插值权重（φ_face = coef·φ_owner + …）
    double f_pressureValue;            // 面上压力值 φ_face

    Face(int id, const std::vector<int>& nodeIDs, const std::vector<Vector>& nodeCoords);
    bool isBoundary() const { return neighborCell == -1; } //判断自己是否属于边界面

    // ---------- 在 owner/neighbor 重心已知时，分解 A_j = E_j + T_j ----------
/**
 * @param Cp      ownerCell 的重心坐标
 * @param Cn      neighborCell 的重心坐标
 * @param method  选择三种校正策略之一
 */
    void computeFaceVectors
    (
        const Vector& Cp,
        const Vector& Cn,
        NormalVectorCorrectionMethod method = NormalVectorCorrectionMethod::MinimumCorrection
    );


private:
    void computeGeometry();  // 根据维度自动计算几何属性（中点、法向量、长度或面积）
};


//#pragma once
//#include <vector>
//#include <algorithm>
//#include "UserDefineVarType.h"
//using namespace std;
//
//class Face 
//{
//public:
//    int id;                              // 面的编号
//    vector<int> nodeIDs;            // 构成该面节点的编号（一般为2个）
//    vector<Vector> nodeCoords;      // 构成面节点的坐标
//    double length;                       // 面的长度（边长）
//    Vector normal;                       // 面的法向量（在二维问题中取 (dy, -dx)）
//    Vector midpoint;                     // 面的中点坐标    
//    int ownerCell;                       // face所属的单元：ownerCell为第一个包含该面单元编号，
//    int neighborCell;                    // 如果面同时属于两个单元，则neighborCell为第二个单元编号，否则为 -1 表示边界面
//
//    // 构造函数中初始化编号及节点信息，并计算几何属性
//    Face(int id, const vector<int>& nodeIDs, const vector<Vector>& nodeCoords)
//        : id(id), nodeIDs(nodeIDs), nodeCoords(nodeCoords), length(0.0), ownerCell(-1), neighborCell(-1) 
//        { 
//            computeGeometry();
//        }
//
//    // 根据组成面节点坐标计算中点、边长和法向量（仅适用于二维直线面）
//    void computeGeometry() 
//    {
//        if (nodeCoords.size() >= 2) 
//        {
//			midpoint = (nodeCoords[0] + nodeCoords[1]) / 2.0; // 计算面的中点坐标
//			Vector diff = nodeCoords[1] - nodeCoords[0];      // 计算面的边向量
//			length = diff.Mag();							 // 计算面的边长  
//			normal = Vector(diff.m_y, -diff.m_x, 0.0);       // 计算面的法向量(二维情形下：法向量 = (dy, -dx, 0))
//            double n_norm = normal.Mag();
//            if (n_norm != 0)
//				normal = normal / n_norm;               // 单位化处理
//        }
//    }
//};

//class MatrixFace
//{
//public:
//	////the Geometry&Mesh information of the face///////
//	int f_ID;					// face ID
//	Vector f_centerPosition;	// face center position
//	Vector f_vectorArea;		// normal area vector of the face
//	vector<int> f_inNodeID;		// node ID included in the face
//	vector<int> f_inCellID;		// cell ID included in the face
//	int f_ownerID;				// owner cell ID
//	int f_neighborID;			// neighbor cell ID
//	Vector f_ownerToNeighbor;	// vector from owner cell to neighbor cell (d)
//
//	//////the Physical information of the face///////
//	double f_pressure;			// face pressure value
//	Vector f_pressureGrad;		// face pressure gradient
//
//	///the Numerical information of the face
//	Vector f_vectorE;			// face vector E
//	Vector f_vectorT;			// face vector T
//	double faceDiscreCoef;		// face discretization coefficient
//	double faceSource;			// face source term
//	double f_linearInterpolationCoef; // linear interpolation coefficient at the face
//
//public:
//	MatrixFace();
//};
#pragma once
#include <vector>
#include <cmath>
#include "UserDefineVarType.h"
#include "AABB.h"
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
    vector<int> FaceNodeIDs;            // 构成该面节点的编号（2D 为2个，3D 为3个）
    vector<Vector> FaceNodeCoords;      // 构成面节点的坐标
    double length;                       // 面的长度（2D边长 or 3D面积）
    Vector normal;                       // 面的法向量（单位向量）
    Vector midpoint;                     // 面的几何中心
    int ownerCell;                       // 所属的单元编号
    int neighborCell;                    // 相邻的单元编号（无为 -1）
	AABB boundingBox;                   // 面的包围盒（用于碰撞检测等）

    Vector ownerToNeighbor;            // 从 ownerCell.center 指向 neighborCell.center 的向量
    Vector vectorE;                    // 分解得到的正交分量（沿 d 方向）
    Vector vectorT;                    // 分解得到的非正交分量
    double faceDiscreCoef = 0.0;       // 面的离散化系数 a_face
    double faceSource = 0.0;        // 面的源项

    // —— 下面这些字段用于压力插值和梯度（待修正） ——  
    Vector f_pressureGrad;             // 面上压力梯度
    double f_linearInterpolationCoef;  // 插值权重（φ_face = coef·φ_owner + …）
    double f_pressureValue;            // 面上压力值 φ_face

    Face(int id, const vector<int>& nodeIDs, const vector<Vector>& nodeCoords);
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
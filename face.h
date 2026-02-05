#pragma once

#include <vector>
#include <cmath>
#include "UserDefineVarType.h"
#include "AABB.h"


// =========================================================
// 非正交校正方法枚举
// =========================================================
enum class NormalVectorCorrectionMethod 
{
    MinimumCorrection,    // 最小修正值法（Minimum Correction Approach）
    OrthogonalCorrection, // 正交修正法（Orthogonal Correction Approach）
    OverRelaxed          // 超松弛修正法（Over-relaxed Approach）
};

class Face 
{
public:
    // =========================================================
    // 1. 成员属性 (Member Properties)
    // =========================================================

    // --- 拓扑结构与标识 (Topology & Identity) ---
    int id;                              // 面的编号
    std::vector<int> FaceNodeIDs;        // 构成该面节点的编号（2D 为2个，3D 为3个）
    int physicalGroupId = -1;            // 默认 -1 表示内部面，正数表示边界 ID
    int ownerCell;                       // 所属的单元编号
    int neighborCell;                    // 相邻的单元编号（无为 -1）

	int ownerCell_index = -1;             // 所属单元在 cells_ 中的索引
	int neighborCell_index = -1;		  // 相邻单元在 cells_ 中的索引

    // --- 几何属性 (Geometry) ---
    std::vector<Vector> FaceNodeCoords;  // 构成面节点的坐标
    double length;                       // 面的长度（2D边长 or 3D面积）
    Vector normal;                       // 面的法向量（单位向量）
    Vector midpoint;                     // 面的几何中心
	AABB boundingBox;                   // 面的包围盒（用于碰撞检测等）

    // --- 有限体积离散化几何参数 (Discretization Data) ---
    Vector ownerToNeighbor;            // 从 ownerCell.center 指向 neighborCell.center 的向量
    Vector vectorE;                    // 分解得到的正交分量（沿 d 方向）
    Vector vectorT;                    // 分解得到的非正交分量
    double f_linearInterpolationCoef;  // 插值权重（φ_face = coef·φ_owner + …）
   
    // =========================================================
    // 2. 成员行为 (Member Behaviors)
    // =========================================================

    // --- 构造与基础状态 ---
    Face(int id, const std::vector<int>& nodeIDs, const std::vector<Vector>& nodeCoords);

    // 判断自己是否属于边界面
    bool isBoundary() const { return neighborCell == -1; } //判断自己是否属于边界面

    // --- 几何与向量计算 ---

    /**
     * @brief 在 owner/neighbor 重心已知时，分解 A_j = E_j + T_j
     * @param Cp      ownerCell 的重心坐标
     * @param Cn      neighborCell 的重心坐标
     * @param method   选择三种校正策略之一
     */
    void computeFaceVectors
    (
        const Vector& Cp,
        const Vector& Cn,
        NormalVectorCorrectionMethod method = NormalVectorCorrectionMethod::MinimumCorrection
    );
    
    // 边界面的向量计算
    void computeFaceVectorsBoundary(
        const Vector& Cp,
        NormalVectorCorrectionMethod method = NormalVectorCorrectionMethod::OrthogonalCorrection
    );

private:
    // --- 内部辅助函数 ---
    void computeGeometry();  // 根据维度自动计算几何属性（中点、法向量、长度或面积）
};
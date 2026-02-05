#pragma once

#include <vector>
#include <cmath>
#include "UserDefineVarType.h" // Vector
#include "Node.h"
#include "Face.h"
#include "2D_FractureMeshElement.h"

// 前向声明
class FractureElement_2D;

// =========================================================
// 裂缝内部边 (Fracture Internal Edge / Face)
// 对应层级: Level 2 (微观拓扑的连接件)
// 物理含义: 连接两个裂缝网格单元的界面，用于计算裂缝内流动的通量
// =========================================================
class FractureEdge_2D
{
public:
    // =========================================================
    // 1. 基础属性
    // =========================================================
    int id;                         ///< 全局唯一编号
    int nodeIndices[2];             ///< 构成边的两个节点 (局部索引)

    int ownerCellID;                ///< 拥有的单元 (Local ID)
    int neighborCellID;             ///< 邻接的单元 (Local ID, -1表示边界)

    int ownerCell_index = -1;       ///< Owner 在 fractureElements 数组中的直接下标 (0-based)
    int neighborCell_index = -1;    ///< Neighbor 在 fractureElements 数组中的直接下标 (-1 表示无)

    double f_linearInterpolationCoef = 0.5; ///< 线性插值权重 w (phi_f = w*O + (1-w)*N)

    Vector ownerToNeighbor;         ///< 两个单元中心的距离向量 (d_ON = x_N - x_O)

    // =========================================================
    // 2. 几何属性
    // =========================================================
    double length;                  ///< 边长 (等同于 3D Face 的 Area)
    Vector midpoint;                ///< 中点
    Vector normal;                  ///< 单位法向量 (位于裂缝平面内，垂直于边，指向 Neighbor)

    // =========================================================
    // 3. FVM 离散化向量 (用于非正交修正)
    // =========================================================
    // 对应 Face.h 中的定义
    Vector vectorE;                 ///< 正交分量 (沿着 d 方向)
    Vector vectorT;                 ///< 非正交修正分量 (T = A - E)
    double interpolationCoef;       ///< 线性插值系数//以替代

    // =========================================================
    // 4. 构造与行为
    // =========================================================
    FractureEdge_2D();
    FractureEdge_2D(int _id, int n1, int n2);

    /**
     * @brief 计算基础几何 (依赖相邻单元)
     * @param allNodes 节点列表
     * @param owner    拥有的单元对象
     * @param neighbor 邻接的单元对象 (若是边界，传入 owner 自身或特殊处理)
     * @details
     * In-plane Normal = Normalize( EdgeVector x SurfaceNormal )
     * 其中 SurfaceNormal = Normalize( owner.normal + neighbor.normal )
     */
    void computeGeometry(const std::vector<Node>& allNodes,
        const FractureElement_2D& owner,
        const FractureElement_2D& neighbor);

    /**
    * @brief 计算 FVM 向量 (E, T) - 支持多种修正策略
    * @param Cp Owner 单元中心
    * @param Cn Neighbor 单元中心
    * @param method 非正交修正方法 (默认 Orthogonal)
    * @details 复用 Face 类中的数学逻辑，分解 A = E + T
    */
    void computeFVMVectors(const Vector& Cp, const Vector& Cn,
        NormalVectorCorrectionMethod method = NormalVectorCorrectionMethod::OrthogonalCorrection);
};
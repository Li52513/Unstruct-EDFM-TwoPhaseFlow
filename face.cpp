#include "Face.h"

Face::Face(int id, const std::vector<int>& nodeIDs, const std::vector<Vector>& nodeCoords)
    : id(id), nodeIDs(nodeIDs), nodeCoords(nodeCoords), length(0.0), ownerCell(-1), neighborCell(-1)
{
    computeGeometry();
}

void Face::computeGeometry()
{
    size_t n = nodeCoords.size(); //记录节点数量
    if (n == 2) 
    {
        // 2D 面（线段）
        midpoint = (nodeCoords[0] + nodeCoords[1]) / 2.0;
        Vector diff = nodeCoords[1] - nodeCoords[0];
        length = diff.Mag();
        normal = Vector(diff.m_y, -diff.m_x, 0.0);
        double normMag = normal.Mag();
        if (normMag > 1e-12)
            normal = normal / normMag;
    }
    else if (n == 3)
    {
        // 3D 面（三角形面片）
        const Vector& p0 = nodeCoords[0];
        const Vector& p1 = nodeCoords[1];
        const Vector& p2 = nodeCoords[2];

        // 面的几何中心（三角形重心）
        midpoint = (p0 + p1 + p2) / 3.0;

        // 计算法向量：normal = (p1 - p0) × (p2 - p0)
        Vector v1 = p1 - p0;
        Vector v2 = p2 - p0;
        normal = v1 & v2;  // 使用已重载的叉乘运算符 &
        double area = 0.5 * normal.Mag(); // 面积为叉积长度的一半
        length = area; // 3D 中将面积赋值给 length（可改为 area 字段）

        // 单位化法向量
        if (normal.Mag() > 1e-12)
            normal = normal / normal.Mag();
    }
}

void Face::computeFaceVectors(const Vector& Cp, const Vector& Cn, NormalVectorCorrectionMethod method)
{
    // 1) 计算从 owner 到 neighbor 的向量 d
    ownerToNeighbor = Cn - Cp;
    double d_norm = ownerToNeighbor.Mag();
    if (d_norm < 1e-12) 
    {
        // 重心重合或太近，直接置零
        vectorE = Vector(0.0, 0.0, 0.0);
        vectorT = Vector(0.0, 0.0, 0.0);
        return;
    }
    Vector ej = ownerToNeighbor / d_norm;  // 单位方向向量

    // 2) 构造面积矢量 A_j = normal * length
    Vector Aj = normal * length;
    
    if ((Aj * ej) < 0)
    {
		Aj = -Aj; //		// 确保法向量与单位向量同向
    }

    // 3) 依据不同算法分解 A_j = E_j + T_j
    double Aj_dot_ej = Aj * ej;  // 点积

    switch (method) 
    {
    case NormalVectorCorrectionMethod::MinimumCorrection:
        // E_j = (A_j·e_j) e_j
        vectorE = ej * Aj_dot_ej;
        break;

    case NormalVectorCorrectionMethod::OrthogonalCorrection:
        // E_j = |A_j| e_j, 其中 |A_j| = length (normal 已单位化)
        vectorE = ej * length;
        break;

    case NormalVectorCorrectionMethod::OverRelaxed:
    {
        // E_j = |A_j|^2 / (A_j·e_j) * e_j
        double Aj_norm = length;
        double factor = (Aj_norm * Aj_norm) / Aj_dot_ej;
        vectorE = ej * factor;
        break;
    }
    }

    // 4) 非正交分量
    vectorT = Aj - vectorE;
}

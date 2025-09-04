#include "Face.h"

Face::Face(int id, const std::vector<int>& nodeIDs, const std::vector<Vector>& nodeCoords)
    : id(id), FaceNodeIDs(nodeIDs), FaceNodeCoords(nodeCoords), length(0.0), ownerCell(-1), neighborCell(-1)
{
    computeGeometry();
}

void Face::computeGeometry()
{
    size_t n = FaceNodeCoords.size(); //记录节点数量
    if (n == 2) 
    {
        // 2D 面（线段）
        midpoint = (FaceNodeCoords[0] + FaceNodeCoords[1]) / 2.0;
        Vector diff = FaceNodeCoords[1] - FaceNodeCoords[0];
        length = diff.Mag();
		//cout << "Face " << id << " 2D 面长度: " << length << endl;
        normal = Vector(diff.m_y, -diff.m_x, 0.0);
		//cout << "Face " << id << " 2D 面法向量: ("
			//<< normal.m_x << ", " << normal.m_y << ", " << normal.m_z << ")" << endl;
        double normMag = normal.Mag();
        if (normMag > 1e-12)
			normal = normal / normMag; // 单位化法向量

        if (!FaceNodeCoords.empty())
        {
			// 计算包围盒
			Vector minCoord = FaceNodeCoords[0];
			Vector maxCoord = FaceNodeCoords[0];
			for (const auto& coord : FaceNodeCoords)
			{
				if (coord.m_x < minCoord.m_x) minCoord.m_x = coord.m_x;
				if (coord.m_y < minCoord.m_y) minCoord.m_y = coord.m_y;
				if (coord.m_z < minCoord.m_z) minCoord.m_z = coord.m_z;
				if (coord.m_x > maxCoord.m_x) maxCoord.m_x = coord.m_x;
				if (coord.m_y > maxCoord.m_y) maxCoord.m_y = coord.m_y;
				if (coord.m_z > maxCoord.m_z) maxCoord.m_z = coord.m_z;
			}
			boundingBox = AABB(minCoord, maxCoord);
        }

    }
    else if (n == 3)
    {
        // 3D 面（三角形面片）
        const Vector& p0 = FaceNodeCoords[0];
        const Vector& p1 = FaceNodeCoords[1];
        const Vector& p2 = FaceNodeCoords[2];

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
   /* cout << "Face " << id << " 的面积矢量 A_j: (" 
         << Aj.m_x << ", " << Aj.m_y << ", " << Aj.m_z << ") 长度: " 
         << length << endl;*/

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


/**
 * @brief 处理边界面时的分解 A_j = E_j + T_j
 *
 * @param Cp      ownerCell 的重心坐标
 * @param method  选择三种校正策略之一（目前未使用）
 * 边界面没邻居中心坐标 Cn，但做 E/T 分解只需面法向和模长即可
 * 用 rPF 保证 normal 始终朝外，避免节点序反向导致的“内法向”。
 */
void Face::computeFaceVectorsBoundary(const Vector& Cp,
    NormalVectorCorrectionMethod /*method*/)
{
    // 1) 依据 owner中心 -> 面心 的向量与当前 normal 的夹角，修正法向指向外侧
    Vector rPF = midpoint - Cp;                 // owner中心 -> 面心
    double sgn = (normal * rPF) >= 0.0 ? 1.0 : -1.0;
    Vector nhat = normal * sgn;                 // 外法向（单位向量）

    // 2) 面积矢量 A = |A| * n_hat；2D 情况下 |A| = 边长 length
    Vector Aj = nhat * length;

    // 3) 边界面分解：取正交主分量，非正交分量置零（T=0）
    vectorE = Aj;                               // = |A| n_hat
    vectorT = Vector(0.0, 0.0, 0.0);

    // 4) 给 ownerToNeighbor 一个有意义的单位方向（外法向），便于后续使用
    ownerToNeighbor = nhat;
}
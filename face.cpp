#include "Face.h"
#include "UserDefineVarType.h"

//构造函数
Face::Face(int id, const std::vector<int>& nodeIDs, const std::vector<Vector>& nodeCoords)
    : id(id), FaceNodeIDs(nodeIDs), FaceNodeCoords(nodeCoords), length(0.0), ownerCell(-1), neighborCell(-1)
{
    computeGeometry();
}

// =========================================================
// 【2D/3D】构建网格面的几何信息：单位法向量、面积/长度、几何中心、包围盒
// =========================================================
void Face::computeGeometry()
{
    size_t n = FaceNodeCoords.size(); //记录节点数量
    if (n < 2) return;

    //适用于2D和3D
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
    
    // --- Case 1: 2D 面 (线段) ---
    if (n == 2) 
    {
        // 2D 面（线段）
        midpoint = (FaceNodeCoords[0] + FaceNodeCoords[1]) / 2.0;   //中点计算
        Vector diff = FaceNodeCoords[1] - FaceNodeCoords[0];        
        length = diff.Mag();                                        //线段长度计算

        // 2D 法向：(dy, -dx)
        normal = Vector(diff.m_y, -diff.m_x, 0.0);
        double normMag = normal.Mag();
        if (normMag > 1e-12)
			normal = normal / normMag; // 单位化法向量   
    }

    // --- Case 2: 3D 面 (三角形/四边形/多边形) ---
    else
    {
        // 1. 计算几何中心 (Simple Average)
        midpoint = Vector(0.0, 0.0, 0.0);
        for (const auto& p : FaceNodeCoords) {
            midpoint = midpoint + p;
        }
        midpoint = midpoint / static_cast<double>(n);

        // 2. 计算面积矢量 (Area Vector) 和 面积 (length)
        // 方法：将多边形三角剖分为扇形 (Fan triangulation)，求叉积之和
        // 假设多边形是凸的或节点顺序良好的简单多边形
        Vector areaVec(0.0, 0.0, 0.0);
        const Vector& p0 = FaceNodeCoords[0];

        for (size_t i = 1; i < n - 1; ++i) {
            Vector v1 = FaceNodeCoords[i] - p0;
            Vector v2 = FaceNodeCoords[i + 1] - p0;
            // 注意：使用 & 作为叉乘运算符
            areaVec = areaVec + (v1 & v2);
        }
        areaVec = areaVec * 0.5;

        // 3. 设置属性
        // 在 3D 模式下，成员变量 'length' 存储的是【面积】
        length = areaVec.Mag();

        if (length > 1e-12) {
            normal = areaVec / length; // 单位化法向
        }
        else {
            // 退化面处理
            normal = Vector(0.0, 0.0, 0.0);
        }

    }

    
}
// =========================================================
// 【2D/3D】计算面法矢量、正交投影矢量、非正交修正矢量
// =========================================================

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
    double Aj_dot_ej = Aj * ej;  // 点积
    
    if ((Aj * ej) < 0)
    {
		Aj = -Aj; //		// 确保法向量与单位向量同向
        Aj_dot_ej = -Aj_dot_ej;
    }
    // —— 关键容错：近正交（Aj·ej≈0）时的保护 —— 
    const double eps = 1e-14;   // 可按你的坐标尺度调整
    if (Aj_dot_ej < eps) Aj_dot_ej = eps;

    // 3) 依据不同算法分解 A_j = E_j + T_j
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
   
    //5） 正交插值权重gamma 投影到 e_j（0..1 之间）
    double D = ownerToNeighbor * ej;                        // 应等于 d_norm
    double s = (midpoint - Cp) * ej;                        // owner 到 面心 在 e_j 上的投影
    double gamma = (std::abs(D) > 1e-14) ? (s / D) : 0.5;   // 容错
    if (gamma < 0.0) gamma = 0.0;
    if (gamma > 1.0) gamma = 1.0;
    f_linearInterpolationCoef = gamma;
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
    Vector rPF = midpoint - Cp;
    double sgn = (normal * rPF) >= 0.0 ? 1.0 : -1.0;

    // ↓↓↓ 把外指法向写回，并确保为单位向量 ↓↓↓
    normal = normal * sgn;
    double nmag2 = normal.m_x * normal.m_x + normal.m_y * normal.m_y + normal.m_z * normal.m_z;
    if (nmag2 > 0.0) {
        double invn = 1.0 / std::sqrt(nmag2);
        normal = normal * invn;
    }

    // 2) 面积矢量 A = |A| * n_hat；2D 情况下 |A| = 边长 length（3D 为面积）
    Vector Aj = normal * length;

    // 3) 边界面分解
    vectorE = Aj;
    vectorT = Vector(0.0, 0.0, 0.0);

    // 4) owner→外部 的方向（与 normal 一致）
    ownerToNeighbor = rPF;
    f_linearInterpolationCoef = 1.0;
}
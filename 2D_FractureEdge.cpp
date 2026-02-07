#include "2D_FractureEdge.h"
#include <iostream>

FractureEdge_2D::FractureEdge_2D()
    : id(-1), ownerCellID(-1), neighborCellID(-1), parentFractureID(-1), ownerCellLocalIndex(-1), neighborCellLocalIndex(-1), length(0.0), interpolationCoef(0.5)
{
    nodeIndices[0] = -1;
    nodeIndices[1] = -1;
}

FractureEdge_2D::FractureEdge_2D(int _id, int n1, int n2)
    : id(_id), ownerCellID(-1), neighborCellID(-1), parentFractureID(-1), ownerCellLocalIndex(-1), neighborCellLocalIndex(-1), length(0.0), interpolationCoef(0.5)
{
    nodeIndices[0] = n1;
    nodeIndices[1] = n2;
}

void FractureEdge_2D::computeGeometry(const std::vector<Node>& allNodes,
    const FractureElement_2D& owner,
    const FractureElement_2D& neighbor)
{
    const Vector& p1 = allNodes[nodeIndices[0]].coord;
    const Vector& p2 = allNodes[nodeIndices[1]].coord;

    // 1. 长度与中点
    Vector edgeVec = p2 - p1;
    length = edgeVec.Mag();
    midpoint = (p1 + p2) * 0.5;

    if (length > 1e-14) {

        // 取 Owner 和 Neighbor 法向的平均值，代表边所在位置的切平面法向
        Vector surfaceNormal = (owner.normal + neighbor.normal);
        double smag = surfaceNormal.Mag();

        if (smag < 1e-14) {
            // 极端情况：两单元法向相反（折叠）？或者都为0
            surfaceNormal = owner.normal;
        }
        else {
            surfaceNormal = surfaceNormal / smag;
        }

        Vector tangent = edgeVec / length;

        // 计算面内法向 (垂直于边，切于表面)
        normal = tangent & surfaceNormal;

        // 再次归一化
        double nMag = normal.Mag();
        if (nMag > 1e-14) normal = normal / nMag;
    }
    else {
        normal = Vector(0, 0, 0);
    }
}

void FractureEdge_2D::computeFVMVectors(const Vector& Cp, const Vector& Cn, NormalVectorCorrectionMethod method)
{
    // 参考 face.cpp 中的 computeFaceVectors 逻辑
    // A = area * normal (这里 area 就是 length * aperture? 不，通常 FVM 中 Face Area 就是截面积)
    // 但在 2D 裂缝流形中，我们将积分降维。
    // 这里我们定义 Area Vector A = length * normal
    // (注意：后续在构建方程时，需要乘以裂缝开度 aperture)

    Vector A = normal * length;

    // 距离向量 d = Cn - Cp
    Vector d = Cn - Cp; // Owner -> Neighbor 向量
    double d_magSq = d.Mag2(); // |d|^2

    // 0. 重合异常处理
    if (d_magSq < 1e-20) {
        vectorE = A;
        vectorT = Vector(0, 0, 0);
        interpolationCoef = 0.5;
        return;
    }

    // 1. 检查并校正法向方向
     // FVM 约定：法向必须从 Owner 指向 Neighbor (即 A 与 d 的夹角应 < 90度)
    if ((A * d) < 0.0) {
        normal = normal * -1.0;
        A = A * -1.0;
    }

    // 2. 准备单位向量 e_d
    double d_mag = std::sqrt(d_magSq);
    Vector e_d = d / d_mag;
    double A_dot_d = A * d;

    // 3. 根据算法分解 A = E + T
    switch (method)
    {
    case NormalVectorCorrectionMethod::MinimumCorrection:
    {
        // E = (A . e_d) * e_d
        // 这种方法使得 |E| 等于 A 在 d 上的投影长度
        // T = A - E 与 d 垂直（正交分解）
        // 优点：当网格极其扭曲时，能保证稳定性
        double projection = A_dot_d / d_mag;
        vectorE = e_d * projection;
    }
    break;

    case NormalVectorCorrectionMethod::OrthogonalCorrection:
    {
        // E = |A| * e_d
        // 这种方法假设理想情况下 A 平行于 d。
        // 当网格正交时，E = A, T = 0。
        // 这是最常用的默认方法。
        vectorE = e_d * length;
    }
    break;

    case NormalVectorCorrectionMethod::OverRelaxed:
    {
        // E = (|A|^2 / (A . e_d)) * e_d
        // 这种方法也被称为 "Over-relaxed approach"
        // 当 A 与 d 夹角较大时，E 的模长会显著增大，增强扩散项的作用
        // 注意防止除零
        if (std::abs(A_dot_d) > 1e-14) {
            double factor = (length * length) / A_dot_d; // |A|^2 = length^2
            vectorE = d * factor; // 注意: d = e_d * |d|, 这里逻辑可能有误，请参照 Face.cpp

            // 修正参照 Face.cpp 的逻辑: vectorE = ej * factor; 其中 ej 是单位向量
            // factor = (Aj_norm * Aj_norm) / Aj_dot_ej; 
            // A . e_d = A_dot_d / d_mag
            double Aj_dot_ej = A_dot_d / d_mag;
            double factor_corrected = (length * length) / Aj_dot_ej;
            vectorE = e_d * factor_corrected;
        }
        else {
            // 退回到正交修正
            vectorE = e_d * length;
        }
    }
    break;
    }

    // 4. 计算非正交分量 T
    vectorT = A - vectorE;

    // 5. 计算插值系数 (Geometrical Interpolation Coefficient)
    // 用于计算面上的值： phi_f = g * phi_N + (1-g) * phi_P
    // 简单投影法：计算面心到 P 点在 d 方向上的投影占比
    Vector Pf = midpoint - Cp;
    double proj_Pf = Pf * e_d; // P->Face 在 d 上的投影长度

    interpolationCoef = proj_Pf / d_mag;

    // 限制范围 [0, 1] 保证数值稳定性
    if (interpolationCoef < 0.0) interpolationCoef = 0.0;
    if (interpolationCoef > 1.0) interpolationCoef = 1.0;
}

void FractureEdge_2D::computeBoundaryGeometry(const std::vector<Node>& allNodes,
    const FractureElement_2D& owner)
{
    const Vector& p1 = allNodes[nodeIndices[0]].coord;
    const Vector& p2 = allNodes[nodeIndices[1]].coord;

    // 1. 长度与中点
    Vector edgeVec = p2 - p1;
    length = edgeVec.Mag();
    midpoint = (p1 + p2) * 0.5;

    // 2. 计算平面内向外法向量 (In-plane Outward Normal)
    if (length > 1e-14) {
        // 归一化边向量
        Vector tangent = edgeVec / length;
        // 初始猜测：假设是逆时针，法向 = 切向 x 单元面法向
        normal = tangent & owner.normal;
        // 归一化
        double nMag = normal.Mag();
        if (nMag > 1e-14) {
            normal = normal / nMag;
        }
        else {
            normal = Vector(0, 0, 0);
        }

        // B. [核心修复] 方向自洽性校验 (Eliminate Winding Order Assumption)
        // 计算从单元质心指向边中点的向量
        Vector centroidToMid = midpoint - owner.centroid;

        // 判断：如果法向与(质心->中点)向量背道而驰 (点积 < 0)，说明指回了内部
        // 必须翻转为指向外部
        if ((normal * centroidToMid) < 0.0) {
            normal = normal * -1.0;
        }
    }
    else {
        normal = Vector(0, 0, 0);
    }
    // 3. FVM 插值系数处理
    // 对于边界，通常设为 1.0 (完全由 Owner 决定) 或 0.5 (视具体边界格式定)
    // 这里保持默认或设为特定值
    interpolationCoef = 0.0; // 边界边通常不需要插值系数，或者设为0表示完全偏向Owner

}




   
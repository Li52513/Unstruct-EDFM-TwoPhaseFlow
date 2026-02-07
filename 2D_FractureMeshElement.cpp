#include "2D_FractureMeshElement.h"

#include "Node.h"

#include <iostream>
#include <cmath>
#include <limits> // 用于 numeric_limits

// =========================================================
// 构造函数实现
// =========================================================

FractureElement_2D::FractureElement_2D()
    : id(-1), localIndex(-1), parentFractureID(-1), solverIndex(-1), area(0.0), centroid(0.0, 0.0, 0.0)
{
}

FractureElement_2D::FractureElement_2D(int _id, const std::vector<int>& _nodes)
    : id(_id), localIndex(-1), parentFractureID(-1), solverIndex(-1), nodeIndices(_nodes), area(0.0), centroid(0.0, 0.0, 0.0)
{
}

// =========================================================
// 成员行为实现
// =========================================================
void FractureElement_2D::computeGeometry(const std::vector<Node>& allFracNodes)
{
    if (nodeIndices.size() < 3) return;

    // 1. 计算质心与包围盒
    Vector sum(0, 0, 0);
    double maxVal = std::numeric_limits<double>::max();
    double minVal = std::numeric_limits<double>::lowest();
    Vector minCoords(maxVal, maxVal, maxVal); // 使用 minCoords 避免混淆
    Vector maxCoords(minVal, minVal, minVal);

    std::vector<Vector> pts;
    pts.reserve(nodeIndices.size());

    for (int idx : nodeIndices) {
        if (idx < 0 || idx >= (int)allFracNodes.size()) continue;
        const Vector& p = allFracNodes[idx].coord;
        pts.push_back(p);
        sum = sum + p;

        if (p.m_x < minCoords.m_x) minCoords.m_x = p.m_x;
        if (p.m_y < minCoords.m_y) minCoords.m_y = p.m_y;
        if (p.m_z < minCoords.m_z) minCoords.m_z = p.m_z;

        if (p.m_x > maxCoords.m_x) maxCoords.m_x = p.m_x;
        if (p.m_y > maxCoords.m_y) maxCoords.m_y = p.m_y;
        if (p.m_z > maxCoords.m_z) maxCoords.m_z = p.m_z;
    }
    centroid = sum / static_cast<double>(pts.size());

    // [修正] 使用 AABB 构造函数
    boundingBox = AABB(minCoords, maxCoords);

    // 2. 计算法向与面积 (针对扭曲四边形)
    if (pts.size() == 3) {
        Vector v1 = pts[1] - pts[0];
        Vector v2 = pts[2] - pts[0];
        Vector n = v1 & v2;
        area = n.Mag() * 0.5;
        normal = (area > 1e-14) ? n / n.Mag() : Vector(0, 0, 1);
    }
    else if (pts.size() >= 4) {
        Vector n1 = (pts[1] - pts[0]) & (pts[2] - pts[0]);
        double a1 = n1.Mag() * 0.5;

        Vector n2 = (pts[2] - pts[0]) & (pts[3] - pts[0]);
        double a2 = n2.Mag() * 0.5;

        area = a1 + a2;
        if (area > 1e-14) {
            // 面积加权平均法向
            Vector n_sum = n1 + n2;

            // 手动计算归一化
            double n_mag = n_sum.Mag();

            if (n_mag > 1e-12) {
                normal.m_x = n_sum.m_x / n_mag;
                normal.m_y = n_sum.m_y / n_mag;
                normal.m_z = n_sum.m_z / n_mag;
            }
            else {
                // 如果和向量的模太小，使用默认法向或第一个三角形的法向
                double n1_mag = n1.Mag();
                if (n1_mag > 1e-12) {
                    normal.m_x = n1.m_x / n1_mag;
                    normal.m_y = n1.m_y / n1_mag;
                    normal.m_z = n1.m_z / n1_mag;
                }
                else {
                    // 最后保底：设为单位Z轴
                    normal = Vector(0, 0, 1);
                }
            }
        }
        else {
            // 面积太小，无法计算有效法向
            normal = Vector(0, 0, 0);
        }
    }
}
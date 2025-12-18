#pragma once
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

#include "Node.h"
#include "UserDefineVarType.h"

// Internal helper functions for fracture geometry operations.
// Kept header-only to avoid changing build configuration.
namespace FractureGeom
{
    // ===== 2D helpers =====
    // ==================== 2D基础几何函数 ====================

    // 功能：计算2D向量的叉积 (仅考虑x,y分量)
    // 参数：a,b - 输入向量
    // 返回：标量叉积值 (a.x*b.y - a.y*b.x)
    inline double cross2D(const Vector& a, const Vector& b)
    {
        return a.m_x * b.m_y - a.m_y * b.m_x;
    }

    // 功能：计算三角形面积
    // 参数：A,B,C - 三角形三个顶点
    // 返回：三角形面积 (取绝对值)
    inline double triArea(const Vector& A, const Vector& B, const Vector& C)
    {
        return 0.5 * std::abs(cross2D(B - A, C - A));
    }

    // 功能：计算三角形质心
    // 参数：A,B,C - 三角形三个顶点
    // 返回：质心坐标 (三个顶点的平均值)
    inline Vector triCentroid(const Vector& A, const Vector& B, const Vector& C)
    {
        return (A + B + C) / 3.0;
    }

    inline double pointToSegmentDistance_local(
        const Vector& p,
        const Vector& a,
        const Vector& b)
    {
        Vector ab = b - a;
        double L2 = ab.m_x * ab.m_x + ab.m_y * ab.m_y + ab.m_z * ab.m_z;
        if (L2 <= 1e-30) {
            double dx = p.m_x - a.m_x, dy = p.m_y - a.m_y, dz = p.m_z - a.m_z;
            return std::sqrt(dx * dx + dy * dy + dz * dz);
        }
        double t = ((p - a) * ab) / L2;
        if (t < 0.0) t = 0.0; else if (t > 1.0) t = 1.0;
        Vector proj = a + ab * t;
        double dx = p.m_x - proj.m_x, dy = p.m_y - proj.m_y, dz = p.m_z - proj.m_z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    // 2D segment–segment intersection test with tolerance.
    inline bool segSegIntersect2D(
        const Vector& p,
        const Vector& p2,
        const Vector& q,
        const Vector& q2,
        Vector& out,
        double eps = 1e-12)
    {
        Vector r = p2 - p;
        Vector s = q2 - q;
        double denom = cross2D(r, s);
        if (std::fabs(denom) < eps) return false; // parallel or nearly parallel

        Vector qp = q - p;
        double t = cross2D(qp, s) / denom;
        double u = cross2D(qp, r) / denom;
        if (t < -eps || t > 1.0 + eps || u < -eps || u > 1.0 + eps) return false;

        t = std::min(1.0, std::max(0.0, t));
        out = p + r * t;
        return true;
    }

    // Count intersections between a segment and polygon edges (deduplicated).
    inline int countIntersectionsWithPolygon(
        const std::vector<Vector>& poly,
        const Vector& s0,
        const Vector& s1,
        std::vector<Vector>& hits,
        double eps = 1e-12)
    {
        hits.clear();
        if (poly.size() < 3) return 0;
        const size_t n = poly.size();
        for (size_t i = 0; i < n; ++i) {
            const Vector& a = poly[i];
            const Vector& b = poly[(i + 1) % n];
            Vector x{};
            if (segSegIntersect2D(s0, s1, a, b, x, eps)) {
                bool dup = false;
                for (const auto& h : hits) {
                    double dx = h.m_x - x.m_x, dy = h.m_y - x.m_y;
                    if (dx * dx + dy * dy < 1e-20) { dup = true; break; }
                }
                if (!dup) hits.push_back(x);
            }
        }
        return (int)hits.size();
    }

    // Split polygon by line s0->s1 (Sutherland–Hodgman style), producing two polygons.
    inline void splitPolygonByLine(
        const std::vector<Vector>& poly,
        const Vector& s0,
        const Vector& s1,
        std::vector<Vector>& pos, // cross>=0 side
        std::vector<Vector>& neg, // cross<=0 side
        double eps = 1e-12)
    {
        pos.clear(); neg.clear();
        const size_t n = poly.size();
        if (n < 3) return;

        Vector nvec = s1 - s0;
        auto sideVal = [&](const Vector& p)->double { return cross2D(nvec, p - s0); };

        for (size_t i = 0; i < n; ++i) {
            const Vector& A = poly[i];
            const Vector& B = poly[(i + 1) % n];
            double sa = sideVal(A);
            double sb = sideVal(B);

            bool AinPos = (sa > eps);
            bool AinNeg = (sa < -eps);
            bool Aon = (!AinPos && !AinNeg);

            bool BinPos = (sb > eps);
            bool BinNeg = (sb < -eps);

            if (AinPos || Aon) pos.push_back(A);
            if (AinNeg || Aon) neg.push_back(A);

            if ((AinPos && BinNeg) || (AinNeg && BinPos)) {
                Vector X{};
                segSegIntersect2D(A, B, s0, s1, X, eps);
                pos.push_back(X);
                neg.push_back(X);
            }
        }

        auto dedup = [](std::vector<Vector>& P) {
            if (P.size() < 3) return;
            std::vector<Vector> out; out.reserve(P.size());
            for (size_t i = 0; i < P.size(); ++i) {
                const auto& cur = P[i];
                const auto& prv = P[(i + P.size() - 1) % P.size()];
                double dx = cur.m_x - prv.m_x, dy = cur.m_y - prv.m_y;
                if (dx * dx + dy * dy > 1e-24) out.push_back(cur);
            }
            P.swap(out);
            };
        dedup(pos);
        dedup(neg);
    }

    // Triangle 3-point Gauss integral for ∫_△ d(x) dS (returns integral, not mean).
    inline double triGauss3Integral(
        const Vector& V1,
        const Vector& V2,
        const Vector& V3,
        const Vector& segStart,
        const Vector& segEnd)
    {
        double area = triArea(V1, V2, V3);
        if (area <= 0.0) return 0.0;

        Vector P1 = V1 * (1.0 / 6.0) + V2 * (1.0 / 6.0) + V3 * (2.0 / 3.0);
        Vector P2 = V1 * (1.0 / 6.0) + V2 * (2.0 / 3.0) + V3 * (1.0 / 6.0);
        Vector P3 = V1 * (2.0 / 3.0) + V2 * (1.0 / 6.0) + V3 * (1.0 / 6.0);

        double d1 = pointToSegmentDistance_local(P1, segStart, segEnd);
        double d2 = pointToSegmentDistance_local(P2, segStart, segEnd);
        double d3 = pointToSegmentDistance_local(P3, segStart, segEnd);

        return (area / 3.0) * (d1 + d2 + d3);
    }

    inline Vector polygonCentroid(const std::vector<Vector>& P, double& areaAbs)
    {
        const size_t n = P.size();
        double A2 = 0.0, Cx = 0.0, Cy = 0.0;
        for (size_t i = 0; i < n; ++i) {
            const auto& p = P[i];
            const auto& q = P[(i + 1) % n];
            double cr = p.m_x * q.m_y - q.m_x * p.m_y;
            A2 += cr;
            Cx += (p.m_x + q.m_x) * cr;
            Cy += (p.m_y + q.m_y) * cr;
        }
        double A = 0.5 * A2;
        areaAbs = std::abs(A);
        if (std::abs(A) < 1e-16) {
            double mx = 0.0, my = 0.0;
            for (auto& p : P) { mx += p.m_x; my += p.m_y; }
            return Vector{ mx / n, my / n, 0.0 };
        }
        return Vector{ Cx / (3.0 * A2), Cy / (3.0 * A2), 0.0 };
    }

    inline double fanGaussAverage(
        const std::vector<Vector>& poly,
        const Vector& segStart,
        const Vector& segEnd)
    {
        if (poly.size() < 3) return 0.0;
        double polyArea = 0.0;
        Vector C = polygonCentroid(poly, polyArea);
        if (polyArea <= 0.0) return 0.0;

        double integral = 0.0;
        const size_t n = poly.size();
        for (size_t i = 0; i < n; ++i) {
            const Vector& A = poly[i];
            const Vector& B = poly[(i + 1) % n];
            integral += triGauss3Integral(C, A, B, segStart, segEnd);
        }
        return integral / polyArea;
    }

    inline bool pointOnSegment2D(const Vector& P, const Vector& A, const Vector& B, double eps = 1e-12)
    {
        Vector AP = P - A, AB = B - A;
        double cross = std::fabs(AP.m_x * AB.m_y - AP.m_y * AB.m_x);
        if (cross > eps) return false;
        double dot = AP.m_x * AB.m_x + AP.m_y * AB.m_y;
        if (dot < -eps) return false;
        double ab2 = AB.m_x * AB.m_x + AB.m_y * AB.m_y;
        if (dot - ab2 > eps) return false;
        return true;
    }

    // 2D polygon point-in test (odd-even rule), including boundary.
    inline bool pointInPolygon2D(
        const Vector& P,
        const std::vector<int>& polyNodeIDs,
        const std::unordered_map<int, Node>& nodes,
        double eps = 1e-12)
    {
        const int n = (int)polyNodeIDs.size();
        if (n < 3) return false;

        for (int i = 0; i < n; ++i) {
            const Vector& A = nodes.at(polyNodeIDs[i]).coord;
            const Vector& B = nodes.at(polyNodeIDs[(i + 1) % n]).coord;
            if (pointOnSegment2D(P, A, B, eps)) return true;
        }

        bool inside = false;
        for (int i = 0, j = n - 1; i < n; j = i++) {
            const Vector& Pi = nodes.at(polyNodeIDs[i]).coord;
            const Vector& Pj = nodes.at(polyNodeIDs[j]).coord;
            bool intersect = ((Pi.m_y > P.m_y) != (Pj.m_y > P.m_y));
            if (intersect) {
                double x_cross = Pj.m_x + (Pi.m_x - Pj.m_x)
                    * ((P.m_y - Pj.m_y) / (Pi.m_y - Pj.m_y + 1e-30));
                if (x_cross > P.m_x - eps) inside = !inside;
            }
        }
        return inside;
    }

 
    inline bool pointToSegmentDistance(const Vector& p, const Vector& s, const Vector& e)
    {
        Vector se = e - s;
        double segLenSq = se.m_x * se.m_x + se.m_y * se.m_y;
        if (segLenSq < 1e-12) return (p - s).Mag();
        double t = ((p - s) * se) / segLenSq;
        t = max(0.0, std::min(1.0, t));
        Vector projection = s + se * t;
        return (p - projection).Mag();
    }

    inline bool pointInTriangle(const Vector& p, const Vector& a, const Vector& b, const Vector& c)
    {
        Vector v0 = c - a, v1 = b - a, v2 = p - a;
        double dot00 = v0 * v0, dot01 = v0 * v1, dot02 = v0 * v2;
        double dot11 = v1 * v1, dot12 = v1 * v2;
        double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
        double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        double v = (dot00 * dot12 - dot01 * dot02) * invDenom;
        return (u > 1e-6) && (v > 1e-6) && ((u + v) < 1 - 1e-6);
    }
} // namespace FractureGeom


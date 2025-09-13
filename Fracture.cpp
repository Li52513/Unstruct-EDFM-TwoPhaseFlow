#include "Fracture.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include "Fluid.h"
#include "Matrix.h"
#include "Mesh.h"
#include <unordered_map>
#include <numeric>


Fracture::Fracture(const Vector& s, const Vector& e)  
   : start(s), end(e), id(0)  
{  

}

// ====== 小工具：2D 叉积z分量、三角形面积与重心、多边形重心 ======
static inline double cross2D(const Vector& a, const Vector& b) {
    return a.m_x * b.m_y - a.m_y * b.m_x;
}

static inline double triArea(const Vector& A, const Vector& B, const Vector& C) {
    return 0.5 * std::abs(cross2D(B - A, C - A));
}

static inline Vector triCentroid(const Vector& A, const Vector& B, const Vector& C) {
    return (A + B + C) / 3.0;
}

static inline Vector lerp(const Vector& a, const Vector& b, double t) {
    return a * (1.0 - t) + b * t;
}

static inline double pointToSegmentDistance_local(const Vector& p,
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


// 线段-线段相交（严格在内部相交，端点相交视作相交；平行/共线返回 false）
static bool segSegIntersect2D(const Vector& p, const Vector& p2,
    const Vector& q, const Vector& q2,
    Vector& out, double eps = 1e-12)
{
    Vector r = p2 - p;
    Vector s = q2 - q;
    double denom = cross2D(r, s);
    if (std::fabs(denom) < eps) return false; // 平行或几乎平行

    Vector qp = q - p;
    double t = cross2D(qp, s) / denom;
    double u = cross2D(qp, r) / denom;
    if (t < -eps || t > 1.0 + eps || u < -eps || u > 1.0 + eps) return false;

    // 夹在两段范围内
    t = std::min(1.0, std::max(0.0, t));
    out = p + r * t;
    return true;
}

// 统计裂缝段与单元多边形边的交点数
static int countIntersectionsWithPolygon(const std::vector<Vector>& poly,
    const Vector& s0, const Vector& s1,
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
            // 去重：和已收集交点相距很近就跳过
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

// 按直线 s0->s1 劈分多边形（Sutherland–Hodgman 风格，生成两侧多边形）
static void splitPolygonByLine(const std::vector<Vector>& poly,
    const Vector& s0, const Vector& s1,
    std::vector<Vector>& pos, // 叉积>=0 一侧
    std::vector<Vector>& neg, // 叉积<=0 另一侧
    double eps = 1e-12)
{
    pos.clear(); neg.clear();
    const size_t n = poly.size();
    if (n < 3) return;

    Vector nvec = s1 - s0; // 直线方向
    auto sideVal = [&](const Vector& p)->double {
        return cross2D(nvec, p - s0); // >0 在“左侧”，<0 在“右侧”
        };

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
        bool Bon = (!BinPos && !BinNeg);

        // 把 A 放到对应侧（边界点放两侧，以保守封闭）
        if (AinPos || Aon) pos.push_back(A);
        if (AinNeg || Aon) neg.push_back(A);

        // 边 AB 是否跨越直线？
        if ((AinPos && BinNeg) || (AinNeg && BinPos)) {
            Vector X{};
            segSegIntersect2D(A, B, s0, s1, X, eps); // 与无限直线的交点
            pos.push_back(X);
            neg.push_back(X);
        }
    }

    // 可能出现重复点/共线退化；去掉连续重复
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
    dedup(pos); dedup(neg);
}


// 三角形三点高斯对 ∫_△ d(x) dS 的数值近似（返回“积分值”，不是均值）
// 顶点：V1,V2,V3
static double triGauss3Integral(const Vector& V1, const Vector& V2, const Vector& V3,
    const Vector& segStart, const Vector& segEnd)
{
    // 三点：barycentric (1/6,1/6,2/3) 的全排列；权重 = Area/3
    auto area = triArea(V1, V2, V3);
    if (area <= 0) return 0.0;

    // 三个采样点
    Vector P1 = V1 * (1.0 / 6.0) + V2 * (1.0 / 6.0) + V3 * (2.0 / 3.0);
    Vector P2 = V1 * (1.0 / 6.0) + V2 * (2.0 / 3.0) + V3 * (1.0 / 6.0);
    Vector P3 = V1 * (2.0 / 3.0) + V2 * (1.0 / 6.0) + V3 * (1.0 / 6.0);

    double d1 = pointToSegmentDistance_local(P1, segStart, segEnd);
    double d2 = pointToSegmentDistance_local(P2, segStart, segEnd);
    double d3 = pointToSegmentDistance_local(P3, segStart, segEnd);

    // ∫ ≈ (Area/3) * (d1 + d2 + d3)
    return (area / 3.0) * (d1 + d2 + d3);
}

//计算多边形的重心
static inline Vector polygonCentroid(const std::vector<Vector>& P, double& areaAbs) {
    const size_t n = P.size();
    double A2 = 0.0, Cx = 0.0, Cy = 0.0; // A2=2*Area(with sign)
    for (size_t i = 0; i < n; ++i) {
        const auto& p = P[i];
        const auto& q = P[(i + 1) % n];
        double cr = p.m_x * q.m_y - q.m_x * p.m_y;
        A2 += cr; Cx += (p.m_x + q.m_x) * cr; Cy += (p.m_y + q.m_y) * cr;
    }
    double A = 0.5 * A2; areaAbs = std::abs(A);
    if (std::abs(A) < 1e-16) { // 退化：用顶点均值
        double mx = 0, my = 0;
        for (auto& p : P) { mx += p.m_x; my += p.m_y; }
        return Vector{ mx / n, my / n, 0.0 };
    }
    return Vector{ Cx / (3.0 * A2), Cy / (3.0 * A2), 0.0 };
}

// 对一个简单多边形做“中心扇形 + 三点高斯”，返回“平均距离”
static double fanGaussAverage(const std::vector<Vector>& poly,
    const Vector& segStart, const Vector& segEnd)
{
    if (poly.size() < 3) return 0.0;
    double polyArea = 0.0;
    Vector C = polygonCentroid(poly, polyArea);
    if (polyArea <= 0) return 0.0;

    double integral = 0.0;
    const size_t n = poly.size();
    for (size_t i = 0; i < n; ++i) {
        const Vector& A = poly[i];
        const Vector& B = poly[(i + 1) % n];
        integral += triGauss3Integral(C, A, B, segStart, segEnd);
    }
    return integral / polyArea;
}


// ====== 新增：穿越感知 + 三点高斯的平均距离 ======
double Fracture::computeCrossAwareAverageDistance(
    const Cell& cell,
    const std::unordered_map<int, Node>& meshNodes,
    const Vector& segStart, const Vector& segEnd)
{
    // 取单元多边形
    std::vector<Vector> poly; poly.reserve(cell.CellNodeIDs.size());
    for (int nid : cell.CellNodeIDs) poly.push_back(meshNodes.at(nid).coord);
    if (poly.size() < 3) {
        // 退化：回退到中心点距离
        return pointToSegmentDistance(cell.center, segStart, segEnd);
    }

    // 判断是否“穿越”
    std::vector<Vector> hits;
    int nHit = countIntersectionsWithPolygon(poly, segStart, segEnd, hits);
    if (nHit >= 2) {
        // 按直线把多边形一分为二，两个子多边形分别 fan + Gauss，再面积加权
        std::vector<Vector> Ppos, Pneg;
        splitPolygonByLine(poly, segStart, segEnd, Ppos, Pneg);

        double areaPos = 0.0, areaNeg = 0.0;
        if (Ppos.size() >= 3) { Vector c = polygonCentroid(Ppos, areaPos); (void)c; }
        if (Pneg.size() >= 3) { Vector c = polygonCentroid(Pneg, areaNeg); (void)c; }
        double A = areaPos + areaNeg;
        if (A <= 0) {
            // 极端退化：回退到原扇形 + 高斯
            return fanGaussAverage(poly, segStart, segEnd);
        }

        double favPos = (areaPos > 0 && Ppos.size() >= 3) ? fanGaussAverage(Ppos, segStart, segEnd) : 0.0;
        double favNeg = (areaNeg > 0 && Pneg.size() >= 3) ? fanGaussAverage(Pneg, segStart, segEnd) : 0.0;

        return (favPos * areaPos + favNeg * areaNeg) / A;
    }
    else {
        // 不穿越：直接对原多边形做扇形 + 三点高斯
        return fanGaussAverage(poly, segStart, segEnd);
    }
}


double Fracture::computeAreaWeightedDistance(const Cell& cell,
    const std::unordered_map<int, Node>& meshNodes,
    const Vector& segStart, const Vector& segEnd)
{
    // 取多边形顶点（按 cell 拓扑顺序）
    std::vector<Vector> poly;
    poly.reserve(cell.CellNodeIDs.size());
    for (int nid : cell.CellNodeIDs) poly.push_back(meshNodes.at(nid).coord);

    if (poly.size() < 3) {
        // 极端/退化：回退到中心距离
        return pointToSegmentDistance(cell.center, segStart, segEnd);
    }

    double polyArea = 0.0;
	Vector C = polygonCentroid(poly, polyArea); // 多边形重心

    double sum = 0.0, sumA = 0.0;
    const size_t n = poly.size();
    for (size_t i = 0; i < n; ++i) {
        const Vector& A = poly[i];
        const Vector& B = poly[(i + 1) % n];
        double Ai = triArea(C, A, B);
        if (Ai <= 0) continue;
        Vector Gi = triCentroid(C, A, B);
        double di = pointToSegmentDistance(Gi, segStart, segEnd);
        sum += di * Ai;
        sumA += Ai;
    }

    if (sumA <= 0) {
        // 再次回退
        return pointToSegmentDistance(cell.center, segStart, segEnd);
    }
    return sum / sumA;
}



void Fracture::computeAABB()
{
    Vector minCoord, maxCoord;
	minCoord.m_x = std::min(start.m_x, end.m_x);
	minCoord.m_y = std::min(start.m_y, end.m_y);
	minCoord.m_z = std::min(start.m_z, end.m_z);

	maxCoord.m_x = std::max(start.m_x, end.m_x);
	maxCoord.m_y = std::max(start.m_y, end.m_y);
	maxCoord.m_z = std::max(start.m_z, end.m_z);
	fractureBox = AABB(minCoord, maxCoord); // 计算包围盒
}


/*====利用基岩网格边界对裂缝进行网格划分并排序---基于lineSegmentIntersection子函数 *注：对于非结构化网格需要进行排序===*/
void Fracture::DetectFracturetoMeshFaceIntersections(const Mesh& mesh,const vector<Cell>& meshCells,const unordered_map<int, Node>& meshNodes,bool useAABBFilter)
{
    /* =========================================================
   0. 基本准备：裂缝自身 AABB 以及要遍历的面列表
   =========================================================*/

	AABB fractureBox(start, end); // 计算裂缝的包围盒
	vector<int> faceIDList; // 用于存储需要遍历的面 ID

    if (useAABBFilter)
    {
		faceIDList = mesh.getCandidateFacesFromBins(fractureBox); // 获取与裂缝包围盒相交的面 ID 列表
        candidateFaceCount_ = static_cast<int>(faceIDList.size());
    }
    else
    {
        const auto& faces = mesh.getFaces();
		faceIDList.resize(faces.size());
		iota(faceIDList.begin(), faceIDList.end(), 1); // 如果不使用 AABB 过滤，则遍历所有面
    }

	////打印候选面 ID 列表
 //     std::cout << "[Detect] face candidates = "
 //              << faceIDList.size() << " / "
 //              << mesh.getFaces().size() << "\n";

    /* =========================================================
       1. 精确检测：裂缝-面交点 & 去重
    =========================================================*/
    intersections.clear();
	const auto&  faces= mesh.getFaces(); // 获取网格面的引用

    for (int fid : faceIDList)
    {
		const Face& face = faces[fid - 1]; // 获取面对象（ID 从 1 开始）
		if (face.FaceNodeCoords.size() < 2) continue; // 跳过无效面
        
        /* —— 二次快速剔除：裂缝盒 vs. 面盒 —— */

		if (useAABBFilter && !fractureBox.overlaps(face.boundingBox)) continue; // 如果裂缝包围盒与面包围盒不相交则跳过

		const Vector& A = face.FaceNodeCoords[0]; // 网格面的两个节点坐标
		const Vector& B = face.FaceNodeCoords[1]; // 网格面的两个节点坐标

		Vector ip; // 存储交点坐标

		if (!lineSegmentIntersection(start, end, A, B, ip)) continue; // 判断裂缝(start→end) 与 面(A→B) 的线段是否相交

        //在裂缝上的归一化位置 param ∈ [0,1]
		double segLen = (end - start).Mag(); // 计算裂缝总长度用于计算归一化参数
		double t = (segLen > 1e-8 ? (ip - start).Mag() / segLen : 0.0); // 计算归一化参数，用于后面排序
		// 如果本 face 上已有一个非常接近的点，就跳过
		bool dup = false;
		for (const auto& e : intersections)
			if (e.edgeID == face.id && (e.point - ip).Mag() < 1e-8)
			{
				dup = true;
				break;
			}
		if (dup) continue; // 如果是重复点则跳过
		// 添加新的交点
		intersections.emplace_back(-1, ip, face.id, t, false, 0, IntersectionOrigin::FracFace);
        
    }

    /* =========================================================
       2. 起点 / 终点 也要算交点（若尚未包含）
       =========================================================*/
    bool startIncluded = false, endIncluded = false;
    for (auto& e : intersections)
    {
        if ((e.point - start).Mag() < 1e-6) startIncluded = true;  //如果其中有一个交点与起点重合那么起点就是交点
        if ((e.point - end).Mag() < 1e-6) endIncluded = true;      //如果其中有一个交点与终点重合那么终点就是交点
    }
    if (!startIncluded)
        intersections.emplace_back
        (
			-1,                     //id(稍后再排)
			start,                  //起点坐标
			-1,                     //关联的 faceID
			0.0,                    // param    
            false,                  // isFF（不是裂缝–裂缝交点）
            0,                      // globalFFID
            IntersectionOrigin::FracStart  // 来源：裂缝起点
        );// 起点
    if (!endIncluded)
        intersections.emplace_back
        (
            -1,
            end,
            -1,
            1.0,
            false,
            0,
            IntersectionOrigin::FracEnd    // 【新增】来源：裂缝终点
        );   // 终点

    /* =========================================================
      3. 排序（按 param）并暂编号
      =========================================================*/
    std::sort(intersections.begin(), intersections.end(),
        [](const FractureIntersectionPointByMatrixMesh& a, const FractureIntersectionPointByMatrixMesh& b)
        {
            return a.param < b.param;
        });
    for (int i = 0; i < intersections.size(); ++i)
        intersections[i].id = i + 1;
    /* =========================================================
   4. 建立 “Cell → 交点” 映射，只保留每个 Cell 最靠前的两个点
   =========================================================*/
    std::unordered_map<int, std::vector<FractureIntersectionPointByMatrixMesh>> perCell;

    for (auto& ip : intersections) {
        int cid = -1;

        if (ip.edgeID >= 1 && ip.edgeID <= (int)faces.size()) {
            // 通过边直接找到宿主 Cell
            const Face& f = faces[ip.edgeID - 1];
            cid = (f.ownerCell >= 0 ? f.ownerCell : f.neighborCell);
        }
        else {
            // 起点/终点：靠几何定位
            cid = findContainingCell(ip.point, meshCells, meshNodes);
        }

        if (cid >= 0)  perCell[cid].push_back(ip);
        else {
            std::cerr << "[Warn] 交点无法匹配 Cell (faceID="
                << ip.edgeID << ")\n";
        }
    }

    std::vector<FractureIntersectionPointByMatrixMesh> filtered;
    for (auto& kv : perCell) {
        auto& v = kv.second;
        if (v.empty()) continue;
        std::sort(v.begin(), v.end(),
            [](auto& a, auto& b) { return a.param < b.param; });
        int take = std::min(2, (int)v.size());
        filtered.insert(filtered.end(), v.begin(), v.begin() + take);
    }

    /* =========================================================
       5. 距离过近(<1e-8) 的点再剔除，最终重编号
       =========================================================*/
    std::sort(filtered.begin(), filtered.end(),
        [](auto& a, auto& b) { return a.param < b.param; });

    std::vector<FractureIntersectionPointByMatrixMesh> cleaned;
    for (const auto& ip : filtered) {
        if (!cleaned.empty() && (ip.point - cleaned.back().point).Mag() < 1e-8)
            continue;
        cleaned.push_back(ip);
    }

    for (int i = 0; i < (int)cleaned.size(); ++i)
        cleaned[i].id = i + 1;

    intersections.swap(cleaned);
}
/*==根据裂缝交点划分裂缝单元，并计算每个单元的长度（computeSegmentLength）、所属网格单元（findContainingCell）以及平均距离（computeDistanceFromCenter&computeAverageDistanceFromNodes）===*/
//void Fracture::subdivide(const vector<Cell>& meshCells, const unordered_map<int, Node>& meshNodes,  bool useCenterDistance)   //useCenterDistance 默认采用computeAverageDistanceFromNodes 如果给定true 则采用computeDistanceFromCenter计算
//{
//    elements.clear();
//    constexpr double tolLen = 1e-10;     // 段长阈值
//    constexpr double tolPoint = 1e-8;    // 判断中点在三角形内部时的容差
//
//    for (int i = 0; i + 1 < (int)intersections.size(); ++i)
//    {
//        const auto& I0 = intersections[i];
//        const auto& I1 = intersections[i + 1];
//       
//        Vector P0 = I0.point, P1 = I1.point;
//        double segLen = computeSegmentLength(P0, P1);
//
//		if (segLen < tolLen)
//		{
//            cerr << "Warning: Zero-length fracture segment detected, "
//                "check duplicate intersection points." <<"point is:"<<"("<< P0.m_x<<","<< P0.m_y<<")" << endl;
//			continue;   // 跳过该段，防止除零
//		}
//
//        Vector mid = computeMidpoint(P0, P1);
//        int associatedCellID = findContainingCell(mid, meshCells, meshNodes);
//        if (associatedCellID < 0) 
//        {
//            std::cerr << "[提示] 裂缝段 param=["
//                << I0.param << "," << I1.param
//                << "] 与任何单元不匹配，可能与网格面平行或擦肩而过，已跳过\n";
//            continue;
//        }
//
//        double avgDist = -1;
//        if (associatedCellID != -1)
//        {
//            for (const auto& cell : meshCells)
//            {
//                if (cell.id == associatedCellID)
//                {
//                    avgDist = useCenterDistance
//                        ? computeDistanceFromCenter(cell, P0, P1)
//                        : computeAverageDistanceFromNodes(cell, meshNodes, P0, P1);
//                    break;
//                }
//            }
//        }
//        else
//        {
//            cerr << "Warning: Midpoint of segment " << i + 1 << " not inside any cell." << std::endl;
//            avgDist = -1.0;
//        }
//        int newElemId = static_cast<int>(elements.size()) + 1;
//        FractureElement elem(newElemId, associatedCellID, segLen, avgDist);
//        elem.param0 = I0.param;
//        elem.param1 = I1.param;
//        elem.isFFatStart = I0.isFF;
//        elem.isFFatEnd = I1.isFF;
//        elem.gIDstart = I0.globalFFID;
//        elem.gIDend = I1.globalFFID;
//        elements.push_back(elem);
//    }
//}

// ====== 三个重载：首选显式策略，其它两个是兼容转发 ======
void Fracture::subdivide(const std::vector<Cell>& meshCells,
    const std::unordered_map<int, Node>& meshNodes,
    DistanceMetric metric)
{
    elements.clear();
    constexpr double tolLen = 1e-10;
    constexpr double tolPoint = 1e-8;

    for (int i = 0; i + 1 < (int)intersections.size(); ++i)
    {
        const auto& I0 = intersections[i];
        const auto& I1 = intersections[i + 1];

        Vector P0 = I0.point, P1 = I1.point;
        double segLen = computeSegmentLength(P0, P1);
        if (segLen < tolLen) {
            std::cerr << "Warning: Zero-length fracture segment @("
                << P0.m_x << "," << P0.m_y << ")\n";
            continue;
        }

        Vector mid = computeMidpoint(P0, P1);
        int associatedCellID = findContainingCell(mid, meshCells, meshNodes);
        if (associatedCellID < 0) {
            std::cerr << "[提示] 裂缝段 param=[" << I0.param << "," << I1.param
                << "] 未匹配到单元，已跳过\n";
            continue;
        }

        double avgDist = -1.0;
        // 找到 cell 对象
        const Cell* host = nullptr;
        for (const auto& c : meshCells) if (c.id == associatedCellID) { host = &c; break; }
        if (!host) {
            std::cerr << "[警告] 找不到 id=" << associatedCellID << " 的单元\n";
            continue;
        }

        // ——— 这里切换三种算法 ———
        switch (metric) {
        case DistanceMetric::CellCenter:
            avgDist = computeDistanceFromCenter(*host, P0, P1);
            break;
        case DistanceMetric::NodeAverage:
            avgDist = computeAverageDistanceFromNodes(*host, meshNodes, P0, P1);
            break;
        case DistanceMetric::AreaWeight:
            avgDist = computeAreaWeightedDistance(*host, meshNodes, P0, P1);
            break;
        case DistanceMetric::CrossAwareGauss:   // ★ 新增
            avgDist = computeCrossAwareAverageDistance(*host, meshNodes, P0, P1);
            break;
        default:
            avgDist = computeAreaWeightedDistance(*host, meshNodes, P0, P1);
            break;
        }

        int newElemId = static_cast<int>(elements.size()) + 1;
        FractureElement elem(newElemId, associatedCellID, segLen, avgDist);
        elem.param0 = I0.param; elem.param1 = I1.param;
        elem.isFFatStart = I0.isFF; elem.isFFatEnd = I1.isFF;
        elem.gIDstart = I0.globalFFID; elem.gIDend = I1.globalFFID;
        elements.push_back(elem);
    }
}

// 兼容旧：bool 转策略
void Fracture::subdivide(const std::vector<Cell>& meshCells,
    const std::unordered_map<int, Node>& meshNodes,
    bool useCenterDistance)
{
    subdivide(meshCells, meshNodes,
        useCenterDistance ? DistanceMetric::CellCenter
        : DistanceMetric::NodeAverage);
}

// 兼容旧：无参版本（默认用面积加权）
void Fracture::subdivide(const std::vector<Cell>& meshCells,
    const std::unordered_map<int, Node>& meshNodes)
{
    subdivide(meshCells, meshNodes, DistanceMetric::AreaWeight);
}




int Fracture::locateSegment(double param) const
{
    for (int i = 0; i < (int)elements.size(); ++i) {
        if (param >= elements[i].param0 && param <= elements[i].param1)
            return i;
    }
    return -1;  // 找不到就返回 -1
}

void Fracture::computeGeometryCouplingCoefficientgeomCIandgeomAlpha()
{
    constexpr double eps_d = 1e-12;
    // 逐段填充纯几何耦合系数
    for (auto& elem : elements)
    {
        //基岩网格各端点到裂缝段的平均距离d 
        double d = max(elem.avgDistance, eps_d);
        // 横截面积 = 段长 * 裂缝厚度（2D 为1）
        double A_fr = elem.length * 1;  //Note: 2D中 此时为1m 
        // geomCI: 纯几何耦合  CI = Af/d  ref:Modeling study of the thermal-hydraulic-mechanical coupling process for EGS based on the framework of EDFM and XFEM
        elem.geomCI = A_fr / d;
		// geomAlpha: 纯几何 α  alpha= 2/ L ref:Modeling study of the thermal-hydraulic-mechanical coupling process for EGS based on the framework of EDFM and XFEM
        elem.geomAlpha = 2.0 / elem.length;
    }
}

// -----------------------------------------------------------------------------

/*====点到线段的距离（pointToSegmentDistance）====*/
double Fracture::pointToSegmentDistance(const Vector& p, const Vector& s, const Vector& e) {
    Vector se = e - s;
    double segLenSq = se.m_x * se.m_x + se.m_y * se.m_y;
    if (segLenSq < 1e-12) return (p - s).Mag();
    double t = ((p - s) * se) / segLenSq;
    t = max(0.0, std::min(1.0, t));
    Vector projection = s + se * t;
    return (p - projection).Mag();
}

/*====判断点是否在三角形内 （pointInTriangle）====*/
bool Fracture::pointInTriangle(const Vector& p, const Vector& a, const Vector& b, const Vector& c) 
{
    Vector v0 = c - a, v1 = b - a, v2 = p - a;
    double dot00 = v0 * v0, dot01 = v0 * v1, dot02 = v0 * v2;
    double dot11 = v1 * v1, dot12 = v1 * v2;
    double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    return (u > 1e-6) && (v > 1e-6) && ((u + v) < 1 - 1e-6);

}

/*====判断两线段是否有交点且在两线段内 （lineSegmentIntersection）====*/
bool Fracture::lineSegmentIntersection(const Vector& p, const Vector& q, const Vector& r, const Vector& s, Vector& ip) 
{
    if (std::max(p.m_x, q.m_x) < std::min(r.m_x, s.m_x) ||
        std::min(p.m_x, q.m_x) > std::max(r.m_x, s.m_x) ||
        std::max(p.m_y, q.m_y) < std::min(r.m_y, s.m_y) ||
        std::min(p.m_y, q.m_y) > std::max(r.m_y, s.m_y)) {
        return false; // bounding box 不相交
    }
    const double tol = 1e-10;
    // 精确几何检测阶段--向量叉积法
    Vector pq = q - p;
    Vector rs = s - r;

    double denom = pq.m_x * rs.m_y - pq.m_y * rs.m_x;

    // 1) 判断是否平行（或共线）
    if (fabs(denom) < tol) 
    {
        return false;  // 平行或共线，不考虑重合情况
    }

    Vector pr = r - p;

    double t = (pr.m_x * rs.m_y - pr.m_y * rs.m_x) / denom;
    double u = (pr.m_x * pq.m_y - pr.m_y * pq.m_x) / denom;

    // 2) 判断 t、u 是否都在 [0,1] 内（带容差）
    if (t < -tol || t > 1.0 + tol || u < -tol || u > 1.0 + tol) 
    {
        return false;  // 交点不在线段范围内
    }

    // 3) 计算交点
    ip = p + pq * t;
    
    return true;
}

/*====计算裂缝段长度（computeSegmentLength）====*/
double Fracture::computeSegmentLength(const Vector& a, const Vector& b) {
    return (b - a).Mag();
}

/*====计算裂缝段中点坐标（computeMidpoint）====*/
Vector Fracture::computeMidpoint(const Vector& a, const Vector& b) 
{
    return (a + b) / 2.0;
}


// —— 工具1：点是否在二维线段上（含端点）——
static inline bool pointOnSegment2D(const Vector& P, const Vector& A, const Vector& B, double eps = 1e-12) {
    Vector AP = P - A, AB = B - A;
    // 面积（叉积模）≈0 表示共线
    double cross = std::fabs(AP.m_x * AB.m_y - AP.m_y * AB.m_x);
    if (cross > eps) return false;
    // 点积范围判断（投影在线段内部）
    double dot = AP.m_x * AB.m_x + AP.m_y * AB.m_y;
    if (dot < -eps) return false;
    double ab2 = AB.m_x * AB.m_x + AB.m_y * AB.m_y;
    if (dot - ab2 > eps) return false;
    return true;
}

// —— 工具2：二维多边形点内测（奇偶规则），包含边界 ——
// 顶点按环序给出（你的 CellNodeIDs 就是环序）
static bool pointInPolygon2D(const Vector& P,
    const std::vector<int>& polyNodeIDs,
    const std::unordered_map<int, Node>& nodes,
    double eps = 1e-12)
{
    const int n = (int)polyNodeIDs.size();
    if (n < 3) return false;

    // 边界命中优先
    for (int i = 0; i < n; ++i) {
        const Vector& A = nodes.at(polyNodeIDs[i]).coord;
        const Vector& B = nodes.at(polyNodeIDs[(i + 1) % n]).coord;
        if (pointOnSegment2D(P, A, B, eps)) return true;
    }

    // 射线法：向 +x 发射
    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        const Vector& Pi = nodes.at(polyNodeIDs[i]).coord;
        const Vector& Pj = nodes.at(polyNodeIDs[j]).coord;

        // 只用 xy；忽略 z
        bool intersect = ((Pi.m_y > P.m_y) != (Pj.m_y > P.m_y));
        if (intersect) {
            double x_cross = Pj.m_x + (Pi.m_x - Pj.m_x) * ((P.m_y - Pj.m_y) / (Pi.m_y - Pj.m_y + 1e-30));
            if (x_cross > P.m_x - eps) inside = !inside;
        }
    }
    return inside;
}


/*====判断交点是否在基岩网格内，并找到包含该点的 Cell（findContainingCell）====*/
//int Fracture::findContainingCell(const Vector& point, const std::vector<Cell>& cells, const std::unordered_map<int, Node>& nodes)
//{
//    for (const auto& cell : cells) 
//    {
//        if (cell.CellNodeIDs.size() != 3) 
//            continue;
//        Vector a = nodes.at(cell.CellNodeIDs[0]).coord;
//        Vector b = nodes.at(cell.CellNodeIDs[1]).coord;
//        Vector c = nodes.at(cell.CellNodeIDs[2]).coord;
//        if (pointInTriangle(point, a, b, c)) 
//            return cell.id;
//        else if (pointInTriangle(point, a, b, c))
//        {
//            std::cerr << "中点落在单元边界或顶点上，建议加密网格或调整裂缝。" << std::endl;
//        }
//    }
//    return -1;
//}

// —— 替换你的函数：支持三角/四边形/任意 2D 多边形 ——
// 约定：工作在 XY 平面（你的裂缝/基岩2D情形正是如此）
int Fracture::findContainingCell(const Vector& point,
    const std::vector<Cell>& cells,
    const std::unordered_map<int, Node>& nodes)
{
    int bestCid = -1;
    double bestD2 = 1e300;

    for (const auto& cell : cells) {
        const auto& ids = cell.CellNodeIDs;
        if ((int)ids.size() < 3) continue; // 至少三点

        // 快速 AABB 剪枝（提速）
        double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
        for (int nid : ids) {
            const Vector& v = nodes.at(nid).coord;
            xmin = std::min(xmin, v.m_x); xmax = std::max(xmax, v.m_x);
            ymin = std::min(ymin, v.m_y); ymax = std::max(ymax, v.m_y);
        }
        if (point.m_x < xmin - 1e-12 || point.m_x > xmax + 1e-12 ||
            point.m_y < ymin - 1e-12 || point.m_y > ymax + 1e-12) continue;

        // 真正点内测：任意多边形（含三角/四边形）
        if (pointInPolygon2D(point, ids, nodes)) 
        {
            // 多个候选时，选重心更近者（更稳）

            double d2 = (point - cell.center).Mag2();
            if (d2 < bestD2) { bestD2 = d2; bestCid = cell.id; }
        }
    }
    return bestCid; // -1 表示未找到
}

/*====计算裂缝段到基岩网格单元的平均距离--方法1利用 Cell 中的各节点计算====*/
double Fracture::computeAverageDistanceFromNodes(const Cell& cell, const std::unordered_map<int, Node>& meshNodes, const Vector& segStart, const Vector& segEnd)
{
    double sumDist = 0.0;
    for (int nid : cell.CellNodeIDs) {
        Vector nodeCoord = meshNodes.at(nid).coord;
        sumDist += pointToSegmentDistance(nodeCoord, segStart, segEnd);
    }
    return sumDist / cell.CellNodeIDs.size();
}

/*====计算裂缝段到基岩网格单元的平均距离--方法2利用 Cell 中点计算====*/
double Fracture::computeDistanceFromCenter(const Cell& cell, const Vector& segStart, const Vector& segEnd) {
    return pointToSegmentDistance(cell.center, segStart, segEnd);
}


void Fracture::sortAndRenumberIntersections()
{
    sort(intersections.begin(), intersections.end(),
        [](const FractureIntersectionPointByMatrixMesh& a,
            const FractureIntersectionPointByMatrixMesh& b)
        { return a.param < b.param; });

    int nid = 1;
    for (auto& I : intersections) I.id = nid++;
}


////  Fracture::computeFractureDiscreteCoefficients
////      • 端面导流系数 aW/aE           : 2ρ kf w / μ(Li+Lj)
////      • 有效渗透率 k_eff             : (Af φf + Ac φm) / (Af φf/kf + Ac φm/km)
////      • 裂缝-基岩连接指数 CI (b_fr)  : (kf/μ) Af / d      ← 文献Tingyu师兄公式 ΘA/d 单位：[m³⋅Pa⁻¹⋅s⁻¹]
////      现在处理的是裂缝在不同位置保持开度相同的情况，如果以后每段开度不同，只需把 aperture 换成 elem.aperture_i 即可。
//// -----------------------------------------------------------------------------
//
//
//
//void Fracture::computeFractureDiscreteCoefficients(const FluidProperties& water,
//    const Matrix& matrix,
//    Mesh& mesh)
//{
//    computeGeometryCouplingCoefficientgeomCIandgeomAlpha();
//    constexpr double eps_d = 1e-6;
//
//    constexpr double eps_mu = 1e-30;
//
//
//    // 物理 θ_phys = k_f / μ
//    double theta_phys = permeability / (fluid.fluid_viscosity + eps_mu);
//	cout << "物理 θ_phys = k_f / μ = " << theta_phys << endl;
//    // 物理流动 coef_phys = 2·ρ·k_f·w / μ （用于 aW/aE）
//    double coef_phys = 2.0 * fluid.fluid_rho * permeability * aperture
//        / (fluid.fluid_viscosity + eps_mu);
//    for (size_t i = 0; i < elements.size(); ++i) {
//        auto& elem = elements[i];
//        // —— 1. 端面导流系数 aW/aE ——  
//        double L = elem.length;
//        double Lp = (i > 0) ? elements[i - 1].length : 0.0;
//        double Ln = (i + 1 < elements.size()) ? elements[i + 1].length : 0.0;
//        if (i == 0) {
//            elem.aW_fr = coef_phys / L;
//            elem.aE_fr = coef_phys / (L + Ln);
//        }
//        else if (i + 1 == elements.size()) {
//            elem.aE_fr = coef_phys / L;
//            elem.aW_fr = coef_phys / (L + Lp);
//        }
//        else {
//            elem.aW_fr = coef_phys / (L + Lp);
//            elem.aE_fr = coef_phys / (L + Ln);
//        }
//        // —— 2. 物性 CI：CI_phys = θ_phys * geomCI ——  
//        elem.b_fr = theta_phys * elem.geomCI;
//        // —— 3. 总系数 aP ——  
//        elem.aP_fr = elem.aW_fr + elem.aE_fr + elem.b_fr;
//        // —— 4. 物理 α_phys（TI 用）——  
//        elem.alpha_fr = 2.0 * permeability * aperture
//            / (fluid.fluid_viscosity * elem.length + eps_mu);
//        // —— 5. 写回宿主单元，只存物性 CI_phys ——  
//        auto it = mesh.cellId2index.find(elem.cellID);
//        if (it != mesh.cellId2index.end()) {
//            Cell& host = mesh.cells[it->second];
//            host.CI.push_back(elem.b_fr);
//            host.CI_belongFraction.push_back(elem.id);
//            host.fracturePressure.push_back(elem.p_fr);
//        }
//    }
// 
//}
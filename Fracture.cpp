#include "Fracture.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include "Fluid.h"
#include "Matrix.h"
#include "Mesh.h"
#include <unordered_map>
#include <numeric>
#include "FractureGeometryUtils.h"


Fracture::Fracture(const Vector& s, const Vector& e)  
   : start(s), end(e), id(0)  
{  

}

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
        return FractureGeom::pointToSegmentDistance(cell.center, segStart, segEnd);
    }

    // 判断是否“穿越”
    std::vector<Vector> hits;
    int nHit = FractureGeom::countIntersectionsWithPolygon(poly, segStart, segEnd, hits);
    if (nHit >= 2) {
        // 按直线把多边形一分为二，两个子多边形分别 fan + Gauss，再面积加权
        std::vector<Vector> Ppos, Pneg;
        FractureGeom::splitPolygonByLine(poly, segStart, segEnd, Ppos, Pneg);

        double areaPos = 0.0, areaNeg = 0.0;
        if (Ppos.size() >= 3) { Vector c = FractureGeom::polygonCentroid(Ppos, areaPos); (void)c; }
        if (Pneg.size() >= 3) { Vector c = FractureGeom::polygonCentroid(Pneg, areaNeg); (void)c; }
        double A = areaPos + areaNeg;
        if (A <= 0) {
            // 极端退化：回退到原扇形 + 高斯
            return FractureGeom::fanGaussAverage(poly, segStart, segEnd);
        }

        double favPos = (areaPos > 0 && Ppos.size() >= 3) ? FractureGeom::fanGaussAverage(Ppos, segStart, segEnd) : 0.0;
        double favNeg = (areaNeg > 0 && Pneg.size() >= 3) ? FractureGeom::fanGaussAverage(Pneg, segStart, segEnd) : 0.0;

        return (favPos * areaPos + favNeg * areaNeg) / A;
    }
    else {
        // 不穿越：直接对原多边形做扇形 + 三点高斯
        return FractureGeom::fanGaussAverage(poly, segStart, segEnd);
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
        return FractureGeom::pointToSegmentDistance(cell.center, segStart, segEnd);
    }

    double polyArea = 0.0;
	Vector C = FractureGeom::polygonCentroid(poly, polyArea); // 多边形重心

    double sum = 0.0, sumA = 0.0;
    const size_t n = poly.size();
    for (size_t i = 0; i < n; ++i) {
        const Vector& A = poly[i];
        const Vector& B = poly[(i + 1) % n];
        double Ai = FractureGeom::triArea(C, A, B);
        if (Ai <= 0) continue;
        Vector Gi = FractureGeom::triCentroid(C, A, B);
        double di = FractureGeom::pointToSegmentDistance(Gi, segStart, segEnd);
        sum += di * Ai;
        sumA += Ai;
    }

    if (sumA <= 0) {
        // 再次回退
        return FractureGeom::pointToSegmentDistance(cell.center, segStart, segEnd);
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

		if (!lineSegmentIntersection(start, end, A, B, ip)) continue; // 判断裂缝(start→end) 与 面(A→B) 的线段是否相交,并输出交点坐标到 ip

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

// ====== 根据裂缝交点对裂缝进行划分 ======
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
        case DistanceMetric::CrossAwareGauss:  
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
        if (FractureGeom::pointInPolygon2D(point, ids, nodes))
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
        sumDist += FractureGeom::pointToSegmentDistance(nodeCoord, segStart, segEnd);
    }
    return sumDist / cell.CellNodeIDs.size();
}

/*====计算裂缝段到基岩网格单元的平均距离--方法2利用 Cell 中点计算====*/
double Fracture::computeDistanceFromCenter(const Cell& cell, const Vector& segStart, const Vector& segEnd) {
    return FractureGeom::pointToSegmentDistance(cell.center, segStart, segEnd);
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
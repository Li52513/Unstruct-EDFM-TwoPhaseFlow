#include "Fracture.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <numeric>

#include "UserDefineVarType.h"

#include "Mesh.h"
#include "Cell.h"
#include "Node.h"
#include "8_DOP.h"

#include "GeometryCalculate.h"
#include "FractureIntersectionPoint.h"
#include "FractureElement.h"

// ==================
// 构造函数实现
// ===================
Fracture::Fracture(const Vector& s, const Vector& e)  
   : start(s), end(e), id(0)  
{  

}

// ==================
// 静态函数实现
// ===================
/**
* 计算裂缝段长度
*/
double Fracture::computeSegmentLength(const Vector& a, const Vector& b) {
    return (b - a).Mag();
}

/**
* 计算裂缝段中点坐标
*/
Vector Fracture::computeMidpoint(const Vector& a, const Vector& b)
{
    return (a + b) / 2.0;
}

/**
* 精确判断两线段是否有交点且在两线段内
*/
bool Fracture::lineSegmentIntersection(const Vector& p, const Vector& q, const Vector& r, const Vector& s, Vector& ip)
{
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
/**
* 判断2D情况点是否在某个单元内，支持三角/四边形/任意 2D 多边形。约定在x,y平面，没有z方向
*/
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
        if (GeomCalculate::pointInPolygon2D(point, ids, nodes))
        {
            // 多个候选时，选重心更近者（更稳）

            double d2 = (point - cell.center).Mag2();
            if (d2 < bestD2) { bestD2 = d2; bestCid = cell.id; }
        }
    }
    return bestCid; // -1 表示未找到
}

/**
* 计算裂缝段到基岩网格单元的平均距离--方法1利用 Cell 中的各节点计算
*/
double Fracture::computeAverageDistanceFromNodes(const Cell& cell, const std::unordered_map<int, Node>& meshNodes, const Vector& segStart, const Vector& segEnd)
{
    double sumDist = 0.0;
    for (int nid : cell.CellNodeIDs) {
        Vector nodeCoord = meshNodes.at(nid).coord;
        sumDist += GeomCalculate::pointToSegmentDistance(nodeCoord, segStart, segEnd);
    }
    return sumDist / cell.CellNodeIDs.size();
}

/**
* 计算裂缝段到基岩网格单元的平均距离--方法2利用 Cell 中点计算
*/
double Fracture::computeDistanceFromCenter(const Cell& cell, const Vector& segStart, const Vector& segEnd) {
    return GeomCalculate::pointToSegmentDistance(cell.center, segStart, segEnd);
}

/**
* 计算裂缝段到基岩网格单元的平均距离--方法3 高精度方法。将单元划分为多个三角形，通过面积加权积分计算平均距离
*/
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
        return GeomCalculate::pointToSegmentDistance(cell.center, segStart, segEnd);
    }

    double polyArea = 0.0;
    Vector C = GeomCalculate::polygonCentroid(poly, polyArea); // 多边形重心

    double sum = 0.0, sumA = 0.0;
    const size_t n = poly.size();
    for (size_t i = 0; i < n; ++i) 
    {
        const Vector& A = poly[i];
        const Vector& B = poly[(i + 1) % n];
        double Ai = GeomCalculate::triArea(C, A, B);
        if (Ai <= 0) continue;
        Vector Gi = GeomCalculate::triCentroid(C, A, B);
        double di = GeomCalculate::pointToSegmentDistance(Gi, segStart, segEnd);
        sum += di * Ai;
        sumA += Ai;
    }

    if (sumA <= 0) {
        // 再次回退
        return GeomCalculate::pointToSegmentDistance(cell.center, segStart, segEnd);
    }
    return sum / sumA;
}

/**
* 计算裂缝段到基岩网格单元的平均距离--方法4  最高精度方法。考虑裂缝是否完全贯穿单元，若贯穿则分割多边形分别积分（高斯积分+扇形分割）
*/
double Fracture::computeCrossAwareAverageDistance(
    const Cell& cell,
    const std::unordered_map<int, Node>& meshNodes,
    const Vector& segStart, const Vector& segEnd)
{
    // 取单元多边形
    std::vector<Vector> poly; poly.reserve(cell.CellNodeIDs.size());
    for (int nid : cell.CellNodeIDs) poly.push_back(meshNodes.at(nid).coord);
    if (poly.size() < 3)
    {
        // 退化：回退到中心点距离
        return GeomCalculate::pointToSegmentDistance(cell.center, segStart, segEnd);
    }
    // 判断是否“穿越”
    std::vector<Vector> hits;
    int nHit = GeomCalculate::countIntersectionsWithPolygon(poly, segStart, segEnd, hits);
    if (nHit >= 2) {
        // 按直线把多边形一分为二，两个子多边形分别 fan + Gauss，再面积加权
        std::vector<Vector> Ppos, Pneg;
        GeomCalculate::splitPolygonByLine(poly, segStart, segEnd, Ppos, Pneg);

        double areaPos = 0.0, areaNeg = 0.0;
        if (Ppos.size() >= 3) { Vector c = GeomCalculate::polygonCentroid(Ppos, areaPos); (void)c; }
        if (Pneg.size() >= 3) { Vector c = GeomCalculate::polygonCentroid(Pneg, areaNeg); (void)c; }
        double A = areaPos + areaNeg;
        if (A <= 0) {
            // 极端退化：回退到原扇形 + 高斯
            return GeomCalculate::fanGaussAverage(poly, segStart, segEnd);
        }

        double favPos = (areaPos > 0 && Ppos.size() >= 3) ? GeomCalculate::fanGaussAverage(Ppos, segStart, segEnd) : 0.0;
        double favNeg = (areaNeg > 0 && Pneg.size() >= 3) ? GeomCalculate::fanGaussAverage(Pneg, segStart, segEnd) : 0.0;

        return (favPos * areaPos + favNeg * areaNeg) / A;
    }
    else
    {
        // 不穿越：直接对原多边形做扇形 + 三点高斯
        return GeomCalculate::fanGaussAverage(poly, segStart, segEnd);
    }
}

// ==================
// 成员行为函数实现
// ===================

/*
* 根据起点和终点，计算fractureBox 属性
*/
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

/*
 * 利用基岩网格边界、以及裂缝的起点和终点记录裂缝的交点集合并排序，不包括裂缝–裂缝交点-重载
 */
void Fracture::DetectFracturetoMeshFaceIntersections_improved
(
    const Mesh& mesh,
    const std::vector<Cell>& meshCells,
    const std::unordered_map<int, Node>& meshNodes,
    IntersectionSearchStrategy_2D strategy,
    IntersectionStatistics_2D& stats
)
{
    // 重置统计器
    stats.reset();

    // 0. 更新裂缝自身基础属性
    computeAABB(); // 更新 fractureBox (AABB)

    // [新增] 预计算 8-DOP (仅当策略需要时构建)
    bool use8DOP = (strategy == IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP ||
        strategy == IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    Box8DOP frac8DOP;
    if (use8DOP) {
        frac8DOP.fromSegment(start, end);
    }

    // 1. 候选面列表生成 (Candidate Generation)
    std::vector<int> faceIDList;
    const auto& faces = mesh.getFaces();

    switch (strategy)
    {
    case IntersectionSearchStrategy_2D::GridIndexing:
    case IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP:
        // 【方法3 & 方法4】：矩形区域 Bin 查询
        faceIDList = mesh.getCandidateFacesFromBins(fractureBox);
        break;

    case IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA:
        // 【方法5】：DDA 射线追踪查询 (最精准的候选集)
        faceIDList = mesh.getCandidateFacesFromBins_DDA(start, end);
        break;

    case IntersectionSearchStrategy_2D::BruteForce:
    case IntersectionSearchStrategy_2D::GlobalAABB:
    default:
        // 【方法1 & 方法2】：暴力/全局遍历 -> 构造全集 [1, 2, ..., N]
        if (!faces.empty()) {
            faceIDList.resize(faces.size());
            std::iota(faceIDList.begin(), faceIDList.end(), 1);
        }
        break;
    }

    /// 记录统计指标1：候选面数量
    stats.candidateCount = static_cast<long long>(faceIDList.size());
    candidateFaceCount_ = static_cast<int>(faceIDList.size()); // 更新成员变量

    // 2. 精确检测循环 (Intersection Testing Loop)
    intersections.clear();

    for (int fid : faceIDList)
    {
        /// 安全性检查
        if (fid < 1 || fid > static_cast<int>(faces.size())) continue;
        int fIdx = mesh.getFaceIndex(fid);
        if (fIdx == -1) {
            // 这是一个异常情况，理论上候选列表里的 ID 应该都存在
            continue;
        }
        const Face& face = faces[fIdx];

        if (face.FaceNodeCoords.size() < 2) continue;

        /// 策略分支：决定是否执行包围盒过滤 (Culling)
        bool proceedToGeometry = false;

        if (strategy == IntersectionSearchStrategy_2D::BruteForce)
        {
            // 【方法1】：纯暴力，跳过所有 culling
            proceedToGeometry = true;
        }
        else if (use8DOP)
        {
            // 【方法4 & 方法5】：8-DOP 精筛
            // Box8DOP::overlaps 内部已包含 AABB 检测，因此无需重复检测 AABB
            stats.aabbCheckCount++;
            if (frac8DOP.overlaps(face.boundingBox)) {
                proceedToGeometry = true;
            }
        }
        else
        {
            // 【方法2 & 方法3】：标准 AABB 粗筛
            stats.aabbCheckCount++;
            if (this->fractureBox.overlaps(face.boundingBox)) {
                proceedToGeometry = true;
            }
        }

        /// 精确几何求交 (开销最大)
        if (proceedToGeometry)
        {
            stats.geometryCheckCount++; // 记录昂贵的计算次数

            const Vector& A = face.FaceNodeCoords[0];
            const Vector& B = face.FaceNodeCoords[1];
            Vector ip;

            // 执行线段-线段求交
            if (!lineSegmentIntersection(start, end, A, B, ip))
                continue;

            // --- 后续逻辑：计算参数、去重、存储 ---
            double segLen = (end - start).Mag();
            double t = (segLen > 1e-12 ? (ip - start).Mag() / segLen : 0.0);

            // 局部去重 (防止同一个面因误差算两次)
            bool dup = false;
            for (const auto& e : intersections) {
                if (e.edgeID == face.id && (e.point - ip).Mag() < 1e-8) {
                    dup = true;
                    break;
                }
            }
            if (dup) continue;

            // 成功添加交点
            intersections.emplace_back(-1, ip, face.id, t, false, 0, IntersectionOrigin::FracFace);
            stats.trueIntersectionCount++; // 记录真实命中
        }
    }

    // 3. 后处理：起点/终点包含性检测 (逻辑保持原有)
    bool startIncluded = false;
    bool endIncluded = false;

    for (const auto& e : intersections) {
        if ((e.point - start).Mag() < 1e-6) startIncluded = true;
        if ((e.point - end).Mag() < 1e-6) endIncluded = true;
    }

    if (!startIncluded) {
        intersections.emplace_back(-1, start, -1, 0.0, false, 0, IntersectionOrigin::FracStart);
    }
    if (!endIncluded) {
        intersections.emplace_back(-1, end, -1, 1.0, false, 0, IntersectionOrigin::FracEnd);
    }

    // 4. 排序与全局编号
    std::sort(intersections.begin(), intersections.end(),
        [](const FractureIntersectionPoint& a, const FractureIntersectionPoint& b) {
            return a.param < b.param;
        });

    for (size_t i = 0; i < intersections.size(); ++i) {
        intersections[i].id = static_cast<int>(i) + 1;
    }

    // 5. Cell 映射与过滤 (每个 Cell 最多保留 2 个交点)
    // 注意：此处使用了临时 map，若性能极其敏感可优化，但对于 2D 通常不是瓶颈
    std::unordered_map<int, std::vector<FractureIntersectionPoint>> perCell;
    for (const auto& ip : intersections)
    {
        int cid = -1;

        if (ip.edgeID >= 1 && ip.edgeID <= static_cast<int>(faces.size()))
        {
            // 通过边直接找到宿主 Cell
            const Face& f = faces[ip.edgeID - 1];
            cid = (f.ownerCell >= 0 ? f.ownerCell : f.neighborCell);
        }
        else
        {
            // 起点/终点：靠几何定位
            cid = findContainingCell(ip.point, meshCells, meshNodes);
        }

        if (cid >= 0)
        {
            perCell[cid].push_back(ip);
        }
        else
        {
            // 这是一个常见警告，当裂缝端点正好落在域外时发生，通常可忽略或作为边界处理
            // std::cerr << "[Warn] 交点无法匹配 Cell (faceID=" << ip.edgeID << ")\n";
        }
    }

    std::vector<FractureIntersectionPoint> finalPoints;
    finalPoints.reserve(intersections.size());

    // 策略：仅保留每个 Cell 内有效的线段端点
    // 这里保持原有逻辑：按 param 顺序，去除极近的重复点
    for (const auto& ip : intersections) {
        if (!finalPoints.empty() && (ip.point - finalPoints.back().point).Mag() < 1e-8) {
            continue; // 跳过重复
        }
        finalPoints.push_back(ip);
    }

    // 更新回去
    intersections = std::move(finalPoints);
    // 再次重编号以保持 ID 连续性
    for (size_t i = 0; i < intersections.size(); ++i) {
        intersections[i].id = static_cast<int>(i) + 1;
    }
}

/*
 * 排序交点并重新编号
 */
void Fracture::sortAndRenumberIntersections()
{
    sort(intersections.begin(), intersections.end(),
        [](const FractureIntersectionPoint& a,
            const FractureIntersectionPoint& b)
        { return a.param < b.param; });

    int nid = 1;
    for (auto& I : intersections) I.id = nid++;
}

/*
 * 根据裂缝交点对裂缝进行划分
 */
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

/*
 * 计算几何耦合系数 geomCI 和 geomAlpha（待修正）
 */
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

/*
 * 给定 param，定位它属于哪一段
 */
int Fracture::locateSegment(double param) const
{
    for (int i = 0; i < (int)elements.size(); ++i) {
        if (param >= elements[i].param0 && param <= elements[i].param1)
            return i;
    }
    return -1;  // 找不到就返回 -1
}

















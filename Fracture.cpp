#include "Fracture.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include "Fluid.h"
#include "Matrix.h"
#include "Mesh.h"
#include <unordered_map>


Fracture::Fracture(const Vector& s, const Vector& e)  
   : start(s), end(e), id(0)  
{  
     
}


/*====利用基岩网格边界对裂缝进行网格划分并排序---基于lineSegmentIntersection子函数 *注：对于非结构化网格需要进行排序===*/
void Fracture::DetectFracturetoMeshFaceIntersections(const std::vector<Face>& meshFaces,const std::vector<Cell>& meshCells,const std::unordered_map<int, Node>& meshNodes)
{
    // 1) 面内交点计算 & 去重
    intersections.clear();
    for (const auto& face : meshFaces)
    {
        if (face.FaceNodeCoords.size() < 2) continue;   //跳过无效面
        Vector A = face.FaceNodeCoords[0];              //网格面的两个节点坐标    
        Vector B = face.FaceNodeCoords[1];              //网格面的两个节点坐标
        Vector ip;                                  //存储交点坐标
        // 判断裂缝(start→end) 与 面(A→B) 的线段相交
        if (lineSegmentIntersection(start, end, A, B, ip))
        {
            double segLen = (end - start).Mag(); // 计算裂缝总长度用于计算归一化参数
            double t = (segLen > 1e-8 ? (ip - start).Mag() / segLen : 0.0); //计算归一化参数，用于后面排序
            // 如果本 face 上已有一个非常接近的点，就跳过
            bool dup = false;
            for (auto& e : intersections)
            {
                if (e.edgeID == face.id && (e.point - ip).Mag() < 1e-8)
                {
                    dup = true; break; 
                }
            }

            if (!dup)
                intersections.emplace_back
                (
                    -1,                // id（稍后再排）
                    ip,                // 交点坐标
                    face.id,           // 关联的 faceID
                    t,                 // param
                    false,             // isFF（不是裂缝–裂缝交点）
                    0,                 // globalFFID
                    IntersectionOrigin::FracFace   // 来源：裂缝–网格面
                );
           
        }

    }
        // 2) 起点和终点作为交点（若未被添加）
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

    // 3) 按 param 全局升序排序 & 赋临时 ID ---
    sort(intersections.begin(), intersections.end(),
        [](const FractureIntersectionPointByMatrixMesh& a, const FractureIntersectionPointByMatrixMesh& b)
        {
            return a.param < b.param;
        });
    for (int i = 0; i < intersections.size(); ++i)
        intersections[i].id = i + 1;

    //4）找到交点对应的Cell并建立映射 ---
	unordered_map<int, vector<FractureIntersectionPointByMatrixMesh>> perCell; //用于建立 Cell 和交点的映射关系
    for (auto& ip : intersections)
    {
        int cid = -1;
        // 情况 A：标准面交点
        if (ip.edgeID >= 1 && ip.edgeID <= (int)meshFaces.size()) 
        {
            const Face& f = meshFaces[ip.edgeID - 1];
            cid = (f.ownerCell >= 0 ? f.ownerCell : f.neighborCell);  //什么情况会取neighbor neighbor是什么
        }
        // 情况 B：起点/终点（edgeID == -1），需要几何定位
        else    
        {
            cid = findContainingCell(ip.point, meshCells, meshNodes);
        }

        if (cid >= 0) 
        {
            perCell[cid].push_back(ip);
        }
        else 
        {
            cerr << "[警告] 裂缝交点无法分配宿主 Cell, 点("
                << ip.point.m_x << "," << ip.point.m_y << ") "
                << "来自 faceID=" << ip.edgeID << "\n";
        }
    }

    // 5) 每个 Cell 内只保留最前面的两个点  ***用来处理贯穿网格顶点的情况****
    vector<FractureIntersectionPointByMatrixMesh> filtered;
    for (auto& kv : perCell)
    {
        auto& vecPoints = kv.second;
        if (vecPoints.empty()) continue;
        sort(vecPoints.begin(), vecPoints.end(),
            [](auto& a, auto& b) { return a.param < b.param; });
        int take = std::min(2, (int)vecPoints.size());
        for (int i = 0; i < take; ++i)
            filtered.push_back(vecPoints[i]);
    }


    // --- 6) 用过滤结果更新 intersections，并重新排序与编号 ---
    intersections = std::move(filtered);
    sort(intersections.begin(), intersections.end(),
        [](auto& a, auto& b) { return a.param < b.param; });

    // 再次剔除距离过近（<1e-8）的交点
    vector<FractureIntersectionPointByMatrixMesh> cleaned;
    for (const auto& ip : intersections) 
    {
        if (!cleaned.empty() &&
            (ip.point - cleaned.back().point).Mag() < 1e-8) 
        {
            continue; // 跳过过近重复点
        }
        cleaned.push_back(ip);
    }

    // 重编号
    for (int i = 0; i < (int)cleaned.size(); ++i)
        cleaned[i].id = i + 1;

    intersections = move(cleaned);
}
/*==根据裂缝交点划分裂缝单元，并计算每个单元的长度（computeSegmentLength）、所属网格单元（findContainingCell）以及平均距离（computeDistanceFromCenter&computeAverageDistanceFromNodes）===*/
void Fracture::subdivide(const vector<Cell>& meshCells, const unordered_map<int, Node>& meshNodes,  bool useCenterDistance)   //useCenterDistance 默认采用computeAverageDistanceFromNodes 如果给定true 则采用computeDistanceFromCenter计算
{
    elements.clear();
    constexpr double tolLen = 1e-10;     // 段长阈值
    constexpr double tolPoint = 1e-8;    // 判断中点在三角形内部时的容差

    for (int i = 0; i + 1 < (int)intersections.size(); ++i)
    {
        const auto& I0 = intersections[i];
        const auto& I1 = intersections[i + 1];
       
        Vector P0 = I0.point, P1 = I1.point;
        double segLen = computeSegmentLength(P0, P1);

		if (segLen < tolLen)
		{
            cerr << "Warning: Zero-length fracture segment detected, "
                "check duplicate intersection points." <<"point is:"<<"("<< P0.m_x<<","<< P0.m_y<<")" << endl;
			continue;   // 跳过该段，防止除零
		}

        Vector mid = computeMidpoint(P0, P1);
        int associatedCellID = findContainingCell(mid, meshCells, meshNodes);
        if (associatedCellID < 0) 
        {
            std::cerr << "[提示] 裂缝段 param=["
                << I0.param << "," << I1.param
                << "] 与任何单元不匹配，可能与网格面平行或擦肩而过，已跳过\n";
            continue;
        }

        double avgDist = -1;
        if (associatedCellID != -1)
        {
            for (const auto& cell : meshCells)
            {
                if (cell.id == associatedCellID)
                {
                    avgDist = useCenterDistance
                        ? computeDistanceFromCenter(cell, P0, P1)
                        : computeAverageDistanceFromNodes(cell, meshNodes, P0, P1);
                    break;
                }
            }
        }
        else
        {
            cerr << "Warning: Midpoint of segment " << i + 1 << " not inside any cell." << std::endl;
            avgDist = -1.0;
        }
        int newElemId = static_cast<int>(elements.size()) + 1;
        FractureElement elem(newElemId, associatedCellID, segLen, avgDist);
        elem.param0 = I0.param;
        elem.param1 = I1.param;
        elem.isFFatStart = I0.isFF;
        elem.isFFatEnd = I1.isFF;
        elem.gIDstart = I0.globalFFID;
        elem.gIDend = I1.globalFFID;
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
        // 横截面积 = 段长 * 裂缝开度
        double A_fr = elem.length * elem.aperture;
        // geomCI: 纯几何耦合  CI = Af/d  ref:Modeling study of the thermal-hydraulic-mechanical coupling process for EGS based on the framework of EDFM and XFEM
        elem.geomCI = A_fr / d;
		// geomAlpha: 纯几何 α  alpha= 2·Af / L ref:Modeling study of the thermal-hydraulic-mechanical coupling process for EGS based on the framework of EDFM and XFEM
        elem.geomAlpha = 2.0 * elem.aperture / elem.length;
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
    const double tol = 1e-10;
    //采用AABB包围盒剔除完全不可能相交的情况
        // 0) === AABB 剔除阶段 ===
    double min_px = std::min(p.m_x, q.m_x), max_px = std::max(p.m_x, q.m_x);
    double min_py = std::min(p.m_y, q.m_y), max_py = std::max(p.m_y, q.m_y);
    double min_rx = std::min(r.m_x, s.m_x), max_rx = std::max(r.m_x, s.m_x);
    double min_ry = std::min(r.m_y, s.m_y), max_ry = std::max(r.m_y, s.m_y);

    // 如果AABB不相交，立即返回false
    if (max_px < min_rx || max_rx < min_px ||
        max_py < min_ry || max_ry < min_py)
    {
        return false;
    }

    // 随后进入精确几何检测阶段--向量叉积法
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

/*====判断交点是否在基岩网格内，并找到包含该点的 Cell（findContainingCell）====*/
int Fracture::findContainingCell(const Vector& point, const std::vector<Cell>& cells, const std::unordered_map<int, Node>& nodes)
{
    for (const auto& cell : cells) 
    {
        if (cell.CellNodeIDs.size() != 3) 
            continue;
        Vector a = nodes.at(cell.CellNodeIDs[0]).coord;
        Vector b = nodes.at(cell.CellNodeIDs[1]).coord;
        Vector c = nodes.at(cell.CellNodeIDs[2]).coord;
        if (pointInTriangle(point, a, b, c)) 
            return cell.id;
        else if (pointInTriangle(point, a, b, c))
        {
            std::cerr << "中点落在单元边界或顶点上，建议加密网格或调整裂缝。" << std::endl;
        }
    }
    return -1;
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
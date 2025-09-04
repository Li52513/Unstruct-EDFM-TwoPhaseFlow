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
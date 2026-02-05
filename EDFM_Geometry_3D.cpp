#include "EDFM_Geometry_3D.h"
#include <iostream>
#include <cmath>
#include <algorithm>

// =========================================================
// Type 0: Coplanar Polygon Overlap Implementation 实现多边形裁剪算法。这里我们采用稳健的“投影-裁剪-反投影”策略
// =========================================================
bool EDFM_Geometry_3D::GetCoplanarPolygonIntersection(
	const std::vector<Vector>& poly1,
	const std::vector<Vector>& poly2,
	const Vector& normal,
	std::vector<Vector>& outPoints)
{
	if (poly1.size() < 3 || poly2.size() < 3) return false;

	int axis = GetDominantAxis(normal);
	auto Project = [&](const Vector& v) -> Vector {
		if (axis == 0) return Vector(v.m_y, v.m_z, 0);
		if (axis == 1) return Vector(v.m_x, v.m_z, 0);
		return Vector(v.m_x, v.m_y, 0);
		};

	std::vector<Vector> p1_2d, p2_2d;
	for (const auto& v : poly1) p1_2d.push_back(Project(v));
	for (const auto& v : poly2) p2_2d.push_back(Project(v));

	std::vector<Vector> intersection2d;

	// [CRITICAL UPDATE] 增强的顶点包含检测 (含边界)
	// A. Poly1 vertices in Poly2
	for (const auto& v : p1_2d) {
		if (IsPointInPolygon2D(v, p2_2d)) intersection2d.push_back(v);
	}
	// B. Poly2 vertices in Poly1
	for (const auto& v : p2_2d) {
		if (IsPointInPolygon2D(v, p1_2d)) intersection2d.push_back(v);
	}

	// C. Edge Intersections (Full Checking)
	for (size_t i = 0; i < p1_2d.size(); ++i) {
		Vector A = p1_2d[i];
		Vector B = p1_2d[(i + 1) % p1_2d.size()];

		for (size_t j = 0; j < p2_2d.size(); ++j) {
			Vector C = p2_2d[j];
			Vector D = p2_2d[(j + 1) % p2_2d.size()];

			// 检测线段相交
			double det = (B.m_x - A.m_x) * (D.m_y - C.m_y) - (D.m_x - C.m_x) * (B.m_y - A.m_y);

			if (std::abs(det) > 1e-10) {
				// 非平行，求交点
				double t1 = ((C.m_x - A.m_x) * (D.m_y - C.m_y) - (D.m_x - C.m_x) * (C.m_y - A.m_y)) / det;
				double t2 = ((C.m_x - A.m_x) * (B.m_y - A.m_y) - (B.m_x - A.m_x) * (C.m_y - A.m_y)) / det;

				// [FIX] 使用宽松判定 (>= -eps, <= 1+eps)
				if (t1 >= -1e-9 && t1 <= 1.0 + 1e-9 && t2 >= -1e-9 && t2 <= 1.0 + 1e-9) {
					Vector inter(A.m_x + t1 * (B.m_x - A.m_x), A.m_y + t1 * (B.m_y - A.m_y), 0);
					intersection2d.push_back(inter);
				}
			}
			else {
				// [NEW] 平行或共线处理：检测端点是否在线段上
				// 解决共线重叠区域漏算问题
				if (IsPointOnSegment(A, C, D)) intersection2d.push_back(A);
				if (IsPointOnSegment(B, C, D)) intersection2d.push_back(B);
				if (IsPointOnSegment(C, A, B)) intersection2d.push_back(C);
				if (IsPointOnSegment(D, A, B)) intersection2d.push_back(D);
			}
		}
	}

	if (intersection2d.empty()) return false;

	// 4. 反投影回 3D
	Vector P0 = poly1[0];
	double planeD = P0 * normal;

	for (const auto& v2 : intersection2d) {
		double x = v2.m_x;
		double y = v2.m_y;
		double z = 0;

		if (axis == 0) { // Drop X
			y = v2.m_x; z = v2.m_y;
			if (std::abs(normal.m_x) > 1e-10)
				x = (planeD - y * normal.m_y - z * normal.m_z) / normal.m_x;
			else x = P0.m_x;
		}
		else if (axis == 1) { // Drop Y
			x = v2.m_x; z = v2.m_y;
			if (std::abs(normal.m_y) > 1e-10)
				y = (planeD - x * normal.m_x - z * normal.m_z) / normal.m_y;
			else y = P0.m_y;
		}
		else { // Drop Z
			x = v2.m_x; y = v2.m_y;
			if (std::abs(normal.m_z) > 1e-10)
				z = (planeD - x * normal.m_x - y * normal.m_y) / normal.m_z;
			else z = P0.m_z;
		}
		outPoints.push_back(Vector(x, y, z));
	}

	return !outPoints.empty();
}

// =========================================================
// Type 1: Point inside Convex Polyhedron
// =========================================================
bool EDFM_Geometry_3D::IsPointInConvexCell(const Vector& point, const Cell& cell, const Mesh& mesh)
{
	// 遍历单元的所有面
	for (int faceID : cell.CellFaceIDs)
	{
		int fIdx = mesh.getFaceIndex(faceID);
		if (fIdx == -1) {
			// 在生产环境中可选择抛出异常或忽略
			continue;
		}
		const Face& face = mesh.getFaces()[fIdx];

		// 1. 获取面上的任意一点 (通常取第一个顶点)
		if (face.FaceNodeCoords.empty()) continue;
		Vector facePoint = face.FaceNodeCoords[0];

		// 2. 确定指向单元内部的法向
		// Face.normal 通常定义为从 Owner 指向 Neighbor
		// 如果当前 cell 是 Owner，则 OutwardNormal = Face.normal
		// 如果当前 cell 是 Neighbor，则 OutwardNormal = -Face.normal
		// 我们需要 InwardNormal = -OutwardNormal
		Vector outwardNormal = face.normal;
		if (face.neighborCell == cell.id) // 当前 Cell 是 Neighbor
		{
			outwardNormal = face.normal * (-1.0);
		}

		// 3. 计算点到面的有向距离 (或点积)
		// 向量 V = P - FaceCenter
		// Distance = V . OutwardNormal
		// 如果点在内部，它必须在所有面的 "内侧" (Distance <= EPSILON)
		// 注意：这里取决于 Normal 指向。标准 FVM 约定 Normal 指向外。
		// 所以点在内部意味着 (P - FacePt) . OutwardNormal <= 0

		Vector v = point - facePoint;
		double dist = v * outwardNormal;

		if (dist > 1e-6)
		{
			// 在某个面的外侧，则必定在单元外
			return false;
		}
	}
	return true; // 在所有面的内侧
}

// =========================================================
// Type 2 & Type 3: Segment-Polygon Intersection
// =========================================================
bool EDFM_Geometry_3D::IntersectSegmentPolygon(
	const Vector& segStart,
	const Vector& segEnd,
	const std::vector<Vector>& facePolygon,
	const Vector& faceNormal,
	Vector& outPoint)
{
	if (facePolygon.size() < 3) return false;

	// 1. 计算线段向量与面的关系
	Vector dir = segEnd - segStart;
	double len = dir.Mag();
	if (len < EPSILON) return false;// 线段太短
	Vector dirNorm = dir / len;

	// 2. 检查线段是否平行于面 (Dot product 接近 0)
	double dot = faceNormal * dirNorm;
	if (std::abs(dot) < EPSILON) return false; // 平行不相交

	// 3. 计算线面交点参数 t
	// Plane equation: (P - V0) . N = 0
	// Line equation: P = Start + t * Dir
	// (Start + t*Dir - V0) . N = 0
	// t = (V0 - Start) . N / (Dir . N)
	Vector v0 = facePolygon[0];
	double t = ((v0 - segStart) * faceNormal) / (dir * faceNormal); // 注意这里用的是未归一化的 dir

	// 4. 检查 t 是否在线段范围内 [0, 1]
	// 允许微小的容差
	if (t < -EPSILON || t > 1.0 + EPSILON) return false;

	// 计算交点坐标
	outPoint = segStart + dir * t;

	// 5. 检查交点是否在多边形内部 (Point in Polygon)
	// 3D 判别比较复杂，通常投影到 2D 主平面进行判断
	int axis = GetDominantAxis(faceNormal);

	// 构造投影后的 2D 多边形和点
	std::vector<Vector> poly2D;
	Vector pt2D;

	auto Project = [&](const Vector& v) -> Vector {
		if (axis == 0) return Vector(v.m_y, v.m_z, 0); // Drop X
		if (axis == 1) return Vector(v.m_x, v.m_z, 0); // Drop Y
		return Vector(v.m_x, v.m_y, 0);                // Drop Z
		};

	for (const auto& v : facePolygon) poly2D.push_back(Project(v));
	pt2D = Project(outPoint);

	return IsPointInPolygon2D(pt2D, poly2D);

}

// Type 3 实现直接复用 Type 2 的逻辑，因为数学本质一样
bool EDFM_Geometry_3D::IntersectPolygonSegment(
	const Vector& segStart,
	const Vector& segEnd,
	const std::vector<Vector>& polyCoords,
	const Vector& polyNormal,
	Vector& outPoint)
{
	return IntersectSegmentPolygon(segStart, segEnd, polyCoords, polyNormal, outPoint);
}


// =========================================================
// Moller-Trumbore 算法实现线段-三角形求交
// =========================================================
bool EDFM_Geometry_3D::IntersectSegmentTriangle(
	const Vector& segStart,
	const Vector& segEnd,
	const Vector& t1, const Vector& t2, const Vector& t3,
	const Vector& triNormal,
	Vector& outPoint)
{
	const double EPS = 1e-9;
	Vector edge1 = t2 - t1;
	Vector edge2 = t3 - t1;

	Vector dir = segEnd - segStart;
	double len = dir.Mag();
	if (len < EPS) return false;
	Vector D = dir / len; // 归一化方向

	Vector h = D & edge2; // Cross product
	double a = edge1 * h; // Dot product

	if (a > -EPS && a < EPS) return false; // 平行

	double f = 1.0 / a;
	Vector s = segStart - t1;
	double u = f * (s * h);

	if (u < -EPS || u > 1.0 + EPS) return false;

	Vector q = s & edge1;
	double v = f * (D * q);

	if (v < -EPS || u + v > 1.0 + EPS) return false;

	// 计算 t (在线段上的比例)
	double t = f * (edge2 * q);

	// 检查 t 是否在线段长度范围内 [0, Length]
	// 注意：上面的 t 是相对于归一化 D 的距离，所以范围是 [0, len]
	if (t >= -EPS && t <= len + EPS)
	{
		outPoint = segStart + D * t;
		return true;
	}

	return false;
}

// =========================================================
// Helper: 2D Point in Polygon (Ray Casting / Winding Number)
// =========================================================
bool EDFM_Geometry_3D::IsPointInPolygon2D(const Vector& pt, const std::vector<Vector>& polygon)
{
	// [NEW] 先检查是否在边界上 (Explicit Boundary Check)
	size_t n = polygon.size();
	for (size_t i = 0; i < n; ++i) {
		size_t next = (i + 1) % n;
		if (IsPointOnSegment(pt, polygon[i], polygon[next])) return true;
	}

	// Ray Casting (Standard)
	bool inside = false;
	size_t j = n - 1;
	for (size_t i = 0; i < n; i++)
	{
		if (((polygon[i].m_y > pt.m_y) != (polygon[j].m_y > pt.m_y)) &&
			(pt.m_x < (polygon[j].m_x - polygon[i].m_x) * (pt.m_y - polygon[i].m_y) /
				(polygon[j].m_y - polygon[i].m_y + 1e-15) + polygon[i].m_x))
		{
			inside = !inside;
		}
		j = i;
	}
	return inside;
}

int EDFM_Geometry_3D::GetDominantAxis(const Vector& normal)
{
	double absX = std::abs(normal.m_x);
	double absY = std::abs(normal.m_y);
	double absZ = std::abs(normal.m_z);

	if (absX >= absY && absX >= absZ) return 0; // X is dominant
	if (absY >= absX && absY >= absZ) return 1; // Y is dominant
	return 2; // Z is dominant
}

// =========================================================
// Helper: 检测 2D 点是否在线段上 (含端点)
// =========================================================
 bool EDFM_Geometry_3D::IsPointOnSegment(const Vector& p, const Vector& a, const Vector& b, double eps)
{
	// 1. 快速包围盒检测
	if (p.m_x < std::min(a.m_x, b.m_x) - eps || p.m_x > std::max(a.m_x, b.m_x) + eps ||
		p.m_y < std::min(a.m_y, b.m_y) - eps || p.m_y > std::max(a.m_y, b.m_y) + eps)
		return false;

	// 2. 共线检测 (叉乘接近0)
	Vector ap = p - a;
	Vector ab = b - a;
	double cross = ap.m_x * ab.m_y - ap.m_y * ab.m_x;
	return std::abs(cross) < eps;
}


 // =========================================================
// [Task 1 实现] Tribox Overlap Algorithm
// Reference: Tomas Akenine-Möller, "Fast 3D Triangle-Box Overlap Testing"
// =========================================================
 bool EDFM_Geometry_3D::TriBoxOverlap(const Vector& boxCenter, const Vector& boxHalfSize,
	 const Vector& t1, const Vector& t2, const Vector& t3)
 {
	 // 1. 将三角形变换到 Box 的局部坐标系 (Box 中心为原点)
	 // 这样做可以简化计算，利用 Box 的对称性
	 Vector v0 = t1 - boxCenter;
	 Vector v1 = t2 - boxCenter;
	 Vector v2 = t3 - boxCenter;

	 // 2. [SAT Test 1] Box Normals (AABB 快速测试)
	 // 检查三角形的 AABB 是否与 Box AABB 分离
	 // X 轴
	 double minX = std::min({ v0.m_x, v1.m_x, v2.m_x });
	 double maxX = std::max({ v0.m_x, v1.m_x, v2.m_x });
	 if (minX > boxHalfSize.m_x || maxX < -boxHalfSize.m_x) return false;

	 // Y 轴
	 double minY = std::min({ v0.m_y, v1.m_y, v2.m_y });
	 double maxY = std::max({ v0.m_y, v1.m_y, v2.m_y });
	 if (minY > boxHalfSize.m_y || maxY < -boxHalfSize.m_y) return false;

	 // Z 轴
	 double minZ = std::min({ v0.m_z, v1.m_z, v2.m_z });
	 double maxZ = std::max({ v0.m_z, v1.m_z, v2.m_z });
	 if (minZ > boxHalfSize.m_z || maxZ < -boxHalfSize.m_z) return false;

	 // 3. [SAT Test 2] Plane/Box Overlap (三角形所在平面)
	 Vector e0 = v1 - v0;
	 Vector e1 = v2 - v1;
	 Vector normal = e0 & e1; // 叉乘得到法向

	 if (!PlaneBoxOverlap(normal, v0, boxHalfSize)) return false;

	 // 4. [SAT Test 3] 9 Cross Products (Edges x Axes)
	 // 这是 Tribox 算法最精髓的部分，用于剔除那些 AABB 相交但实际从“缝隙”中穿过的情况
	 // 我们需要测试 9 个轴：Box的三轴(x,y,z) 与 Triangle的三边(e0,e1,e2) 的叉积

	 double p0, p1, p2, rad, min, max;
	 // e0 = v1 - v0, e1 = v2 - v1, e2 = v0 - v2
	 Vector e2 = v0 - v2;

	 // 宏定义辅助测试逻辑，避免代码冗余
	 // p_i 是顶点在测试轴上的投影
	 // rad 是 Box 在测试轴上的投影半径
	 // 如果 [min, max] 与 [-rad, rad] 不重叠，则分离

	 // --- Test 1: Axis X0 = (1,0,0) x e0 = (0, -e0.z, e0.y) ---
	 // p = v dot axis. axis_x=0, axis_y=-ez, axis_z=ey
	 p0 = v0.m_z * v1.m_y - v0.m_y * v1.m_z; // v0 跟 v1 在此轴投影相同(因为e0=v1-v0)
	 p2 = v2.m_z * e0.m_y - v2.m_y * e0.m_z;
	 rad = boxHalfSize.m_y * std::abs(e0.m_z) + boxHalfSize.m_z * std::abs(e0.m_y);
	 min = std::min(p0, p2); max = std::max(p0, p2);
	 if (min > rad || max < -rad) return false;

	 // --- Test 2: Axis X1 = (1,0,0) x e1 ---
	 p0 = v0.m_z * e1.m_y - v0.m_y * e1.m_z;
	 p1 = v1.m_z * e1.m_y - v1.m_y * e1.m_z;
	 rad = boxHalfSize.m_y * std::abs(e1.m_z) + boxHalfSize.m_z * std::abs(e1.m_y);
	 min = std::min(p0, p1); max = std::max(p0, p1);
	 if (min > rad || max < -rad) return false;

	 // --- Test 3: Axis X2 = (1,0,0) x e2 ---
	 p0 = v0.m_z * e2.m_y - v0.m_y * e2.m_z;
	 p1 = v1.m_z * e2.m_y - v1.m_y * e2.m_z;
	 rad = boxHalfSize.m_y * std::abs(e2.m_z) + boxHalfSize.m_z * std::abs(e2.m_y);
	 min = std::min(p0, p1); max = std::max(p0, p1);
	 if (min > rad || max < -rad) return false;

	 // --- Test 4: Axis Y0 = (0,1,0) x e0 = (e0.z, 0, -e0.x) ---
	 p0 = -v0.m_z * v1.m_x + v0.m_x * v1.m_z;
	 p2 = -v2.m_z * e0.m_x + v2.m_x * e0.m_z;
	 rad = boxHalfSize.m_x * std::abs(e0.m_z) + boxHalfSize.m_z * std::abs(e0.m_x);
	 min = std::min(p0, p2); max = std::max(p0, p2);
	 if (min > rad || max < -rad) return false;

	 // --- Test 5: Axis Y1 = (0,1,0) x e1 ---
	 p0 = -v0.m_z * e1.m_x + v0.m_x * e1.m_z;
	 p1 = -v1.m_z * e1.m_x + v1.m_x * e1.m_z;
	 rad = boxHalfSize.m_x * std::abs(e1.m_z) + boxHalfSize.m_z * std::abs(e1.m_x);
	 min = std::min(p0, p1); max = std::max(p0, p1);
	 if (min > rad || max < -rad) return false;

	 // --- Test 6: Axis Y2 = (0,1,0) x e2 ---
	 p0 = -v0.m_z * e2.m_x + v0.m_x * e2.m_z;
	 p1 = -v1.m_z * e2.m_x + v1.m_x * e2.m_z;
	 rad = boxHalfSize.m_x * std::abs(e2.m_z) + boxHalfSize.m_z * std::abs(e2.m_x);
	 min = std::min(p0, p1); max = std::max(p0, p1);
	 if (min > rad || max < -rad) return false;

	 // --- Test 7: Axis Z0 = (0,0,1) x e0 = (-e0.y, e0.x, 0) ---
	 p0 = -v0.m_y * v1.m_x + v0.m_x * v1.m_y;
	 p2 = -v2.m_y * e0.m_x + v2.m_x * e0.m_y;
	 rad = boxHalfSize.m_x * std::abs(e0.m_y) + boxHalfSize.m_y * std::abs(e0.m_x);
	 min = std::min(p0, p2); max = std::max(p0, p2);
	 if (min > rad || max < -rad) return false;

	 // --- Test 8: Axis Z1 = (0,0,1) x e1 ---
	 p0 = -v0.m_y * e1.m_x + v0.m_x * e1.m_y;
	 p1 = -v1.m_y * e1.m_x + v1.m_x * e1.m_y;
	 rad = boxHalfSize.m_x * std::abs(e1.m_y) + boxHalfSize.m_y * std::abs(e1.m_x);
	 min = std::min(p0, p1); max = std::max(p0, p1);
	 if (min > rad || max < -rad) return false;

	 // --- Test 9: Axis Z2 = (0,0,1) x e2 ---
	 p0 = -v0.m_y * e2.m_x + v0.m_x * e2.m_y;
	 p1 = -v1.m_y * e2.m_x + v1.m_x * e2.m_y;
	 rad = boxHalfSize.m_x * std::abs(e2.m_y) + boxHalfSize.m_y * std::abs(e2.m_x);
	 min = std::min(p0, p1); max = std::max(p0, p1);
	 if (min > rad || max < -rad) return false;

	 return true; // 通过所有测试，确认相交
 }

 // 辅助: 平面与 Box 相交测试
 bool EDFM_Geometry_3D::PlaneBoxOverlap(const Vector& normal, const Vector& vert, const Vector& maxbox)
 {
	 Vector vmin, vmax;

	 // 找到 Box 上距离平面法向最近和最远的点 (p-vertex and n-vertex)
	 if (normal.m_x > 0.0) { vmin.m_x = -maxbox.m_x; vmax.m_x = maxbox.m_x; }
	 else { vmin.m_x = maxbox.m_x;  vmax.m_x = -maxbox.m_x; }

	 if (normal.m_y > 0.0) { vmin.m_y = -maxbox.m_y; vmax.m_y = maxbox.m_y; }
	 else { vmin.m_y = maxbox.m_y;  vmax.m_y = -maxbox.m_y; }

	 if (normal.m_z > 0.0) { vmin.m_z = -maxbox.m_z; vmax.m_z = maxbox.m_z; }
	 else { vmin.m_z = maxbox.m_z;  vmax.m_z = -maxbox.m_z; }

	 // 检查平面位置
	 // Plane eq: normal * (x - vert) = 0  =>  normal * x = normal * vert
	 double d = normal * vert; // Plane distance constant (unnormalized)

	 // 投影间隔 [minD, maxD]
	 double minD = normal * vmin;
	 double maxD = normal * vmax;

	 if (minD > d || maxD < d) return false; // 平面完全在 Box 之外

	 return true;
 }
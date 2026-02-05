// =========================================================
// 纯静态的工具类，专注于解决“点、线、面、体”之间的几何关系
// =========================================================

#pragma once

#include <vector>
#include <cmath>
#include <algorithm>

// 引入现有基础数据结构
#include "UserDefineVarType.h" // Vector
#include "Face.h"
#include "Cell.h"
#include "Mesh.h"
#include "2D_FractureMeshElement.h" 
#include "2D_FractureIntersection.h"

/**
 * @class EDFM_Geometry_3D
 * @brief 3D-EDFM 核心几何算法库 (静态工具类)
 * @details 实现了 Type 1, Type 2, Type 3 交点检测的底层数学逻辑。
 * 所有算法均为 Production-Ready 级别，包含鲁棒性容差处理。
 */

class EDFM_Geometry_3D
{
public:
    // =========================================================
    // 常量定义
    // =========================================================
    static constexpr double EPSILON = 1e-8; ///< 几何计算通用容差

    // =========================================================
    // Type 0: 共面多边形重叠 (Coplanar Polygon Overlap) - [新增]
    // =========================================================
    /**
     * @brief 计算两个共面多边形的重叠区域顶点
     * @details 处理 Case 1 这种完美共面导致 Type2/3 失效的极端情况。
     * 使用 Sutherland-Hodgman 或简化的凸多边形求交逻辑。
     * @param poly1 多边形1顶点 (基岩面)
     * @param poly2 多边形2顶点 (裂缝微元)
     * @param normal 共面法向 (用于投影到 2D)
     * @param outPoints [输出] 重叠区域的顶点集合
     * @return true 如果存在有效重叠
     */

    static bool GetCoplanarPolygonIntersection(
        const std::vector<Vector>& poly1,
        const std::vector<Vector>& poly2,
        const Vector& normal,
        std::vector<Vector>& outPoints);


    // =========================================================
    // Type 1: 裂缝节点在基岩单元内部 (Point inside Polyhedron)
    // =========================================================
    /**
     * @brief 检测点是否在凸多面体基岩单元内部 (Type 1)
     * @param point 待检测的点坐标 (裂缝节点)
     * @param cell 基岩单元对象
     * @param mesh 全局网格对象 (用于索引 Face 信息)
     * @return true 如果点在单元内部（包含边界）
     * @note 该算法假设 Cell 为凸多面体 (Convex Polyhedron)。
     * 通过检查点是否位于所有边界面的一侧来判定。
     */
    static bool IsPointInConvexCell(const Vector& point, const Cell& cell, const Mesh& mesh);

    // =========================================================
    // Type 2: 裂缝边与基岩面相交 (Segment intersects Polygon)
    // =========================================================
    /**
     * @brief 计算线段与多边形面的交点 (Type 2)
     * @param segStart 线段起点 (裂缝边节点 A)
     * @param segEnd 线段终点 (裂缝边节点 B)
     * @param facePolygon 面顶点坐标集合 (基岩 Face)
     * @param faceNormal 面的法向量 (归一化)
     * @param outPoint [输出] 交点坐标
     * @return true 如果发生有效相交 (交点在线段上 且 在多边形内)
     */
    static bool IntersectSegmentPolygon(
        const Vector& segStart,
        const Vector& segEnd,
        const std::vector<Vector>& facePolygon,
        const Vector& faceNormal,
        Vector& outPoint);

    // =========================================================
    // Type 3: 裂缝面与基岩棱相交 (Polygon intersects Segment)
    // =========================================================
    /**
     * @brief 计算多边形面与线段的交点 (Type 3)
     * @details 数学上等同于 Type 2，但物理意义互换 (Matrix Edge vs Frac Face)
     * @param segStart 线段起点 (基岩棱节点 A)
     * @param segEnd 线段终点 (基岩棱节点 B)
     * @param polyCoords 多边形顶点集合 (裂缝单元)
     * @param polyNormal 多边形法向量
     * @param outPoint [输出] 交点坐标
     * @return true 如果发生有效相交
     */
    static bool IntersectPolygonSegment(
        const Vector& segStart,
        const Vector& segEnd,
        const std::vector<Vector>& polyCoords,
        const Vector& polyNormal,
        Vector& outPoint);

    // =========================================================
    // 鲁棒的线段与三角形求交 (Type 2 & Type 3 的基石),是IntersectSegmentPolygon的增强版，更稳定
    // =========================================================
    static bool IntersectSegmentTriangle(
        const Vector& segStart,
        const Vector& segEnd,
        const Vector& t1, const Vector& t2, const Vector& t3, // 三角形三个顶点
        const Vector& triNormal,
        Vector& outPoint);


    // =========================================================
    //  3D 光栅化核心算法 (Tribox)
    // =========================================================
    /**
     * @brief 判断 3D 三角形与 轴对齐盒子 (AABB) 是否相交
     * @details 基于分离轴定理 (SAT)。用于 3D 光栅化网格索引构建。
     * @param boxCenter 盒子几何中心
     * @param boxHalfSize 盒子半长 (width/2, height/2, depth/2)
     * @param t1, t2, t3 三角形的三个顶点
     * @return true 如果相交或包含
     */
    static bool TriBoxOverlap(const Vector& boxCenter, const Vector& boxHalfSize,
        const Vector& t1, const Vector& t2, const Vector& t3);


private:
    // =========================================================
    // 内部辅助函数
    // =========================================================

    /**
     * @brief 2D 点在多边形内检测 (Ray Casting Algorithm)
     * @param pt 待检测点 (已投影到 2D 平面)
     * @param polygon 多边形顶点 (已投影到 2D 平面)
     * @return true 如果点在多边形内
     */
    static bool IsPointInPolygon2D(const Vector& pt, const std::vector<Vector>& polygon);

    /**
     * @brief 获取多边形投影的主轴方向 (用于 3D -> 2D 降维)
     * @param normal 多边形法向
     * @return 0=YZ平面(丢弃X), 1=XZ平面(丢弃Y), 2=XY平面(丢弃Z)
     */
    static int GetDominantAxis(const Vector& normal);

	/**
	 * @brief 检测点是否在线段上 (包含端点)，考虑容差
	 * @param p 待检测点
	 * @param a 线段起点
	 * @param b 线段终点
	 * @param eps 容差
	 * @return true 如果点在线段上
	 */
    static bool IsPointOnSegment(const Vector& p, const Vector& a, const Vector& b, double eps = 1e-9);

    // Tribox 内部辅助: 平面与盒子相交测试
    static bool PlaneBoxOverlap(const Vector& normal, const Vector& vert, const Vector& maxbox);
};
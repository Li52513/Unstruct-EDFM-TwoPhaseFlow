#pragma once
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>
#include "Node.h"
#include "UserDefineVarType.h"

namespace GeomCalculate
{
	// 功能：计算2D向量的叉积 (仅考虑x,y分量)// 参数：a,b - 输入向量// 返回：标量叉积值 (a.x*b.y - a.y*b.x)
	double cross2D(const Vector& a, const Vector& b);

	// 功能：计算三角形面积// 参数：A,B,C - 三角形三个顶点// 返回：三角形面积 (取绝对值)
	double triArea(const Vector& A, const Vector& B, const Vector& C);

	// 功能：计算三角形质心// 参数：A,B,C - 三角形三个顶点// 返回：质心坐标 (三个顶点的平均值)
	Vector triCentroid(const Vector& A, const Vector& B, const Vector& C);

	// 功能：计算点 p 到空间线段 ab 的最短距离// 参数：p - 点坐标, a - 线段起点, b - 线段终点// 返回：最短距离值
	double pointToSegmentDistance_local(const Vector& p, const Vector& a, const Vector& b);
	
	// 功能：计算两条二维线段的交点（如果存在）// 参数：p,p2 - 第一条线段端点，q,q2 - 第二条线段端点，out - 输出交点坐标，eps - 容差// 返回：是否存在交点（true/false）
	bool segSegIntersect2D(const Vector& p, const Vector& p2, const Vector& q, const Vector& q2, Vector& out, double eps = 1e-12);
	
	// 功能：计算线段 s0->s1 与多边形 poly 的所有交点// 参数：poly - 多边形顶点列表，s0,s1 - 线段端点，hits - 输出交点列表，eps - 容差// 返回：交点数量
	int countIntersectionsWithPolygon(const std::vector<Vector>& poly, const Vector& s0, const Vector& s1, std::vector<Vector>& hits, double eps = 1e-12);

	// 功能：用线段 s0->s1 将多边形 poly 分割成两部分// 参数：poly - 多边形顶点列表，s0,s1 - 线段端点，pos - 输出正侧多边形顶点列表，neg - 输出负侧多边形顶点列表，eps - 容差// 返回：无
	void splitPolygonByLine(const std::vector<Vector>& poly, const Vector& s0, const Vector& s1, std::vector<Vector>& pos, std::vector<Vector>& neg, double eps = 1e-12);

	// 功能：计算三角形 C-A-B 上的三点高斯积分，积分核为点到线段 segStart-segEnd 的距离// 参数：V1,V2,V3 - 三角形顶点，segStart,segEnd - 线段端点// 返回：积分结果
	double triGauss3Integral(const Vector& V1, const Vector& V2, const Vector& V3, const Vector& segStart, const Vector& segEnd);

	// 功能：计算二维多边形的质心和面积// 参数：P - 多边形顶点列表，areaAbs - 输出多边形面积（绝对值）// 返回：多边形质心坐标
	Vector polygonCentroid(const std::vector<Vector>& P, double& areaAbs);

	// 功能：对二维多边形 poly 关于线段 segStart-segEnd 做扇形划分并进行三点高斯积分平均// 参数：poly - 多边形顶点列表，segStart,segEnd - 线段端点// 返回：积分平均结果
	double fanGaussAverage(const std::vector<Vector>& poly, const Vector& segStart, const Vector& segEnd);

	// 功能：判断点 P 是否在二维线段 AB 上（含端点）// 参数：P - 待测点，A,B - 线段端点，eps - 容差// 返回：是否在线段上（true/false）
	bool pointOnSegment2D(const Vector& P, const Vector& A, const Vector& B, double eps = 1e-12);

	// 功能：判断点 P 是否在二维多边形内;参数：P - 待测点，polyNodeIDs - 多边形顶点对应的节点ID列表，nodes - 节点ID到节点对象的映射，eps - 容差;返回：是否在多边形内（true/false）
	bool pointInPolygon2D(const Vector& P, const std::vector<int>& polyNodeIDs, const std::unordered_map<int, Node>& nodes, double eps = 1e-12);

	// 功能：计算点 p 到线段 s-e 的最短距离// 参数：p - 点坐标, s - 线段起点, e - 线段终点// 返回：最短距离值
	double pointToSegmentDistance(const Vector& p, const Vector& s, const Vector& e);

	// 功能：判断点 p 是否在三角形 abc 内（严格内，不含边界）// 参数：p - 待测点，a,b,c - 三角形顶点// 返回：是否在三角形内（true/false）
	double pointInTriangle(const Vector& p, const Vector& a, const Vector& b, const Vector& c);
}
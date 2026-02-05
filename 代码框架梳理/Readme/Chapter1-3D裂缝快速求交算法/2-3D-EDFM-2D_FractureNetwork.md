```C++
#pragma once

#include <vector>
#include <string>
#include <iostream>

#include "UserDefineVarType.h"
#include "2D_Fracture.h"
#include "2D_FracIndex.h" // 引入新建立的索引头文件

// =========================================================
// 辅助结构定义
// =========================================================

// =========================================================
// F-F 交互微观线段 (Micro Segment)
// =========================================================
struct IntersectionSegment
{
    Vector start;
    Vector end;
    double length;
    int cellID_1; // 属于 Frac1 的局部单元 ID
    int cellID_2; // 属于 Frac2 的局部单元 ID

    // 构造
    IntersectionSegment(const Vector& s, const Vector& e, int c1, int c2)
        : start(s), end(e), cellID_1(c1), cellID_2(c2)
    {
        length = (e - s).Mag();
    }
};

// =========================================================
// F-F 交互对象 (Macro Intersection Object)
// =========================================================
struct FracFracIntersectionObject
{
    int id;                 ///< 全局交互 ID
    int fracID_1;           ///< 裂缝 1 全局 ID
    int fracID_2;           ///< 裂缝 2 全局 ID

    std::vector<IntersectionSegment> segments; ///< [升级] 包含多段微观线段
    double totalLength;     ///< 总相交长度

    FracFracIntersectionObject() : id(-1), fracID_1(-1), fracID_2(-1), totalLength(0.0) {}
};

/**
 * @brief F-F 求交搜索策略
 */
enum class FFIntersectionStrategy { BruteForce, Octree_Optimized };


class FractureNetwork_2D
{
public:
    // =========================================================
    // 1. 成员属性
    // =========================================================
    std::vector<Fracture_2D> fractures;                  ///< 宏观裂缝集合
    std::vector<FracFracIntersectionObject> ffIntersections; // F-F相交信息

    // 全局索引缓存
    FracElemIndex_2D fracElemIndex;                      ///< 单元全局索引
    bool isIndexValid = false;                           ///< 索引是否有效标记

    // =========================================================
    // 2. 管理行为
    // =========================================================

    FractureNetwork_2D() = default;

    /**
     * @brief 添加裂缝
     * @param frac 2D 裂缝对象
     */
    void addFracture(const Fracture_2D& frac);



    /**
         * @brief 一键划分所有裂缝的网格
         * @param nU 第一方向划分数
         * @param nV 第二方向划分数
         * @param method  非正交修正方法
         */
    void meshAllFractures(int nU, int nV, NormalVectorCorrectionMethod method = NormalVectorCorrectionMethod::OrthogonalCorrection);

    /**
     * @brief 重建全局单元索引
     */
    void rebuildGlobalIndex();

    // =========================================================
    // 3. F-F 求交核心
    // =========================================================
    /**
     * @brief [重构] 检测裂缝-裂缝相交 (Element-wise)
     * @details 支持非共面裂缝，生成微观线段集合
     * @param strategy 策略开关
     */
    void DetectFractureFractureIntersections(FFIntersectionStrategy strategy = FFIntersectionStrategy::Octree_Optimized);

    /**
    * @brief 获取裂缝列表 (const 访问)
    */
    const std::vector<Fracture_2D>& getFractures() const { return fractures; }

    // =========================================================
    // 4. I/O
    // =========================================================
    /**
     * @brief 导出网络信息 (包含每条裂缝的网格和 F-F 交线)
     * @param prefix 文件前缀
     */
    void exportNetworkToTxt(const std::string& prefix) const;

    // 检查接口
    void inspectFractureEdges(int fracID, const std::string& prefix) const;

    /**
     * @brief [新增] 检查并导出 F-F 交线及微观网格对应关系至 CSV
     * @details 输出交线分段 ID、几何坐标，以及该线穿过的 (FracA_Cell, FracB_Cell) ID 对
     * 用于验证微观拓扑连接 (NNC) 的正确性
     * @param prefix 文件名前缀
     */
    void inspectIntersections(const std::string& prefix) const;

private:


    // =========================================================
     // 核心求交算法 (Core Intersection Algorithms)
     // =========================================================

     /**
      * @brief [核心] 两个空间三角形求交
      * @details 计算 Tri1 和 Tri2 的交线段。基于 Moller 算法与区间重叠法。
      * @param p1, p2, p3 三角形 1 的顶点
      * @param q1, q2, q3 三角形 2 的顶点
      * @param outStart 输出交线段起点
      * @param outEnd   输出交线段终点
      * @return true 如果存在有效线段 (Point contact 返回 false)
      */
    bool _intersectTriangleTriangle(const Vector& p1, const Vector& p2, const Vector& p3,
        const Vector& q1, const Vector& q2, const Vector& q3,
        Vector& outStart, Vector& outEnd);

    /**
     * @brief [辅助] 计算三角形与平面的截取线段
     * @param p1, p2, p3 三角形顶点
     * @param planeN 平面法向
     * @param planeD 平面参数 d (N*X + d = 0)
     * @param outI1 输出截取点 1
     * @param outI2 输出截取点 2
     * @return true 如果存在截取线段 (2个交点)
     */
    bool _getTriPlaneIntersection(const Vector& p1, const Vector& p2, const Vector& p3,
        const Vector& planeN, double planeD,
        Vector& outI1, Vector& outI2);



    // [新增] 两个单元 (四边形) 求交
    void _intersectElementElement(const FractureElement_2D& e1, const std::vector<Node>& nodes1,
        const FractureElement_2D& e2, const std::vector<Node>& nodes2,
        std::vector<IntersectionSegment>& outSegs);

};

```

```C++
#include "2D_FractureNetwork.h"
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <limits>
#include <chrono>

// =========================================================
// 内部辅助类：本地化八叉树 (Local Octree for FractureElements)
// =========================================================
namespace {
    struct OctNode {
        AABB box;
        std::vector<int> cellIndices; // 存储单元索引
        OctNode* children[8];
        bool isLeaf;

        OctNode(const AABB& b) : box(b), isLeaf(true) {
            for (int i = 0; i < 8; ++i) children[i] = nullptr;
        }
        ~OctNode() {
            for (int i = 0; i < 8; ++i) if (children[i]) delete children[i];
        }
    };

    class FracElementOctree {
    public:
        OctNode* root;
        const std::vector<FractureElement_2D>& elements;
        int maxDepth;
        int maxItems;

        FracElementOctree(const std::vector<FractureElement_2D>& elems, const AABB& totalBox)
            : elements(elems), maxDepth(8), maxItems(20) {
            root = new OctNode(totalBox);
            for (size_t i = 0; i < elems.size(); ++i) {
                insert(root, (int)i, 0);
            }
        }
        ~FracElementOctree() { delete root; }

        void query(const AABB& queryBox, std::vector<int>& results) {
            queryRec(root, queryBox, results);
        }

    private:
        void insert(OctNode* node, int idx, int depth) {
            // [修正] 使用 overlaps 替代 intersects
            if (!node->box.overlaps(elements[idx].boundingBox)) return;

            if (node->isLeaf) {
                node->cellIndices.push_back(idx);
                if (depth < maxDepth && node->cellIndices.size() >(size_t)maxItems) {
                    split(node);
                    // Re-distribute
                    for (int oldIdx : node->cellIndices) {
                        for (int k = 0; k < 8; ++k) insert(node->children[k], oldIdx, depth + 1);
                    }
                    node->cellIndices.clear();
                }
            }
            else {
                for (int k = 0; k < 8; ++k) insert(node->children[k], idx, depth + 1);
            }
        }

        void split(OctNode* node) {
            node->isLeaf = false;
            // [修正] 手动计算中心点，使用 min/max 替代 minP/maxP
            Vector minP = node->box.min;
            Vector maxP = node->box.max;
            Vector c = (minP + maxP) * 0.5;

            // Generate 8 children boxes
            for (int i = 0; i < 8; ++i) {
                Vector childMin, childMax;
                childMin.m_x = (i & 1) ? c.m_x : minP.m_x;
                childMax.m_x = (i & 1) ? maxP.m_x : c.m_x;
                childMin.m_y = (i & 2) ? c.m_y : minP.m_y;
                childMax.m_y = (i & 2) ? maxP.m_y : c.m_y;
                childMin.m_z = (i & 4) ? c.m_z : minP.m_z;
                childMax.m_z = (i & 4) ? maxP.m_z : c.m_z;

                // 使用构造函数 AABB(min, max) - 假设 AABB 有此构造函数，或 min/max 是 public
                // 根据 AABB.h: AABB(Vector p1, Vector p2) 会自动处理 min/max
                node->children[i] = new OctNode(AABB(childMin, childMax));
            }
        }

        void queryRec(OctNode* node, const AABB& qBox, std::vector<int>& results) {
            // [修正] 使用 overlaps
            if (!node->box.overlaps(qBox)) return;

            if (node->isLeaf) {
                for (int idx : node->cellIndices) {
                    // Double check AABB
                    if (elements[idx].boundingBox.overlaps(qBox)) {
                        results.push_back(idx);
                    }
                }
            }
            else {
                for (int k = 0; k < 8; ++k) queryRec(node->children[k], qBox, results);
            }
        }
    };
} // namespace
// =========================================================
// 管理行为实现
// =========================================================

void FractureNetwork_2D::addFracture(const Fracture_2D& frac)
{
    fractures.push_back(frac);
    isIndexValid = false; // 结构改变，索引失效
}

void FractureNetwork_2D::meshAllFractures(int nU, int nV, NormalVectorCorrectionMethod method)
{
    std::cout << "[FracNet] Meshing all " << fractures.size() << " fractures (" << nU << "x" << nV << ")..." << std::endl;
    for (auto& frac : fractures) {
        frac.MeshFractureSurface(nU, nV,method);
    }
    // 网格改变后，建议重建索引
    rebuildGlobalIndex();

}

void FractureNetwork_2D::rebuildGlobalIndex()
{
    fracElemIndex = buildFracElemIndex_2D(*this);
    isIndexValid = true;
    std::cout << "[FracNet] Global index rebuilt. Total micro elements: " << fracElemIndex.total << std::endl;
}

// 检查接口
void FractureNetwork_2D::inspectFractureEdges(int fracID, const std::string & prefix) const {
    for (const auto& frac : fractures) {
        if (frac.id == fracID) {
            frac.exportEdgesDetailToCSV(prefix);
            return;
        }
    }
}


// =========================================================
// F-F 求交核心实现
// =========================================================
void FractureNetwork_2D::DetectFractureFractureIntersections(FFIntersectionStrategy strategy)
{
    std::cout << "[FracNet] Detecting F-F Intersections (Element-wise)..." << std::endl;
    ffIntersections.clear();

    auto startTime = std::chrono::high_resolution_clock::now();
    long long pairsChecked = 0;

    size_t N = fractures.size();
    int intersectID = 0;

    // 1. 预计算宏观 AABB
    std::vector<AABB> macroBoxes(N);
    for (size_t i = 0; i < N; ++i) {
        if (fractures[i].fracNodes.empty()) continue;
        Vector minP(1e30, 1e30, 1e30), maxP(-1e30, -1e30, -1e30);
        for (const auto& n : fractures[i].fracNodes) {
            if (n.coord.m_x < minP.m_x) minP.m_x = n.coord.m_x;
            if (n.coord.m_y < minP.m_y) minP.m_y = n.coord.m_y;
            if (n.coord.m_z < minP.m_z) minP.m_z = n.coord.m_z;
            if (n.coord.m_x > maxP.m_x) maxP.m_x = n.coord.m_x;
            if (n.coord.m_y > maxP.m_y) maxP.m_y = n.coord.m_y;
            if (n.coord.m_z > maxP.m_z) maxP.m_z = n.coord.m_z;
        }
        macroBoxes[i] = AABB(minP, maxP);
    }

    for (size_t i = 0; i < N; ++i) {
        // Octree 构建 (仅在 Octree 模式下，且针对 Fracture J 构建)
        // 这里为了简化循环逻辑，我们每次对内层循环构建一次，或者只在必要时构建
        // 最佳实践：外层 I，内层 J。如果 J 很大，对 J 建树。
        // 但这里我们对每一对 (I, J) 进行判断。

        for (size_t j = i + 1; j < N; ++j) {

            // 2. 宏观 AABB 粗筛
            if (!macroBoxes[i].overlaps(macroBoxes[j])) continue;

            FracFracIntersectionObject ffObj;
            ffObj.id = intersectID;
            ffObj.fracID_1 = fractures[i].id;
            ffObj.fracID_2 = fractures[j].id;

            const auto& F1 = fractures[i];
            const auto& F2 = fractures[j];

            if (strategy == FFIntersectionStrategy::BruteForce) {
                for (const auto& e1 : F1.fracCells) {
                    if (!e1.boundingBox.overlaps(macroBoxes[j])) continue; // 微元 vs 宏观 预筛
                    for (const auto& e2 : F2.fracCells) {
                        if (e1.boundingBox.overlaps(e2.boundingBox)) {
                            pairsChecked++;
                            _intersectElementElement(e1, F1.fracNodes, e2, F2.fracNodes, ffObj.segments);
                        }
                    }
                }
            }
            else { // Octree_Optimized
                // 为 F2 构建临时八叉树 (仅当 F2 单元数较多时)
                // 注意：这会有构建开销，但对于 N*N 比较，收益巨大
                FracElementOctree octreeF2(F2.fracCells, macroBoxes[j]);
                std::vector<int> candidates;
                candidates.reserve(50);

                for (const auto& e1 : F1.fracCells) {
                    if (!e1.boundingBox.overlaps(macroBoxes[j])) continue;

                    candidates.clear();
                    octreeF2.query(e1.boundingBox, candidates);

                    // =========================================================
                    // 增加去重逻辑，防止跨节点单元被重复计算
                    // =========================================================
                    std::sort(candidates.begin(), candidates.end());
                    auto last = std::unique(candidates.begin(), candidates.end());
                    candidates.erase(last, candidates.end());
                    // =========================================================

                    for (int idx2 : candidates) {
                        const auto& e2 = F2.fracCells[idx2];
                        pairsChecked++;
                        _intersectElementElement(e1, F1.fracNodes, e2, F2.fracNodes, ffObj.segments);
                    }
                }
            }

            if (!ffObj.segments.empty()) {
                double totalLen = 0.0;
                for (auto& s : ffObj.segments) totalLen += s.length;
                ffObj.totalLength = totalLen;
                ffIntersections.push_back(ffObj);
                intersectID++;
            }
        }
    }
    std::cout << "  -> Pairs Checked: " << pairsChecked << std::endl;

    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = endTime - startTime;

    std::cout << "  -> Strategy: " << (strategy == FFIntersectionStrategy::BruteForce ? "Brute Force" : "Octree/AABB Optimized") << std::endl;
    std::cout << "  -> Time Cost: " << elapsed.count() << " s" << std::endl;
    std::cout << "  -> Pairs Checked: " << pairsChecked << std::endl;
    std::cout << "  -> Intersections: " << ffIntersections.size() << " macro pairs." << std::endl;
}

// =========================================================
// 微元求交逻辑
// =========================================================
// 两个单元求交 (4 对三角形测试)
void FractureNetwork_2D::_intersectElementElement(const FractureElement_2D& e1, const std::vector<Node>& nodes1,
    const FractureElement_2D& e2, const std::vector<Node>& nodes2,
    std::vector<IntersectionSegment>& outSegs)
{
    // 定义简单三角形结构
    struct Tri { Vector a, b, c; };
    std::vector<Tri> tris1, tris2;

    auto getP1 = [&](int idx) { return nodes1[e1.nodeIndices[idx]].coord; };
    auto getP2 = [&](int idx) { return nodes2[e2.nodeIndices[idx]].coord; };

    // 拆分 E1 -> 2 Tris
    if (e1.nodeIndices.size() >= 3) tris1.push_back({ getP1(0), getP1(1), getP1(2) });
    if (e1.nodeIndices.size() >= 4) tris1.push_back({ getP1(0), getP1(2), getP1(3) });

    // 拆分 E2 -> 2 Tris
    if (e2.nodeIndices.size() >= 3) tris2.push_back({ getP2(0), getP2(1), getP2(2) });
    if (e2.nodeIndices.size() >= 4) tris2.push_back({ getP2(0), getP2(2), getP2(3) });

    Vector s, e;
    for (const auto& t1 : tris1) {
        for (const auto& t2 : tris2) {
            if (_intersectTriangleTriangle(t1.a, t1.b, t1.c, t2.a, t2.b, t2.c, s, e)) {
                if ((e - s).Mag() > 1e-9) {
                    outSegs.emplace_back(s, e, e1.id, e2.id);
                }
            }
        }
    }
}

// =========================================================
// 三角形求交核心实现 (Robust Moller Variant)
// =========================================================

// Moller 算法实现 (完整版)
bool FractureNetwork_2D::_intersectTriangleTriangle(const Vector& p1, const Vector& p2, const Vector& p3,
    const Vector& q1, const Vector& q2, const Vector& q3,
    Vector& outStart, Vector& outEnd)
{
    const double EPS = 1e-9;

    // 1. Plane 2 方程
    Vector q21 = q2 - q1; Vector q31 = q3 - q1;
    Vector N2 = q21 & q31; double n2Mag = N2.Mag();
    if (n2Mag < EPS) return false;
    N2 = N2 / n2Mag; double d2 = -(N2 * q1);

    // 距离
    double dp1 = (N2 * p1) + d2; double dp2 = (N2 * p2) + d2; double dp3 = (N2 * p3) + d2;
    if ((dp1 > EPS && dp2 > EPS && dp3 > EPS) || (dp1 < -EPS && dp2 < -EPS && dp3 < -EPS)) return false;

    // 2. Plane 1 方程
    Vector p21 = p2 - p1; Vector p31 = p3 - p1;
    Vector N1 = p21 & p31; double n1Mag = N1.Mag();
    if (n1Mag < EPS) return false;
    N1 = N1 / n1Mag; double d1 = -(N1 * p1);

    // 距离
    double dq1 = (N1 * q1) + d1; double dq2 = (N1 * q2) + d1; double dq3 = (N1 * q3) + d1;
    if ((dq1 > EPS && dq2 > EPS && dq3 > EPS) || (dq1 < -EPS && dq2 < -EPS && dq3 < -EPS)) return false;

    // 3. 交线方向
    Vector D = N1 & N2;
    double lenD = D.Mag();
    if (lenD < EPS) return false; // 共面或平行
    Vector Dir = D / lenD;

    // 4. 计算 Tri1 在交线上的区间
    Vector t1s, t1e;
    if (!_getTriPlaneIntersection(p1, p2, p3, N2, d2, t1s, t1e)) return false;

    // 5. 计算 Tri2 在交线上的区间
    Vector t2s, t2e;
    if (!_getTriPlaneIntersection(q1, q2, q3, N1, d1, t2s, t2e)) return false;

    // 6. 区间求交 (投影到 Dir)
    double s1 = t1s * Dir; double e1 = t1e * Dir;
    if (s1 > e1) { std::swap(s1, e1); std::swap(t1s, t1e); }

    double s2 = t2s * Dir; double e2 = t2e * Dir;
    if (s2 > e2) { std::swap(s2, e2); std::swap(t2s, t2e); }

    double start = std::max(s1, s2);
    double end = std::min(e1, e2);

    if (end - start > EPS) {
        // 反投影回 3D
        double ratioS = (start - s1) / (e1 - s1);
        double ratioE = (end - s1) / (e1 - s1);
        outStart = t1s + (t1e - t1s) * ratioS;
        outEnd = t1s + (t1e - t1s) * ratioE;
        return true;
    }

    return false;
}

bool FractureNetwork_2D::_getTriPlaneIntersection(const Vector& p1, const Vector& p2, const Vector& p3,
    const Vector& planeN, double planeD,
    Vector& outI1, Vector& outI2)
{
    const double EPS = 1e-10;
    double d1 = (planeN * p1) + planeD;
    double d2 = (planeN * p2) + planeD;
    double d3 = (planeN * p3) + planeD;

    std::vector<Vector> pts;
    auto checkEdge = [&](const Vector& A, const Vector& B, double dA, double dB) {
        if ((dA > EPS && dB < -EPS) || (dA < -EPS && dB > EPS)) {
            double t = dA / (dA - dB);
            pts.push_back(A + (B - A) * t);
        }
        };
    checkEdge(p1, p2, d1, d2);
    checkEdge(p2, p3, d2, d3);
    checkEdge(p3, p1, d3, d1);

    // 边界点处理
    if (std::abs(d1) <= EPS) pts.push_back(p1);
    if (std::abs(d2) <= EPS) pts.push_back(p2);
    if (std::abs(d3) <= EPS) pts.push_back(p3);

    // 去重
    if (pts.size() > 2) {
        // 简单去重逻辑
        for (size_t i = 0; i < pts.size(); ++i) {
            for (size_t j = i + 1; j < pts.size(); ) {
                if ((pts[i] - pts[j]).Mag() < EPS) pts.erase(pts.begin() + j);
                else ++j;
            }
        }
    }

    if (pts.size() >= 2) {
        outI1 = pts[0];
        outI2 = pts[1];
        return true;
    }
    return false;
}

// =========================================================
// 导出更新
// =========================================================
void FractureNetwork_2D::exportNetworkToTxt(const std::string& prefix) const
{
    if (fractures.empty()) {
        std::cout << "[Warning] FractureNetwork_2D::exportNetworkToTxt - No fractures to export." << std::endl;
    }

    for (const auto& frac : fractures) {
        // 调用 Fracture_2D 类已有的 exportToTxt 成员函数
        frac.exportToTxt(prefix);
    }
    // 导出 F-F 详细信息 (适应新结构)
    std::string ffFileName = prefix + "_FF_Intersections.txt";
    std::ofstream ff(ffFileName);

    if (ff.is_open()) {
        ff << std::fixed << std::setprecision(6);
        // Header: LineID, StartX,Y,Z, EndX,Y,Z, Frac1, Frac2, Cell1, Cell2, SegLength
        for (const auto& obj : ffIntersections) {
            for (const auto& seg : obj.segments) {
                ff << obj.id << " "
                    << seg.start.m_x << " " << seg.start.m_y << " " << seg.start.m_z << " "
                    << seg.end.m_x << " " << seg.end.m_y << " " << seg.end.m_z << " "
                    << obj.fracID_1 << " " << obj.fracID_2 << " "
                    << seg.cellID_1 << " " << seg.cellID_2 << " "
                    << seg.length << "\n";
            }
        }
    }
}

void FractureNetwork_2D::inspectIntersections(const std::string& prefix) const
{
    std::string fileName = prefix + "_Intersection_MicroMap.csv";
    std::ofstream fs(fileName);

    if (!fs.is_open()) {
        std::cerr << "[Error] FractureNetwork_2D::inspectIntersections - Cannot open file: " << fileName << std::endl;
        return;
    }

    // 写入 CSV 表头
    // SegmentID: 每一小段交线的唯一ID (此处使用全局累加计数模拟)
    // MacroID: 属于哪一大组宏观交线
    fs << "MacroIntersectID,Frac1_ID,Frac2_ID,SegLength,"
        << "StartX,StartY,StartZ,EndX,EndY,EndZ,"
        << "CellID_on_Frac1,CellID_on_Frac2\n";

    fs << std::fixed << std::setprecision(8);

    long long totalSegments = 0;

    // 遍历所有宏观交线对象
    for (const auto& obj : ffIntersections) {
        // 遍历该宏观交线下的所有微观线段
        for (const auto& seg : obj.segments) {
            fs << obj.id << ","
                << obj.fracID_1 << "," << obj.fracID_2 << "," << seg.length << ","
                << seg.start.m_x << "," << seg.start.m_y << "," << seg.start.m_z << ","
                << seg.end.m_x << "," << seg.end.m_y << "," << seg.end.m_z << ","
                << seg.cellID_1 << "," << seg.cellID_2 << "\n";

            totalSegments++;
        }
    }

    fs.close();
    std::cout << "[Inspect] Intersection micro-topology map exported to: " << fileName << std::endl;
    std::cout << "          Total micro-segments exported: " << totalSegments << std::endl;
}
```


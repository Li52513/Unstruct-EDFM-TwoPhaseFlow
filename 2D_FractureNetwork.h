#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <iomanip>
#include <unordered_map>

#include "UserDefineVarType.h"
#include "2D_Fracture.h"
#include "2D_FracIndex.h" // 引入新建立的索引头文件
#include "2D_FractureMeshElement.h"
#include "2D_FractureEdge.h"

#include "Node.h"

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

    std::vector<IntersectionSegment> segments; ///< 包含多段微观线段
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
     * @brief 重建全局单元索引 (旧接口，仅用于 FracElemIndex_2D 统计)
     */
    void rebuildGlobalIndex();

    // =========================================================
    // [New] 索引分配核心接口
    // =========================================================
    /**
     * @brief 分配全局求解器索引 (Solver Indices)
     * @details 遍历所有裂缝的所有微观单元，为其分配全局唯一的矩阵行号。
     * 索引范围从 startOffset 开始递增。
     * @param startOffset 起始偏移量 (通常为基岩网格总数 nMatrix)
     * @return 分配后的下一个可用索引 (即 total_DOF)
     */
    int distributeSolverIndices(int startOffset);


    /**
     * @brief [FVM核心] 重建所有裂缝边的 FVM 离散化属性 (Index & Geometry)
     * @details
     * [修改说明] 弃用旧的 Offset 推导逻辑，改为直接读取单元的 solverIndex。
     * 1. 将局部边 (Owner=LocalID) 转换为 全局边 (Owner=SolverIndex)。
     * 2. 补全 d_ON, f_linearInterpolationCoef 等 FVM 专用字段。
     * 结果存储在 globalEdges_ 中。
     * @param nodesMap 全局节点映射 (用于计算 d_ON 的几何距离)
     */
    void rebuildEdgeProperties(const std::unordered_map<int, Node>& nodesMap);

    // =========================================================
    // 3. F-F 求交核心
    // =========================================================
    /**
     * @brief   检测裂缝-裂缝相交 (Element-wise)
     * @details 支持非共面裂缝，生成微观线段集合
     * @param strategy 策略开关
     */
    void DetectFractureFractureIntersections(FFIntersectionStrategy strategy = FFIntersectionStrategy::Octree_Optimized);

    /**
     * @brief 构建 裂缝ID -> Fracture对象指针 的快速查找表
     * @details 消除对 fractures 向量存储顺序的依赖，提供 O(1) 的裂缝对象查询。
     * 用于 TransmissibilitySolver 中快速获取裂缝开度(aperture)等属性。
     * @return 映射表 map<FractureID, const Fracture_2D*>
     */
    std::unordered_map<int, const Fracture_2D*> buildFractureIDMap() const;

    /**
    * @brief 获取裂缝列表 (const 访问)
    */
    const std::vector<Fracture_2D>& getFractures() const { return fractures; }

    /**
     * @brief 获取整个网络中所有的微观裂缝边
     * @return 包含所有边的列表 (用于全局组装)
     */
    std::vector<FractureEdge_2D> getAllFractureEdges() const;

    /**
     * @brief 获取 FVM 求解器专用的全局边列表
     * @return 包含所有已处理好 Index 的裂缝边
     */
    std::vector<FractureEdge_2D>& getGlobalEdges() { return globalEdges_; }
    const std::vector<FractureEdge_2D>& getGlobalEdges() const { return globalEdges_; }

    // =========================================================
    // 4. I/O
    // =========================================================
    /**
     * @brief 导出裂缝网络信息 (包含每条裂缝的网格和 F-F 交线)
     * @param prefix 文件前缀
     */
    void exportNetworkToTxt(const std::string& prefix) const;

    /**
     * @brief 导出特定宏观裂缝裂缝微观裂缝的网格面非正交信息至.csv
     * @param fracID 特定裂缝的ID
     * @param prefix 文件前缀
     */
    void inspectFractureEdges(int fracID, const std::string& prefix) const;

    /**
     * @brief 检查并导出 F-F 交线及微观网格对应关系至 CSV
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

    //  两个单元 (四边形) 求交
    void _intersectElementElement(const FractureElement_2D& e1, const std::vector<Node>& nodes1,
        const FractureElement_2D& e2, const std::vector<Node>& nodes2,
        std::vector<IntersectionSegment>& outSegs);
    
    // FVM 求解器使用的全局边列表 (Flat List)
    std::vector<FractureEdge_2D> globalEdges_;

};

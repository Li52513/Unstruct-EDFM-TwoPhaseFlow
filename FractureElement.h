#pragma once
#include <vector>

#include "FractureCommon.h"


struct FractureElement  ///// 描述裂缝单元（裂缝段）的结构体
{

    // ————— 编号信息 —————  （在哪里）
    int id;                  ///< 裂缝单元序号
    int cellID;              ///< 裂缝单元所在基岩网格单元编号

    /** * @brief 所属宏观裂缝的编号 (Parent Fracture ID)
     * @details 对应 FractureNetwork::fractures 向量中的下标 (0-based)
     * 由 MeshManager::BuildGlobalSystemIndexing() 统一分配
     */
    int parentFractureID = -1;

    /** * @brief 全局线性方程组中的行号 (Row Index in AX=B)
     * @details
     * - 范围: [n_matrix_dofs, n_matrix_dofs + n_fracture_dofs - 1]
     * - 初始化为 -1，由 MeshManager::BuildGlobalSystemIndexing() 统一分配
     * - 作用: 在组装 Jacobian 矩阵时，直接作为 global_row_index 使用，避免运行时计算
     */
	int solverIndex = -1;	 ///< 裂缝单元在求解器中的索引（连续）
    // —————— 几何信息 —————— （有多大）
    double aperture = 1.0;          ///< 裂缝单元开度 w
    double length = 0.0;            ///< 裂缝单元长度 L
    double avgDistance = 0.0;       ///< 平均距离 d
    double geomCI = 0.0;            ///< 纯几何耦合系数  CI_geom = (L·w) / d
    double geomAlpha = 0.0;         ///< 纯几何 α_geom = 2·w / L  可以进一步修正
    //——————裂缝段类型———————
    FractureElementType type = FractureElementType::Conductive;
    ///< 裂缝段类型（阻塞/导流）

    // —————— 物性耦合 & 离散系数 —————— （储存离散系数）
    double alpha_fr = 0.0;   ///< TI 计算用 α_phys = 2·w·k / (μ·L)
    double k_eff = 0.0;   ///< 并联系统有效渗透率
    double aW_fr = 0.0;   ///< 左邻接段导流系数
    double aE_fr = 0.0;   ///< 右邻接段导流系数
    double b_fr = 0.0;   ///< 裂缝–基岩源项系数 CI_phys = θ·A / d
    double aP_fr = 0.0;   ///< 本段总离散系数 aW+aE+b
    // —————— 归一化参数区间 ——————
    double param0 = 0.0;   ///< 本段起始归一化位置
    double param1 = 0.0;   ///< 本段终止归一化位置

	// ———— 裂缝–裂缝 连接信息 ————
    struct FractureNeighbor
    {
        int fracID;         ///< 邻居所属宏观裂缝索引 (Index in FractureNetwork)
        int elemID;         ///< 邻居单元局部 Local ID 1-based
        int solverIndex;    ///< 邻居单元的全局求解器索引 (Row in AX=B)
        int intersectionID; ///< 发生交汇的全局交点 ID (GlobalFFPoint::id)

        // 构造函数
        FractureNeighbor(int fID, int eID, int sIdx, int interID)
            : fracID(fID), elemID(eID), solverIndex(sIdx), intersectionID(interID) {
        }
    };

    // 切向几何邻居 (Intra-Fracture Siblings)
    // 存储同一条裂缝上相邻微元的 Global Solver Index
    // 0: Prev Sibling (Towards Start), 1: Next Sibling (Towards End)
    std::vector<int> tangentNeighbors;

    /// 该列表存储所有与本单元产生物理连接（相交）的其他裂缝单元
    /// 注意：即使是多条裂缝交于一点，所有相关连接都会被平铺记录在这里
    std::vector<FractureNeighbor> neighbors;

    // —————— 辅助标记 ——————
	bool isFFatStart = false; ///< 裂缝段起点是否为裂缝–裂缝交点
    bool isFFatEnd = false;   ///< 终点是否是裂缝交点
    int  gIDstart = -1;       ///< 起点对应的 GlobalFFPoint ID
    int  gIDend = -1;         ///< 终点对应的 GlobalFFPoint ID


    // —— 裂缝–裂缝 TI 交换信息 ——
    struct FFF_Exchange
    {
        int peerFracID{ -1 }; // 对端裂缝编号
        int peerSegID{ -1 };  // 对端裂缝段编号
        int peerGlobalSeg{ -1 };   // 对端全局段索引（便于快速定位）
        int atGlobalFF{ -1 };	  // 所在全局 FF 交点编号
        double TIw{ 0.0 };         // 水相交汇导纳（Star–Delta 后的“边”）
        double TIg{ 0.0 };         // 气相交汇导纳
    };
    std::vector<FFF_Exchange> ffEx;


    FractureElement(int id_, int cellID_, double len = 0.0, double avgDist = 0.0);
};
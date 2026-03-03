#pragma once

#include "MeshManager.h"
#include "2D_FieldManager.h"

/**
 * @class TransmissibilitySolver_2D
 * @brief 2D-EDFM 静态传导率求解器
 * @details
 * 负责计算基于 2D 几何与静态物性的 NNC (基岩-裂缝) 和 F-F (裂缝-裂缝) 传导率。
 * 架构与 3D 版本保持一致，适配 MeshManager (2D) 和 FieldManager_2D。
 *
 * * [Physics Formula]:
 * - NNC Flow (Matrix-Fracture):
 * T = (L_seg * h) / ( d_m / K_mn + w_f / (2 * K_f) )
 * where:
 * L_seg = fracture segment length in cell
 * h     = thickness (assumed 1.0 for 2D)
 * d_m   = distance from matrix center to fracture line
 * K_mn  = matrix permeability projected to fracture normal
 *
 * - FF Flow (Fracture-Fracture Intersection Point):
 * T = h / ( d_1 / (w_1 * K_1) + d_2 / (w_2 * K_2) )
 * where:
 * Intersection is a point in 2D (Contact length = thickness h)
 * d_1, d_2 = distance from segment centroid to intersection point
 * K_1, K_2 = tangential permeability of fractures
 */
class TransmissibilitySolver_2D
{
public:

    /**
     * @brief 计算基岩-基岩 (Matrix-Matrix) 静态面传导率 (Flow & Heat)
     * @details
     * 遍历所有基岩网格的内部面，利用非正交分解提取的正交分量模长 |E| 作为有效面积，
     * 计算相邻基岩网格间流体和热量的稳态传导系数，结果存入面心场。
     * @param meshMgr 网格管理器 (提供拓扑关系与几何距离)
     * @param fieldMgr 场管理器 (提供 K, Lam_m 等物性场及输出目标)
     */
    static void Calculate_Transmissibility_Matrix(const MeshManager& meshMgr, FieldManager_2D& fieldMgr);

    /**
 * @brief 计算宏观裂缝内部相邻单元之间的静态传导率 (FI: Fracture Internal)
 * @details 处理同一根宏观裂缝剖分后的线段单元间的连通性。采用串联阻力模型。
 * @param meshMgr 2D 网格管理器
 * @param fieldMgr 2D 场管理器
 */
    static void Calculate_Transmissibility_FractureInternal(const MeshManager& meshMgr, FieldManager_2D& fieldMgr);

    /**
     * @brief 计算基岩-裂缝 (NNC) 静态传导率 (Flow & Heat)
     * @details
     * 遍历所有基岩单元及其包含的裂缝段，计算传导率并存入 nncFields。
     * 索引顺序遵循 Matrix Cell Global ID 顺序，确保确定性。
     * * @param meshMgr 网格管理器 (提供拓扑关系)
     * @param fieldMgr 场管理器 (提供 K, Phi, Lambda 等物性场及输出目标)
     */
    static void Calculate_Transmissibility_NNC(const MeshManager& meshMgr, FieldManager_2D& fieldMgr);    /**
     * @brief 计算裂缝-裂缝 (FF) 静态传导率 (Flow & Heat)
     * @details
     * 路线 B：以交点 ID 为核心提取 Junction 关联裂缝段，并执行 Star-Delta 变换，
     * 生成扁平化的 FF 传导率对 (Pairwise FF Connections)。
     * 输出场长度等于所有 Junction 组合对数量总和，而非 globalFFPts 数量。
     * @param meshMgr 网格管理器 (提供 FractureNetwork、裂缝段与交点拓扑)
     * @param fieldMgr 场管理器 (提供 K, Phi, Lambda 等物性场及输出目标)
     */
    static void Calculate_Transmissibility_FF(const MeshManager& meshMgr, FieldManager_2D& fieldMgr);

};
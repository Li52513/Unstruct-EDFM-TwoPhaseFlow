#pragma once

#include "3D_MeshManager.h"
#include "3D_FieldManager.h"

/**
 * @class TransmissibilitySolver_3D
 * @brief 3D-EDFM 静态传导率求解器
 * @details 负责计算基于几何与静态物性 (K, Lambda) 的 NNC 和 F-F 传导率。
 * 计算结果为 "几何因子 * 静态物性"，不包含随压力/饱和度变化的流度 (Mobility)。
 * * 物理公式:
 * - T_Flow = (K_eff * Area) / Distance
 * - T_Heat = (Lambda_eff * Area) / Distance
 */

class TransmissibilitySolver_3D
{
public:
    /**
     * @brief 计算基岩-裂缝 (NNC) 静态传导率 (Flow & Heat)
     * @details
     * 1. 渗流: T_NNC_Flow = Area / ( d_m/K_m + w_f/(4*K_f) )  [调和平均]
     * 2. 传热: T_NNC_Heat = Area / ( d_m/Lam_m + w_f/(4*Lam_f) ) [调和平均]
     * 结果存入 FieldManager 的 nncFields 中。
     * * @param meshMgr 网格管理器 (提供 NNC 拓扑与几何参数)
     * @param fieldMgr 场管理器 (提供 K, Phi, Lambda 场数据，及输出目标)
     */
    static void Calculate_Transmissibility_NNC(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr);

    /**
     * @brief 计算裂缝-裂缝 (F-F) 静态传导率 (Flow & Heat)
     * @details
     * 遍历所有裂缝交线段 (Intersection Segments)，计算微元间的传导率。
     * 1. 渗流: T_FF_Flow = (K_avg * Area) / Distance
     * 2. 传热: T_FF_Heat = (Lam_avg * Area) / Distance
     * 结果存入 FieldManager 的 ffFields (实际上是 nncFields 的 FF 部分) 中。
     * * @param meshMgr 网格管理器 (提供 F-F 交线拓扑)
     * @param fieldMgr 场管理器 (提供裂缝物性场，及输出目标)
     */
    static void Calculate_Transmissibility_FF(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr);
};
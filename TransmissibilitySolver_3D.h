#pragma once

#include "3D_MeshManager.h"
#include "3D_FieldManager.h"

/**
 * @class TransmissibilitySolver_3D
 * @brief 3D-EDFM 静态传导率求解器 (High-Performance Optimized)
 * @details
 * 负责计算基于几何与静态物性的 NNC 和 F-F 传导率。
 * * [Optimization Log]:
 * - 采用 SolverIndex 直接访问 FieldManager 内部扁平化数组，消除 Map/Offset 查找开销 (O(1) Access)。
 * - 支持基岩各向异性渗透率张量 (Kxx, Kyy, Kzz) 在裂缝法向上的投影。
 * - 全面采用 OpenMP 并行化加速大规模交互对的计算。
 * * [Physics Formula]:
 * - NNC Flow (Matrix-Fracture):
 * T = A_int / ( d_m / K_mn + w_f / (2 * K_f) )
 * where d_m = distance from matrix center to fracture plane
 * K_mn = matrix permeability projected to fracture normal
 * w_f/2 = distance from fracture surface to center (half-aperture)
 * - FF Flow (Fracture-Fracture intersection):
 * T = L_int / ( d_1 / (w_1 * K_1) + d_2 / (w_2 * K_2) )
 * Strict series resistance model through the intersection line.
 */
class TransmissibilitySolver_3D
{
public:
    /**
     * @brief 计算基岩-裂缝 (NNC) 静态传导率 (Flow & Heat)
     * @details
     * 使用标准的半宽距离模型 (Factor=2) 计算裂缝侧阻力
     * 结果存入 FieldManager 的 nncFields 中。
     * * @param meshMgr 网格管理器 (提供 NNC 拓扑与几何参数)
     * @param fieldMgr 场管理器 (提供 K, Phi, Lambda 场数据，及输出目标)
     */
    static void Calculate_Transmissibility_NNC(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr);

    /**
     * @brief 计算裂缝-裂缝 (F-F) 静态传导率 (Flow & Heat)
     * @details
     * 遍历所有裂缝交线段 (Intersection Segments)，使用改进的串联电阻模型，分别考虑交线两侧裂缝的开度和渗透率。
     * 1. 渗流: T_FF_Flow = (K_avg * Area) / Distance
     * 2. 传热: T_FF_Heat = (Lam_avg * Area) / Distance
     * 结果存入 FieldManager 的 ffFields (实际上是 nncFields 的 FF 部分) 中。
     * * @param meshMgr 网格管理器 (提供 F-F 交线拓扑)
     * @param fieldMgr 场管理器 (提供裂缝物性场，及输出目标)
     */
    static void Calculate_Transmissibility_FF(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr);
};
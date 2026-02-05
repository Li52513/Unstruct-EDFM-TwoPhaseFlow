#pragma once

#include <vector>
#include <string>
#include <memory>

#include "3D_MeshManager.h"
#include "3D_FieldManager.h"

/**
 * @struct FractureGlobalParams
 * @brief 裂缝通用物性参数 (除导流能力外的其他属性)
 * @details 用于补充 Fracture_2D 对象中未存储的热/流体属性。
 */
struct FractureGlobalParams
{
    double phi = 1.0;           ///< 裂缝孔隙度 [-] (通常为 1.0，部分填充裂缝可小于 1.0)
    double rho = 0.0;           ///< 裂缝内岩石/支撑剂密度 [kg/m^3] (若为空裂缝则忽略)
    double cp = 0.0;            ///< 裂缝比热容 [J/(kg*K)]
    double lambda = 0.0;        ///< 裂缝导热系数 [W/(m*K)] (支撑剂或流体背景值)

    double kn_to_kt_ratio = 1.0;///< 法向/切向渗透率比值 (<=1.0). 用于模拟裂缝面伤害(Skin).
    // Kn = Kt * ratio. NNC计算使用 Kn, FF计算使用 Kt.

// 默认构造
    FractureGlobalParams() = default;

    // 全参构造
    FractureGlobalParams(double p, double r, double c, double l, double ratio = 1.0)
        : phi(p), rho(r), cp(c), lambda(l), kn_to_kt_ratio(ratio) {
    }
};

/**
 * @class FractureProperties_3D
 * @brief 3D 裂缝静态物性管理器
 * @details 负责解析宏观裂缝对象 (Fracture_2D)，计算微观渗透率 (K = Fcd/w)，
 * 并将所有物性参数 (K_n, K_t, w, Phi, etc.) 填充到 FieldManager 的裂缝域场中。
 */
class FractureProperties_3D
{
public:
    /**
     * @brief 构造函数
     * @param meshMgr 3D网格管理器 (提供裂缝网络与拓扑索引)
     * @param fieldMgr 3D场管理器 (用于写入裂缝物性场)
     */
    FractureProperties_3D(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr);

    /**
     * @brief 设置裂缝全局通用参数
     * @details 设置那些未存储在 Fracture_2D 对象中的属性 (如热容、孔隙度等)
     */
    void setGlobalProperties(const FractureGlobalParams& params);

    /**
     * @brief 执行裂缝物性初始化
     * @details
     * 1. 遍历 FractureNetwork 中的每一条宏观裂缝。
     * 2. 读取 aperture 和 conductivity。
     * 3. 计算 Kt = conductivity / aperture, Kn = Kt * ratio。
     * 4. 转换单位 (mD -> m^2)。
     * 5. 利用全局索引 offset 将属性写入对应的微元场。
     */
    void InitializeFractureProperties();

private:
    const MeshManager_3D& meshMgr_;
    FieldManager_3D& fieldMgr_;

    FractureGlobalParams globalParams_; ///< 全局通用参数

    // 单位转换常数: 1 mD = 9.869233e-16 m^2
    static constexpr double MD_TO_M2 = 9.869233e-16;
};
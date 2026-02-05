#pragma once

#include <vector>
#include <string>
#include <memory>

#include "3D_MeshManager.h"
#include "3D_FieldManager.h"
#include "SolverContrlStrName_op.h"
#include "CapRelPerm.h" // [新增] 引入相渗/毛管力模型

/**
 * @struct LinearInitParams
 * @brief 线性分布初始化参数包
 * @details 用于描述物理场在空间中的线性分布:
 * Val(x) = refVal + grad_x * (x - x_ref) + grad_y * (y - y_ref) + grad_z * (z - z_ref)
 */
struct LinearInitParams
{
    double refVal = 0.0;        ///< 参考点的值
    double x_ref = 0.0;         ///< 参考点 X 坐标
    double y_ref = 0.0;         ///< 参考点 Y 坐标
    double z_ref = 0.0;         ///< 参考点 Z 坐标

    double grad_x = 0.0;        ///< X 方向梯度
    double grad_y = 0.0;        ///< Y 方向梯度
    double grad_z = 0.0;        ///< Z 方向梯度 (如 rho * g)

    LinearInitParams() = default;
    explicit LinearInitParams(double val) : refVal(val) {}
    LinearInitParams(double val, double zRef, double gradZ)
        : refVal(val), z_ref(zRef), grad_z(gradZ) {
    }
};

/**
 * @class VariableInitializer_3D
 * @brief 3D 主变量场构建与初始化器
 * @details 负责创建基岩和裂缝中的代求变量场 (Pressure, Temperature, Saturation)，
 * 并提供基于几何位置的线性初始化功能。
 * * 核心职责:
 * 1. Create: 使用 createMatrixScalar/createFractureScalar 强制分配内存。
 * 2. Initialize: 计算 P, T, S 初始分布。
 * 3. Derived: 针对两相流，基于初始饱和度计算并填充 Pc 和 Kr 场。
 */
class VariableInitializer_3D
{
public:
    VariableInitializer_3D(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr);

    // =========================================================
    // 单相流初始化 (Single Phase)
    // =========================================================

    /**
     * @brief 初始化单相流主变量
     */
    void InitSinglePhaseState(const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::TemperatureEquation_String& tConfig,
        const LinearInitParams& pInit,
        const LinearInitParams& tInit);

    // =========================================================
    // IMPES 两相流初始化 (Two Phase)
    // =========================================================

    /**
     * @brief 初始化 IMPES 两相流主变量及辅助变量 (Pc, Kr)
     * @param pConfig 压力方程配置
     * @param tConfig 温度方程配置
     * @param sConfig 饱和度方程配置
     * @param pInit 压力初始化参数
     * @param tInit 温度初始化参数
     * @param sInit 饱和度初始化参数
     * @param vgParams VG模型参数 (用于计算初始 Pc)
     * @param rpParams 相渗参数 (用于计算初始 Kr)
     */
    void InitIMPESState(const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::TemperatureEquation_String& tConfig,
        const PhysicalProperties_string_op::SaturationEquation_String& sConfig,
        const LinearInitParams& pInit,
        const LinearInitParams& tInit,
        const LinearInitParams& sInit,
        const VGParams& vgParams,        // [新增]
        const RelPermParams& rpParams);  // [新增]

private:
    // =========================================================
    // 内部通用实现
    // =========================================================

    /**
     * @brief 通用标量场初始化 (Matrix) - 强制创建
     */
    void CreateAndInitMatrixScalar(const std::string& tagName,
        const std::string& oldTagName,
        const std::string& prevTagName,
        const LinearInitParams& params);

    /**
     * @brief 通用标量场初始化 (Fracture) - 强制创建
     */
    void CreateAndInitFractureScalar(const std::string& tagName,
        const std::string& oldTagName,
        const std::string& prevTagName,
        const LinearInitParams& params);

    /**
     * @brief 计算并填充两相流衍生变量 (Pc, Kr)
     * @details 基于已初始化的 Sw 场，计算 Pc, krw, krg 并填入场管理器
     */
    void CalculateDerivedTwoPhaseFields(const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::SaturationEquation_String& sConfig,
        const VGParams& vg,
        const RelPermParams& rp);

    const MeshManager_3D& meshMgr_;
    FieldManager_3D& fieldMgr_;
};
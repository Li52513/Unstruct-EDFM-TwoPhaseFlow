#pragma once

#include <vector>
#include <string>
#include <memory>
#include <cmath>

#include "3D_MeshManager.h"
#include "3D_FieldManager.h"
#include "SolverContrlStrName_op.h" // 使用优化后的配置头文件
#include "WaterPropertyTable.h"
#include "PropertiesSummary.h"

/**
 * @class WaterProperties_3D
 * @brief 3D 水物性管理器 (动态更新 & IAPWS-95)
 * @details 负责管理基岩和裂缝中的水相属性。
 * * 核心功能：
 * - Constant: 全场赋常数
 * - IAPWS-95: 基于 P, T 场查表插值 (WaterPropertyTable)
 * - 压缩系数计算 (dRho/dP, 数值差分)
 * - 有效热物性计算
 */
class WaterProperties_3D
{
public:
    /**
     * @brief 构造函数 (依赖注入模式)
     * @param fieldMgr 3D场管理器
     * @param pConfig 压力方程配置 (提供 pressure_field 名称)
     * @param tConfig 温度方程配置 (提供 temperatue_field 名称)
     */
    WaterProperties_3D(FieldManager_3D& fieldMgr,
        const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::TemperatureEquation_String& tConfig);

    // =========================================================
    // 1. Matrix Domain (基岩) 操作
    // =========================================================

    /**
     * @brief [基岩] 常量赋值模式
     */
    void UpdateMatrix_Constant(const WaterProperties& params);

    /**
     * @brief [基岩] IAPWS-95 插值模式
     * @details 根据配置中的 P 和 T 场，动态更新 rho, mu, cp, lambda, enthalpy
     */
    void UpdateMatrix_IAPWS();

    /**
     * @brief [基岩] 计算水压缩系数 (dRho/dP)
     * @details 采用数值中心差分法计算
     */
    void CalculateCompressibilityMatrix();

    /**
     * @brief [基岩] 计算有效热物性 (Effective Thermal Properties)
     * @details 若为单相水流，计算饱和水状态下的岩石混合热物性
     */
    void UpdateEffectiveThermalPropertiesMatrix();

    // =========================================================
    // 2. Fracture Domain (裂缝) 操作
    // =========================================================

    /**
     * @brief [裂缝] 常量赋值模式
     */
    void UpdateFracture_Constant(const WaterProperties& params);

    /**
     * @brief [裂缝] IAPWS-95 插值模式
     */
    void UpdateFracture_IAPWS();

    /**
     * @brief [裂缝] 计算水压缩系数
     */
    void CalculateCompressibilityFracture();

    /**
     * @brief [裂缝] 计算有效热物性
     */
    void UpdateEffectiveThermalPropertiesFracture();

private:
    FieldManager_3D& fieldMgr_;

    // 注入的配置对象
    PhysicalProperties_string_op::PressureEquation_String pConfig_;
    PhysicalProperties_string_op::TemperatureEquation_String tConfig_;

    // 数值差分步长 [Pa]
    const double deltaP_ = 100.0;
};
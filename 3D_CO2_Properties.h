#pragma once

#include <vector>
#include <string>
#include <memory>
#include <cmath>

#include "3D_MeshManager.h"
#include "3D_FieldManager.h"
#include "SolverContrlStrName_op.h" // 使用优化后的配置头文件
#include "CO2PropertyTable.h"

/**
 * @struct CO2PropertyParams
 * @brief CO2 常物性参数包 (用于 Constant 模式或 Fallback)
 */
struct CO2PropertyParams
{
    double rho = 0.0;       ///< 密度 [kg/m^3]
    double mu = 0.0;        ///< 粘度 [Pa.s]
    double cp = 0.0;        ///< 比热容 [J/(kg.K)]
    double lambda = 0.0;    ///< 导热系数 [W/(m.K)]
    double enthalpy = 0.0;  ///< 比焓 [J/kg]

    CO2PropertyParams() = default;
    CO2PropertyParams(double r, double m, double c, double l, double h)
        : rho(r), mu(m), cp(c), lambda(l), enthalpy(h) {
    }
};

/**
 * @class CO2Properties_3D
 * @brief 3D CO2 物性管理器 (动态更新 & 状态方程)
 * @details 负责管理基岩和裂缝中的 CO2 相属性。
 * 设计采用了 "求解器无关" (Solver Agnostic) 模式，通过构造函数注入 P/T 场名称配置。
 * * 功能包括：
 * - Constant: 全场赋常数
 * - Span-Wagner: 基于 P, T 场查表插值
 * - 压缩系数计算 (dRho/dP, 数值差分)
 * - 有效热物性计算 (岩石-流体混合)
 */
class CO2Properties_3D
{
public:
    /**
     * @brief 构造函数 (依赖注入模式)
     * @param fieldMgr 3D场管理器
     * @param pConfig 压力方程配置 (提供 pressure_field 名称)
     * @param tConfig 温度方程配置 (提供 temperatue_field 名称)
     */
    CO2Properties_3D(FieldManager_3D& fieldMgr,
        const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::TemperatureEquation_String& tConfig);

    // =========================================================
    // 1. Matrix Domain (基岩) 操作
    // =========================================================

    /**
     * @brief [基岩] 常量赋值模式
     * @param params CO2 常数参数
     */
    void UpdateMatrix_Constant(const CO2PropertyParams& params);

    /**
     * @brief [基岩] Span-Wagner 插值模式
     * @details 根据配置中的 P 和 T 场，动态更新 rho, mu, cp, lambda, enthalpy
     */
    void UpdateMatrix_SpanWagner();

    /**
     * @brief [基岩] 计算 CO2 压缩系数 (dRho/dP)
     * @details 采用数值中心差分法计算，结果存入 drho_g_dp 场
     */
    void CalculateCompressibilityMatrix();

    /**
     * @brief [基岩] 计算有效热物性 (Effective Thermal Properties)
     * @details
     * C_eff = phi*rho_g*cp_g + (1-phi)*rho_r*cp_r
     * L_eff = phi*lam_g + (1-phi)*lam_r
     */
    void UpdateEffectiveThermalPropertiesMatrix();

    // =========================================================
    // 2. Fracture Domain (裂缝) 操作
    // =========================================================

    /**
     * @brief [裂缝] 常量赋值模式
     */
    void UpdateFracture_Constant(const CO2PropertyParams& params);

    /**
     * @brief [裂缝] Span-Wagner 插值模式
     */
    void UpdateFracture_SpanWagner();

    /**
     * @brief [裂缝] 计算 CO2 压缩系数
     */
    void CalculateCompressibilityFracture();

    /**
     * @brief [裂缝] 计算有效热物性
     */
    void UpdateEffectiveThermalPropertiesFracture();

private:
    FieldManager_3D& fieldMgr_;

    // 注入的配置对象 (存储场名称)
    PhysicalProperties_string_op::PressureEquation_String pConfig_;
    PhysicalProperties_string_op::TemperatureEquation_String tConfig_;

    // 数值差分步长 [Pa]
    const double deltaP_ = 100.0;
};
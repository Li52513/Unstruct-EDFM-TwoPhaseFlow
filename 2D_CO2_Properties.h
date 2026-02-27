/**
 * @file 2D_CO2_Properties.h
 * @brief 2D CO2流体物理属性管理子模块
 * @details 负责 2D-EDFM 框架下 CO2 相的标量物性动态更新。
 * 提供常物性模式、基于 Span-Wagner 状态方程的查表模式，及压缩系数和有效热物性的计算。
 */

#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <omp.h>
#include <cmath>
#include <algorithm>

#include "2D_FieldManager.h"
#include "SolverContrlStrName_op.h"
#include "PropertiesSummary.h"
#include "CO2PropertyTable.h"

 /**
  * @class CO2Properties_2D
  * @brief 2D CO2相物性管理器 (IMPES 标量模式专用)
  */
class CO2Properties_2D
{
public:
    /**
     * @brief 构造函数
     * @param fieldMgr 2D场管理器引用
     * @param pConfig 压力方程配置 (提供系统字典)
     * @param tConfig 温度方程配置
     */
    CO2Properties_2D(FieldManager_2D& fieldMgr,
        const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::TemperatureEquation_String& tConfig);

    ~CO2Properties_2D() = default;

    // =========================================================
    // 基础物性更新
    // =========================================================

    /** @brief 基岩域：常物性更新 */
    void UpdateMatrix_Constant(const CO2Properties& params);

    /** @brief 裂缝域：常物性更新 */
    void UpdateFracture_Constant(const CO2Properties& params);

    /** @brief 基岩域：基于 Span-Wagner 查表更新 (使用 p_g 作为计算压力) */
    void UpdateMatrix_SpanWagner();

    /** @brief 裂缝域：基于 Span-Wagner 查表更新 (使用 p_g 作为计算压力) */
    void UpdateFracture_SpanWagner();

    // =========================================================
    // 衍生参数与耦合计算
    // =========================================================

    /** @brief 基岩域：计算 CO2 等温压缩系数 C_g */
    void CalculateCompressibilityMatrix();

    /** @brief 裂缝域：计算 CO2 等温压缩系数 C_g */
    void CalculateCompressibilityFracture();

    /** @brief 基岩域：计算饱和 CO2 状态下的有效热物性 */
    void UpdateEffectiveThermalPropertiesMatrix();

    /** @brief 裂缝域：计算饱和 CO2 状态下的有效热物性 */
    void UpdateEffectiveThermalPropertiesFracture();

private:
    FieldManager_2D& fieldMgr_;
    PhysicalProperties_string_op::TemperatureEquation_String tConfig_;
    std::string pg_name_; ///< CO2相专用压力场名称 (p_g = p_w + P_c)
};
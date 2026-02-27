/**
 * @file 2D_Water_Properties.h
 * @brief 2D 水相流体物理属性管理子模块
 * @details 负责 2D-EDFM 框架下水相 (Water) 的标量物性动态更新。
 * 提供常物性模式、基于 IAPWS-95 状态方程的查表模式，以及水相压缩系数和有效热物性的计算。
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
#include "WaterPropertyTable.h"

 /**
  * @class WaterProperties_2D
  * @brief 2D 水相物性管理器 (IMPES 标量模式专用)
  */
class WaterProperties_2D
{
public:
    /**
     * @brief 构造函数
     * @param fieldMgr 2D场管理器引用
     * @param pConfig 压力方程配置 (用于获取水相压力 p_w)
     * @param tConfig 温度方程配置 (用于获取系统温度 T)
     */
    WaterProperties_2D(FieldManager_2D& fieldMgr,
        const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::TemperatureEquation_String& tConfig);

    ~WaterProperties_2D() = default;

    // =========================================================
    // 基础物性更新 (密度、粘度、比热、比焓、导热系数)
    // =========================================================

    /**
     * @brief 基岩域：常物性更新
     * @param params 水相恒定物理属性包
     */
    void UpdateMatrix_Constant(const WaterProperties& params);

    /**
     * @brief 裂缝域：常物性更新
     * @param params 水相恒定物理属性包
     */
    void UpdateFracture_Constant(const WaterProperties& params);

    /**
     * @brief 基岩域：基于 IAPWS-95 状态方程查表更新
     */
    void UpdateMatrix_IAPWS();

    /**
     * @brief 裂缝域：基于 IAPWS-95 状态方程查表更新
     */
    void UpdateFracture_IAPWS();

    // =========================================================
    // 衍生参数计算 (等温压缩系数)
    // =========================================================

    /**
     * @brief 基岩域：计算水相流体等温压缩系数 C_w = (1/rho) * (d_rho / d_P)
     * @details 采用高精度自适应中心差分数值求导。
     */
    void CalculateCompressibilityMatrix();

    /**
     * @brief 裂缝域：计算水相流体等温压缩系数 C_w
     */
    void CalculateCompressibilityFracture();

    // =========================================================
    // 耦合参数计算 (岩石-流体有效热物性)
    // =========================================================

    /**
     * @brief 基岩域：计算饱和水状态下的有效体积热容与有效导热系数
     */
    void UpdateEffectiveThermalPropertiesMatrix();

    /**
     * @brief 裂缝域：计算饱和水状态下的有效体积热容与有效导热系数
     */
    void UpdateEffectiveThermalPropertiesFracture();

private:
    FieldManager_2D& fieldMgr_;
    PhysicalProperties_string_op::PressureEquation_String pConfig_;
    PhysicalProperties_string_op::TemperatureEquation_String tConfig_;
};
/**
 * @file 2D_PhysicalPropertiesManager.h
 * @brief 2D 物性综合管理器 (Facade & Lifecycle Manager) - 工业级模块化版
 * @details 负责统一管理 2D-EDFM 框架下的基岩、裂缝、CO2、水相的物性子模块。
 * 核心职责：
 * 1. 依赖注入：持久化 2D 网格与场管理器，并向下分发给 4 个子模块。
 * 2. 静态初始化：统一调度基岩与裂缝物性子模块的构建。
 * 3. IMPES 动态更新：调度双精度流体状态方程查表及耦合热物性计算。
 * 4. FIM 动态更新：提供基于 ADVar 的模板化全场物性、雅可比梯度计算与安全组装。
 */

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <omp.h>

 // 引入核心依赖
#include "MeshManager.h"
#include "2D_FieldManager.h"
#include "SolverContrlStrName_op.h"
#include "PropertiesSummary.h"
#include "AABB.h"

// 引入分离后的子模块
#include "2D_RockSolidProperties.h"
#include "2D_FractureProperties.h"
#include "2D_CO2_Properties.h"
#include "2D_Water_Properties.h"

// 引入 AD 评估器
#include "AD_FluidEvaluator.h"

/**
 * @class PhysicalPropertiesManager_2D
 * @brief 2D 物性综合管理器
 */
class PhysicalPropertiesManager_2D
{
public:
    /**
     * @brief 构造函数 (依赖注入模式)
     * @param meshMgr 2D网格管理器
     * @param fieldMgr 2D场管理器
     * @param pConfig 压力方程配置 (用于流体查表字段名称索引)
     * @param tConfig 温度方程配置 (用于流体查表字段名称索引)
     */
    PhysicalPropertiesManager_2D(const MeshManager& meshMgr,
        FieldManager_2D& fieldMgr,
        const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::TemperatureEquation_String& tConfig);

    ~PhysicalPropertiesManager_2D();

    // =========================================================
    // 阶段 1: 静态属性初始化 (Initialization Phase)
    // =========================================================

    /**
     * @brief 初始化基岩静态物性 (背景值 + 局部区域非均质)
     * @param bg 背景岩石参数
     * @param regions 局部区域列表 <区域名称, <包围盒(AABB), 局部岩石参数>>
     */
    void InitRockProperties(const SolidProperties_RockMatrix& bg,
        const std::vector<std::pair<std::string, std::pair<AABB, SolidProperties_RockMatrix>>>& regions);

    /**
     * @brief 初始化裂缝静态物性
     * @param params 裂缝全局参数
     */
    void InitFractureProperties(const SolidProperties_Frac& params);

    // =========================================================
    // 阶段 2: 动态属性更新 (IMPES / Double 标量模式)
    // =========================================================

    /**
     * @brief 更新 CO2 相的标量物性状态 (支持常物性)
     * @param isConstant 若为 true，则全场赋予 CO2 标准物理常数；若为 false，调用 Span-Wagner 查表
     */
    void UpdateFluidEOS_CO2(bool isConstant = false);

    /**
     * @brief 更新水相的标量物性状态 (支持常物性)
     * @param isConstant 若为 true，则全场赋予水标准物理常数；若为 false，调用 IAPWS-95 查表
     */
    void UpdateFluidEOS_Water(bool isConstant = false);

    /**
     * @brief 更新所有流体的标量物性状态 (一键更新封装)
     * @param isConstant 是否统一采用常物性模式
     */
    void UpdateAllFluidEOS(bool isConstant = false);

    /**
     * @brief 更新基于纯 CO2 饱和状态下的岩石-流体有效热物性 (C_eff, Lambda_eff)
     */
    void UpdateEffectiveThermal_CO2();

    /**
     * @brief 更新基于纯水饱和状态下的岩石-流体有效热物性 (C_eff, Lambda_eff)
     */
    void UpdateEffectiveThermal_Water();

    // =========================================================
    // 阶段 3: FIM 全隐式自动微分属性更新 (ADVar 模板模式)
    // =========================================================

    /**
     * @brief FIM 框架下的水相状态方程 AD 场更新 (包含基岩与裂缝)
     * @tparam N 独立自变量数量
     * @param p_field_name 压力场名称
     * @param t_field_name 温度场名称
     * @param isConstant 若为 true，仅赋予恒定初值且梯度严格为 0
     */
    template<int N>
    void UpdateFluidEOS_Water_AD(const std::string& p_field_name, const std::string& t_field_name, bool isConstant = false)
    {
        using namespace PhysicalProperties_string_op;
        Water w_str;

        // ---------------------------------------------------------
        // 1. 基岩域 (Matrix Domain) 更新
        // ---------------------------------------------------------
        auto pF_m = fieldMgr_.getMatrixADScalar<N>(p_field_name);
        auto tF_m = fieldMgr_.getMatrixADScalar<N>(t_field_name);
        if (!pF_m || !tF_m) {
            std::cerr << "[PPM_2D Error] Missing Matrix AD fields for Water EOS update." << std::endl;
            return;
        }

        auto rho_m = fieldMgr_.getOrCreateMatrixADScalar<N>(w_str.rho_tag);
        auto mu_m = fieldMgr_.getOrCreateMatrixADScalar<N>(w_str.mu_tag);
        auto cp_m = fieldMgr_.getOrCreateMatrixADScalar<N>(w_str.cp_tag);
        auto cv_m = fieldMgr_.getOrCreateMatrixADScalar<N>(w_str.cv_tag);
        auto h_m = fieldMgr_.getOrCreateMatrixADScalar<N>(w_str.h_tag);
        auto k_m = fieldMgr_.getOrCreateMatrixADScalar<N>(w_str.k_tag);

        int nMatrixCells = static_cast<int>(meshMgr_.mesh().getCells().size());

#pragma omp parallel for
        for (int i = 0; i < nMatrixCells; ++i) {
            if (isConstant) {
                rho_m->data[i] = ADVar<N>(1000.0);
                mu_m->data[i] = ADVar<N>(1e-3);
                cp_m->data[i] = ADVar<N>(4200.0);
                cv_m->data[i] = ADVar<N>(4182.0);
                h_m->data[i] = ADVar<N>(1.0e5);
                k_m->data[i] = ADVar<N>(0.6);
            }
            else {
                auto props = AD_Fluid::Evaluator::evaluateWater<N>(pF_m->data[i], tF_m->data[i]);
                rho_m->data[i] = props.rho;
                mu_m->data[i] = props.mu;
                cp_m->data[i] = props.cp;
                cv_m->data[i] = props.cv;
                h_m->data[i] = props.h;
                k_m->data[i] = props.k;
            }
        }

        // ---------------------------------------------------------
        // 2. 裂缝域 (Fracture Domain) 更新
        // ---------------------------------------------------------
        auto pF_f = fieldMgr_.getFractureADScalar<N>(p_field_name);
        auto tF_f = fieldMgr_.getFractureADScalar<N>(t_field_name);
        if (!pF_f || !tF_f) return;

        auto rho_f = fieldMgr_.getOrCreateFractureADScalar<N>(w_str.rho_tag);
        auto mu_f = fieldMgr_.getOrCreateFractureADScalar<N>(w_str.mu_tag);
        auto cp_f = fieldMgr_.getOrCreateFractureADScalar<N>(w_str.cp_tag);
        auto cv_f = fieldMgr_.getOrCreateFractureADScalar<N>(w_str.cv_tag);
        auto h_f = fieldMgr_.getOrCreateFractureADScalar<N>(w_str.h_tag);
        auto k_f = fieldMgr_.getOrCreateFractureADScalar<N>(w_str.k_tag);

        const auto& fractures = meshMgr_.fracture_network().fractures;
        size_t totalFracCells = rho_f->data.size();

        // 预计算局部偏移量，绝对安全并行化
        std::vector<size_t> localOffsets(fractures.size() + 1, 0);
        for (size_t i = 0; i < fractures.size(); ++i) {
            localOffsets[i + 1] = localOffsets[i] + fractures[i].elements.size();
        }

        for (size_t i = 0; i < fractures.size(); ++i) {
            size_t baseOffset = localOffsets[i];
            const auto& frac = fractures[i];

#pragma omp parallel for
            for (int localIdx = 0; localIdx < static_cast<int>(frac.elements.size()); ++localIdx) {
                size_t arrayIdx = baseOffset + localIdx;
                if (arrayIdx < totalFracCells) {
                    if (isConstant) {
                        rho_f->data[arrayIdx] = ADVar<N>(1000.0);
                        mu_f->data[arrayIdx] = ADVar<N>(1e-3);
                        cp_f->data[arrayIdx] = ADVar<N>(4200.0);
                        cv_f->data[arrayIdx] = ADVar<N>(4182.0);
                        h_f->data[arrayIdx] = ADVar<N>(1.0e5);
                        k_f->data[arrayIdx] = ADVar<N>(0.6);
                    }
                    else {
                        auto props = AD_Fluid::Evaluator::evaluateWater<N>(pF_f->data[arrayIdx], tF_f->data[arrayIdx]);
                        rho_f->data[arrayIdx] = props.rho;
                        mu_f->data[arrayIdx] = props.mu;
                        cp_f->data[arrayIdx] = props.cp;
                        cv_f->data[arrayIdx] = props.cv;
                        h_f->data[arrayIdx] = props.h;
                        k_f->data[arrayIdx] = props.k;
                    }
                }
            }
        }
    }

    /**
     * @brief FIM 框架下的 CO2 相状态方程 AD 场更新 (包含基岩与裂缝)
     * @tparam N 独立自变量数量
     * @param p_field_name 压力场名称 (如 "p_g")
     * @param t_field_name 温度场名称
     * @param isConstant 若为 true，仅赋予恒定初值且梯度严格为 0
     */
    template<int N>
    void UpdateFluidEOS_CO2_AD(const std::string& p_field_name, const std::string& t_field_name, bool isConstant = false)
    {
        using namespace PhysicalProperties_string_op;
        CO2 co2_str;

        // ---------------------------------------------------------
        // 1. 基岩域 (Matrix Domain) 更新
        // ---------------------------------------------------------
        auto pF_m = fieldMgr_.getMatrixADScalar<N>(p_field_name);
        auto tF_m = fieldMgr_.getMatrixADScalar<N>(t_field_name);
        if (!pF_m || !tF_m) {
            std::cerr << "[PPM_2D Error] Missing Matrix AD fields for CO2 EOS update." << std::endl;
            return;
        }

        auto rho_m = fieldMgr_.getOrCreateMatrixADScalar<N>(co2_str.rho_tag);
        auto mu_m = fieldMgr_.getOrCreateMatrixADScalar<N>(co2_str.mu_tag);
        auto cp_m = fieldMgr_.getOrCreateMatrixADScalar<N>(co2_str.cp_tag);
        auto cv_m = fieldMgr_.getOrCreateMatrixADScalar<N>(co2_str.cv_tag);
        auto h_m = fieldMgr_.getOrCreateMatrixADScalar<N>(co2_str.h_tag);
        auto k_m = fieldMgr_.getOrCreateMatrixADScalar<N>(co2_str.k_tag);

        int nMatrixCells = static_cast<int>(meshMgr_.mesh().getCells().size());

#pragma omp parallel for
        for (int i = 0; i < nMatrixCells; ++i) {
            if (isConstant) {
                rho_m->data[i] = ADVar<N>(800.0);
                mu_m->data[i] = ADVar<N>(1.48e-5);
                cp_m->data[i] = ADVar<N>(1100.0);
                cv_m->data[i] = ADVar<N>(850.0);
                h_m->data[i] = ADVar<N>(3.0e5);
                k_m->data[i] = ADVar<N>(0.03);
            }
            else {
                auto props = AD_Fluid::Evaluator::evaluateCO2<N>(pF_m->data[i], tF_m->data[i]);
                rho_m->data[i] = props.rho;
                mu_m->data[i] = props.mu;
                cp_m->data[i] = props.cp;
                cv_m->data[i] = props.cv;
                h_m->data[i] = props.h;
                k_m->data[i] = props.k;
            }
        }

        // ---------------------------------------------------------
        // 2. 裂缝域 (Fracture Domain) 更新
        // ---------------------------------------------------------
        auto pF_f = fieldMgr_.getFractureADScalar<N>(p_field_name);
        auto tF_f = fieldMgr_.getFractureADScalar<N>(t_field_name);
        if (!pF_f || !tF_f) return;

        auto rho_f = fieldMgr_.getOrCreateFractureADScalar<N>(co2_str.rho_tag);
        auto mu_f = fieldMgr_.getOrCreateFractureADScalar<N>(co2_str.mu_tag);
        auto cp_f = fieldMgr_.getOrCreateFractureADScalar<N>(co2_str.cp_tag);
        auto cv_f = fieldMgr_.getOrCreateFractureADScalar<N>(co2_str.cv_tag);
        auto h_f = fieldMgr_.getOrCreateFractureADScalar<N>(co2_str.h_tag);
        auto k_f = fieldMgr_.getOrCreateFractureADScalar<N>(co2_str.k_tag);

        const auto& fractures = meshMgr_.fracture_network().fractures;
        size_t totalFracCells = rho_f->data.size();

        std::vector<size_t> localOffsets(fractures.size() + 1, 0);
        for (size_t i = 0; i < fractures.size(); ++i) {
            localOffsets[i + 1] = localOffsets[i] + fractures[i].elements.size();
        }

        for (size_t i = 0; i < fractures.size(); ++i) {
            size_t baseOffset = localOffsets[i];
            const auto& frac = fractures[i];

#pragma omp parallel for
            for (int localIdx = 0; localIdx < static_cast<int>(frac.elements.size()); ++localIdx) {
                size_t arrayIdx = baseOffset + localIdx;
                if (arrayIdx < totalFracCells) {
                    if (isConstant) {
                        rho_f->data[arrayIdx] = ADVar<N>(800.0);
                        mu_f->data[arrayIdx] = ADVar<N>(1.48e-5);
                        cp_f->data[arrayIdx] = ADVar<N>(1100.0);
                        cv_f->data[arrayIdx] = ADVar<N>(850.0);
                        h_f->data[arrayIdx] = ADVar<N>(3.0e5);
                        k_f->data[arrayIdx] = ADVar<N>(0.03);
                    }
                    else {
                        auto props = AD_Fluid::Evaluator::evaluateCO2<N>(pF_f->data[arrayIdx], tF_f->data[arrayIdx]);
                        rho_f->data[arrayIdx] = props.rho;
                        mu_f->data[arrayIdx] = props.mu;
                        cp_f->data[arrayIdx] = props.cp;
                        cv_f->data[arrayIdx] = props.cv;
                        h_f->data[arrayIdx] = props.h;
                        k_f->data[arrayIdx] = props.k;
                    }
                }
            }
        }
    }

    // =========================================================
    // 访问器 (Accessors) - 提供底层子模块的精细控制能力
    // =========================================================
    RockSolidProperties_2D& getRockMgr() const { return *rockMgr_; }
    FractureProperties_2D& getFracMgr() const { return *fracMgr_; }
    CO2Properties_2D& getCO2Mgr() const { return *co2Mgr_; }
    WaterProperties_2D& getWaterMgr() const { return *waterMgr_; }

private:
    const MeshManager& meshMgr_;
    FieldManager_2D& fieldMgr_;

    // 子模块组件 (拥有绝对生命周期所有权)
    std::unique_ptr<RockSolidProperties_2D> rockMgr_;
    std::unique_ptr<FractureProperties_2D> fracMgr_;
    std::unique_ptr<CO2Properties_2D> co2Mgr_;
    std::unique_ptr<WaterProperties_2D> waterMgr_;
};
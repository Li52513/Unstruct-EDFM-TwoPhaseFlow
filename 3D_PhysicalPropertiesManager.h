#pragma once

#include <memory>
#include <string>
#include <vector>

// 引入核心依赖
#include "3D_MeshManager.h"
#include "3D_FieldManager.h"
#include "SolverContrlStrName_op.h"

// 引入子模块头文件
#include "3D_RockSolidProperties.h"
#include "3D_FractureProperties.h"
#include "3D_CO2_Properties.h"
#include "3D_Water_Properties.h"


/**
 * @class PhysicalPropertiesManager_3D
 * @brief 3D 物性综合管理器 (Facade & Lifecycle Manager)
 * @details 负责统一管理基岩、裂缝、CO2、水相的物性模块。
 * * 核心职责：
 * 1. 依赖注入：将 MeshManager, FieldManager 和 Config 分发给子模块。
 * 2. 初始化管控：提供静态物性（Rock, Fracture）的初始化入口。
 * 3. 运行时更新：提供流体状态方程 (EOS) 的独立更新接口，支持单相/多相灵活切换。
 */
class PhysicalPropertiesManager_3D
{
public:
    /**
     * @brief 构造函数
     * @param meshMgr 3D网格管理器
     * @param fieldMgr 3D场管理器
     * @param pConfig 压力方程配置 (用于流体查表)
     * @param tConfig 温度方程配置 (用于流体查表)
     */
    PhysicalPropertiesManager_3D(const MeshManager_3D& meshMgr,
        FieldManager_3D& fieldMgr,
        const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::TemperatureEquation_String& tConfig);

    // 析构函数 (PIMPL 模式需在 cpp 实现)
    ~PhysicalPropertiesManager_3D();

    // =========================================================
    // 阶段 1: 静态属性初始化 (Initialization Phase)
    // =========================================================

    /**
     * @brief 初始化基岩静态物性 (背景值 + 局部区域)
     * @param bg 背景岩石参数
     * @param regions 局部区域列表 <名称, <包围盒, 参数>>
     */
    void InitRockProperties(const RockPropertyParams& bg,
        const std::vector<std::pair<std::string, std::pair<BoundingBox3D, RockPropertyParams>>>& regions);

    /**
     * @brief 初始化裂缝静态物性 (导流能力 -> 渗透率)
     * @param params 裂缝全局参数 (比热、密度、各向异性比率等)
     */
    void InitFractureProperties(const FractureGlobalParams& params);

    // =========================================================
    // 阶段 2: 动态属性更新 (Runtime Loop)
    // =========================================================

    /**
     * @brief 更新 CO2 的 EOS 状态 (Span-Wagner)
     * @details 更新基岩与裂缝中 CO2 的密度、粘度、焓、导热系数及压缩系数。
     * 适用于: 单相 CO2 模拟，或两相流中更新气相物性。
     */
    void UpdateFluidEOS_CO2();

    /**
     * @brief 更新水相的 EOS 状态 (IAPWS-95)
     * @details 更新基岩与裂缝中水的密度、粘度、焓、导热系数及压缩系数。
     * 适用于: 单相水模拟，或两相流中更新水相物性。
     */
    void UpdateFluidEOS_Water();

    /**
     * @brief [便利函数] 更新所有流体的 EOS 状态
     * @details 依次调用 UpdateFluidEOS_CO2 和 UpdateFluidEOS_Water。
     */
    void UpdateAllFluidEOS();

    /**
     * @brief 更新 CO2 单相有效热物性 (C_eff, Lambda_eff)
     * @details 假设孔隙完全被 CO2 填充 (S_g = 1.0)。
     * 注意：这将覆盖 C_eff 场。
     */
    void UpdateEffectiveThermal_CO2();

    /**
     * @brief 更新水相单相有效热物性 (C_eff, Lambda_eff)
     * @details 假设孔隙完全被水填充 (S_w = 1.0)。
     * 注意：这将覆盖 C_eff 场。
     */
    void UpdateEffectiveThermal_Water();

    // =========================================================
    // 访问器 (Accessors)
    // =========================================================
    RockSolidProperties_3D& getRockMgr() const { return *rockMgr_; }
    FractureProperties_3D& getFracMgr() const { return *fracMgr_; }
    CO2Properties_3D& getCO2Mgr() const { return *co2Mgr_; }
    WaterProperties_3D& getWaterMgr() const { return *waterMgr_; }

private:
    // 子模块实例 (拥有所有权)
    std::unique_ptr<RockSolidProperties_3D> rockMgr_;
    std::unique_ptr<FractureProperties_3D> fracMgr_;
    std::unique_ptr<CO2Properties_3D> co2Mgr_;
    std::unique_ptr<WaterProperties_3D> waterMgr_;
};
/**
 * @file 2D_PhysicalPropertiesManager.cpp
 * @brief 2D 物性综合管理器具体实现 (统一分发流)
 */

#include "2D_PhysicalPropertiesManager.h"

using namespace PhysicalProperties_string_op;

// =========================================================
// 构造与析构
// =========================================================

PhysicalPropertiesManager_2D::PhysicalPropertiesManager_2D(
    const MeshManager& meshMgr,
    FieldManager_2D& fieldMgr,
    const PressureEquation_String& pConfig,
    const TemperatureEquation_String& tConfig)
    : meshMgr_(meshMgr)
    , fieldMgr_(fieldMgr)
    , rockMgr_(std::make_unique<RockSolidProperties_2D>(meshMgr, fieldMgr))
    , fracMgr_(std::make_unique<FractureProperties_2D>(meshMgr, fieldMgr))
    , co2Mgr_(std::make_unique<CO2Properties_2D>(fieldMgr, pConfig, tConfig))
    , waterMgr_(std::make_unique<WaterProperties_2D>(fieldMgr, pConfig, tConfig))
{
    std::cout << "[PropMgr_2D] Submodule architecture successfully initialized." << std::endl;
}

PhysicalPropertiesManager_2D::~PhysicalPropertiesManager_2D() = default;

// =========================================================
// 1. 基岩与裂缝静态属性初始化 (调度)
// =========================================================

void PhysicalPropertiesManager_2D::InitRockProperties(
    const SolidProperties_RockMatrix& bg,
    const std::vector<std::pair<std::string, std::pair<AABB, SolidProperties_RockMatrix>>>& regions)
{
    rockMgr_->setBackgroundProperties(bg);
    for (const auto& region : regions) {
        rockMgr_->addRegion(region.first, region.second.first, region.second.second);
    }
    rockMgr_->InitializeRockProperties();
}

void PhysicalPropertiesManager_2D::InitFractureProperties(const SolidProperties_Frac& params)
{
    fracMgr_->setGlobalProperties(params);
    fracMgr_->InitializeFractureProperties();
}

// =========================================================
// 2. 纯 Double (IMPES) 动态属性更新 (调度)
// =========================================================

void PhysicalPropertiesManager_2D::UpdateFluidEOS_CO2(bool isConstant)
{
    if (isConstant) {
        // 使用结构体实例聚合初始化，赋予安全的物理兜底参数
        CO2Properties constProps;
        constProps.rho = 800.0;
        constProps.mu = 1.48e-5;
        constProps.cp = 1100.0;
        constProps.cv = 850.0;
        constProps.k = 0.03;
        constProps.h = 3.0e5;

        co2Mgr_->UpdateMatrix_Constant(constProps);
        co2Mgr_->UpdateFracture_Constant(constProps);
    }
    else {
        co2Mgr_->UpdateMatrix_SpanWagner();
        co2Mgr_->UpdateFracture_SpanWagner();
    }

    // 联动触发压缩系数计算
    co2Mgr_->CalculateCompressibilityMatrix();
    co2Mgr_->CalculateCompressibilityFracture();
}

void PhysicalPropertiesManager_2D::UpdateFluidEOS_Water(bool isConstant)
{
    if (isConstant) {
        // 使用结构体实例聚合初始化，赋予安全的物理兜底参数
        WaterProperties constProps;
        constProps.rho = 1000.0;
        constProps.mu = 1e-3;
        constProps.cp = 4200.0;
        constProps.cv = 4182.0;
        constProps.k = 0.6;
        constProps.h = 1.0e5;

        waterMgr_->UpdateMatrix_Constant(constProps);
        waterMgr_->UpdateFracture_Constant(constProps);
    }
    else {
        waterMgr_->UpdateMatrix_IAPWS();
        waterMgr_->UpdateFracture_IAPWS();
    }

    // 联动触发压缩系数计算
    waterMgr_->CalculateCompressibilityMatrix();
    waterMgr_->CalculateCompressibilityFracture();
}

void PhysicalPropertiesManager_2D::UpdateAllFluidEOS(bool isConstant)
{
    UpdateFluidEOS_Water(isConstant);
    UpdateFluidEOS_CO2(isConstant);
}

// =========================================================
// 3. 有效热物理参数更新 (C_eff, Lambda_eff)
// =========================================================

void PhysicalPropertiesManager_2D::UpdateEffectiveThermal_Water()
{
    waterMgr_->UpdateEffectiveThermalPropertiesMatrix();
    waterMgr_->UpdateEffectiveThermalPropertiesFracture();
}

void PhysicalPropertiesManager_2D::UpdateEffectiveThermal_CO2()
{
    co2Mgr_->UpdateEffectiveThermalPropertiesMatrix();
    co2Mgr_->UpdateEffectiveThermalPropertiesFracture();
}
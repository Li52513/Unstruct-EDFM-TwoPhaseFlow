#include "3D_PhysicalPropertiesManager.h"
#include <iostream>

using namespace PhysicalProperties_string_op;

// =========================================================
// 构造与析构
// =========================================================

PhysicalPropertiesManager_3D::PhysicalPropertiesManager_3D(
    const MeshManager_3D& meshMgr,
    FieldManager_3D& fieldMgr,
    const PressureEquation_String& pConfig,
    const TemperatureEquation_String& tConfig)
    : rockMgr_(std::make_unique<RockSolidProperties_3D>(meshMgr, fieldMgr))
    , fracMgr_(std::make_unique<FractureProperties_3D>(meshMgr, fieldMgr))
    , co2Mgr_(std::make_unique<CO2Properties_3D>(fieldMgr, pConfig, tConfig))
    , waterMgr_(std::make_unique<WaterProperties_3D>(fieldMgr, pConfig, tConfig))
{
    // 构造函数体留空，所有初始化已在列表中完成
    std::cout << "[PropMgr] Manager initialized with dependency injection." << std::endl;
}

PhysicalPropertiesManager_3D::~PhysicalPropertiesManager_3D() = default;

// =========================================================
// 静态初始化实现
// =========================================================

void PhysicalPropertiesManager_3D::InitRockProperties(const RockPropertyParams& bg,
    const std::vector<std::pair<std::string, std::pair<BoundingBox3D, RockPropertyParams>>>& regions)
{
    std::cout << "[PropMgr] Initializing Rock Properties..." << std::endl;

    // 1. 设置背景
    rockMgr_->setBackgroundProperties(bg);

    // 2. 添加所有区域
    for (const auto& region : regions) {
        rockMgr_->addRegion(region.first, region.second.first, region.second.second);
    }

    // 3. 执行初始化 (写入场数据)
    rockMgr_->InitializeRockProperties();
}

void PhysicalPropertiesManager_3D::InitFractureProperties(const FractureGlobalParams& params)
{
    std::cout << "[PropMgr] Initializing Fracture Properties..." << std::endl;

    // 1. 设置全局参数
    fracMgr_->setGlobalProperties(params);

    // 2. 执行初始化 (写入场数据)
    fracMgr_->InitializeFractureProperties();
}

// =========================================================
// 动态更新实现
// =========================================================

void PhysicalPropertiesManager_3D::UpdateFluidEOS_CO2()
{
    // CO2 更新 (Matrix & Fracture)
    // Span-Wagner 插值
    co2Mgr_->UpdateMatrix_SpanWagner();
    co2Mgr_->UpdateFracture_SpanWagner();
    // 压缩系数计算
    co2Mgr_->CalculateCompressibilityMatrix();
    co2Mgr_->CalculateCompressibilityFracture();
}

void PhysicalPropertiesManager_3D::UpdateFluidEOS_Water()
{
    // Water 更新 (Matrix & Fracture)
    // IAPWS-95 插值
    waterMgr_->UpdateMatrix_IAPWS();
    waterMgr_->UpdateFracture_IAPWS();
    // 压缩系数计算
    waterMgr_->CalculateCompressibilityMatrix();
    waterMgr_->CalculateCompressibilityFracture();
}

void PhysicalPropertiesManager_3D::UpdateAllFluidEOS()
{
    // 依次调用，保持原有的一键更新功能
    UpdateFluidEOS_CO2();
    UpdateFluidEOS_Water();
}

void PhysicalPropertiesManager_3D::UpdateEffectiveThermal_CO2()
{
    // 计算基于 CO2 饱和 (Sg=1) 的有效热物性
    // 注意：这将覆写 C_eff 和 lambda_eff 场
    co2Mgr_->UpdateEffectiveThermalPropertiesMatrix();
    co2Mgr_->UpdateEffectiveThermalPropertiesFracture();
}

void PhysicalPropertiesManager_3D::UpdateEffectiveThermal_Water()
{
    // 计算基于水饱和 (Sw=1) 的有效热物性
    // 注意：这将覆写 C_eff 和 lambda_eff 场
    waterMgr_->UpdateEffectiveThermalPropertiesMatrix();
    waterMgr_->UpdateEffectiveThermalPropertiesFracture();
}
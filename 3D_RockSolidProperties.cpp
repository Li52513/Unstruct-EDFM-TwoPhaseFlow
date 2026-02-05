#include "3D_RockSolidProperties.h"
#include "SolverContrlStrName_op.h"
#include <iostream>
#include <iomanip>

using namespace PhysicalProperties_string_op;

RockSolidProperties_3D::RockSolidProperties_3D(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr)
    : meshMgr_(meshMgr), fieldMgr_(fieldMgr)
{
}

void RockSolidProperties_3D::setBackgroundProperties(const RockPropertyParams& params)
{
    backgroundProps_ = params;
}

void RockSolidProperties_3D::addRegion(const std::string& name, const BoundingBox3D& box, const RockPropertyParams& params)
{
    RockRegion3D region;
    region.name = name;
    region.bounds = box;
    region.props = params;
    regions_.push_back(region);
}

void RockSolidProperties_3D::InitializeRockProperties()
{
    std::cout << "[RockProp] Initializing 3D matrix rock properties..." << std::endl;

    // 1. 获取场引用 (若不存在则创建，初始值为 0.0)
    Rock rockStr;
    // 渗透率 (Permeability)
	auto& Kxx = fieldMgr_.getOrCreateMatrixScalar(rockStr.k_xx_tag)->data;
    auto& Kyy = fieldMgr_.getOrCreateMatrixScalar(rockStr.k_yy_tag)->data;
    auto& Kzz = fieldMgr_.getOrCreateMatrixScalar(rockStr.k_zz_tag)->data;
    // 标量属性
    auto& Phi = fieldMgr_.getOrCreateMatrixScalar(rockStr.phi_tag)->data;
    auto& Rho = fieldMgr_.getOrCreateMatrixScalar(rockStr.rho_tag)->data;
    auto& Cp = fieldMgr_.getOrCreateMatrixScalar(rockStr.cp_tag)->data;
    auto& Lam = fieldMgr_.getOrCreateMatrixScalar(rockStr.lambda_tag)->data;
    auto& Cr = fieldMgr_.getOrCreateMatrixScalar(rockStr.c_r_tag)->data;

    // 2. 获取基岩单元列表
    const auto& cells = meshMgr_.mesh().getCells();
    size_t numCells = cells.size();

    // 统计各区域命中次数 (用于调试)
    std::vector<int> regionCounts(regions_.size(), 0);
    int backgroundCount = 0;

    // 3. 遍历所有单元进行赋值
#pragma omp parallel for
    for (int i = 0; i < numCells; ++i)
    {
        const auto& cell = cells[i];
        Vector center = cell.center;

        // 默认使用背景值
        RockPropertyParams currentP = backgroundProps_;
        bool isRegionFound = false;

        // 检查局部区域 (反向遍历，实现"后盖前"的图层逻辑)
        for (int r = static_cast<int>(regions_.size()) - 1; r >= 0; --r)
        {
            if (regions_[r].bounds.contains(center))
            {
                currentP = regions_[r].props;

                // 非并行下统计
#pragma omp atomic 
                regionCounts[r]++;

                isRegionFound = true;
                break; // 找到最高优先级的区域后立即停止
            }
        }

        if (!isRegionFound) {
            backgroundCount++;
        }

        // 4. 赋值并进行单位转换 (mD -> m^2)
        Kxx[i] = currentP.k_xx * MD_TO_M2;
        Kyy[i] = currentP.k_yy * MD_TO_M2;
        Kzz[i] = currentP.k_zz * MD_TO_M2;

        // 直接赋值其他标量
        Phi[i] = currentP.phi;
        Rho[i] = currentP.rho;
        Cp[i] = currentP.cp;
        Lam[i] = currentP.lambda;
        Cr[i] = currentP.c_r;
    }

    // 5. 输出统计报告
    std::cout << "  -> Total Matrix Cells: " << numCells << std::endl;
    std::cout << "  -> Background Region: " << backgroundCount << " cells" << std::endl;
    for (size_t r = 0; r < regions_.size(); ++r)
    {
        std::cout << "  -> Region [" << regions_[r].name << "]: " << regionCounts[r] << " cells" << std::endl;
    }
    std::cout << "[RockProp] Initialization completed." << std::endl;
}
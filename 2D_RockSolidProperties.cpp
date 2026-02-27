/**
 * @file 2D_RockSolidProperties.cpp
 * @brief 2D 基岩静态物理属性管理子模块实现
 */

#include "2D_RockSolidProperties.h"

using namespace PhysicalProperties_string_op;

RockSolidProperties_2D::RockSolidProperties_2D(const MeshManager& meshMgr, FieldManager_2D& fieldMgr)
    : meshMgr_(meshMgr), fieldMgr_(fieldMgr)
{
}

void RockSolidProperties_2D::setBackgroundProperties(const SolidProperties_RockMatrix& bg)
{
    bg_ = bg;
}

void RockSolidProperties_2D::addRegion(const std::string& name, const AABB& box, const SolidProperties_RockMatrix& props)
{
    regions_.push_back({ name, {box, props} });
}

void RockSolidProperties_2D::InitializeRockProperties()
{
    std::cout << "[RockProps_2D] Initializing Matrix Rock Properties..." << std::endl;
    Rock rockTags;

    // 1. 在 FieldManager 中分配标量场内存，并用背景值赋初值
    auto k_xx = fieldMgr_.createMatrixScalar(rockTags.k_xx_tag, bg_.kxx);
    auto k_yy = fieldMgr_.createMatrixScalar(rockTags.k_yy_tag, bg_.kyy);
    auto phi = fieldMgr_.createMatrixScalar(rockTags.phi_tag, bg_.phi_r);
    auto rho_r = fieldMgr_.createMatrixScalar(rockTags.rho_tag, bg_.rho_r);
    auto cp_r = fieldMgr_.createMatrixScalar(rockTags.cp_tag, bg_.cp_r);
    auto lambda = fieldMgr_.createMatrixScalar(rockTags.lambda_tag, bg_.k_r);
    auto comp = fieldMgr_.createMatrixScalar(rockTags.c_r_tag, bg_.compressibility);

    const auto& cells = meshMgr_.mesh().getCells();
    int nCells = static_cast<int>(cells.size());

    // 2. 内联几何判定函数
    auto isInsideAABB = [](const Vector& p, const AABB& box) {
        return (p.m_x >= box.min.m_x && p.m_x <= box.max.m_x &&
            p.m_y >= box.min.m_y && p.m_y <= box.max.m_y);
        };

    // 3. 多线程遍历网格并应用局部非均质覆盖
#pragma omp parallel for
    for (int i = 0; i < nCells; ++i)
    {
        const Vector& pos = cells[i].center;
        bool inRegion = false;

        // 逆序遍历区域，实现后定义者优先覆盖
        for (auto it = regions_.rbegin(); it != regions_.rend(); ++it)
        {
            if (isInsideAABB(pos, it->second.first)) {
                const auto& localProps = it->second.second;
                k_xx->data[i] = localProps.kxx;
                k_yy->data[i] = localProps.kyy;
                phi->data[i] = localProps.phi_r;
                rho_r->data[i] = localProps.rho_r;
                cp_r->data[i] = localProps.cp_r;
                lambda->data[i] = localProps.k_r;
                comp->data[i] = localProps.compressibility;
                inRegion = true;
                break;
            }
        }

        // 兜底策略：若未落入任何区域，确保写入背景值（防御性重写）
        if (!inRegion) {
            k_xx->data[i] = bg_.kxx;
            k_yy->data[i] = bg_.kyy;
            phi->data[i] = bg_.phi_r;
            rho_r->data[i] = bg_.rho_r;
            cp_r->data[i] = bg_.cp_r;
            lambda->data[i] = bg_.k_r;
            comp->data[i] = bg_.compressibility;
        }
    }
    std::cout << "[RockProps_2D] Matrix Rock Properties initialization completed." << std::endl;
}
/**
 * @file 2D_FractureProperties.cpp
 * @brief 2D 裂缝静态物理属性管理子模块实现
 */

#include "2D_FractureProperties.h"

using namespace PhysicalProperties_string_op;

FractureProperties_2D::FractureProperties_2D(const MeshManager& meshMgr, FieldManager_2D& fieldMgr)
    : meshMgr_(meshMgr), fieldMgr_(fieldMgr)
{
}

void FractureProperties_2D::setGlobalProperties(const SolidProperties_Frac& params)
{
    params_ = params;
}

void FractureProperties_2D::InitializeFractureProperties()
{
    std::cout << "[FracProps_2D] Initializing Fracture Properties..." << std::endl;
    Fracture_string fracTags;

    // 1. 在 FieldManager 中分配裂缝专属标量场内存，赋全局初值
    auto k_n = fieldMgr_.createFractureScalar(fracTags.k_n_tag, params_.permeability);
    auto k_t = fieldMgr_.createFractureScalar(fracTags.k_t_tag, params_.permeability);
    auto phi = fieldMgr_.createFractureScalar(fracTags.phi_tag, params_.phi_f);
    auto ap = fieldMgr_.createFractureScalar(fracTags.aperture_tag, params_.aperture);
    auto rho = fieldMgr_.createFractureScalar(fracTags.rho_tag, params_.rho_f);
    auto cp = fieldMgr_.createFractureScalar(fracTags.cp_tag, params_.cp_f);
    auto lambda = fieldMgr_.createFractureScalar(fracTags.lambda_tag, params_.k_f);
    auto comp = fieldMgr_.createFractureScalar(fracTags.c_r_tag, params_.compressibility);

    const auto& fractures = meshMgr_.fracture_network().fractures;

    // 2. 预计算 O(N) 安全局部偏移量，支撑严格无锁的并行化
    std::vector<size_t> localOffsets(fractures.size() + 1, 0);
    for (size_t i = 0; i < fractures.size(); ++i) {
        localOffsets[i + 1] = localOffsets[i] + fractures[i].elements.size();
    }

    size_t totalFracCells = k_n->data.size();

    // 3. 多线程遍历网格并写入属性
    for (size_t i = 0; i < fractures.size(); ++i) {
        size_t baseOffset = localOffsets[i];
        const auto& frac = fractures[i];

#pragma omp parallel for
        for (int localIdx = 0; localIdx < static_cast<int>(frac.elements.size()); ++localIdx) {
            size_t arrayIdx = baseOffset + localIdx;

            // 边界越界防护
            if (arrayIdx < totalFracCells) {
                k_n->data[arrayIdx] = params_.permeability;
                k_t->data[arrayIdx] = params_.permeability;
                phi->data[arrayIdx] = params_.phi_f;

                // 孔径逻辑：优先从解析好的几何网格对象读取，若无则使用全局参数
                ap->data[arrayIdx] = (frac.elements[localIdx].aperture > 0.0) ?
                    frac.elements[localIdx].aperture : params_.aperture;

                rho->data[arrayIdx] = params_.rho_f;
                cp->data[arrayIdx] = params_.cp_f;
                lambda->data[arrayIdx] = params_.k_f;
                comp->data[arrayIdx] = params_.compressibility;
            }
        }
    }
    std::cout << "[FracProps_2D] Fracture Properties initialization completed." << std::endl;
}
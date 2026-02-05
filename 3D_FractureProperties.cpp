#include "3D_FractureProperties.h"
#include "SolverContrlStrName_op.h"
#include <iostream>
#include <iomanip>

using namespace PhysicalProperties_string_op;

FractureProperties_3D::FractureProperties_3D(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr)
    : meshMgr_(meshMgr), fieldMgr_(fieldMgr)
{
    // 设置默认参数
    globalParams_ = FractureGlobalParams(1.0, 0.0, 0.0, 0.0, 1.0);
}

void FractureProperties_3D::setGlobalProperties(const FractureGlobalParams& params)
{
    globalParams_ = params;
}

void FractureProperties_3D::InitializeFractureProperties()
{
    std::cout << "[FracProp] Initializing 3D fracture properties..." << std::endl;

    // 1. 准备字符串标签
    Fracture_string fracStr;

    // 2. 获取/创建场引用 (Fracture Domain)
    // 渗透率 & 开度
    auto& Kt = fieldMgr_.getOrCreateFractureScalar(fracStr.k_t_tag, 1e-12)->data;
    auto& Kn = fieldMgr_.getOrCreateFractureScalar(fracStr.k_n_tag, 1e-12)->data;
    auto& Ap = fieldMgr_.getOrCreateFractureScalar(fracStr.aperture_tag, 1e-3)->data;

    // 标量属性
    auto& Phi = fieldMgr_.getOrCreateFractureScalar(fracStr.phi_tag, 1.0)->data;
    auto& Rho = fieldMgr_.getOrCreateFractureScalar(fracStr.rho_tag, 0.0)->data;
    auto& Cp = fieldMgr_.getOrCreateFractureScalar(fracStr.cp_tag, 0.0)->data;
    auto& Lam = fieldMgr_.getOrCreateFractureScalar(fracStr.lambda_tag, 0.0)->data;

    // 3. 获取裂缝网络数据
    const auto& frNet = meshMgr_.fracture_network();
    const auto& fractures = frNet.getFractures();
    const auto& indexer = frNet.fracElemIndex;

    // 校验索引有效性
    if (indexer.offset.size() < fractures.size()) {
        std::cerr << "[Error] Fracture Indexer offset size mismatch! Call rebuildGlobalIndex first." << std::endl;
        return;
    }

    size_t totalMicroElements = 0;

    // 4. 遍历所有宏观裂缝
    // #pragma omp parallel for // 裂缝数量通常不多，且内部赋值简单，并行收益取决于微元总数
    for (size_t i = 0; i < fractures.size(); ++i)
    {
        const auto& frac = fractures[i];

        // 4.1 读取宏观属性
        double w = frac.aperture; // [m]
        double Fcd = frac.conductivity; // [mD*m]

        // 4.2 计算渗透率 [mD] -> [m^2]
        // Kt = Fcd / w
        double Kt_val_mD = (w > 1e-15) ? (Fcd / w) : 0.0;
        double Kt_val_m2 = Kt_val_mD * MD_TO_M2;

        // Kn = Kt * ratio
        double Kn_val_m2 = Kt_val_m2 * globalParams_.kn_to_kt_ratio;

        // 4.3 获取全局索引偏移量
        size_t startOffset = indexer.offset[i];
        size_t numElems = frac.fracCells.size();

        // 4.4 批量赋值给该裂缝下的所有微观单元
        for (size_t localIdx = 0; localIdx < numElems; ++localIdx)
        {
            size_t globalIdx = startOffset + localIdx;

            // 几何与渗透率
            Ap[globalIdx] = w;
            Kt[globalIdx] = Kt_val_m2;
            Kn[globalIdx] = Kn_val_m2;

            // 通用标量属性
            Phi[globalIdx] = globalParams_.phi;
            Rho[globalIdx] = globalParams_.rho;
            Cp[globalIdx] = globalParams_.cp;
            Lam[globalIdx] = globalParams_.lambda;
        }

        totalMicroElements += numElems;
    }

    // 5. 输出统计
    std::cout << "  -> Processed " << fractures.size() << " macro fractures." << std::endl;
    std::cout << "  -> Total " << totalMicroElements << " micro elements assigned." << std::endl;
    std::cout << "  -> Kn/Kt Ratio: " << globalParams_.kn_to_kt_ratio << std::endl;
    std::cout << "[FracProp] Initialization completed." << std::endl;
}
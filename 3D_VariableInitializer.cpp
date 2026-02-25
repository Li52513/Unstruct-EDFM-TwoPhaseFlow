#include "3D_VariableInitializer.h"
#include <iostream>
#include <omp.h>

using namespace PhysicalProperties_string_op;

// =========================================================
// 辅助函数: 线性分布计算
// =========================================================
inline double ComputeLinearValue(const Vector& pos, const LinearInitParams& p)
{
    return p.refVal +
        p.grad_x * (pos.m_x - p.x_ref) +
        p.grad_y * (pos.m_y - p.y_ref) +
        p.grad_z * (pos.m_z - p.z_ref);
}

// =========================================================
// 构造函数
// =========================================================
VariableInitializer_3D::VariableInitializer_3D(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr)
    : meshMgr_(meshMgr), fieldMgr_(fieldMgr)
{
}

// =========================================================
// 外部接口实现
// =========================================================

void VariableInitializer_3D::InitSinglePhaseState(const PressureEquation_String& pConfig,
    const TemperatureEquation_String& tConfig,
    const LinearInitParams& pInit,
    const LinearInitParams& tInit)
{
    std::cout << "[VarInit] Initializing Single-Phase State..." << std::endl;

    // 1. 初始化压力场 (P, P_old, P_prev)
    CreateAndInitMatrixScalar(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit);
    CreateAndInitFractureScalar(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit);

    // 2. 初始化温度场 (T, T_old, T_prev)
    CreateAndInitMatrixScalar(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit);
    CreateAndInitFractureScalar(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit);

    std::cout << "[VarInit] Single-Phase Initialization Completed." << std::endl;
}

void VariableInitializer_3D::InitIMPESState(const PressureEquation_String& pConfig,
    const TemperatureEquation_String& tConfig,
    const SaturationEquation_String& sConfig,
    const LinearInitParams& pInit,
    const LinearInitParams& tInit,
    const LinearInitParams& sInit,
    const CapRelPerm::VGParams& vgParams,
    const CapRelPerm::RelPermParams& rpParams)
{
    std::cout << "[VarInit] Initializing IMPES Two-Phase State..." << std::endl;

    // 1. 初始化基础变量 (P, T, Sw)
    // ---------------------------------------------------------
    // 压力
    CreateAndInitMatrixScalar(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit);
    CreateAndInitFractureScalar(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit);

    // 温度
    CreateAndInitMatrixScalar(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit);
    CreateAndInitFractureScalar(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit);

    // 饱和度
    CreateAndInitMatrixScalar(sConfig.saturation, sConfig.saturation_old, sConfig.saturation_prev, sInit);
    CreateAndInitFractureScalar(sConfig.saturation, sConfig.saturation_old, sConfig.saturation_prev, sInit);

    // 2. 初始化 IMPES 衍生变量 (Pc, Kr)
    // ---------------------------------------------------------
    // 这部分逻辑被提取到单独函数，因为需要同时处理 Matrix 和 Fracture 的逻辑，
    // 且需要访问刚初始化的 Sw 场
    CalculateDerivedTwoPhaseFields(pConfig, sConfig, vgParams, rpParams);

    // 3. 气相压力初始化 P_g = P_w + Pc (这里假设 pInit 是 P_w)
    // 注意：CalculateDerivedTwoPhaseFields 已经计算了 Pc，这里可以进一步修正 P_g
    // 为简化，这里先给 P_g 赋初值 P_w，后续在求解循环开始时通常会再次 Update Pc 和 Pg
    CreateAndInitMatrixScalar(pConfig.pressure_g, pConfig.pressure_g_old, pConfig.pressure_g_prev, pInit);
    CreateAndInitFractureScalar(pConfig.pressure_g, pConfig.pressure_g_old, pConfig.pressure_g_prev, pInit);

    std::cout << "[VarInit] IMPES Initialization Completed." << std::endl;
}

// =========================================================
// 内部核心实现
// =========================================================

void VariableInitializer_3D::CreateAndInitMatrixScalar(const std::string& tagName,
    const std::string& oldTagName,
    const std::string& prevTagName,
    const LinearInitParams& params)
{
    // [Updated] 使用 createMatrixScalar 强制分配/重置内存
    auto& field = fieldMgr_.createMatrixScalar(tagName, params.refVal)->data;
    auto& fieldOld = fieldMgr_.createMatrixScalar(oldTagName, params.refVal)->data;
    auto& fieldPrev = fieldMgr_.createMatrixScalar(prevTagName, params.refVal)->data;

    const auto& cells = meshMgr_.mesh().getCells();
    int nCells = static_cast<int>(cells.size());

    if (field.size() != nCells) {
        std::cerr << "[Error] Matrix field size mismatch during init!" << std::endl;
        return;
    }

#pragma omp parallel for
    for (int i = 0; i < nCells; ++i)
    {
        const Vector& center = cells[i].center;
        double val = ComputeLinearValue(center, params);

        field[i] = val;
        fieldOld[i] = val;
        fieldPrev[i] = val;
    }
}

void VariableInitializer_3D::CreateAndInitFractureScalar(const std::string& tagName,
    const std::string& oldTagName,
    const std::string& prevTagName,
    const LinearInitParams& params)
{
    // [Updated] 使用 createFractureScalar 强制分配/重置内存
    auto& field = fieldMgr_.createFractureScalar(tagName, params.refVal)->data;
    auto& fieldOld = fieldMgr_.createFractureScalar(oldTagName, params.refVal)->data;
    auto& fieldPrev = fieldMgr_.createFractureScalar(prevTagName, params.refVal)->data;

    const auto& frNet = meshMgr_.fracture_network();
    const auto& fractures = frNet.getFractures();
    const auto& indexer = frNet.fracElemIndex;

    if (indexer.offset.size() < fractures.size()) return;

    size_t totalFracCells = field.size();

    // 遍历宏观裂缝 -> 微观单元
    for (size_t i = 0; i < fractures.size(); ++i)
    {
        const auto& frac = fractures[i];
        size_t startOffset = indexer.offset[i];
        size_t numElems = frac.fracCells.size();

#pragma omp parallel for
        for (int localIdx = 0; localIdx < static_cast<int>(numElems); ++localIdx)
        {
            size_t globalIdx = startOffset + localIdx;
            const Vector& center = frac.fracCells[localIdx].centroid;
            double val = ComputeLinearValue(center, params);

            if (globalIdx < totalFracCells) {
                field[globalIdx] = val;
                fieldOld[globalIdx] = val;
                fieldPrev[globalIdx] = val;
            }
        }
    }
}

void VariableInitializer_3D::CalculateDerivedTwoPhaseFields(const PressureEquation_String& pConfig,
    const SaturationEquation_String& sConfig,
    const CapRelPerm::VGParams& vg,
    const CapRelPerm::RelPermParams& rp)
{
    // 准备 Tag 名称 (从 PhysicalProperties_string_op 中获取)
    TwoPhaseAux auxTags;
    Water watTags;
    PhysicalProperties_string_op::CO2 co2Tags;

    // ----------------------------
    // 1. Matrix Domain
    // ----------------------------
    {
        // 输入场: Sw
        const auto& Sw_vec = fieldMgr_.getMatrixScalar(sConfig.saturation)->data;

        // 输出场: Pc, krw, krg (强制创建)
        // 注意: pConfig.Pc_field 可能与 auxTags.Pc_tag 相同，以 config 为准
        auto& Pc_vec = fieldMgr_.createMatrixScalar(pConfig.Pc_field, 0.0)->data;
        auto& krw_vec = fieldMgr_.createMatrixScalar(watTags.k_rw_tag, 0.0)->data;
        auto& krg_vec = fieldMgr_.createMatrixScalar(co2Tags.k_rg_tag, 0.0)->data;

        int n = static_cast<int>(Sw_vec.size());
#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            double sw = Sw_vec[i];

            // 计算 Pc
            Pc_vec[i] = pc_vG(sw, vg);

            // 计算 kr
            double krw, krg;
            kr_Mualem_vG(sw, vg, rp, krw, krg);
            krw_vec[i] = krw;
            krg_vec[i] = krg;
        }
    }

    // ----------------------------
    // 2. Fracture Domain
    // ----------------------------
    {
        const auto& Sw_vec = fieldMgr_.getFractureScalar(sConfig.saturation)->data;

        // 裂缝可能使用不同的相对渗透率曲线 (通常裂缝是直线的 kr=Sw)，
        // 但此处遵循 "adaptable to parameters" 的要求，暂时使用传入的 VG/RP 参数。
        // 若需区分，可在接口中增加 fractureVGParams。
        auto& Pc_vec = fieldMgr_.createFractureScalar(pConfig.Pc_field, 0.0)->data;
        auto& krw_vec = fieldMgr_.createFractureScalar(watTags.k_rw_tag, 0.0)->data;
        auto& krg_vec = fieldMgr_.createFractureScalar(co2Tags.k_rg_tag, 0.0)->data;

        int n = static_cast<int>(Sw_vec.size());
#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            double sw = Sw_vec[i];
            Pc_vec[i] = pc_vG(sw, vg);

            double krw, krg;
            kr_Mualem_vG(sw, vg, rp, krw, krg);
            krw_vec[i] = krw;
            krg_vec[i] = krg;
        }
    }
}
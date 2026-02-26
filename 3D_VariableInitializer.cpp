#include "3D_VariableInitializer.h"
#include <iostream>
#include <omp.h>

using namespace PhysicalProperties_string_op;



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
    std::cout << "[VarInit_3D] Initializing Single-Phase State..." << std::endl;

    CreateAndInitMatrixScalar(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit);
    CreateAndInitFractureScalar(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit);

    CreateAndInitMatrixScalar(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit);
    CreateAndInitFractureScalar(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit);

    std::cout << "[VarInit_3D] Single-Phase Initialization Completed." << std::endl;
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
    std::cout << "[VarInit_3D] Initializing IMPES Two-Phase State..." << std::endl;

    // 1. 初始化基础变量 (P, T, Sw)
    CreateAndInitMatrixScalar(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit);
    CreateAndInitFractureScalar(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit);

    CreateAndInitMatrixScalar(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit);
    CreateAndInitFractureScalar(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit);

    CreateAndInitMatrixScalar(sConfig.saturation, sConfig.saturation_old, sConfig.saturation_prev, sInit);
    CreateAndInitFractureScalar(sConfig.saturation, sConfig.saturation_old, sConfig.saturation_prev, sInit);

    // 2. 初始化 IMPES 衍生变量 (Pc, Kr)
    CalculateDerivedTwoPhaseFields(sConfig, vgParams, rpParams);

    // 3. 气相压力初始化 P_g = P_w + Pc (这里为简化，给予同样的基准场初始化，由后续物性更新修正)
    TwoPhaseState_String tpState;
    CreateAndInitMatrixScalar(tpState.p_g_field, tpState.p_g_old_field, tpState.p_g_prev_field, pInit);
    CreateAndInitFractureScalar(tpState.p_g_field, tpState.p_g_old_field, tpState.p_g_prev_field, pInit);

    std::cout << "[VarInit_3D] IMPES Initialization Completed." << std::endl;
}

// =========================================================
// 内部核心实现
// =========================================================

void VariableInitializer_3D::CreateAndInitMatrixScalar(const std::string& tagName,
    const std::string& oldTagName,
    const std::string& prevTagName,
    const LinearInitParams& params)
{
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
    auto& field = fieldMgr_.createFractureScalar(tagName, params.refVal)->data;
    auto& fieldOld = fieldMgr_.createFractureScalar(oldTagName, params.refVal)->data;
    auto& fieldPrev = fieldMgr_.createFractureScalar(prevTagName, params.refVal)->data;

    const auto& frNet = meshMgr_.fracture_network();
    const auto& fractures = frNet.getFractures();
    const auto& indexer = frNet.fracElemIndex;

    if (indexer.offset.size() < fractures.size()) return;

    size_t totalFracCells = field.size();

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

void VariableInitializer_3D::CalculateDerivedTwoPhaseFields(const SaturationEquation_String& sConfig,
    const CapRelPerm::VGParams& vg,
    const CapRelPerm::RelPermParams& rp)
{
    TwoPhaseState_String tpState;
    Water watTags;
    PhysicalProperties_string_op::CO2 co2Tags;

    // ----------------------------
    // 1. Matrix Domain
    // ----------------------------
    {
        const auto& Sw_vec = fieldMgr_.getMatrixScalar(sConfig.saturation)->data;

        auto& Pc_vec = fieldMgr_.createMatrixScalar(tpState.Pc_field, 0.0)->data;
        auto& krw_vec = fieldMgr_.createMatrixScalar(watTags.k_rw_tag, 0.0)->data;
        auto& krg_vec = fieldMgr_.createMatrixScalar(co2Tags.k_rg_tag, 0.0)->data;

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

    // ----------------------------
    // 2. Fracture Domain
    // ----------------------------
    {
        const auto& Sw_vec = fieldMgr_.getFractureScalar(sConfig.saturation)->data;

        auto& Pc_vec = fieldMgr_.createFractureScalar(tpState.Pc_field, 0.0)->data;
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
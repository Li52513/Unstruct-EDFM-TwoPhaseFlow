#include "2D_VariableInitializer.h"

using namespace PhysicalProperties_string_op;

VariableInitializer_2D::VariableInitializer_2D(const MeshManager& meshMgr, FieldManager_2D& fieldMgr)
    : meshMgr_(meshMgr), fieldMgr_(fieldMgr) {
}

void VariableInitializer_2D::InitSinglePhaseState(const PressureEquation_String& pConfig,
    const TemperatureEquation_String& tConfig,
    const LinearInitParams& pInit,
    const LinearInitParams& tInit)
{
    std::cout << "[VarInit_2D] Initializing Single-Phase State..." << std::endl;
    CreateAndInitMatrixScalar(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit);
    CreateAndInitFractureScalar(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit);

    CreateAndInitMatrixScalar(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit);
    CreateAndInitFractureScalar(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit);
    std::cout << "[VarInit_2D] Single-Phase Initialization Completed." << std::endl;
}

void VariableInitializer_2D::InitIMPESState(const PressureEquation_String& pConfig,
    const TemperatureEquation_String& tConfig,
    const SaturationEquation_String& sConfig,
    const LinearInitParams& pInit,
    const LinearInitParams& tInit,
    const LinearInitParams& sInit,
    const CapRelPerm::VGParams& vgParams,
    const CapRelPerm::RelPermParams& rpParams)
{
    std::cout << "[VarInit_2D] Initializing IMPES Two-Phase State..." << std::endl;

    CreateAndInitMatrixScalar(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit);
    CreateAndInitFractureScalar(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit);

    CreateAndInitMatrixScalar(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit);
    CreateAndInitFractureScalar(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit);

    CreateAndInitMatrixScalar(sConfig.saturation, sConfig.saturation_old, sConfig.saturation_prev, sInit);
    CreateAndInitFractureScalar(sConfig.saturation, sConfig.saturation_old, sConfig.saturation_prev, sInit);

    CalculateDerivedTwoPhaseFields(sConfig, vgParams, rpParams);

    // IMPES: 初始化衍生量 p_g
    TwoPhaseState_String tpState;
    CreateAndInitMatrixScalar(tpState.p_g_field, tpState.p_g_old_field, tpState.p_g_prev_field, pInit);
    CreateAndInitFractureScalar(tpState.p_g_field, tpState.p_g_old_field, tpState.p_g_prev_field, pInit);

    std::cout << "[VarInit_2D] IMPES Initialization Completed." << std::endl;
}

void VariableInitializer_2D::CreateAndInitMatrixScalar(const std::string& tagName, const std::string& oldTagName, const std::string& prevTagName, const LinearInitParams& params)
{
    auto& field = fieldMgr_.createMatrixScalar(tagName, params.refVal)->data;
    auto& fieldOld = fieldMgr_.createMatrixScalar(oldTagName, params.refVal)->data;
    auto& fieldPrev = fieldMgr_.createMatrixScalar(prevTagName, params.refVal)->data;

    const auto& cells = meshMgr_.mesh().getCells();
    int nCells = static_cast<int>(cells.size());

#pragma omp parallel for
    for (int i = 0; i < nCells; ++i) {
        double val = ComputeLinearValue(cells[i].center, params);
        field[i] = val; fieldOld[i] = val; fieldPrev[i] = val;
    }
}

void VariableInitializer_2D::CreateAndInitFractureScalar(const std::string& tagName, const std::string& oldTagName, const std::string& prevTagName, const LinearInitParams& params)
{
    auto& field = fieldMgr_.createFractureScalar(tagName, params.refVal)->data;
    auto& fieldOld = fieldMgr_.createFractureScalar(oldTagName, params.refVal)->data;
    auto& fieldPrev = fieldMgr_.createFractureScalar(prevTagName, params.refVal)->data;

    const auto& fractures = meshMgr_.fracture_network().fractures;
    size_t totalFracCells = field.size();
    const auto& cells = meshMgr_.mesh().getCells();

    // 局部偏移量预计算
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
                int hostCellID = frac.elements[localIdx].cellID;
                double val = ComputeLinearValue(cells[hostCellID].center, params);
                field[arrayIdx] = val; fieldOld[arrayIdx] = val; fieldPrev[arrayIdx] = val;
            }
        }
    }
}

void VariableInitializer_2D::CalculateDerivedTwoPhaseFields(const SaturationEquation_String& sConfig, const CapRelPerm::VGParams& vg, const CapRelPerm::RelPermParams& rp)
{
    TwoPhaseState_String tpState;
    Water watTags;
    PhysicalProperties_string_op::CO2 co2Tags;

    // 1. Matrix
    {
        const auto& Sw_vec = fieldMgr_.getMatrixScalar(sConfig.saturation)->data;
        auto& Pc_vec = fieldMgr_.createMatrixScalar(tpState.Pc_field, 0.0)->data;
        auto& krw_vec = fieldMgr_.createMatrixScalar(watTags.k_rw_tag, 0.0)->data;
        auto& krg_vec = fieldMgr_.createMatrixScalar(co2Tags.k_rg_tag, 0.0)->data;

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(Sw_vec.size()); ++i) {
            double sw = Sw_vec[i];
            Pc_vec[i] = pc_vG(sw, vg);
            kr_Mualem_vG(sw, vg, rp, krw_vec[i], krg_vec[i]);
        }
    }

    // 2. Fracture
    {
        const auto& Sw_vec = fieldMgr_.getFractureScalar(sConfig.saturation)->data;
        auto& Pc_vec = fieldMgr_.createFractureScalar(tpState.Pc_field, 0.0)->data;
        auto& krw_vec = fieldMgr_.createFractureScalar(watTags.k_rw_tag, 0.0)->data;
        auto& krg_vec = fieldMgr_.createFractureScalar(co2Tags.k_rg_tag, 0.0)->data;

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(Sw_vec.size()); ++i) {
            double sw = Sw_vec[i];
            Pc_vec[i] = pc_vG(sw, vg);
            kr_Mualem_vG(sw, vg, rp, krw_vec[i], krg_vec[i]);
        }
    }
}
#pragma once

#include <vector>
#include <string>
#include <memory>
#include <omp.h>
#include <iostream>

#include "MeshManager.h" 
#include "2D_FieldManager.h"
#include "SolverContrlStrName_op.h"
#include "CapRelPerm_HD.h" 
#include "InitParams.h" // 引入公共独立参数

/**
 * @class VariableInitializer_2D
 * @brief 2D-EDFM 主变量场构建与初始化器 (纯净版，与 solverIndex 彻底解耦)
 */
class VariableInitializer_2D
{
public:
    VariableInitializer_2D(const MeshManager& meshMgr, FieldManager_2D& fieldMgr);

    void InitSinglePhaseState(const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::TemperatureEquation_String& tConfig,
        const LinearInitParams& pInit,
        const LinearInitParams& tInit);

    void InitIMPESState(const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::TemperatureEquation_String& tConfig,
        const PhysicalProperties_string_op::SaturationEquation_String& sConfig,
        const LinearInitParams& pInit,
        const LinearInitParams& tInit,
        const LinearInitParams& sInit,
        const CapRelPerm::VGParams& vgParams,
        const CapRelPerm::RelPermParams& rpParams);

    template<int N>
    void InitFIMState(const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::TemperatureEquation_String& tConfig,
        const PhysicalProperties_string_op::SaturationEquation_String& sConfig,
        const LinearInitParams& pInit,
        const LinearInitParams& tInit,
        const LinearInitParams& sInit,
        int pIdx = 0, int sIdx = 1, int tIdx = 2)
    {
        std::cout << "[VarInit_2D] Initializing FIM State (N=" << N << ") with AD Seeding..." << std::endl;

        CreateAndInitMatrixADScalar<N>(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit, pIdx);
        CreateAndInitFractureADScalar<N>(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit, pIdx);

        CreateAndInitMatrixADScalar<N>(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit, tIdx);
        CreateAndInitFractureADScalar<N>(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit, tIdx);

        if (sIdx >= 0) {
            CreateAndInitMatrixADScalar<N>(sConfig.saturation, sConfig.saturation_old, sConfig.saturation_prev, sInit, sIdx);
            CreateAndInitFractureADScalar<N>(sConfig.saturation, sConfig.saturation_old, sConfig.saturation_prev, sInit, sIdx);
        }
        std::cout << "[VarInit_2D] FIM State Seeding Completed." << std::endl;
    }

private:
    void CreateAndInitMatrixScalar(const std::string& tagName, const std::string& oldTagName, const std::string& prevTagName, const LinearInitParams& params);
    void CreateAndInitFractureScalar(const std::string& tagName, const std::string& oldTagName, const std::string& prevTagName, const LinearInitParams& params);
    void CalculateDerivedTwoPhaseFields(const PhysicalProperties_string_op::SaturationEquation_String& sConfig, const CapRelPerm::VGParams& vg, const CapRelPerm::RelPermParams& rp);

    template<int N>
    void CreateAndInitMatrixADScalar(const std::string& tagName, const std::string& oldTagName, const std::string& prevTagName, const LinearInitParams& params, int seedIdx)
    {
        auto& field = fieldMgr_.createMatrixADScalar<N>(tagName)->data;
        auto& fieldOld = fieldMgr_.createMatrixADScalar<N>(oldTagName)->data;
        auto& fieldPrev = fieldMgr_.createMatrixADScalar<N>(prevTagName)->data;

        const auto& cells = meshMgr_.mesh().getCells();
        int nCells = static_cast<int>(cells.size());

#pragma omp parallel for
        for (int i = 0; i < nCells; ++i) {
            double val = ComputeLinearValue(cells[i].center, params);
            ADVar<N> adVal(val, seedIdx);
            field[i] = adVal; fieldOld[i] = adVal; fieldPrev[i] = adVal;
        }
    }

    template<int N>
    void CreateAndInitFractureADScalar(const std::string& tagName, const std::string& oldTagName, const std::string& prevTagName, const LinearInitParams& params, int seedIdx)
    {
        auto& field = fieldMgr_.createFractureADScalar<N>(tagName)->data;
        auto& fieldOld = fieldMgr_.createFractureADScalar<N>(oldTagName)->data;
        auto& fieldPrev = fieldMgr_.createFractureADScalar<N>(prevTagName)->data;

        const auto& fractures = meshMgr_.fracture_network().fractures;
        size_t totalFracCells = field.size();
        const auto& cells = meshMgr_.mesh().getCells();

        // [局部偏移量计算] 无假设、绝对安全的 O(N) 预处理，完美支持 OpenMP
        std::vector<size_t> localOffsets(fractures.size() + 1, 0);
        for (size_t i = 0; i < fractures.size(); ++i) {
            localOffsets[i + 1] = localOffsets[i] + fractures[i].elements.size();
        }

        for (size_t i = 0; i < fractures.size(); ++i) {
            size_t baseOffset = localOffsets[i];
            const auto& frac = fractures[i];

#pragma omp parallel for
            for (int localIdx = 0; localIdx < static_cast<int>(frac.elements.size()); ++localIdx) {
                // 纯正的局部物理场数组下标
                size_t arrayIdx = baseOffset + localIdx;

                if (arrayIdx < totalFracCells) {
                    // 2D 裂缝线段代理质心: 利用其所在的基岩单元中心 (高度精确且无假设)
                    int hostCellID = frac.elements[localIdx].cellID;
                    double val = ComputeLinearValue(cells[hostCellID].center, params);

                    ADVar<N> adVal(val, seedIdx);
                    field[arrayIdx] = adVal;
                    fieldOld[arrayIdx] = adVal;
                    fieldPrev[arrayIdx] = adVal;
                }
            }
        }
    }

    const MeshManager& meshMgr_;
    FieldManager_2D& fieldMgr_;
};
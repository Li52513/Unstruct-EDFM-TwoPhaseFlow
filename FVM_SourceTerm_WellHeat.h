#pragma once
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FVM_WellDOF.h"

namespace FVM {
    namespace SourceTerm {

        // --- new: 注入井 Dirichlet 约束配置 -----------------------------------------
        struct InjectorDirichletPinOptions {
            bool   enable = false;
            double relativeWeight = 0.0;  ///< 乘在 cp*q_in*Tin 上的倍数
            double minWeight = 0.0;  ///< 绝对下限 (W/K)
        };
        // -----------------------------------------------------------------------------

        bool add_temperature_source_from_qm_single(
            MeshManager& mgr,
            FieldRegistry& reg,
            WellDOF::Role role,
            const char* PI_name,
            const char* mask_name,
            const char* pw_field,
            double Tin,
            const char* p_cell,
            const char* cp_field,
            double thickness,
            bool accumulate,
            bool verbose,
            const InjectorDirichletPinOptions& pin = InjectorDirichletPinOptions{});   // --- new ---
    } // namespace SourceTerm
} // namespace FVM
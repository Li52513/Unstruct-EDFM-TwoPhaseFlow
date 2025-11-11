#pragma once
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FVM_WellDOF.h"

namespace FVM {
    namespace SourceTerm {


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
            bool verbose);   // --- new ---
    } // namespace SourceTerm
} // namespace FVM
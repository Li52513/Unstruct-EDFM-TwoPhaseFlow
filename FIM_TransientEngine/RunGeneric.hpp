#pragma once

#include "Types.hpp"
#include "Diagnostics.hpp"
#include "StateSync.hpp"
#include "StepKernels.hpp"

namespace FIM_Engine {

    template <int N, typename MeshMgrType, typename FieldMgrType>
    inline void RunGenericFIMTransient(
        const std::string& caseName,
        MeshMgrType& mgr,
        FieldMgrType& fm,
        const InitialConditions& ic,
        const std::vector<WellScheduleStep>& wells,
        const TransientSolverParams& params,
        SolverRoute route = SolverRoute::FIM,
        const TransientOptionalModules<MeshMgrType, FieldMgrType>& modules =
        TransientOptionalModules<MeshMgrType, FieldMgrType>());

} // namespace FIM_Engine

#include "RunGeneric_impl.hpp"

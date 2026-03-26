#pragma once

#include <string>
#include <vector>

#include "BoundaryConditionManager.h"

enum class VTKBCTransportKind {
    Pressure,
    Temperature,
    Saturation
};

enum class VTKBCPrimaryFluidModel {
    Water,
    CO2
};

enum class VTKBCWaterFamilyDerivePolicy {
    LegacyWaterEvaluator,
    FollowPrimaryFluid
};

struct VTKBCVariableBinding {
    std::string field_name;
    const BoundarySetting::BoundaryConditionManager* bc = nullptr;
    VTKBCTransportKind transport_kind = VTKBCTransportKind::Pressure;
};

struct VTKBoundaryVisualizationContext {
    std::vector<VTKBCVariableBinding> bindings;

    // Keep legacy compatibility by default.
    // In single-phase runs, *_w_* tags can represent the configured primary fluid.
    VTKBCWaterFamilyDerivePolicy water_family_policy = VTKBCWaterFamilyDerivePolicy::LegacyWaterEvaluator;
    VTKBCPrimaryFluidModel primary_fluid_model = VTKBCPrimaryFluidModel::Water;

    // Export controls / numerical guards
    bool export_bc_mask = true;
    double geom_floor = 1.0e-12;
    double coeff_floor = 1.0e-20;
    double denom_eps = 1.0e-14;

    const VTKBCVariableBinding* findBinding(const std::string& field_name) const {
        for (const auto& b : bindings) {
            if (b.field_name == field_name && b.bc) return &b;
        }
        return nullptr;
    }

    const VTKBCVariableBinding* findFirstByKind(VTKBCTransportKind kind) const {
        for (const auto& b : bindings) {
            if (b.transport_kind == kind && b.bc) return &b;
        }
        return nullptr;
    }

    bool usePrimaryForWaterFamily() const {
        return water_family_policy == VTKBCWaterFamilyDerivePolicy::FollowPrimaryFluid;
    }
};

/**
 * @file Well_WellControlTypes.h
 * @brief Well control types and schedule step definitions.
 */
#ifndef WELL_WELLCONTROLTYPES_H
#define WELL_WELLCONTROLTYPES_H

#include <string>
#include <vector>

enum class WellControlMode {
    BHP,    ///< Bottom hole pressure control
    Rate    ///< Rate control
};

enum class WellTargetDomain {
    Matrix,     ///< Matrix domain
    Fracture    ///< Fracture domain
};

enum class WellComponentMode {
    Water,  ///< Water only
    Gas,    ///< Gas only (CO2)
    Total   ///< Total flux
};

/**
 * @brief Rate target type.
 */
enum class WellRateTargetType {
    MassRate,       ///< Target in mass rate [kg/s]
    StdVolumeRate   ///< Target in standard volume rate [m^3/s]
};

enum class WellAxis {
    X,
    Y,
    Z,
    None    ///< For 2D cases or explicit user override
};

/**
 * @brief Completion id semantic space.
 */
enum class CompletionIdSpace {
    SolverIndex,        ///< completion_id is global solver block index
    FractureLocalIndex, ///< completion_id is fracture-local index
    AutoLegacy          ///< Backward-compatible legacy inference
};

/**
 * @brief One completion specification.
 */
struct WellCompletionSpec {
    int completion_id = -1;
    CompletionIdSpace completion_id_space = CompletionIdSpace::AutoLegacy;
    int completion_solver_index = -1;   ///< Normalized global solver block index
};

struct WellScheduleStep {
    double t_start = 0.0;
    double t_end = 0.0;
    std::string well_name = "";
    WellTargetDomain domain = WellTargetDomain::Matrix;
    WellControlMode control_mode = WellControlMode::BHP;
    double target_value = 0.0;
    WellComponentMode component_mode = WellComponentMode::Total;
    double rw = 0.1;
    double skin = 0.0;
    WellAxis well_axis = WellAxis::None;

    // Legacy single-completion fields:
    // Matrix: cell local index | Fracture: solverIndex or fracture-local index.
    int completion_id = -1;
    CompletionIdSpace completion_id_space = CompletionIdSpace::AutoLegacy;
    int completion_solver_index = -1;

    // Multi-completion model (single completion is completions.size()==1).
    std::vector<WellCompletionSpec> completions;

    double wi_override = -1.0;
    double L_override = -1.0;

    // For WellComponentMode::Total:
    // - If frac_w + frac_g > 0, assembler normalizes and uses them.
    // - If both are 0, assembler falls back to mobility-weighted auto split.
    double frac_w = 0.0;
    double frac_g = 0.0;

    // Optional injection stream temperature (K).
    // Effective for injection branch in energy equation when > 0.
    double injection_temperature = -1.0;
    bool injection_is_co2 = false;

    /**
     * @brief Rate target type (default keeps current mass-rate behavior).
     */
    WellRateTargetType rate_target_type = WellRateTargetType::MassRate;

    /**
     * @brief Standard-condition reference pressure [Pa] for StdVolumeRate.
     */
    double std_pressure = 101325.0;

    /**
     * @brief Standard-condition reference temperature [K] for StdVolumeRate.
     */
    double std_temperature = 288.15;
};

#endif // WELL_WELLCONTROLTYPES_H

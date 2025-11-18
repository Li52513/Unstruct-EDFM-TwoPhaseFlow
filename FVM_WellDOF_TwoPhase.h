#pragma once
#include <string>
#include <vector>

/**
 * @struct WellDOF_TwoPhase
 * @brief Extends the single-phase WellDOF to handle two-phase injection.
 *
 * This structure defines a well's degree of freedom and its control settings for
 * a two-phase (Water-CO2) system. It adds parameters to specify the phase behavior
 * of the injected fluid.
 */

struct WellDOF_TwoPhase {

    // --- Enums from the original WellDOF ---
    enum class Mode { Pressure, Rate };
    enum class Role { Injector, Producer };

    // --- Core well properties ---
    std::string name;                         // Well name, used for logging and field naming
    Mode   mode = Mode::Pressure;             // Control mode: constant bottom-hole pressure or constant mass rate
    Role   role = Role::Injector;             // Well role: Injector or Producer

    // --- Control target value ---
    double target = 0.0;                      // Target value: BHP in [Pa] or total mass flow rate in [kg/s]

    // === Two-Phase Specific Injection Parameters (for Injectors only) ===
    double Tin = 0.0;                           // Injection temperature T_inj [K]
    double s_w_bh = 1.0;                        // Water saturation of the injected fluid at the wellbore, S_w_bh [-]
                                                // s_w_bh = 1.0 for pure water injection
                                                // s_w_bh = 0.0 for pure CO2 injection
                                                // 0 < s_w_bh < 1 for mixed-phase injection

     // --- Fluid properties of the injected fluid ---
    double mu_w_inj = 3e-4;                   // Viscosity of injected water at wellbore conditions [Pa，s]
    double mu_g_inj = 1.5e-5;                 // Viscosity of injected CO2 at wellbore conditions [Pa，s]
    double rho_w_inj = 980.0;                 // Density of injected water [kg/m^3]
    double rho_g_inj = 700.0;                 // Density of injected CO2 [kg/m^3]
    double cp_w_inj = 4200.0;                 // Heat capacity of injected water [J/(kg，K)]
    double cp_g_inj = 1200.0;                 // Heat capacity of injected CO2 [J/(kg，K)]

    // --- Internal fields for solver linkage ---
    std::string mask_field;                   // Field name for the well's perforation mask (e.g., "mask_INJ1")
    std::string PI_field_w;                   // Field name for the WATER phase productivity index (e.g., "PI_w_INJ1")
    std::string PI_field_g;                   // Field name for the GAS (CO2) phase productivity index (e.g., "PI_g_INJ1")
    int    lid = -1;                          // Index in the global unknown vector if solved simultaneously
};

/**
 * @brief Registers well unknowns for a system with two-phase wells.
 * @param Nc Number of cell unknowns.
 * @param wells A vector of two-phase well degrees of freedom.
 * @return The total number of unknowns (cells + wells).
 */

inline int register_well_unknowns_TwoPhase(int Nc, std::vector<WellDOF_TwoPhase>& wells)
{
    int lid = Nc;
    for (auto& w : wells)
    {
        w.lid = lid++;
    }
    return lid;
}
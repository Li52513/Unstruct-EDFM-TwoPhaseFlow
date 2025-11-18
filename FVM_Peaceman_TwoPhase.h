#pragma once
#include <string>
#include <vector>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "UserDefineVarType.h"

/**
 * @namespace PeacemanTwoPhase
 * @brief Handles pre-computation for the Peaceman well model in a two-phase context.
 *        This module is fully independent of the single-phase well model.
 */
namespace PeacemanTwoPhase {

    /**
     * @struct WellSpec
     * @brief Defines the geometric and physical specifications of a well.
     *
     * This is a self-contained definition for the two-phase module, ensuring
     * no dependency on the single-phase implementation.
     */
    struct WellSpec
    {
        std::string name = "WELL_2P"; // Well name for logging
        Vector pos;                   // Well's projected 2D coordinates
        double rw = 0.1;              // Wellbore radius [m]
        double skin = 0.0;            // Skin factor [-]
        double H = 1.0;               // Effective formation thickness [m]
        double perfRadius = 0.0;      // Perforation radius for selecting cells (>0), or 0 to use nearest N
        int    maxHitCells = 1;       // Number of nearest cells to select if perfRadius is 0
    };

    /**
     * @struct PeacemanParams_TwoPhase
     * @brief Parameters for calculating the geometric Well Index (WI).
     */
    struct PeacemanParams_TwoPhase
    {
        double reFactor = 0.28;       // Equivalent radius factor: r_e = reFactor * sqrt(A_cell/pi)
        double fallbackKh = 1e-14;    // Default horizontal permeability if fields not found [m^2]
    };


    /**
     * @brief Builds the perforation mask and the geometric Well Index (WI) field for a well.
     *
     * @param mgr       MeshManager, providing access to the mesh.
     * @param reg       FieldRegistry, where the mask and WI fields will be stored.
     * @param spec      The well's geometric specification (using the local PeacemanTwoPhase::WellSpec).
     * @param prm       Parameters for the WI calculation.
     * @param mask_name The name for the output perforation mask field.
     * @param WI_name   The name for the output geometric Well Index field.
     */
    void build_well_mask_and_WI(
        MeshManager& mgr, FieldRegistry& reg,
        const WellSpec& spec, // Now uses the WellSpec defined within this file's namespace
        const PeacemanParams_TwoPhase& prm,
        const std::string& mask_name, const std::string& WI_name
    );

} // namespace PeacemanTwoPhase
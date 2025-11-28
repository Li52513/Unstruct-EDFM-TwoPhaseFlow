#pragma once

#include <string>

/**
 * @brief Shared default field names produced by the pressure assembly phase.
 *
 * This lightweight header allows other modules (e.g. flux splitting) to query
 * the canonical names without pulling the entire pressure assembly interface
 * and risking circular includes.
 */
namespace PressureAssemblyFieldNames
{
    inline const char* total_mass_flux()            { return "mf_total"; }
    inline const char* total_vol_flux()             { return "Qf_total"; }
    inline const char* total_velocity()             { return "ufn_total"; }
    inline const char* capillary_correction_flux()  { return "mf_capillary_corr"; }
    inline const char* gravity_correction_flux()    { return "mf_gravity_corr"; }
}


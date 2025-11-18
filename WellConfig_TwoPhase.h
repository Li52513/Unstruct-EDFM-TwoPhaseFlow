#pragma once
#include <string>
#include <vector>
#include "FVM_WellDOF_TwoPhase.h"
#include "FVM_Peaceman_TwoPhase.h" // Now includes the self-contained two-phase header

/**
 * @struct WellConfig_TwoPhase
 * @brief Unified configuration for a two-phase well, fully decoupled from single-phase modules.
 */
struct WellConfig_TwoPhase {
    // --- Basic Identification ---
    std::string name;
    WellDOF_TwoPhase::Role role;

    // --- Geometry and Peaceman Model ---
    // This now uses the WellSpec defined inside the PeacemanTwoPhase namespace
    PeacemanTwoPhase::WellSpec geom;
    PeacemanTwoPhase::PeacemanParams_TwoPhase pm_2p;

    // ... (rest of the struct is identical to the previous version) ...
    WellDOF_TwoPhase::Mode mode;
    double target;
    double Tin;
    double s_w_bh;
    double mu_w_inj;
    double mu_g_inj;
    double rho_w_inj;  
    double rho_g_inj;  
    double cp_w_inj;   
    double cp_g_inj; 

    std::string mask_name;
    std::string WI_name;
    std::string pw_name;
    int lid;

    // Default constructor
    WellConfig_TwoPhase() :
        role(WellDOF_TwoPhase::Role::Injector),
        mode(WellDOF_TwoPhase::Mode::Pressure),
        target(0.0), Tin(0.0), s_w_bh(1.0), mu_w_inj(3e-4), mu_g_inj(1.5e-5), rho_w_inj(980.0), rho_g_inj(700.0), cp_w_inj(4200.0), cp_g_inj(1200.0), lid(-1) {
    }

    void derive_names_if_empty() {
        if (mask_name.empty()) mask_name = "mask_" + name;
        if (WI_name.empty())   WI_name = "WI_" + name;
        if (pw_name.empty())   pw_name = "p_w_" + name;
    }

    WellDOF_TwoPhase toWellDOF() const {
        WellDOF_TwoPhase w;
        w.name = pw_name;
        w.role = role;
        w.mode = mode;
        w.target = target;
        w.Tin = Tin;
        w.s_w_bh = s_w_bh;
        w.mu_w_inj = mu_w_inj;
        w.mu_g_inj = mu_g_inj;
        w.rho_w_inj = rho_w_inj;
        w.rho_g_inj = rho_g_inj;
        w.cp_w_inj = cp_w_inj;
        w.cp_g_inj = cp_g_inj;
        w.mask_field = mask_name;
        w.PI_field_w = WI_name; // Still using this to pass the WI field name
        w.PI_field_g = WI_name;
        w.lid = -1;
        return w;
    }
};

// ... (The helper functions build_masks_and_WI_for_all and 
//      register_well_dofs_for_all_TwoPhase remain IDENTICAL,
//      as they already operate on WellConfig_TwoPhase) ...

// For completeness, here is the build_masks_and_WI_for_all function again.
// Note that no changes are needed inside it.
inline void build_masks_and_WI_for_all(
    MeshManager& mgr, FieldRegistry& reg,
    std::vector<WellConfig_TwoPhase>& wells)
{
    std::cout << "--- Building Masks and Well Indices for all Two-Phase Wells ---\n";
    for (auto& well : wells) {
        well.derive_names_if_empty();
        PeacemanTwoPhase::build_well_mask_and_WI(
            mgr, reg, well.geom, well.pm_2p,
            well.mask_name, well.WI_name
        );
    }
}

/**
 * @brief Converts all well configurations to WellDOF_TwoPhase objects and registers them as unknowns.
 *
 * This function populates the `wells_dof` vector and registers the additional degrees of
 * freedom for the wells with the solver system. It also back-fills the `lid` (local index)
 * into the original `wells_cfg` vector for easy lookup.
 *
 * @param Nc          Number of cell-based unknowns.
 * @param wells_cfg   The input vector of well configurations.
 * @param wells_dof   The output vector of WellDOF_TwoPhase objects for the solver.
 * @return The total number of unknowns (Nc + number of wells).
 */
inline int register_well_dofs_for_all_TwoPhase(
    int Nc,
    std::vector<WellConfig_TwoPhase>& wells_cfg,
    std::vector<WellDOF_TwoPhase>& wells_dof)
{
    wells_dof.clear();
    wells_dof.reserve(wells_cfg.size());
    for (const auto& cfg : wells_cfg) {
        wells_dof.push_back(cfg.toWellDOF());
    }

    const int N_total = register_well_unknowns_TwoPhase(Nc, wells_dof);

    // Back-fill the lid into the config structure for convenience
    for (size_t i = 0; i < wells_cfg.size(); ++i) {
        wells_cfg[i].lid = wells_dof[i].lid;
    }

    return N_total;
}
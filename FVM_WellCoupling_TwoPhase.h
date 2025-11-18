#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include "WellConfig_TwoPhase.h" // Includes DOF_TwoPhase, Peaceman_TwoPhase, etc.
#include "Solver_AssemblerCOO.h"
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FieldAcessForDiscre.h" // For cellScalar
#include "MultiPhaseProperties.h"  // For kr_Mualem_vG, etc.

namespace FVM {
	namespace TwoPhaseWellCoupling {

        // Helper to get cell ID by its index in the mesh's cell vector
        inline int cell_id_by_index(const Mesh& mesh, int cidx) {
            return mesh.getCells()[cidx].id;
        }

        /**
         * @brief Couples two-phase wells to the IMPLICIT PRESSURE equation.
         *
         * This function modifies the pressure equation's sparse system (A, b) for all
         * cells connected to a well. It uses total mobility from the previous time step (n)
         * to calculate the total productivity index.
         *
         * For a well cell 'i' and a well with bottom-hole pressure 'P_bh':
         *   - Contribution to the equation is: PI_total_i * (P_bh - P_i)
         *   - This modifies A(i,i) and the right-hand side b(i).
         *
         * @param sys       The sparse system (A, b) for the pressure equation.
         * @param mgr       The MeshManager.
         * @param reg       The FieldRegistry, containing all property fields.
         * @param well      The two-phase well configuration.
         * @param vg_params Van Genuchten parameters to calculate rel-perm if needed.
         * @param rp_params Mualem parameters for rel-perm.
         */


        inline void couple_well_to_pressure_equation(
            SparseSystemCOO& sys,
            MeshManager& mgr,
            const FieldRegistry& reg,
            const WellDOF_TwoPhase& well,
            const std::vector<int>& lid_of_cell
        )
        {
            auto& mesh = mgr.mesh();
            const int Nc = (int)mesh.getCells().size();
            const int well_lid = well.lid;

            if (well_lid < 0) {
                std::cerr << "Warning [2P-PressureCoupling]: Invalid well LID for " << well.name << ". Skipping." << std::endl;
                return;
            }

            // --- Part 1: Handle the well's own constraint row ---
            if (well.mode == WellDOF_TwoPhase::Mode::Pressure) {
                // For a pressure-controlled well, the constraint is simple: 1 * P_bh = P_target
                sys.addA(well_lid, well_lid, 1.0);
                sys.addb(well_lid, well.target);
                // The P_bh is now a known value for the cell equations.
            }

            // --- Part 2: Modify perforated cell rows and assemble the rate constraint if needed ---
            auto mask_field = reg.get<volScalarField>(well.mask_field);
            auto WI_field = reg.get<volScalarField>(well.PI_field_w); // WI_name is passed here
            auto lambda_t_field = reg.get<volScalarField>("lambda_t"); // Total mobility at time 'n'

            if (!mask_field || !WI_field || !lambda_t_field) {
                std::cerr << "ERROR [2P-PressureCoupling]: Missing required fields for well " << well.name << std::endl;
                return;
            }

            double rate_diag_sum = 0.0; // Only used for rate-controlled wells

            for (int cidx = 0; cidx < Nc; ++cidx) {
                if ((*mask_field)[cidx] <= 0.0) continue;

                const double WI_i = (*WI_field)[cidx];
                const double lambda_t_i = (*lambda_t_field)[cidx];

                if (WI_i <= 0.0 || lambda_t_i <= 0.0) continue;

                const double PI_total_i = WI_i * lambda_t_i;
                const int cell_lid = lid_of_cell[cidx];
                if (cell_lid < 0) continue;

                // Modify the cell's row. The equation includes a source term Q_i = PI_total_i * (P_bh - P_i)
                // This term moves to the LHS: -PI_total_i * P_bh + PI_total_i * P_i

                // Add PI_total_i to the diagonal A[cell_lid, cell_lid]
                sys.addA(cell_lid, cell_lid, PI_total_i);

                if (well.mode == WellDOF_TwoPhase::Mode::Pressure) {
                    // P_bh is known (well.target), so -PI_total_i * P_bh is a source term
                    // that moves to the RHS.
                    sys.addb(cell_lid, PI_total_i * well.target);
                }
                else { // Rate mode
                    // P_bh is an unknown, so -PI_total_i * P_bh is a matrix coupling term
                    sys.addA(cell_lid, well_lid, -PI_total_i);

                    // For the well's rate equation ¦²(PI_i * (P_bh - P_i)) = Q_target,
                    // we assemble the terms related to this cell.
                    rate_diag_sum += PI_total_i;
                    sys.addA(well_lid, cell_lid, -PI_total_i);
                }
            }

            // --- Part 3: Finalize the well's row for rate-controlled mode ---
            if (well.mode == WellDOF_TwoPhase::Mode::Rate) {
                if (rate_diag_sum != 0.0) {
                    // Add the diagonal term ¦²(PI_i) to A[well_lid, well_lid]
                    sys.addA(well_lid, well_lid, rate_diag_sum);
                }
                // Add the target rate Q_target to the RHS b[well_lid]
                sys.addb(well_lid, well.target);
            }
        }

        /**
         * @brief Adds the EXPLICIT source/sink terms from a well to the SATURATION equation's RHS.
         *
         * This function calculates the water phase mass flow rate (Q_w) for each perforated cell
         * and adds it to the right-hand side vector `b_sat`.
         *
         * - For PRODUCERS, Q_w is determined by the reservoir's fractional flow.
         * - For INJECTORS, Q_w is determined by the well's specified injection saturation `s_w_bh`.
         *
         * @param b_sat     The right-hand side vector for the explicit saturation update.
         * @param mgr       The MeshManager.
         * @param reg       The FieldRegistry.
         * @param well      The two-phase well configuration.
         * @param p_bh      The well's bottom-hole pressure (known after pressure solve).
         * @param p_w_new   The newly solved water pressure field at time n+1.
         * @param vg_params Van Genuchten parameters for rel-perm calculations.
         * @param rp_params Mualem parameters for rel-perm.
         */
        inline void add_well_source_to_saturation_rhs
        (
            std::vector<double>& b_sat, // This is the explicit source vector, not the full system's b
            MeshManager& mgr,
            const FieldRegistry& reg,
            const WellDOF_TwoPhase& well,
            double p_bh,
            const volScalarField& p_w_new,
            const VGParams& vg_params,
            const RelPermParams& rp_params
        ) 
        {
            auto& mesh = mgr.mesh();
            const int Nc = (int)mesh.getCells().size();

            auto mask_field = reg.get<volScalarField>(well.mask_field);
            auto WI_field = reg.get<volScalarField>(well.PI_field_w); // WI_name

            if (!mask_field || !WI_field) return;

            const size_t limit = std::min(
                std::min(static_cast<size_t>(Nc), b_sat.size()),
                std::min(mask_field->data.size(), WI_field->data.size()));

            // --- Get all required fields from the PREVIOUS time step 'n' for producers ---
            auto rho_w_old = reg.get<volScalarField>("rho_w"); 
            auto lambda_w_old = reg.get<volScalarField>("lambda_w");
            auto lambda_g_old = reg.get<volScalarField>("lambda_g");
            
            // Safety check for required fields
            if (well.role == WellDOF_TwoPhase::Role::Producer && (!rho_w_old || !lambda_w_old || !lambda_g_old)) {
                std::cerr << "ERROR [SaturationCoupling]: Missing rho_w, lambda_w, or lambda_g field for producer well " << well.name << std::endl;
                return;
            }

            for (size_t cidx = 0; cidx < limit; ++cidx)
            {
                if ((*mask_field)[cidx] <= 0.0) continue;

                const double WI_i = (*WI_field)[cidx];
                const double p_w_cell_new = p_w_new[cidx];
                double Q_w_i_mass = 0.0;

                if (well.role == WellDOF_TwoPhase::Role::Injector)
                {
                    // --- INJECTOR LOGIC (Correct from previous version) ---
                    double s_w_up = well.s_w_bh;
                    double krw_up, krg_up;
                    kr_Mualem_vG(s_w_up, vg_params, rp_params, krw_up, krg_up);

                    double lambda_w_inj = krw_up / std::max(well.mu_w_inj, 1e-12);
                    double lambda_g_inj = krg_up / std::max(well.mu_g_inj, 1e-12);
                    double total_mobility_inj = lambda_w_inj + lambda_g_inj;
                    if (total_mobility_inj > 0)
                    {
                        double Q_total_vol = WI_i * total_mobility_inj * (p_bh - p_w_cell_new);
                        if (Q_total_vol > 0) 
                        {
                            double f_w_inj = lambda_w_inj / total_mobility_inj;
                            double Q_w_vol = Q_total_vol * f_w_inj;
                            Q_w_i_mass = Q_w_vol * well.rho_w_inj;
                        }
                    }
                }
                else 
                {
                    // 1. Get mobilities from the reservoir at the previous time step 'n'
                    const double lambda_w_res = (*lambda_w_old)[cidx];
                    const double lambda_g_res = (*lambda_g_old)[cidx];
                    const double lambda_t_res = lambda_w_res + lambda_g_res;

                    if (lambda_t_res > 0) {
                        // 2. Calculate TOTAL VOLUME rate using reservoir mobilities
                        double Q_total_vol = WI_i * lambda_t_res * (p_w_cell_new - p_bh);

                        if (Q_total_vol > 0) { // Ensure production
                            // 3. Calculate water fractional flow in the reservoir
                            double f_w_res = lambda_w_res / lambda_t_res;

                            // 4. Calculate water VOLUME rate
                            double Q_w_vol = Q_total_vol * f_w_res;

                            // 5. **CORRECTION**: Convert to MASS rate using DENSITY FROM THE RESERVOIR
                            const double rho_w_res = (*rho_w_old)[cidx];
                            Q_w_i_mass = -Q_w_vol * rho_w_res; // Negative sign for withdrawal
                        }
                    }
                }

                // Add the final mass source term to the vector
                b_sat[cidx] += Q_w_i_mass;

            }
        }


        /**
         * @brief Couples two-phase wells to the IMPLICIT TEMPERATURE equation.
         *
         * - For PRODUCERS, heat is withdrawn at the reservoir temperature, adding a term
         *   proportional to the unknown T_i to the diagonal of matrix A.
         * - For INJECTORS, heat is added at a fixed injection temperature T_inj, which
         *   contributes a constant source term to the right-hand side vector b.
         *
         * @param sys       The sparse system (A, b) for the temperature equation.
         * @param mgr       The MeshManager.
         * @param reg       The FieldRegistry.
         * @param well      The two-phase well configuration.
         * @param Qw_field  A temporary field holding the water mass rates for the well (computed during saturation step).
         * @param Qg_field  A temporary field holding the gas mass rates for the well.
         */
         inline void couple_well_to_temperature_equation(
             SparseSystemCOO& sys,
             MeshManager& mgr,
             const FieldRegistry& reg,
             const WellDOF_TwoPhase& well,
             const volScalarField& Qw_well, // Pass in the per-cell flow rates
             const volScalarField& Qg_well
         ) {
             auto& mesh = mgr.mesh();
             const int Nc = (int)mesh.getCells().size();

             // Get heat capacity fields
             auto cp_w_field = reg.get<volScalarField>("cp_w");
             auto cp_g_field = reg.get<volScalarField>("cp_g");

             for (int cidx = 0; cidx < Nc; ++cidx) {
                 const double Q_w_i = Qw_well[cidx];
                 const double Q_g_i = Qg_well[cidx];

                 if (Q_w_i == 0.0 && Q_g_i == 0.0) continue;

                 const int row = cidx; // Assuming cell index is the matrix row

                 if (well.role == WellDOF_TwoPhase::Role::Injector) {
                     // --- INJECTOR ---
                // Heat source: Q_h = Q_w * h_w_inj + Q_g * h_g_inj
                // Assuming h = cp * T
                     double Q_h = Qw_well[cidx] * well.cp_w_inj * well.Tin +
                         Qg_well[cidx] * well.cp_g_inj * well.Tin;
                     sys.addb(row, Q_h); // Add to RHS
                 }
                 else {
                     // --- PRODUCER ---
                     // Heat sink depends on the unknown reservoir temperature T_i:
                     // Q_h = Q_w * cp_w * T_i + Q_g * cp_g * T_i
                     // This contributes to the diagonal of matrix A.
                     // Note: Q_w and Q_g for producers are negative.
                     const double cp_w_i = (*cp_w_field)[cidx];
                     const double cp_g_i = (*cp_g_field)[cidx];

                     // The term is -(Q_w*cp_w + Q_g*cp_g) * T_i.
                     // In the equation A*T = b, this becomes a diagonal term in A.
                     // Since Q_w, Q_g are negative, the contribution is positive.
                     double diag_contribution = -(Q_w_i * cp_w_i + Q_g_i * cp_g_i);
                     sys.addA(row, row, diag_contribution);
                 }
             }
         }


	}
}
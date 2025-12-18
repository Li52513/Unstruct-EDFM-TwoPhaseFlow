#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "FaceFieldRegistry.h"
#include "FieldRegistry.h"
#include "FluxSplitterandSolver.h"
#include "MeshManager.h"
#include "SolverContrlStrName.h"

namespace IMPES_Iteration
{
    struct PhaseMassBalanceDiagnostics
    {
        // --- totals (kg) ---
        double Mw_old = 0.0;
        double Mw_new = 0.0;
        double Mg_old = 0.0;
        double Mg_new = 0.0;

        // --- totals (kg/s) over the step ---
        double Qw_sum = 0.0;
        double Qg_sum = 0.0;
        double Fw_out_sum = 0.0; // sum(div(Fw)) = net boundary outflow (sign included)
        double Fg_out_sum = 0.0;

        // --- residual norms (kg/s) ---
        double global_res_w_rate = 0.0; // sum_i r_w,i
        double global_res_g_rate = 0.0;
        double max_abs_res_w_rate = 0.0;
        double max_abs_res_g_rate = 0.0;
        double l1_abs_res_w_rate = 0.0;
        double l1_abs_res_g_rate = 0.0;

        int max_abs_res_w_cellId = -1;
        int max_abs_res_g_cellId = -1;

        // --- “physical” mass balance mismatch using Mw_new-Mw_old (kg/s) ---
        double mass_balance_err_w_rate = 0.0;
        double mass_balance_err_g_rate = 0.0;
    };

    inline bool computePhaseMassBalanceDiagnostics(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const std::string& sw_field,
        const std::string& sw_old_field,
        const FluxSplitConfig& fluxCfg,
        const std::string& qw_field,          // optional, may be empty
        const std::string& qg_field,          // optional, may be empty
        double dt,
        const std::string& out_res_w_field,   // created/overwritten in reg, unit: kg/s
        const std::string& out_res_g_field,   // created/overwritten in reg, unit: kg/s
        PhaseMassBalanceDiagnostics& diag)
    {
        diag = PhaseMassBalanceDiagnostics{};

        if (!(dt > 0.0))
        {
            std::cerr << "[MassDiag] invalid dt.\n";
            return false;
        }

        Mesh& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& faces = mesh.getFaces();
        const auto& id2idx = mesh.getCellId2Index();
        const std::size_t nCells = cells.size();

        auto sw = reg.get<volScalarField>(sw_field.c_str());
        auto sw_old = reg.get<volScalarField>(sw_old_field.c_str());
        if (!sw || !sw_old)
        {
            std::cerr << "[MassDiag] missing saturation fields '" << sw_field
                      << "' or '" << sw_old_field << "'.\n";
            return false;
        }

        auto phi = reg.get<volScalarField>(PhysicalProperties_string::Rock().phi_tag);
        if (!phi)
        {
            std::cerr << "[MassDiag] missing porosity field '" << PhysicalProperties_string::Rock().phi_tag << "'.\n";
            return false;
        }

        auto rho_w = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_tag);
        auto rho_w_old = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_old_tag);
        auto rho_g = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_tag);
        auto rho_g_old = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_old_tag);
        if (!rho_w || !rho_w_old || !rho_g || !rho_g_old)
        {
            std::cerr << "[MassDiag] missing density fields (current/old).\n";
            return false;
        }

        auto mf_w = freg.get<faceScalarField>(fluxCfg.water_mass_flux.c_str());
        auto mf_g = freg.get<faceScalarField>(fluxCfg.gas_mass_flux.c_str());
        if (!mf_w || !mf_g)
        {
            std::cerr << "[MassDiag] missing face mass flux fields '"
                      << fluxCfg.water_mass_flux << "' or '" << fluxCfg.gas_mass_flux << "'.\n";
            return false;
        }

        std::shared_ptr<volScalarField> q_w = nullptr;
        if (!qw_field.empty())
        {
            q_w = reg.get<volScalarField>(qw_field.c_str());
            if (!q_w)
            {
                std::cerr << "[MassDiag] requested qw field '" << qw_field << "' not found.\n";
                return false;
            }
        }

        std::shared_ptr<volScalarField> q_g = nullptr;
        if (!qg_field.empty())
        {
            q_g = reg.get<volScalarField>(qg_field.c_str());
            if (!q_g)
            {
                std::cerr << "[MassDiag] requested qg field '" << qg_field << "' not found.\n";
                return false;
            }
        }

        // --- compute per-cell divergence of face mass fluxes (owner +, neighbor -) ---
        std::vector<double> div_w(nCells, 0.0);
        std::vector<double> div_g(nCells, 0.0);

        for (const auto& F : faces)
        {
            const int iF = F.id - 1;
            if (iF < 0) continue;
            if ((std::size_t)iF >= mf_w->data.size() || (std::size_t)iF >= mf_g->data.size()) continue;

            const double Fw = (*mf_w)[(std::size_t)iF];
            const double Fg = (*mf_g)[(std::size_t)iF];

            const int ownerId = F.ownerCell;
            const int neighId = F.neighborCell;

            if (ownerId >= 0)
            {
                const std::size_t iOwner = id2idx.at(ownerId);
                div_w[iOwner] += Fw;
                div_g[iOwner] += Fg;
            }
            if (neighId >= 0)
            {
                const std::size_t iNeigh = id2idx.at(neighId);
                div_w[iNeigh] -= Fw;
                div_g[iNeigh] -= Fg;
            }
        }

        // --- allocate output residual fields (kg/s) ---
        auto res_w = reg.getOrCreate<volScalarField>(out_res_w_field, nCells, 0.0);
        auto res_g = reg.getOrCreate<volScalarField>(out_res_g_field, nCells, 0.0);
        if (!res_w || !res_g)
        {
            std::cerr << "[MassDiag] failed to allocate residual fields.\n";
            return false;
        }
        res_w->data.assign(nCells, 0.0);
        res_g->data.assign(nCells, 0.0);

        // --- main loop over cells ---
        const double tiny = 1e-30;
        for (const auto& c : cells)
        {
            if (c.id < 0) continue; // skip ghost
            const std::size_t i = id2idx.at(c.id);

            const double V = c.volume;
            if (V <= tiny) continue;

            const double ph = (*phi)[i];

            const double sw_n = (*sw_old)[i];
            const double sw_np1 = (*sw)[i];
            const double sg_n = 1.0 - sw_n;
            const double sg_np1 = 1.0 - sw_np1;

            const double rw = (*rho_w)[i];
            const double rg = (*rho_g)[i];

            const double rw_old = (*rho_w_old)[i];
            const double rg_old = (*rho_g_old)[i];

            const double Qw = q_w ? (*q_w)[i] : 0.0;
            const double Qg = q_g ? (*q_g)[i] : 0.0;

            // mass totals (kg)
            const double Mw_n = ph * V * rw_old * sw_n;
            const double Mw_np1 = ph * V * rw * sw_np1;
            const double Mg_n = ph * V * rg_old * sg_n;
            const double Mg_np1 = ph * V * rg * sg_np1;

            diag.Mw_old += Mw_n;
            diag.Mw_new += Mw_np1;
            diag.Mg_old += Mg_n;
            diag.Mg_new += Mg_np1;

            // totals (kg/s)
            diag.Qw_sum += Qw;
            diag.Qg_sum += Qg;
            diag.Fw_out_sum += div_w[i];
            diag.Fg_out_sum += div_g[i];

            // discrete residuals (kg/s):
            //   (M^{n+1}-M^n)/dt + div(F) - Q = 0,  where M = φ V ρ S
            const double acc_w = (Mw_np1 - Mw_n) / dt;
            const double acc_g = (Mg_np1 - Mg_n) / dt;

            const double r_w = acc_w + div_w[i] - Qw;
            const double r_g = acc_g + div_g[i] - Qg;

            (*res_w)[i] = r_w;
            (*res_g)[i] = r_g;

            diag.global_res_w_rate += r_w;
            diag.global_res_g_rate += r_g;

            const double aw = std::abs(r_w);
            const double ag = std::abs(r_g);
            diag.l1_abs_res_w_rate += aw;
            diag.l1_abs_res_g_rate += ag;

            if (aw >= diag.max_abs_res_w_rate)
            {
                diag.max_abs_res_w_rate = aw;
                diag.max_abs_res_w_cellId = c.id;
            }
            if (ag >= diag.max_abs_res_g_rate)
            {
                diag.max_abs_res_g_rate = ag;
                diag.max_abs_res_g_cellId = c.id;
            }
        }

        // physical mass-balance mismatch using Mw_new-Mw_old (kg/s)
        diag.mass_balance_err_w_rate = (diag.Mw_new - diag.Mw_old) / dt + diag.Fw_out_sum - diag.Qw_sum;
        diag.mass_balance_err_g_rate = (diag.Mg_new - diag.Mg_old) / dt + diag.Fg_out_sum - diag.Qg_sum;

        return true;
    }

    inline bool appendPhaseMassBalanceCSVRow(
        const std::string& file,
        int step,
        double simTime,
        double dt,
        const PhaseMassBalanceDiagnostics& d)
    {
        if (file.empty()) return true;

        bool needHeader = true;
        {
            std::ifstream ifs(file, std::ios::binary);
            if (ifs)
            {
                ifs.seekg(0, std::ios::end);
                needHeader = (ifs.tellg() == 0);
            }
        }

        std::ofstream ofs(file, std::ios::app);
        if (!ofs)
        {
            std::cerr << "[MassDiag] cannot open '" << file << "' for append.\n";
            return false;
        }

        if (needHeader)
        {
            ofs << "step,time,dt,"
                << "Mw_old,Mw_new,Mg_old,Mg_new,"
                << "Qw_sum,Qg_sum,Fw_out_sum,Fg_out_sum,"
                << "global_res_w_rate,global_res_g_rate,"
                << "mass_balance_err_w_rate,mass_balance_err_g_rate,"
                << "max_abs_res_w_rate,max_abs_res_g_rate,"
                << "l1_abs_res_w_rate,l1_abs_res_g_rate,"
                << "max_abs_res_w_cellId,max_abs_res_g_cellId\n";
        }

        ofs << step << ','
            << std::setprecision(12) << simTime << ','
            << dt << ','
            << d.Mw_old << ','
            << d.Mw_new << ','
            << d.Mg_old << ','
            << d.Mg_new << ','
            << d.Qw_sum << ','
            << d.Qg_sum << ','
            << d.Fw_out_sum << ','
            << d.Fg_out_sum << ','
            << d.global_res_w_rate << ','
            << d.global_res_g_rate << ','
            << d.mass_balance_err_w_rate << ','
            << d.mass_balance_err_g_rate << ','
            << d.max_abs_res_w_rate << ','
            << d.max_abs_res_g_rate << ','
            << d.l1_abs_res_w_rate << ','
            << d.l1_abs_res_g_rate << ','
            << d.max_abs_res_w_cellId << ','
            << d.max_abs_res_g_cellId << '\n';

        return true;
    }
} // namespace IMPES_Iteration

#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "TemperatureBCAdapter.h"
#include "FVM_WellDOF_TwoPhase.h"
#include "FVM_WellCoupling_TwoPhase.h"
#include "DiffusionCentral.h"
#include "ConvectionUpwind.h"
#include "Timeterm_BDF.h"
#include "Solver_AssemblerCOO.h"
#include "LinearSolver_Eigen.h"
#include "Solver_TimeLoopUtils.h"
#include "CapRelPerm.h"
#include "IMPES_CommonUtils.h"

namespace IMPES
{
    inline std::vector<double> gatherFieldVector(
        const FieldRegistry& reg,
        Mesh& mesh,
        const std::string& name,
        const std::vector<int>& lid_of_cell,
        int nUnknowns)
    {
        std::vector<double> out(nUnknowns, 0.0);
        auto fld = reg.get<volScalarField>(name);
        if (!fld) return out;

        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const int r = lid_of_cell[i];
            if (r >= 0 && r < nUnknowns)
            {
                out[r] = (*fld)[i];
            }
        }
        return out;
    }

    inline void scatterVectorToField(
        FieldRegistry& reg,
        Mesh& mesh,
        const std::string& name,
        const std::vector<int>& lid_of_cell,
        const std::vector<double>& vec)
    {
        auto fld = reg.getOrCreate<volScalarField>(name, mesh.getCells().size(), 0.0);
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const int r = lid_of_cell[i];
            if (r >= 0 && static_cast<size_t>(r) < vec.size())
            {
                (*fld)[i] = vec[r];
            }
        }
    }

    struct TemperatureAssemblyConfig
    {
        std::string operator_tag = "T_impes";
        std::string temperature_field = "T";
        std::string temperature_old_field = "T_old";
        std::string temperature_prev_field = "T_prev";
        std::string C_eff_field = "C_eff";
        std::string lambda_eff_field = "lambda_eff";
        std::string mass_flux_water = "mf_w";
        std::string mass_flux_gas = "mf_g";
        std::string cp_water_field = "cp_w";
        std::string cp_gas_field = "cp_g";
        std::string lambda_water_field = "lambda_w";
        std::string lambda_gas_field = "lambda_g";
        std::string rho_water_field = "rho_w";
        std::string rho_gas_field = "rho_g";
        std::string pressure_field = "p_w";
        Vector gravity = { 0.0, -9.81, 0.0 };
        bool enable_diffusion = true;
        bool enable_convection = true;
        bool include_well_heat = true;
        VGParams vg_params;
        RelPermParams rp_params;
    };

    struct TemperatureSolveControls
    {
        TemperatureAssemblyConfig assembly;
        LinearSolverOptions linear;
        double under_relax = 1.0;
        double tol_abs = 1e-6;
        double tol_rel = 1e-6;
        double urf_min = 0.15;
        double urf_max = 0.7;
        double urf_step = 0.05;
        int max_outer_iterations = 10;
        bool report_outer_iterations = false;
    };

    struct TemperatureSolveReport
    {
        double lin_residual = 0.0;
        int lin_iterations = 0;
        double dT_inf = 0.0;
    };

    namespace detail
    {
        inline bool buildConvectionContribution(
            MeshManager& mgr,
            FieldRegistry& reg,
            FaceFieldRegistry& freg,
            const TemperatureBCAdapter& Tbc,
            const std::string& temperature_field,
            const std::string& flux_field,
            const std::vector<std::string>& scalars,
            std::vector<double>& acc_aPP,
            std::vector<double>& acc_aPN,
            std::vector<double>& acc_bP)
        {
            if (flux_field.empty()) return true;

            OperatorFieldNames tmpNames = makeNames("tmp_conv_" + flux_field);
            if (!FVM::Convection::build_FaceCoeffs_Upwind(
                mgr, reg, freg,
                temperature_field,
                flux_field,
                scalars,
                tmpNames,
                Tbc))
            {
                return false;
            }

            auto tmp_aPP = freg.get<faceScalarField>(tmpNames.aPP_conv.c_str());
            auto tmp_aPN = freg.get<faceScalarField>(tmpNames.aPN_conv.c_str());
            auto tmp_bP = freg.get<faceScalarField>(tmpNames.bP_conv.c_str());
            if (!tmp_aPP || !tmp_aPN || !tmp_bP) return false;

            for (size_t i = 0; i < acc_aPP.size(); ++i)
            {
                acc_aPP[i] += (*tmp_aPP)[i];
                acc_aPN[i] += (*tmp_aPN)[i];
                acc_bP[i] += (*tmp_bP)[i];
            }
            return true;
        }

        inline void storeConvectionFaceCoeffs(
            FaceFieldRegistry& freg,
            const OperatorFieldNames& nm,
            const std::vector<double>& aPP,
            const std::vector<double>& aPN,
            const std::vector<double>& bP)
        {
            auto aPPf = freg.getOrCreate<faceScalarField>(nm.aPP_conv.c_str(), aPP.size(), 0.0);
            auto aPNf = freg.getOrCreate<faceScalarField>(nm.aPN_conv.c_str(), aPN.size(), 0.0);
            auto bPf = freg.getOrCreate<faceScalarField>(nm.bP_conv.c_str(), bP.size(), 0.0);
            for (size_t i = 0; i < aPP.size(); ++i)
            {
                (*aPPf)[i] = aPP[i];
                (*aPNf)[i] = aPN[i];
                (*bPf)[i] = bP[i];
            }
        }

        inline void computeWellPhaseRates(
            MeshManager& mgr,
            FieldRegistry& reg,
            const TemperatureAssemblyConfig& cfg,
            const WellDOF_TwoPhase& well,
            const volScalarField& pressure,
            volScalarField& Qw,
            volScalarField& Qg)
        {
            auto mask = reg.get<volScalarField>(well.mask_field);
            auto WI = reg.get<volScalarField>(well.PI_field_w);
            auto lambda_w = reg.get<volScalarField>(cfg.lambda_water_field.c_str());
            auto lambda_g = reg.get<volScalarField>(cfg.lambda_gas_field.c_str());
            auto rho_w = reg.get<volScalarField>(cfg.rho_water_field.c_str());
            auto rho_g = reg.get<volScalarField>(cfg.rho_gas_field.c_str());
            if (!mask || !WI || !lambda_w || !lambda_g || !rho_w || !rho_g) return;

            auto& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const auto& id2idx = mesh.getCellId2Index();

            for (const auto& c : cells)
            {
                if (c.id < 0) continue;
                const size_t i = id2idx.at(c.id);
                if ((*mask)[i] <= 0.0) continue;

                const double WI_i = (*WI)[i];
                const double pw = pressure[i];
                const double p_bh = well.target; // pressure-controlled assumption
                double Qw_mass = 0.0;
                double Qg_mass = 0.0;

                if (well.role == WellDOF_TwoPhase::Role::Injector)
                {
                    double sw = clampValue(well.s_w_bh, 0.0, 1.0);
                    double krw, krg;
                    kr_Mualem_vG(sw, cfg.vg_params, cfg.rp_params, krw, krg);
                    const double lam_w = krw / std::max(well.mu_w_inj, 1e-12);
                    const double lam_g = krg / std::max(well.mu_g_inj, 1e-12);
                    const double lam_t = lam_w + lam_g;
                    if (lam_t > 0.0)
                    {
                        const double Q_total = WI_i * lam_t * (p_bh - pw);
                        if (Q_total > 0.0)
                        {
                            const double fw = lam_w / lam_t;
                            const double fg = lam_g / lam_t;
                            Qw_mass = Q_total * fw * well.rho_w_inj;
                            Qg_mass = Q_total * fg * well.rho_g_inj;
                        }
                    }
                }
                else
                {
                    const double lam_w_res = std::max((*lambda_w)[i], 0.0);
                    const double lam_g_res = std::max((*lambda_g)[i], 0.0);
                    const double lam_t = lam_w_res + lam_g_res;
                    if (lam_t > 0.0)
                    {
                        const double Q_total = WI_i * lam_t * (pw - p_bh);
                        if (Q_total > 0.0)
                        {
                            const double fw = lam_w_res / lam_t;
                            const double fg = lam_g_res / lam_t;
                            Qw_mass = -Q_total * fw * std::max((*rho_w)[i], 0.0);
                            Qg_mass = -Q_total * fg * std::max((*rho_g)[i], 0.0);
                        }
                    }
                }

                Qw[i] = Qw_mass;
                Qg[i] = Qg_mass;
            }
        }

        inline void applyWellHeatSources(
            SparseSystemCOO& sys,
            MeshManager& mgr,
            FieldRegistry& reg,
            const TemperatureAssemblyConfig& cfg,
            const std::vector<WellDOF_TwoPhase>& wells,
            const volScalarField& pressure)
        {
            if (wells.empty()) return;
            const int Nc = static_cast<int>(mgr.mesh().getCells().size());

            for (const auto& well : wells)
            {
                volScalarField Qw_tmp("Qw_tmp", Nc, 0.0);
                volScalarField Qg_tmp("Qg_tmp", Nc, 0.0);
                computeWellPhaseRates(mgr, reg, cfg, well, pressure, Qw_tmp, Qg_tmp);
                FVM::TwoPhaseWellCoupling::couple_well_to_temperature_equation(
                    sys, mgr, reg, well, Qw_tmp, Qg_tmp);
            }
        }
    } // namespace detail

    inline bool assembleAndSolveTemperature_IMPES(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const TemperatureBCAdapter& Tbc,
        const std::vector<WellDOF_TwoPhase>& wells,
        double dt,
        const TemperatureSolveControls& ctrl,
        TemperatureSolveReport* report = nullptr)
    {
        const auto& cfg = ctrl.assembly;
        if (!(dt > 0.0))
        {
            std::cerr << "[IMPES][Temp] invalid dt.\n";
            return false;
        }

        if (!reg.has(cfg.temperature_field) || !reg.has(cfg.temperature_old_field) || !reg.has(cfg.C_eff_field))
        {
            std::cerr << "[IMPES][Temp] missing required temperature fields.\n";
            return false;
        }

        const OperatorFieldNames nm = makeNames(cfg.operator_tag);

        if (!FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Temperature(
            mgr, reg, dt,
            cfg.C_eff_field,
            cfg.temperature_old_field,
            nm.a_time, nm.b_time))
        {
            std::cerr << "[IMPES][Temp] time term failed.\n";
            return false;
        }

        if (cfg.enable_diffusion)
        {
            if (!FVM::Diffusion::build_FaceCoeffs_Central(
                mgr, reg, freg,
                nm.a_f_diff, nm.s_f_diff,
                cfg.temperature_field,
                { "iso:" + cfg.lambda_eff_field },
                "",
                FVM::Diffusion::RhoFaceMethod::Linear,
                cfg.gravity,
                Tbc,
                false,
                0))
            {
                std::cerr << "[IMPES][Temp] diffusion build failed.\n";
                return false;
            }
        }

        const auto& faces = mgr.mesh().getFaces();
        std::vector<double> acc_aPP(faces.size(), 0.0);
        std::vector<double> acc_aPN(faces.size(), 0.0);
        std::vector<double> acc_bP(faces.size(), 0.0);

        if (cfg.enable_convection)
        {
            if (!detail::buildConvectionContribution(
                mgr, reg, freg, Tbc,
                cfg.temperature_field,
                cfg.mass_flux_water,
                { cfg.cp_water_field },
                acc_aPP, acc_aPN, acc_bP))
            {
                std::cerr << "[IMPES][Temp] convection (water) failed.\n";
                return false;
            }
            if (!detail::buildConvectionContribution(
                mgr, reg, freg, Tbc,
                cfg.temperature_field,
                cfg.mass_flux_gas,
                { cfg.cp_gas_field },
                acc_aPP, acc_aPN, acc_bP))
            {
                std::cerr << "[IMPES][Temp] convection (gas) failed.\n";
                return false;
            }
            detail::storeConvectionFaceCoeffs(freg, nm, acc_aPP, acc_aPN, acc_bP);
        }
        else
        {
            auto aPP = freg.getOrCreate<faceScalarField>(nm.aPP_conv.c_str(), faces.size(), 0.0);
            auto aPN = freg.getOrCreate<faceScalarField>(nm.aPN_conv.c_str(), faces.size(), 0.0);
            auto bP = freg.getOrCreate<faceScalarField>(nm.bP_conv.c_str(), faces.size(), 0.0);
            std::fill(aPP->data.begin(), aPP->data.end(), 0.0);
            std::fill(aPN->data.begin(), aPN->data.end(), 0.0);
            std::fill(bP->data.begin(), bP->data.end(), 0.0);
        }

        std::string opExpr = "ddt";
        if (cfg.enable_diffusion) opExpr += "+diffusion";
        if (cfg.enable_convection) opExpr += "+convection";

        SparseSystemCOO sys;
        if (!assemble_COO(mgr, reg, freg, opExpr, nm, &sys))
        {
            std::cerr << "[IMPES][Temp] assemble_COO failed.\n";
            return false;
        }

        auto pressure = reg.get<volScalarField>(cfg.pressure_field.c_str());
        if (!pressure)
        {
            std::cerr << "[IMPES][Temp] missing pressure field '" << cfg.pressure_field << "' for well heat coupling.\n";
            return false;
        }
        if (cfg.include_well_heat)
        {
            detail::applyWellHeatSources(sys, mgr, reg, cfg, wells, *pressure);
        }

        int N = 0;
        auto lid_of_cell = buildUnknownMap(mgr.mesh(), N);
        auto Tvec = gatherFieldVector(reg, mgr.mesh(), cfg.temperature_field, lid_of_cell, N);

        double res = 0.0;
        int it = 0;
        if (!solveCOO_Eigen(sys, Tvec, ctrl.linear, &it, &res))
        {
            return false;
        }

        scatterVectorToField(reg, mgr.mesh(), cfg.temperature_field, lid_of_cell, Tvec);

        double dT = 0.0;
        if (!cfg.temperature_prev_field.empty())
        {
            const double urf = clampValue(ctrl.under_relax, 0.0, 1.0);
            if (urf < 1.0)
            {
                underRelaxInPlace(reg, cfg.temperature_field, cfg.temperature_prev_field, urf);
            }
            dT = maxAbsDiff(reg, cfg.temperature_field, cfg.temperature_prev_field);
            updatePrevIterates(reg, { { cfg.temperature_field, cfg.temperature_prev_field } });
        }

        if (report)
        {
            report->lin_residual = res;
            report->lin_iterations = it;
            report->dT_inf = dT;
        }
        return true;
    }
} // namespace IMPES

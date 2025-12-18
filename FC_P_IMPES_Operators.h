#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "BCAdapter.h"
#include "ConvectionUpwind_Flux.h"
#include "DiffusionCentral.h"
#include "FaceFieldRegistry.h"
#include "FieldRegistry.h"
#include "MeshManager.h"
#include "Solver_AssemblerCOO.h"
#include "SolverContrlStrName.h"

// FC-IMPES-I: phase-by-phase discretization (build operators per phase, then add discretizations).
// This header only provides reusable building blocks; it does NOT change the legacy IMPES path.

namespace FC_P_IMPES_I
{
    struct PhaseOperatorConfig
    {
        // Operator tags used to name face fields via makeNames(tag).
        std::string tag_total = "fc_impes_total";
        std::string tag_w = "fc_impes_w";
        std::string tag_g_p = "fc_impes_g_p";   // gas pressure-driven part (uses p_w)
        std::string tag_g_pc = "fc_impes_g_pc"; // gas capillary part (uses Pc)

        // Primary fields
       
        std::string pressure_field = "p_w_FC_IMPES"; // p_w in current codebase
        std::string Pc_field = "Pc_FC_IMPES";

        // Derived per-cell coefficients (mass mobility)
        std::string rho_coeff_total_field = "fc_rho_coeff_total"; // (λwρw + λgρg)
        std::string rho_coeff_w_field = "fc_rho_coeff_w";         // (λwρw)
        std::string rho_coeff_g_field = "fc_rho_coeff_g";         // (λgρg)

        // Discretization knobs (match current IMPES default first)
        std::vector<std::string> mobility_tokens = { "kxx:kxx", "kyy:kyy", "kzz:kzz" };
        Vector gravity = { 0.0, 0.0, 0.0 };
        int gradient_smoothing = 0;
        FVM::Diffusion::RhoFaceMethod rho_face_method = FVM::Diffusion::RhoFaceMethod::Linear;

        // Flux output names (face fields)
        std::string mf_total = "mf_total_fc";
        std::string mf_w = "mf_w_fc";
        std::string mf_g = "mf_g_fc";
        std::string mf_g_p = "mf_g_p_fc";
        std::string mf_g_pc = "mf_g_pc_fc";
    };

    struct FaceCoeffDiffSummary
    {
        double max_abs_a = 0.0;
        double max_abs_s = 0.0;
        int face_id_max_a = -1;
        int face_id_max_s = -1;
    };

    inline bool computePhaseMassMobilityFields(MeshManager& mgr, FieldRegistry& reg, const PhaseOperatorConfig& cfg)
    {
        const auto lambda_w = reg.get<volScalarField>(PhysicalProperties_string::Water().lambda_w_tag);
        const auto lambda_g = reg.get<volScalarField>(PhysicalProperties_string::CO2().lambda_g_tag);
        const auto rho_w = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_tag);
        const auto rho_g = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_tag);
        if (!lambda_w || !lambda_g || !rho_w || !rho_g)
        {
            std::cerr << "[FC-IMPES-I][Ops] missing lambda/rho fields.\n";
            return false;
        }

        const size_t nCells = mgr.mesh().getCells().size();
        auto rho_coeff_total = reg.getOrCreate<volScalarField>(cfg.rho_coeff_total_field, nCells, 0.0);
        auto rho_coeff_w = reg.getOrCreate<volScalarField>(cfg.rho_coeff_w_field, nCells, 0.0);
        auto rho_coeff_g = reg.getOrCreate<volScalarField>(cfg.rho_coeff_g_field, nCells, 0.0);
        if (!rho_coeff_total || !rho_coeff_w || !rho_coeff_g)
        {
            std::cerr << "[FC-IMPES-I][Ops] failed to allocate rho_coeff fields.\n";
            return false;
        }

        const auto& cells = mgr.mesh().getCells();
        const auto& id2idx = mgr.mesh().getCellId2Index();
        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const double lw = std::max((*lambda_w)[i], 0.0);
            const double lg = std::max((*lambda_g)[i], 0.0);
            const double rw = std::max((*rho_w)[i], 0.0);
            const double rg = std::max((*rho_g)[i], 0.0);
            const double w = lw * rw;
            const double g = lg * rg;
            (*rho_coeff_w)[i] = w;
            (*rho_coeff_g)[i] = g;
            (*rho_coeff_total)[i] = w + g;
        }
        return true;
    }

    inline bool buildFaceCoeffsCentral(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& pbc,
        const std::string& operator_tag,
        const std::string& x_field,
        const std::string& rho_coeff_field,
        const PhaseOperatorConfig& cfg)
    {
        const OperatorFieldNames nm = makeNames(operator_tag);
        const bool ok = FVM::Diffusion::TwoPhaseDiffusionTemplate::build_FaceCoeffs_Central(
            mgr,
            reg,
            freg,
            nm.a_f_diff,
            nm.s_f_diff,
            x_field,
            cfg.mobility_tokens,
            rho_coeff_field,
            "",
            cfg.rho_face_method,
            cfg.gravity,
            pbc,
            false,
            cfg.gradient_smoothing);
        if (!ok)
        {
            std::cerr << "[FC-IMPES-I][Ops] build_FaceCoeffs_Central failed for tag=" << operator_tag << ".\n";
        }
        return ok;
    }

    inline bool buildPhaseMassFluxes(
        MeshManager& mgr,
        const FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& pbc,
        const PhaseOperatorConfig& cfg)
    {
        const OperatorFieldNames nmW = makeNames(cfg.tag_w);
        const OperatorFieldNames nmG_p = makeNames(cfg.tag_g_p);
        const OperatorFieldNames nmG_pc = makeNames(cfg.tag_g_pc);

        if (!FVM::Convection::buildFlux_Darcy_Mass(
            mgr, reg, freg,
            nmW.a_f_diff, nmW.s_f_diff,
            cfg.pressure_field,
            PhysicalProperties_string::Water().rho_tag,
            cfg.mf_w,
            cfg.mf_w + "_Q",
            cfg.mf_w + "_ufn",
            &pbc))
        {
            std::cerr << "[FC-IMPES-I][Flux] buildFlux_Darcy_Mass failed for water.\n";
            return false;
        }

        if (!FVM::Convection::buildFlux_Darcy_Mass(
            mgr, reg, freg,
            nmG_p.a_f_diff, nmG_p.s_f_diff,
            cfg.pressure_field,
            PhysicalProperties_string::CO2().rho_tag,
            cfg.mf_g_p,
            cfg.mf_g_p + "_Q",
            cfg.mf_g_p + "_ufn",
            &pbc))
        {
            std::cerr << "[FC-IMPES-I][Flux] buildFlux_Darcy_Mass failed for gas pressure-part.\n";
            return false;
        }

        if (!FVM::Convection::buildFlux_Darcy_Mass(
            mgr, reg, freg,
            nmG_pc.a_f_diff, nmG_pc.s_f_diff,
            cfg.Pc_field,
            PhysicalProperties_string::CO2().rho_tag,
            cfg.mf_g_pc,
            cfg.mf_g_pc + "_Q",
            cfg.mf_g_pc + "_ufn",
            nullptr))
        {
            std::cerr << "[FC-IMPES-I][Flux] buildFlux_Darcy_Mass failed for gas capillary-part.\n";
            return false;
        }

        const auto mf_w = freg.get<faceScalarField>(cfg.mf_w.c_str());
        const auto mf_g_p = freg.get<faceScalarField>(cfg.mf_g_p.c_str());
        const auto mf_g_pc = freg.get<faceScalarField>(cfg.mf_g_pc.c_str());
        if (!mf_w || !mf_g_p || !mf_g_pc)
        {
            std::cerr << "[FC-IMPES-I][Flux] missing intermediate mf fields.\n";
            return false;
        }

        auto mf_g = freg.getOrCreate<faceScalarField>(cfg.mf_g.c_str(), mf_g_p->data.size(), 0.0);
        auto mf_total = freg.getOrCreate<faceScalarField>(cfg.mf_total.c_str(), mf_g_p->data.size(), 0.0);
        if (!mf_g || !mf_total)
        {
            std::cerr << "[FC-IMPES-I][Flux] failed to allocate mf_g/mf_total.\n";
            return false;
        }

        for (size_t i = 0; i < mf_total->data.size(); ++i)
        {
            (*mf_g)[i] = (*mf_g_p)[i] + (*mf_g_pc)[i];
            (*mf_total)[i] = (*mf_w)[i] + (*mf_g)[i];
        }
        return true;
    }

    inline bool compareFaceCoeffsCombinedVsSummed(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& pbc,
        const PhaseOperatorConfig& cfg,
        FaceCoeffDiffSummary* summary_out = nullptr)
    {
        if (!computePhaseMassMobilityFields(mgr, reg, cfg))
        {
            return false;
        }

        if (!buildFaceCoeffsCentral(mgr, reg, freg, pbc, cfg.tag_total, cfg.pressure_field, cfg.rho_coeff_total_field, cfg))
        {
            return false;
        }
        if (!buildFaceCoeffsCentral(mgr, reg, freg, pbc, cfg.tag_w, cfg.pressure_field, cfg.rho_coeff_w_field, cfg))
        {
            return false;
        }
        if (!buildFaceCoeffsCentral(mgr, reg, freg, pbc, cfg.tag_g_p, cfg.pressure_field, cfg.rho_coeff_g_field, cfg))
        {
            return false;
        }

        const OperatorFieldNames nmT = makeNames(cfg.tag_total);
        const OperatorFieldNames nmW = makeNames(cfg.tag_w);
        const OperatorFieldNames nmG = makeNames(cfg.tag_g_p);

        const auto aT = freg.get<faceScalarField>(nmT.a_f_diff.c_str());
        const auto sT = freg.get<faceScalarField>(nmT.s_f_diff.c_str());
        const auto aW = freg.get<faceScalarField>(nmW.a_f_diff.c_str());
        const auto sW = freg.get<faceScalarField>(nmW.s_f_diff.c_str());
        const auto aG = freg.get<faceScalarField>(nmG.a_f_diff.c_str());
        const auto sG = freg.get<faceScalarField>(nmG.s_f_diff.c_str());
        if (!aT || !sT || !aW || !sW || !aG || !sG)
        {
            std::cerr << "[FC-IMPES-I][Ops] missing face coeff fields for comparison.\n";
            return false;
        }

        FaceCoeffDiffSummary sum;
        const size_t nF = aT->data.size();
        for (size_t i = 0; i < nF; ++i)
        {
            const double da = (*aT)[i] - ((*aW)[i] + (*aG)[i]);
            const double ds = (*sT)[i] - ((*sW)[i] + (*sG)[i]);
            const double ada = std::abs(da);
            const double ads = std::abs(ds);
            if (ada > sum.max_abs_a)
            {
                sum.max_abs_a = ada;
                sum.face_id_max_a = static_cast<int>(i) + 1;
            }
            if (ads > sum.max_abs_s)
            {
                sum.max_abs_s = ads;
                sum.face_id_max_s = static_cast<int>(i) + 1;
            }
        }

        if (summary_out) *summary_out = sum;
        return true;
    }
}


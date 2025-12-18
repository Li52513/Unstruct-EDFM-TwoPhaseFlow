#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "FC_P_IMPES_Operators.h"
#include "PressureEqAssemblerandSolver.h" // reuse SparseSystemCOO, PressureAssemblyConfig/Result, well coupling, time-term

namespace FC_P_IMPES_I
{
    namespace detailed
    {
        inline void applyCOO_cellPart(const SparseSystemCOO& sys, int nCells, const std::vector<double>& x, std::vector<double>& y)
        {
            const int n = std::min(sys.n, nCells);
            y.assign(n, 0.0);
            for (const auto& t : sys.A)
            {
                if (t.r < 0 || t.c < 0) continue;
                if (t.r >= n || t.c >= n) continue;
                if (t.c >= static_cast<int>(x.size())) continue;
                y[t.r] += t.v * x[t.c];
            }
            for (int i = 0; i < n && i < static_cast<int>(sys.b.size()); ++i)
            {
                y[i] -= sys.b[i];
            }
        }

        inline void addSystemScaled(SparseSystemCOO& dst, const SparseSystemCOO& src, double scale)
        {
            if (scale == 0.0) return;
            dst.expandTo(std::max(dst.n, src.n));
            if (dst.b.size() < static_cast<size_t>(dst.n)) dst.b.resize(static_cast<size_t>(dst.n), 0.0);

            if (src.b.size() > dst.b.size()) dst.b.resize(src.b.size(), 0.0);
            for (size_t i = 0; i < src.b.size(); ++i)
            {
                dst.b[i] += scale * src.b[i];
            }

            dst.A.reserve(dst.A.size() + src.A.size());
            for (const auto& t : src.A)
            {
                dst.A.push_back({ t.r, t.c, scale * t.v });
            }
        }
    }

    inline bool assemblePressureTwoPhase_FC(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& pbc,
        const std::vector<WellDOF_TwoPhase>& wells,
        double dt,
        const IMPES_Iteration::PressureAssemblyConfig& legacyCfg,
        IMPES_Iteration::PressureAssemblyResult& result,
        PhaseOperatorConfig cfg = PhaseOperatorConfig{})
    {
        if (dt <= 0.0)
        {
            std::cerr << "[FC-IMPES-I][Pressure] invalid dt.\n";
            return false;
        }

        cfg.pressure_field = legacyCfg.pressure_field;
        cfg.Pc_field = legacyCfg.Pc_field;
        cfg.gravity = legacyCfg.gravity;
        cfg.gradient_smoothing = legacyCfg.gradient_smoothing;

        if (!computePhaseMassMobilityFields(mgr, reg, cfg))
        {
            return false;
        }

        // ---- time term (reuse legacy time-term assembler) ----
        const OperatorFieldNames nmT = makeNames(cfg.tag_total);
        if (!IMPES_Iteration::TimeTerm_IMPES_Pressure(
            mgr,
            reg,
            dt,
            legacyCfg.pressure_old_field,
            legacyCfg.pressure_field,
            nmT.a_time,
            nmT.b_time))
        {
            std::cerr << "[FC-IMPES-I][Pressure] time-term assembly failed.\n";
            return false;
        }

        // ---- diffusion operators: Aw(pw) and Ag(pw) ----
        if (!buildFaceCoeffsCentral(mgr, reg, freg, pbc, cfg.tag_w, legacyCfg.pressure_field, cfg.rho_coeff_w_field, cfg))
            return false;
        if (!buildFaceCoeffsCentral(mgr, reg, freg, pbc, cfg.tag_g_p, legacyCfg.pressure_field, cfg.rho_coeff_g_field, cfg))
            return false;

        SparseSystemCOO sys_time, sys_w, sys_g;
        if (!assemble_COO(mgr, reg, freg, "ddt", nmT, &sys_time))
        {
            std::cerr << "[FC-IMPES-I][Pressure] assemble_COO(ddt) failed.\n";
            return false;
        }
        if (!assemble_COO(mgr, reg, freg, "diffusion", makeNames(cfg.tag_w), &sys_w))
        {
            std::cerr << "[FC-IMPES-I][Pressure] assemble_COO(diffusion, w) failed.\n";
            return false;
        }
        if (!assemble_COO(mgr, reg, freg, "diffusion", makeNames(cfg.tag_g_p), &sys_g))
        {
            std::cerr << "[FC-IMPES-I][Pressure] assemble_COO(diffusion, g) failed.\n";
            return false;
        }

        SparseSystemCOO sys;
        detailed::addSystemScaled(sys, sys_time, 1.0);
        detailed::addSystemScaled(sys, sys_w, 1.0);
        detailed::addSystemScaled(sys, sys_g, 1.0);

        // ---- build unknown map (cells first) ----
        int nCellUnknowns = 0;
        result.cell_lid = buildUnknownMap(mgr.mesh(), nCellUnknowns);
        sys.expandTo(nCellUnknowns);

        // ---- capillary contribution: move Ag(Pc) to RHS ----
        PressureBC::Registry nullReg;
        PressureBCAdapter nullPbc{ nullReg };

        if (!buildFaceCoeffsCentral(mgr, reg, freg, nullPbc, cfg.tag_g_pc, legacyCfg.Pc_field, cfg.rho_coeff_g_field, cfg))
            return false;

        SparseSystemCOO sys_cap;
        if (!assemble_COO(mgr, reg, freg, "diffusion", makeNames(cfg.tag_g_pc), &sys_cap))
        {
            std::cerr << "[FC-IMPES-I][Pressure] assemble_COO(diffusion, g_pc) failed.\n";
            return false;
        }
        sys_cap.expandTo(nCellUnknowns);

        std::vector<double> Pc_vec(nCellUnknowns, 0.0);
        if (!IMPES_Iteration::detailed::gatherFieldToVector(mgr.mesh(), reg, legacyCfg.Pc_field, result.cell_lid, nCellUnknowns, Pc_vec))
        {
            std::cerr << "[FC-IMPES-I][Pressure] gather Pc failed.\n";
            return false;
        }
        std::vector<double> capContribution;
        detailed::applyCOO_cellPart(sys_cap, nCellUnknowns, Pc_vec, capContribution); // = A*Pc - b
        for (int r = 0; r < nCellUnknowns; ++r)
        {
            sys.addb(r, -capContribution[static_cast<size_t>(r)]);
        }

        // ---- wells: keep identical coupling as legacy ----
        int desired_n = sys.n;
        for (const auto& well : wells)
        {
            desired_n = std::max(desired_n, well.lid + 1);
        }
        sys.expandTo(desired_n);

        for (const auto& well : wells)
        {
            FVM::TwoPhaseWellsStrict::couple_well_to_pressure_equation_strict_rate(
                sys,
                mgr,
                reg,
                well,
                result.cell_lid,
                legacyCfg.VG_Parameter.vg_params,
                legacyCfg.VG_Parameter.relperm_params);
        }

        result.system = std::move(sys);
        return true;
    }
}


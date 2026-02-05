#pragma once
#include <algorithm>
#include <cmath>
#include "InitConfig.h"
#include "CapRelPerm.h"

// A lightweight capillary/rel-perm model for debugging:
// - Capillary pressure: linear with effective saturation, Pc = pc_entry * (1 - Se)
// - Rel perms: Corey-type krw = Se^nw, krg = (1-Se)^ng
// Use this to reduce nonlinearity compared to Mualem–van Genuchten.

namespace SimpleCapRelPerm
{
    struct PcConfig
    {
        double pc_entry = 1e5;   // Pc at Se = 0
        double pc_max   = 1e7;   // clamp to avoid overflow
        double Se_tol   = 1e-8;  // tolerance for turning off slope near endpoints
    };

    struct RelPermConfig
    {
        double nw = 2.0;   // Corey exponent for water
        double ng = 2.0;   // Corey exponent for gas
        double krw_end = 1.0; // krw at Se=1
        double krg_end = 1.0; // krg at Se=0
    };

    inline bool params_valid(const VGParams& vg)
    {
        const double denom = 1.0 - vg.Swr - vg.Sgr;
        return (vg.Swr >= 0.0) && (vg.Sgr >= 0.0) && (denom > 1e-12);
    }

    inline double clamp01(double v) { return std::max(0.0, std::min(1.0, v)); }

    inline double Se_from_Sw_simple(double Sw, const VGParams& vg)
    {
        const double denom = 1.0 - vg.Swr - vg.Sgr;
        if (denom <= 1e-12) return 0.0;
        return clamp01((Sw - vg.Swr) / denom);
      
    }

    inline double Sw_from_Se_simple(double Se, const VGParams& vg)
    {
        Se = clamp01(Se);
        const double denom = 1.0 - vg.Swr - vg.Sgr;
        return vg.Swr + Se * denom;
    }

    // Linear Pc: Pc = pc_entry * (1 - Se)
    inline double pc_linear(double Sew, const VGParams& vg, const PcConfig& cfg = {})
    {
        if (!params_valid(vg) || cfg.pc_entry <= 0.0) return 0.0;
        const double pc = cfg.pc_entry * (1.0 - Sew);
        return std::max(0.0, std::min(cfg.pc_max, pc));
    }

    inline double dpc_dSw_linear(double Sew, const VGParams& vg, const PcConfig& cfg = {})
    {
        if (!params_valid(vg) || cfg.pc_entry <= 0.0) return 0.0;
        const double denom = 1.0 - vg.Swr - vg.Sgr;
        if (Sew <= cfg.Se_tol || Sew >= 1.0 - cfg.Se_tol) return 0.0; // flatten near endpoints
        return -cfg.pc_entry / denom;
    }

    inline double Sw_from_pc_linear(double pc, const VGParams& vg, const PcConfig& cfg = {})
    {
        if (!params_valid(vg) || cfg.pc_entry <= 0.0) return 1.0 - vg.Sgr;
        if (pc <= 0.0) return 1.0 - vg.Sgr;
        const double Se = clamp01(1.0 - pc / cfg.pc_entry);
        return Sw_from_Se_simple(Se, vg);
    }

    // Corey rel perms
    inline void kr_corey(double Sew, const VGParams& vg, const RelPermConfig& cfg, double& krw, double& krg)
    {
        krw = cfg.krw_end * std::pow(Sew, cfg.nw);
        krg = cfg.krg_end * std::pow(1.0 - Sew, cfg.ng);
        krw = clamp01(krw);
        krg = clamp01(krg);
    }

    inline void dkr_dSw_corey(double Sw, const VGParams& vg, const RelPermConfig& cfg, double& dkrw_dSw, double& dkrg_dSw)
    {
        dkrw_dSw = dkrg_dSw = 0.0;
        if (!params_valid(vg)) return;

        const double denom = 1.0 - vg.Swr - vg.Sgr;
        const double Se = Se_from_Sw_simple(Sw, vg);

        if (cfg.nw > 0.0 && Se > 0.0)
        {
            dkrw_dSw = cfg.krw_end * cfg.nw * std::pow(Se, cfg.nw - 1.0) * (1.0 / denom);
        }
        if (cfg.ng > 0.0 && Se < 1.0)
        {
            dkrg_dSw = -cfg.krg_end * cfg.ng * std::pow(1.0 - Se, cfg.ng - 1.0) * (1.0 / denom);
        }
    }
} // namespace SimpleCapRelPerm

#pragma once
// CapRelPerm.h
#pragma once
#include "InitConfig.h"
#include <algorithm>
#include <cmath>

// ===== common epsilon & helpers =====//
constexpr double kTiny = 1e-12; // 避免除零等数值问题的小量
constexpr double kPcMax = 1e12; // 最大毛细压力，单位 Pa

inline bool vg_params_valid(const VGParams& vg) 
{
    if (!(vg.alpha > 0.0)) return false;
    if (!(vg.n > 1.0)) return false;
    if (vg.Swr < 0.0 || vg.Sgr < 0.0) return false;
    if (vg.Swr + vg.Sgr >= 1.0 - kTiny) return false;
    return true;
}

// 限幅函数

inline double clamp(double v, double lo, double hi) { return std::max(lo, std::min(hi, v)); } 



//有效饱和度Se = (Sw-Swr)/(1-Swr-Sgr)计算

inline double Se_from_Sw(double Sw, const VGParams& vg)
{
    const double denom = 1.0 - vg.Swr - vg.Sgr;
    if (denom <= kTiny) return 0.0;               // 异常参数直接回退
    return clamp((Sw - vg.Swr) / denom, 0.0, 1.0);
}

//毛细压力 van Genuchten毛细压力模型，pc = (1/alpha)*((Se)^(-1/m)-1)^(1/n) 计算 

inline double pc_vG(double Sw, const VGParams& vg, double pc_max = kPcMax) // Pa
{
    if (!vg_params_valid(vg)) return pc_max;

    double Se = Se_from_Sw(Sw, vg);
    if (Se <= kTiny) return pc_max;               // 极干：返回上限
    if (Se >= 1.0 - kTiny) return 0.0;            // 近饱和：Pc≈0，更稳

    const double m = vg.m();
    const double inv_m = 1.0 / m;
    const double core = std::max(std::pow(Se, -inv_m) - 1.0, 0.0);

    double pc = (1.0 / vg.alpha) * std::pow(core, 1.0 / vg.n); // 单位：Pa
    if (!std::isfinite(pc)) pc = pc_max;
    return clamp(pc, 0.0, pc_max);
}

// 反解：Pc(=Pa) → Sw
inline double Sw_from_pc(double pc, const VGParams& vg)
{
    if (!vg_params_valid(vg)) return 1.0;
    if (pc <= 0.0) return 1.0 - vg.Sgr;  // Pc=0 ⇒ Se=1 ⇒ Sw=Swr+(1-Swr-Sgr)

    const double m = vg.m();
    const double aPn = std::pow(vg.alpha * pc, vg.n);
    double Se = std::pow(1.0 + aPn, -m);
    Se = clamp(Se, 0.0, 1.0);

    const double denom = 1.0 - vg.Swr - vg.Sgr;
    return vg.Swr + Se * denom;
}

// dPc/dSw（把 Pc 当未知量时用）
inline double dpc_dSw_vG(double Sw, const VGParams& vg, double cap = kPcMax)
{
    if (!vg_params_valid(vg)) return 0.0;

    const double denom = std::max(1.0 - vg.Swr - vg.Sgr, kTiny);
    const double Se = Se_from_Sw(Sw, vg);
    if (Se <= kTiny || Se >= 1.0 - kTiny) return 0.0;

    const double m = vg.m();
    const double inv_m = 1.0 / m;
    const double core = std::max(std::pow(Se, -inv_m) - 1.0, kTiny); // >0
    const double core_p = std::pow(core, 1.0 / vg.n - 1.0);           // core^(1/n - 1)

    // d(core)/dSe = (-1/m) * Se^{-1/m - 1}
    double dpc_dSe = (1.0 / vg.alpha) * (1.0 / vg.n) * core_p * (-inv_m) * std::pow(Se, -inv_m - 1.0);
    double dpc_dSw = dpc_dSe / denom;

    if (!std::isfinite(dpc_dSw)) return 0.0;
    if (std::abs(dpc_dSw) > cap) dpc_dSw = (dpc_dSw > 0 ? cap : -cap);
    return dpc_dSw; // 物理上应≤0
}

// dSw/dPc（把 Sw 当未知量时用）
inline double dSw_dpc_vG(double pc, const VGParams& vg)
{
    if (!vg_params_valid(vg)) return 0.0;
    if (pc <= 0.0) return 0.0;

    const double denom = 1.0 - vg.Swr - vg.Sgr;
    const double m = vg.m();
    const double aP = vg.alpha * pc;
    const double aPn = std::pow(aP, vg.n);

    // dSe/dPc = -m [1+(αPc)^n]^(-m-1) * n (αPc)^(n-1) * α
    double dSe_dpc = -m * std::pow(1.0 + aPn, -m - 1.0) * vg.n * std::pow(aP, vg.n - 1.0) * vg.alpha;
    double dSw_dpc = dSe_dpc * denom;

    return std::isfinite(dSw_dpc) ? dSw_dpc : 0.0;
}

// krw = Se^L * [1 - (1 - Se^{1/m})^m]^2
// krg = (1-Se)^L * (1 - Se^{1/m})^{2m}

inline void kr_Mualem_vG(double Sw, const VGParams& vg, const RelPermParams& rp,
    double& krw, double& krg)
{
    const double Se = clamp(Se_from_Sw(Sw, vg), kTiny, 1.0 - kTiny);
    const double m = vg.m();
    const double L = rp.L;

    const double Se_1m = std::pow(Se, 1.0 / m);
    double A = 1.0 - Se_1m;               // A = 1 - Se^{1/m}
    if (A < kTiny) A = kTiny;
    const double krw_base = 1.0 - std::pow(A, m);  // = 1 - (1 - Se^{1/m})^m

    // Mualem–vG
    krw = std::pow(Se, L) * std::pow(krw_base, 2.0);
    krg = std::pow(1.0 - Se, L) * std::pow(1.0 - krw_base, 2.0); // = (1-Se)^L * (A^m)^2

    krw = clamp(krw, 0.0, 1.0);
    krg = clamp(krg, 0.0, 1.0);
}


inline void dkr_dSw_Mualem_vG(double Sw, const VGParams& vg, const RelPermParams& rp,
    double& dkrw_dSw, double& dkrg_dSw)
{
    dkrw_dSw = dkrg_dSw = 0.0;
    if (!vg_params_valid(vg)) return;

    const double L = rp.L; // 连通性参数
    const double denom = std::max(1.0 - vg.Swr - vg.Sgr, kTiny);
    const double Se = clamp(Se_from_Sw(Sw, vg), kTiny, 1.0 - kTiny);
    const double m = vg.m();
    const double inv_m = 1.0 / m;
    const double dSe_dSw = 1.0 / denom;

    // 便捷量
    const double Se_p_1m = std::pow(Se, inv_m);         // Se^{1/m}
    const double A = 1.0 - Se_p_1m;               // A = 1 - Se^{1/m}
    const double Apos = std::max(A, kTiny);
    const double B = 1.0 - std::pow(Apos, m);     // B = 1 - (1 - Se^{1/m})^m = krw_base

    // d(Se^L)/dSe = L * Se^{L-1}
    const double dSeL_dSe = (L == 0.0 ? 0.0 : L * std::pow(Se, L - 1.0));
    // d(1-Se)^L/dSe = -L * (1-Se)^{L-1}
    const double d1mSeL_dSe = (L == 0.0 ? 0.0 : -L * std::pow(1.0 - Se, L - 1.0));

    // dA/dSe = -(1/m)*Se^{1/m - 1}
    const double dA_dSe = -(inv_m)*std::pow(Se, inv_m - 1.0);
    // dB/dSe = - m * A^{m-1} * dA/dSe = (A^{m-1}/m) * Se^{1/m - 1}
    const double dB_dSe = -m * std::pow(Apos, m - 1.0) * dA_dSe;

    // krw = Se^L * B^2
    const double krw = std::pow(Se, L) * B * B;
    // dkrw/dSe = (d(Se^L)/dSe)*B^2 + Se^L * 2*B*dB/dSe
    double dkrw_dSe = dSeL_dSe * (B * B) + std::pow(Se, L) * 2.0 * B * dB_dSe;

    // krg = (1-Se)^L * (1 - B)^2, 其中 (1 - B) = A^m
    const double one_minus_B = 1.0 - B; // = A^m
    // d(1-B)/dSe = -dB/dSe
    double dkrg_dSe = d1mSeL_dSe * (one_minus_B * one_minus_B)
        + std::pow(1.0 - Se, L) * 2.0 * (one_minus_B) * (-dB_dSe);

    dkrw_dSw = dkrw_dSe * dSe_dSw;
    dkrg_dSw = dkrg_dSe * dSe_dSw;

    if (!std::isfinite(dkrw_dSw)) dkrw_dSw = 0.0;
    if (!std::isfinite(dkrg_dSw)) dkrg_dSw = 0.0;
}


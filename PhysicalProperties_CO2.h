#pragma once
#pragma once
#include <algorithm>

inline double __clamp(double x, double a, double b) {
    return x < a ? a : (x > b ? b : x);
}

namespace CO2 {

    // k(T) [W/(m·K)]
    inline double k_W_mK_low(double T) {
        const double a0 = -1.13009e-3, a1 = 2.701846e-5, a2 = 1.571292e-7, a3 = -1.790745e-10, a4 = 6.671466e-14;
        return (((a4 * T + a3) * T + a2) * T + a1) * T + a0;
    }

    inline double k_W_mK_high(double T) {
        const double b0 = -1.0189761e-2, b1 = 9.648924e-5, b2 = -1.750928e-8, b3 = 1.867531e-12;
        return (((b3 * T + b2) * T + b1) * T + b0);
    }

    inline double k_W_mK(double T) {
        const double Tmin = 220.0, Tmid = 1000.0, Tmax = 3273.16;
        T = __clamp(T, Tmin, Tmax);
        return (T <= Tmid ? k_W_mK_low(T) : k_W_mK_high(T));
    }

    // cp_mass(T) [J/(kg·K)]
    inline double cp_mass_low(double T) {
        const double a0 = 448.3863096, a1 = 1.677003266, a2 = -1.291503313e-3, a3 = 4.008266979e-7;
        return ((a3 * T + a2) * T + a1) * T + a0;
    }
    inline double cp_mass_high(double T) {
        const double b0 = 618.3146915, b1 = 1.091034066, b2 = -6.308463656e-4, b3 = 1.73304167e-7, b4 = -1.826643607e-11;
        return ((((b4 * T + b3) * T + b2) * T + b1) * T + b0);
    }
    inline double cp_mass_J_kgK(double T) {
        const double Tmin = 293.0, Tmid = 800.0, Tmax = 3000.0;
        T = __clamp(T, Tmin, Tmax);
        return (T <= Tmid ? cp_mass_low(T) : cp_mass_high(T));
    }

    // cp_molar(T) [J/(mol·K)]（如需导出/诊断用）
    inline double cp_molar_low(double T) {
        const double a0 = 19.72900348, a1 = 0.07378814536, a2 = -5.682612568e-5, a3 = 1.763637588e-8;
        return ((a3 * T + a2) * T + a1) * T + a0;
    }
    inline double cp_molar_high(double T) {
        const double b0 = 27.20587016, b1 = 0.04800550056, b2 = -2.775724594e-5, b3 = 7.62538184e-9, b4 = -8.037229696e-13;
        return ((((b4 * T + b3) * T + b2) * T + b1) * T + b0);
    }
    inline double cp_molar_J_molK(double T) {
        const double Tmin = 293.0, Tmid = 800.0, Tmax = 3000.0;
        T = __clamp(T, Tmin, Tmax);
        return (T <= Tmid ? cp_molar_low(T) : cp_molar_high(T));
    }

    // rho(T) [kg/m^3], T∈[270,2000]
    inline double rho_CO2_kg_m3(double T) {
        const double Tmin = 270.0, Tmax = 2000.0;
        T = __clamp(T, Tmin, Tmax);
        return 539.7 / T;
    }

    // mu(T) [Pa·s], T∈[220,1000]
    inline double mu_CO2_Pa_s(double T) {
        const double Tmin = 220.0, Tmax = 1000.0;
        T = __clamp(T, Tmin, Tmax);
        const double c0 = -9.578337e-7, c1 = 5.774162e-8, c2 = -1.31496e-11, c3 = -7.474451e-15, c4 = 5.109362e-18;
        return (((c4 * T + c3) * T + c2) * T + c1) * T + c0;
    }

}
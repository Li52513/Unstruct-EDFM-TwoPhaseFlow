/**
 * @file AD_FluidEvaluator.h
 * @brief Automatic Differentiation Unified Fluid State Evaluator (AD Fluid Evaluator)
 * @details
 * Bridges double-valued state lookup tables (Span-Wagner / IAPWS lookup tables,
 * or CoolProp analytical EOS) to ADVar gradient images.
 * Provides 6 fluid properties: rho, mu, cp, cv, h, k with full AD gradients.
 * Features 4-tier robustness fallback and phase-mismatch detection.
 *
 * Compile-time switch:
 *   #define USE_COOLPROP_EOS  -> uses CoolProp analytical first_partial_deriv
 *   (default)               -> uses PropertyTable + robust numerical diff
 */

#pragma once

#include "ADVar.hpp"
#include "WaterPropertyTable.h"
#include "CO2PropertyTable.h"
#include "PropertiesSummary.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <exception>

#ifdef USE_COOLPROP_EOS
#include <AbstractState.h>
#include <memory>
#endif

namespace AD_Fluid
{
    // =========================================================
    // 1. Common type definitions
    // =========================================================

    /**
     * @brief Phase mismatch exception (numerical perturbation crosses phase boundary)
     */
    class PhaseMismatchException : public std::exception {
    public:
        const char* what() const noexcept override {
            return "Fluid Phase Mismatch detected across numerical perturbation.";
        }
    };

    /**
     * @brief Differentiation mode enum (V3 integration: tracks whether perturbation crossed a phase boundary)
     */
    enum class DiffMode { Central, Forward, Backward, ZeroGrad };

    /**
     * @struct ADFluidProperties
     * @brief Carries full AD-enabled fluid properties with gradients
     * @tparam N number of primary variables
     */
    template<int N>
    struct ADFluidProperties
    {
        ADVar<N> rho; ///< Density [kg/m^3]
        ADVar<N> mu;  ///< Dynamic viscosity [Pa*s]
        ADVar<N> cp;  ///< Isobaric specific heat [J/(kg*K)]
        ADVar<N> cv;  ///< Isochoric specific heat [J/(kg*K)]
        ADVar<N> h;   ///< Specific enthalpy [J/kg]
        ADVar<N> k;   ///< Thermal conductivity [W/(m*K)]

        bool isFallback = false;                  ///< True if deepest fallback was triggered (derivatives may be degraded)
        DiffMode diff_mode_P = DiffMode::Central; ///< Differentiation mode used for pressure derivative
        DiffMode diff_mode_T = DiffMode::Central; ///< Differentiation mode used for temperature derivative
        bool near_bound = false;                  ///< True if state is near a phase boundary
    };

    // =========================================================
    // 2. Main evaluator class
    // =========================================================

    class Evaluator
    {
    public:
        /**
         * @brief Evaluate water (H2O) fluid properties with AD gradients
         * @tparam N number of primary variables
         * @param P Water pressure (ADVar, Pa)
         * @param T Temperature (ADVar, K)
         * @return ADFluidProperties<N>
         */
        template<int N>
        static ADFluidProperties<N> evaluateWater(const ADVar<N>& P, const ADVar<N>& T)
        {
#ifdef USE_COOLPROP_EOS
            return evaluateWater_CoolProp(P, T);
#else
            double p_val = P.val;
            double t_val = T.val;
            auto& wt = WaterPropertyTable::instance();
            const double p_clamped = wt.clampPressure(p_val);
            const double t_clamped = wt.clampTemperature(t_val);

            WaterProperties base_props;
            ADFluidProperties<N> res;
            res.isFallback = false;
            if (p_clamped != p_val || t_clamped != t_val) {
                res.near_bound = true;
            }

            try {
                base_props = wt.getProperties(p_clamped, t_clamped);
            }
            catch (...) {
                base_props.rho = 1000.0;
                base_props.mu = 1e-3;
                base_props.cp = 4200.0;
                base_props.cv = 4182.0;
                base_props.h = 1.0e5;
                base_props.k = 0.6;
                res.isFallback = true;
            }

            double dP = std::max(10.0, std::min(std::abs(p_clamped) * 1e-6, 1000.0));
            double dT = std::max(0.001, std::min(std::abs(t_clamped) * 1e-6, 0.1));

            double dRho_dP = 0.0, dMu_dP = 0.0, dCp_dP = 0.0, dCv_dP = 0.0, dH_dP = 0.0, dK_dP = 0.0;
            double dRho_dT = 0.0, dMu_dT = 0.0, dCp_dT = 0.0, dCv_dT = 0.0, dH_dT = 0.0, dK_dT = 0.0;

            auto evalP = [&](double p_eval) { return wt.getProperties(p_eval, t_clamped); };
            res.diff_mode_P = computeRobustDerivatives(evalP, p_clamped, dP, base_props, dRho_dP, dMu_dP, dCp_dP, dCv_dP, dH_dP, dK_dP);

            auto evalT = [&](double t_eval) { return wt.getProperties(p_clamped, t_eval); };
            res.diff_mode_T = computeRobustDerivatives(evalT, t_clamped, dT, base_props, dRho_dT, dMu_dT, dCp_dT, dCv_dT, dH_dT, dK_dT);

            if (res.diff_mode_P != DiffMode::Central || res.diff_mode_T != DiffMode::Central) {
                res.near_bound = true;
            }

            AssembleADVar(res.rho, base_props.rho, dRho_dP, dRho_dT, P, T);
            AssembleADVar(res.mu,  base_props.mu,  dMu_dP,  dMu_dT,  P, T);
            AssembleADVar(res.cp,  base_props.cp,  dCp_dP,  dCp_dT,  P, T);
            AssembleADVar(res.cv,  base_props.cv,  dCv_dP,  dCv_dT,  P, T);
            AssembleADVar(res.h,   base_props.h,   dH_dP,   dH_dT,   P, T);
            AssembleADVar(res.k,   base_props.k,   dK_dP,   dK_dT,   P, T);

            return res;
#endif
        }

        /**
         * @brief Evaluate CO2 fluid properties with AD gradients
         */
        template<int N>
        static ADFluidProperties<N> evaluateCO2(const ADVar<N>& P, const ADVar<N>& T)
        {
#ifdef USE_COOLPROP_EOS
            return evaluateCO2_CoolProp(P, T);
#else
            double p_val = P.val;
            double t_val = T.val;
            auto& gt = CO2PropertyTable::instance();
            const double p_clamped = gt.clampPressure(p_val);
            const double t_clamped = gt.clampTemperature(t_val);

            CO2Properties base_props;
            ADFluidProperties<N> res;
            res.isFallback = false;
            if (p_clamped != p_val || t_clamped != t_val) {
                res.near_bound = true;
            }

            try {
                base_props = gt.getProperties(p_clamped, t_clamped);
            }
            catch (...) {
                base_props.rho = 800.0;
                base_props.mu = 1.48e-5;
                base_props.cp = 1100.0;
                base_props.cv = 850.0;
                base_props.h = 3.0e5;
                base_props.k = 0.03;
                res.isFallback = true;
            }

            double dP = std::max(10.0, std::min(std::abs(p_clamped) * 1e-6, 1000.0));
            double dT = std::max(0.001, std::min(std::abs(t_clamped) * 1e-6, 0.1));

            double dRho_dP = 0.0, dMu_dP = 0.0, dCp_dP = 0.0, dCv_dP = 0.0, dH_dP = 0.0, dK_dP = 0.0;
            double dRho_dT = 0.0, dMu_dT = 0.0, dCp_dT = 0.0, dCv_dT = 0.0, dH_dT = 0.0, dK_dT = 0.0;

            auto evalP = [&](double p_eval) { return gt.getProperties(p_eval, t_clamped); };
            res.diff_mode_P = computeRobustDerivatives(evalP, p_clamped, dP, base_props, dRho_dP, dMu_dP, dCp_dP, dCv_dP, dH_dP, dK_dP);

            auto evalT = [&](double t_eval) { return gt.getProperties(p_clamped, t_eval); };
            res.diff_mode_T = computeRobustDerivatives(evalT, t_clamped, dT, base_props, dRho_dT, dMu_dT, dCp_dT, dCv_dT, dH_dT, dK_dT);

            if (res.diff_mode_P != DiffMode::Central || res.diff_mode_T != DiffMode::Central) {
                res.near_bound = true;
            }

            AssembleADVar(res.rho, base_props.rho, dRho_dP, dRho_dT, P, T);
            AssembleADVar(res.mu,  base_props.mu,  dMu_dP,  dMu_dT,  P, T);
            AssembleADVar(res.cp,  base_props.cp,  dCp_dP,  dCp_dT,  P, T);
            AssembleADVar(res.cv,  base_props.cv,  dCv_dP,  dCv_dT,  P, T);
            AssembleADVar(res.h,   base_props.h,   dH_dP,   dH_dT,   P, T);
            AssembleADVar(res.k,   base_props.k,   dK_dP,   dK_dT,   P, T);

            return res;
#endif
        }

        /**
         * @brief Apply a small derivative floor to avoid zero-Jacobian fallback states.
         * @details
         * This routine is used only in the deepest fallback path when robust differencing
         * fails in all directions. It injects physically tiny but non-zero sensitivities,
         * preventing fully decoupled zero-gradient rows in the global Jacobian.
         */
        template<typename TProps>
        static void applyDerivativeFloor(
            const TProps& base_props,
            double reference_value,
            double& drho, double& dmu, double& dcp,
            double& dcv, double& dh,  double& dk)
        {
            const double denom = std::max(std::abs(reference_value), 1.0);
            const double floor_scale = 1.0e-12;
            auto floorMag = [&](double base_val) -> double {
                const double scaled = floor_scale * std::max(std::abs(base_val), 1.0) / denom;
                return std::max(scaled, 1.0e-18);
            };
            auto applyFloorOne = [&](double& dval, double base_val) {
                const double floor_mag = floorMag(base_val);
                if (!std::isfinite(dval)) { dval = floor_mag; return; }
                if (std::abs(dval) < floor_mag) {
                    const double sign = (dval < 0.0) ? -1.0 : 1.0;
                    dval = sign * floor_mag;
                }
            };
            applyFloorOne(drho, base_props.rho);
            applyFloorOne(dmu,  base_props.mu);
            applyFloorOne(dcp,  base_props.cp);
            applyFloorOne(dcv,  base_props.cv);
            applyFloorOne(dh,   base_props.h);
            applyFloorOne(dk,   base_props.k);
        }

        /**
         * @brief Exact full-safe ADVar assembly (prevents uninitialized gradient contamination)
         */
        template<int N>
        static void AssembleADVar(ADVar<N>& target, double val, double d_dP, double d_dT,
            const ADVar<N>& P, const ADVar<N>& T)
        {
            target = ADVar<N>(val);        // constructor ensures grad.setZero()
            target.grad += d_dP * P.grad;  // Eigen += accumulates correctly
            target.grad += d_dT * T.grad;
        }

        /**
         * @brief Check whether two fluid property evaluations are in the same thermodynamic phase
         */
        template<typename TProps>
        static bool isSamePhase(const TProps& p1, const TProps& p2)
        {
            double diff = std::abs(p1.rho - p2.rho);
            double avg  = 0.5 * (p1.rho + p2.rho);
            if (avg > 1e-6 && (diff / avg) > 0.01) return false;
            return true;
        }

        /**
         * @brief Robust numerical derivative computation (4-tier fallback cascade)
         * @return DiffMode indicating which tier succeeded
         */
        template<typename Func, typename TProps>
        static DiffMode computeRobustDerivatives(Func evalFunc, double val, double delta,
            const TProps& base_props,
            double& drho, double& dmu, double& dcp,
            double& dcv,  double& dh,  double& dk)
        {
            try {
                // Tier 1: Central difference
                auto prop_plus  = evalFunc(val + delta);
                auto prop_minus = evalFunc(val - delta);
                if (!isSamePhase(prop_plus, prop_minus)) throw PhaseMismatchException();

                double inv_2d = 1.0 / (2.0 * delta);
                drho = (prop_plus.rho - prop_minus.rho) * inv_2d;
                dmu  = (prop_plus.mu  - prop_minus.mu)  * inv_2d;
                dcp  = (prop_plus.cp  - prop_minus.cp)  * inv_2d;
                dcv  = (prop_plus.cv  - prop_minus.cv)  * inv_2d;
                dh   = (prop_plus.h   - prop_minus.h)   * inv_2d;
                dk   = (prop_plus.k   - prop_minus.k)   * inv_2d;
                return DiffMode::Central;
            }
            catch (...) {
                try {
                    // Tier 2: Forward difference
                    auto prop_plus = evalFunc(val + delta);
                    double inv_d = 1.0 / delta;
                    drho = (prop_plus.rho - base_props.rho) * inv_d;
                    dmu  = (prop_plus.mu  - base_props.mu)  * inv_d;
                    dcp  = (prop_plus.cp  - base_props.cp)  * inv_d;
                    dcv  = (prop_plus.cv  - base_props.cv)  * inv_d;
                    dh   = (prop_plus.h   - base_props.h)   * inv_d;
                    dk   = (prop_plus.k   - base_props.k)   * inv_d;
                    return DiffMode::Forward;
                }
                catch (...) {
                    try {
                        // Tier 3: Backward difference
                        auto prop_minus = evalFunc(val - delta);
                        double inv_d = 1.0 / delta;
                        drho = (base_props.rho - prop_minus.rho) * inv_d;
                        dmu  = (base_props.mu  - prop_minus.mu)  * inv_d;
                        dcp  = (base_props.cp  - prop_minus.cp)  * inv_d;
                        dcv  = (base_props.cv  - prop_minus.cv)  * inv_d;
                        dh   = (base_props.h   - prop_minus.h)   * inv_d;
                        dk   = (base_props.k   - prop_minus.k)   * inv_d;
                        return DiffMode::Backward;
                    }
                    catch (...) {
                        // Tier 4: derivative floor (prevents zero-Jacobian rows)
                        drho = dmu = dcp = dcv = dh = dk = 0.0;
                        applyDerivativeFloor(base_props, val, drho, dmu, dcp, dcv, dh, dk);
                        return DiffMode::ZeroGrad;
                    }
                }
            }
        }

        // =========================================================
        // CoolProp EOS backend (activated by #define USE_COOLPROP_EOS)
        // =========================================================
#ifdef USE_COOLPROP_EOS
    private:
        /**
         * @brief Lightweight POD for CoolProp query results.
         * Field names intentionally match WaterProperties/CO2Properties so that
         * computeRobustDerivatives and applyDerivativeFloor can accept this type
         * without any changes (both are templated on TProps).
         */
        struct CPProps {
            double rho = 0.0;
            double mu  = 0.0;
            double cp  = 0.0;
            double cv  = 0.0;
            double h   = 0.0;
            double k   = 0.0;
        };

        /**
         * @brief Per-thread CoolProp AbstractState for CO2 (HEOS = Span-Wagner EOS)
         */
        static CoolProp::AbstractState& co2_state_tl() {
            thread_local static std::unique_ptr<CoolProp::AbstractState> s;
            if (!s) {
                s.reset(CoolProp::AbstractState::factory("HEOS", "CO2"));
            }
            return *s;
        }

        /**
         * @brief Per-thread CoolProp AbstractState for Water (HEOS = IAPWS-95)
         */
        static CoolProp::AbstractState& water_state_tl() {
            thread_local static std::unique_ptr<CoolProp::AbstractState> s;
            if (!s) {
                s.reset(CoolProp::AbstractState::factory("HEOS", "Water"));
            }
            return *s;
        }

        /**
         * @brief Query CO2 properties at (P, T) using CoolProp.
         */
        static CPProps queryCO2_CP(double P, double T) {
            auto& st = co2_state_tl();
            st.update(CoolProp::PT_INPUTS, P, T);
            return { st.rhomass(), st.viscosity(), st.cpmass(),
                     st.cvmass(), st.hmass(), st.conductivity() };
        }

        /**
         * @brief Query Water properties at (P, T) using CoolProp.
         */
        static CPProps queryWater_CP(double P, double T) {
            auto& st = water_state_tl();
            st.update(CoolProp::PT_INPUTS, P, T);
            return { st.rhomass(), st.viscosity(), st.cpmass(),
                     st.cvmass(), st.hmass(), st.conductivity() };
        }

        /**
         * @brief Try to obtain analytical first partial derivatives from a CoolProp state
         *        that has already been updated to the target (P,T).
         * @return true if all derivatives were obtained successfully;
         *         false if the state is two-phase or if any derivative call threw.
         * @note Viscosity and conductivity derivatives are not reliable in the two-phase
         *       dome; CoolProp throws or returns NaN for those. We detect this and fall
         *       back to the numerical tier in the caller.
         */
        static bool getAnalyticalDerivs_CP(
            CoolProp::AbstractState& st,
            double& dRho_dP, double& dMu_dP, double& dCp_dP,
            double& dCv_dP,  double& dH_dP,  double& dK_dP,
            double& dRho_dT, double& dMu_dT, double& dCp_dT,
            double& dCv_dT,  double& dH_dT,  double& dK_dT)
        {
            try {
                using namespace CoolProp;
                if (st.phase() == iphase_twophase) return false;

                // Pressure derivatives (constant T)
                dRho_dP = st.first_partial_deriv(iDmass,       iP, iT);
                dMu_dP  = st.first_partial_deriv(iviscosity,   iP, iT);
                dCp_dP  = st.first_partial_deriv(iCpmass,      iP, iT);
                dCv_dP  = st.first_partial_deriv(iCvmass,      iP, iT);
                dH_dP   = st.first_partial_deriv(iHmass,       iP, iT);
                dK_dP   = st.first_partial_deriv(iconductivity,iP, iT);

                // Temperature derivatives (constant P)
                dRho_dT = st.first_partial_deriv(iDmass,       iT, iP);
                dMu_dT  = st.first_partial_deriv(iviscosity,   iT, iP);
                dCp_dT  = st.first_partial_deriv(iCpmass,      iT, iP);
                dCv_dT  = st.first_partial_deriv(iCvmass,      iT, iP);
                dH_dT   = st.first_partial_deriv(iHmass,       iT, iP);
                dK_dT   = st.first_partial_deriv(iconductivity,iT, iP);

                // Sanity check: reject non-finite derivatives
                auto finite6 = [](double a, double b, double c, double d, double e, double f) {
                    return std::isfinite(a) && std::isfinite(b) && std::isfinite(c)
                        && std::isfinite(d) && std::isfinite(e) && std::isfinite(f);
                };
                if (!finite6(dRho_dP, dMu_dP, dCp_dP, dCv_dP, dH_dP, dK_dP)) return false;
                if (!finite6(dRho_dT, dMu_dT, dCp_dT, dCv_dT, dH_dT, dK_dT)) return false;

                return true;
            }
            catch (...) {
                return false;
            }
        }

        /**
         * @brief CoolProp-backed CO2 evaluator.
         * Priority: (1) CoolProp analytical derivatives
         *           (2) numerical diff via CoolProp (two-phase / near-critical)
         *           (3) fallback constant values if CoolProp itself fails
         */
        template<int N>
        static ADFluidProperties<N> evaluateCO2_CoolProp(const ADVar<N>& P, const ADVar<N>& T)
        {
            const double p_val = P.val;
            const double t_val = T.val;

            ADFluidProperties<N> res;
            CPProps base;

            // Stage 1: base property query
            try {
                base = queryCO2_CP(p_val, t_val);
            }
            catch (...) {
                base = { 800.0, 1.48e-5, 1100.0, 850.0, 3.0e5, 0.03 };
                res.isFallback = true;
            }

            double dRho_dP=0, dMu_dP=0, dCp_dP=0, dCv_dP=0, dH_dP=0, dK_dP=0;
            double dRho_dT=0, dMu_dT=0, dCp_dT=0, dCv_dT=0, dH_dT=0, dK_dT=0;

            if (!res.isFallback) {
                // Stage 2: analytical derivatives
                bool ok = getAnalyticalDerivs_CP(co2_state_tl(),
                    dRho_dP, dMu_dP, dCp_dP, dCv_dP, dH_dP, dK_dP,
                    dRho_dT, dMu_dT, dCp_dT, dCv_dT, dH_dT, dK_dT);

                if (!ok) {
                    // Stage 3: numerical diff using CoolProp as oracle
                    res.near_bound = true;
                    const double dP = std::max(10.0, std::min(std::abs(p_val) * 1e-6, 1000.0));
                    const double dT = std::max(0.001, std::min(std::abs(t_val) * 1e-6, 0.1));
                    auto evalP = [&](double p) { return queryCO2_CP(p, t_val); };
                    auto evalT = [&](double t) { return queryCO2_CP(p_val, t); };
                    res.diff_mode_P = computeRobustDerivatives(evalP, p_val, dP, base,
                        dRho_dP, dMu_dP, dCp_dP, dCv_dP, dH_dP, dK_dP);
                    res.diff_mode_T = computeRobustDerivatives(evalT, t_val, dT, base,
                        dRho_dT, dMu_dT, dCp_dT, dCv_dT, dH_dT, dK_dT);
                }
            }

            AssembleADVar(res.rho, base.rho, dRho_dP, dRho_dT, P, T);
            AssembleADVar(res.mu,  base.mu,  dMu_dP,  dMu_dT,  P, T);
            AssembleADVar(res.cp,  base.cp,  dCp_dP,  dCp_dT,  P, T);
            AssembleADVar(res.cv,  base.cv,  dCv_dP,  dCv_dT,  P, T);
            AssembleADVar(res.h,   base.h,   dH_dP,   dH_dT,   P, T);
            AssembleADVar(res.k,   base.k,   dK_dP,   dK_dT,   P, T);
            return res;
        }

        /**
         * @brief CoolProp-backed Water evaluator (IAPWS-95).
         * Same priority cascade as evaluateCO2_CoolProp.
         */
        template<int N>
        static ADFluidProperties<N> evaluateWater_CoolProp(const ADVar<N>& P, const ADVar<N>& T)
        {
            const double p_val = P.val;
            const double t_val = T.val;

            ADFluidProperties<N> res;
            CPProps base;

            try {
                base = queryWater_CP(p_val, t_val);
            }
            catch (...) {
                base = { 1000.0, 1e-3, 4200.0, 4182.0, 1.0e5, 0.6 };
                res.isFallback = true;
            }

            double dRho_dP=0, dMu_dP=0, dCp_dP=0, dCv_dP=0, dH_dP=0, dK_dP=0;
            double dRho_dT=0, dMu_dT=0, dCp_dT=0, dCv_dT=0, dH_dT=0, dK_dT=0;

            if (!res.isFallback) {
                bool ok = getAnalyticalDerivs_CP(water_state_tl(),
                    dRho_dP, dMu_dP, dCp_dP, dCv_dP, dH_dP, dK_dP,
                    dRho_dT, dMu_dT, dCp_dT, dCv_dT, dH_dT, dK_dT);

                if (!ok) {
                    res.near_bound = true;
                    const double dP = std::max(10.0, std::min(std::abs(p_val) * 1e-6, 1000.0));
                    const double dT = std::max(0.001, std::min(std::abs(t_val) * 1e-6, 0.1));
                    auto evalP = [&](double p) { return queryWater_CP(p, t_val); };
                    auto evalT = [&](double t) { return queryWater_CP(p_val, t); };
                    res.diff_mode_P = computeRobustDerivatives(evalP, p_val, dP, base,
                        dRho_dP, dMu_dP, dCp_dP, dCv_dP, dH_dP, dK_dP);
                    res.diff_mode_T = computeRobustDerivatives(evalT, t_val, dT, base,
                        dRho_dT, dMu_dT, dCp_dT, dCv_dT, dH_dT, dK_dT);
                }
            }

            AssembleADVar(res.rho, base.rho, dRho_dP, dRho_dT, P, T);
            AssembleADVar(res.mu,  base.mu,  dMu_dP,  dMu_dT,  P, T);
            AssembleADVar(res.cp,  base.cp,  dCp_dP,  dCp_dT,  P, T);
            AssembleADVar(res.cv,  base.cv,  dCv_dP,  dCv_dT,  P, T);
            AssembleADVar(res.h,   base.h,   dH_dP,   dH_dT,   P, T);
            AssembleADVar(res.k,   base.k,   dK_dP,   dK_dT,   P, T);
            return res;
        }
#endif  // USE_COOLPROP_EOS
    };  // class Evaluator

}  // namespace AD_Fluid

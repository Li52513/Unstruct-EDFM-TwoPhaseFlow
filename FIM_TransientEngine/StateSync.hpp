#pragma once

#include "Types.hpp"
#include "../FIM_StateMap.h"
#include "../AD_FluidEvaluator.h"
#include "../SolverContrlStrName_op.h"

#include <algorithm>
#include <filesystem>
#include <string>

#ifdef _WIN32
#include <direct.h>
#define MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define MKDIR(path) mkdir(path, 0777)
#endif

namespace FIM_Engine {

    inline int MatrixBlockCount(const MeshManager& mgr) { return mgr.getMatrixDOFCount(); }
    inline int MatrixBlockCount(const MeshManager_3D& mgr) { return mgr.fracture_network().getSolverIndexOffset(); }

    /**
     * @brief Clamp scalar water saturation into constitutive-safe interval.
     */
    inline double ClampSwForConstitutive(
        double sw,
        const CapRelPerm::VGParams& vg_params,
        double safety_eps = 1.0e-8)
    {
        const double eps = std::max(1.0e-12, std::min(1.0e-2, safety_eps));
        const double lower = std::max(0.0, std::min(1.0, vg_params.Swr + eps));
        const double upper_raw = std::max(0.0, std::min(1.0, 1.0 - vg_params.Sgr - eps));
        const double upper = std::max(lower, upper_raw);
        if (sw < lower) return lower;
        if (sw > upper) return upper;
        return sw;
    }

    /**
     * @brief Clamp AD water saturation into constitutive-safe interval.
     */
    template<int N>
    inline ADVar<N> ClampSwForConstitutive(
        const ADVar<N>& sw,
        const CapRelPerm::VGParams& vg_params,
        double safety_eps = 1.0e-8)
    {
        const double sw_clamped = ClampSwForConstitutive(sw.val, vg_params, safety_eps);
        if (sw_clamped == sw.val) {
            return sw;
        }
        return ADVar<N>(sw_clamped);
    }

    template<int N>
    inline AD_Fluid::ADFluidProperties<N> BuildConstantFluidProps(
        const FluidConstantProperties& cfg,
        const ADVar<N>& T,
        bool no_convection)
    {
        AD_Fluid::ADFluidProperties<N> props;
        props.rho = ADVar<N>(cfg.rho);
        props.mu = ADVar<N>(cfg.mu);
        if (no_convection) {
            props.cp = ADVar<N>(0.0);
            props.cv = ADVar<N>(0.0);
            props.h = ADVar<N>(0.0);
        }
        else {
            props.cp = ADVar<N>(cfg.cp);
            props.cv = ADVar<N>(cfg.cv);
            props.h = ADVar<N>(cfg.cp) * T;
        }
        props.k = ADVar<N>(cfg.k);
        props.isFallback = false;
        props.near_bound = false;
        return props;
    }

    struct FluidPropertyEvalContext {
        SinglePhaseFluidModel single_phase_model = SinglePhaseFluidModel::Water;
        FluidPropertyEvalConfig config = FluidPropertyEvalConfig();
    };

    inline FluidPropertyEvalConfig BuildEffectiveFluidPropertyEvalConfig(
        SinglePhaseFluidModel model,
        const FluidPropertyEvalConfig& requested)
    {
        FluidPropertyEvalConfig cfg = requested;
        if (model == SinglePhaseFluidModel::ConstantWater ||
            model == SinglePhaseFluidModel::ConstantWaterNoConvection) {
            cfg.enable_single_phase_constant = true;
            cfg.single_phase_is_co2 = false;
            cfg.single_phase_no_convection = (model == SinglePhaseFluidModel::ConstantWaterNoConvection);
        }
        else if (!cfg.enable_single_phase_constant &&
                 !cfg.enable_two_phase_constant &&
                 model == SinglePhaseFluidModel::CO2) {
            // Backward-compatible: N=2 single-phase CO2 should remain CO2 EOS by default.
            cfg.single_phase_is_co2 = true;
        }
        return cfg;
    }

    inline FluidPropertyEvalContext BuildFluidPropertyEvalContext(
        SinglePhaseFluidModel model,
        const FluidPropertyEvalConfig& requested = FluidPropertyEvalConfig())
    {
        FluidPropertyEvalContext ctx;
        ctx.single_phase_model = model;
        ctx.config = BuildEffectiveFluidPropertyEvalConfig(model, requested);
        return ctx;
    }

    inline FluidPropertyEvalContext BuildFluidPropertyEvalContext(
        const FluidPropertyEvalConfig& requested)
    {
        FluidPropertyEvalContext ctx;
        ctx.single_phase_model = requested.single_phase_is_co2
            ? SinglePhaseFluidModel::CO2
            : SinglePhaseFluidModel::Water;
        ctx.config = requested;
        return ctx;
    }

    inline bool IsSinglePhaseNoConvectionActive(const FluidPropertyEvalContext& ctx) {
        return ctx.config.enable_single_phase_constant && ctx.config.single_phase_no_convection;
    }

    template<int N>
    inline AD_Fluid::ADFluidProperties<N> EvalSinglePhaseFluid(
        const FluidPropertyEvalContext& ctx,
        const ADVar<N>& P,
        const ADVar<N>& T)
    {
        if (ctx.config.enable_single_phase_constant) {
            const auto& c = ctx.config.single_phase_is_co2 ? ctx.config.gas : ctx.config.water;
            return BuildConstantFluidProps<N>(c, T, ctx.config.single_phase_no_convection);
        }
        if (ctx.config.single_phase_is_co2 || ctx.single_phase_model == SinglePhaseFluidModel::CO2) {
            return AD_Fluid::Evaluator::evaluateCO2<N>(P, T);
        }
        return AD_Fluid::Evaluator::evaluateWater<N>(P, T);
    }

    template<int N>
    inline AD_Fluid::ADFluidProperties<N> EvalTwoPhaseWaterFluid(
        const FluidPropertyEvalContext& ctx,
        const ADVar<N>& P,
        const ADVar<N>& T)
    {
        if (ctx.config.enable_two_phase_constant) {
            return BuildConstantFluidProps<N>(ctx.config.water, T, false);
        }
        return AD_Fluid::Evaluator::evaluateWater<N>(P, T);
    }

    template<int N>
    inline AD_Fluid::ADFluidProperties<N> EvalTwoPhaseGasFluid(
        const FluidPropertyEvalContext& ctx,
        const ADVar<N>& P,
        const ADVar<N>& T)
    {
        if (ctx.config.enable_two_phase_constant) {
            return BuildConstantFluidProps<N>(ctx.config.gas, T, false);
        }
        return AD_Fluid::Evaluator::evaluateCO2<N>(P, T);
    }

    // Legacy helper retained for compatibility with old call sites.
    template<int N>
    inline AD_Fluid::ADFluidProperties<N> EvalPrimaryFluid(
        SinglePhaseFluidModel model,
        const ADVar<N>& P,
        const ADVar<N>& T)
    {
        return EvalSinglePhaseFluid<N>(BuildFluidPropertyEvalContext(model), P, T);
    }

    template<int N>
    inline AD_Fluid::ADFluidProperties<N> EvalPrimaryFluid(
        SinglePhaseFluidModel model,
        const FluidPropertyEvalConfig& fluid_cfg,
        const ADVar<N>& P,
        const ADVar<N>& T)
    {
        return EvalSinglePhaseFluid<N>(BuildFluidPropertyEvalContext(model, fluid_cfg), P, T);
    }

    inline std::string NormalizeOutputRoot(const std::string& output_root) {
        return output_root.empty() ? TransientSolverParams().output_root_dir : output_root;
    }

    inline std::filesystem::path BuildCaseOutputDirPath(
        const std::string& output_root,
        const std::string& caseName)
    {
        return std::filesystem::path(NormalizeOutputRoot(output_root)) / caseName;
    }

    inline std::string BuildCaseOutputDir(
        const std::string& output_root,
        const std::string& caseName)
    {
        return BuildCaseOutputDirPath(output_root, caseName).generic_string();
    }

    inline std::string BuildCaseOutputFile(
        const std::string& output_root,
        const std::string& caseName,
        const std::string& relative_name)
    {
        return (BuildCaseOutputDirPath(output_root, caseName) / relative_name).generic_string();
    }

    inline void MakePath(const std::string& output_root, const std::string& caseName) {
        std::filesystem::create_directories(BuildCaseOutputDirPath(output_root, caseName));
    }

    template <typename FieldMgrType, typename MeshMgrType, int N>
    inline void SyncStateToFieldManager(
        const FIM_StateMap<N>& state,
        FieldMgrType& fm,
        const MeshMgrType& mgr,
        SinglePhaseFluidModel sp_model = SinglePhaseFluidModel::Water,
        const FluidPropertyEvalConfig& fluid_cfg = FluidPropertyEvalConfig(),
        const CapRelPerm::VGParams& vg_params = CapRelPerm::VGParams(),
        const CapRelPerm::RelPermParams& rp_params = CapRelPerm::RelPermParams()) {
        int nMat = MatrixBlockCount(mgr);
        int nTotal = mgr.getTotalDOFCount();
        const FluidPropertyEvalContext fluid_ctx = BuildFluidPropertyEvalContext(sp_model, fluid_cfg);

        const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        const auto satCfg = PhysicalProperties_string_op::SaturationEquation_String::FIM();
        const PhysicalProperties_string_op::Water wCfg;
        const PhysicalProperties_string_op::CO2 gCfg;
        const PhysicalProperties_string_op::TransientInternalFieldNames internalNames;

        auto f_pw = fm.getOrCreateMatrixScalar(pCfg.pressure_field, 0.0);
        auto f_T = fm.getOrCreateMatrixScalar(tCfg.temperatue_field, 0.0);
        auto f_rhow = fm.getOrCreateMatrixScalar(wCfg.rho_tag, 0.0);
        auto f_hw = fm.getOrCreateMatrixScalar(wCfg.h_tag, 0.0);
        auto f_lamw_mob = fm.getOrCreateMatrixScalar(wCfg.lambda_w_tag, 0.0);

        auto frac_pw = fm.getOrCreateFractureScalar(pCfg.pressure_field, 0.0);
        auto frac_T = fm.getOrCreateFractureScalar(tCfg.temperatue_field, 0.0);
        auto frac_rhow = fm.getOrCreateFractureScalar(wCfg.rho_tag, 0.0);
        auto frac_hw = fm.getOrCreateFractureScalar(wCfg.h_tag, 0.0);
        auto frac_lamw_mob = fm.getOrCreateFractureScalar(wCfg.lambda_w_tag, 0.0);

        auto f_P_viz = fm.getOrCreateMatrixScalar(internalNames.pressure_viz, 0.0);
        auto frac_P_viz = fm.getOrCreateFractureScalar(internalNames.pressure_viz, 0.0);

        std::shared_ptr<volScalarField> f_sw, f_Sw_viz, f_rhog, f_hg, f_lamg_mob;
        std::shared_ptr<volScalarField> frac_sw, frac_Sw_viz, frac_rhog, frac_hg, frac_lamg_mob;

        if constexpr (N == 3) {
            f_sw = fm.getOrCreateMatrixScalar(satCfg.saturation, 0.0);
            f_Sw_viz = fm.getOrCreateMatrixScalar(internalNames.saturation_viz, 0.0);
            f_rhog = fm.getOrCreateMatrixScalar(gCfg.rho_tag, 0.0);
            f_hg = fm.getOrCreateMatrixScalar(gCfg.h_tag, 0.0);
            f_lamg_mob = fm.getOrCreateMatrixScalar(gCfg.lambda_g_tag, 0.0);

            frac_sw = fm.getOrCreateFractureScalar(satCfg.saturation, 0.0);
            frac_Sw_viz = fm.getOrCreateFractureScalar(internalNames.saturation_viz, 0.0);
            frac_rhog = fm.getOrCreateFractureScalar(gCfg.rho_tag, 0.0);
            frac_hg = fm.getOrCreateFractureScalar(gCfg.h_tag, 0.0);
            frac_lamg_mob = fm.getOrCreateFractureScalar(gCfg.lambda_g_tag, 0.0);
        }

        for (int i = 0; i < nTotal; ++i) {
            double p = state.P[i];
            double t = state.T[i];
            ADVar<N> P_ad(p), T_ad(t);
            AD_Fluid::ADFluidProperties<N> propsW{};
            if constexpr (N == 3) {
                propsW = EvalTwoPhaseWaterFluid<N>(fluid_ctx, P_ad, T_ad);
            }
            else {
                propsW = EvalSinglePhaseFluid<N>(fluid_ctx, P_ad, T_ad);
            }

            double sw = (N == 3) ? state.Sw[i] : 1.0;
            double sw_constitutive = sw;
            double rho_w = propsW.rho.val;
            double mu_w = propsW.mu.val;
            double krw = 1.0, krg = 0.0;
            AD_Fluid::ADFluidProperties<N> propsG{};

            if constexpr (N == 3) {
                sw_constitutive = ClampSwForConstitutive(sw, vg_params);
                const ADVar<N> Sw_const_ad(sw_constitutive);
                if (i >= nMat) {
                    krw = sw_constitutive;
                    krg = 1.0 - sw_constitutive;
                    propsG = EvalTwoPhaseGasFluid<N>(fluid_ctx, P_ad, T_ad);  // Pc=0
                }
                else {
                    const ADVar<N> Pc_const_ad = CapRelPerm::pc_vG<N>(Sw_const_ad, vg_params);
                    const ADVar<N> P_gas_ad = P_ad + Pc_const_ad;
                    propsG = EvalTwoPhaseGasFluid<N>(fluid_ctx, P_gas_ad, T_ad);
                    ADVar<N> krw_ad, krg_ad;
                    CapRelPerm::kr_Mualem_vG<N>(Sw_const_ad, vg_params, rp_params, krw_ad, krg_ad);
                    krw = krw_ad.val;
                    krg = krg_ad.val;
                }
            }
            double lambda_w_mob = krw / std::max(mu_w, 1e-18);

            if (i < nMat) {
                (*f_pw)[i] = p; (*f_T)[i] = t; (*f_rhow)[i] = rho_w; (*f_hw)[i] = propsW.h.val; (*f_lamw_mob)[i] = lambda_w_mob; (*f_P_viz)[i] = p;
                if constexpr (N == 3) {
                    (*f_sw)[i] = sw_constitutive; (*f_Sw_viz)[i] = sw_constitutive; (*f_rhog)[i] = propsG.rho.val; (*f_hg)[i] = propsG.h.val; (*f_lamg_mob)[i] = krg / std::max(propsG.mu.val, 1e-18);
                }
            }
            else {
                int fi = i - nMat;
                (*frac_pw)[fi] = p; (*frac_T)[fi] = t; (*frac_rhow)[fi] = rho_w; (*frac_hw)[fi] = propsW.h.val; (*frac_lamw_mob)[fi] = lambda_w_mob; (*frac_P_viz)[fi] = p;
                if constexpr (N == 3) {
                    (*frac_sw)[fi] = sw_constitutive; (*frac_Sw_viz)[fi] = sw_constitutive; (*frac_rhog)[fi] = propsG.rho.val; (*frac_hg)[fi] = propsG.h.val; (*frac_lamg_mob)[fi] = krg / std::max(propsG.mu.val, 1e-18);
                }
            }
        }
    }

    template <typename FieldMgrType>
    inline void InjectStaticProperties(FieldMgrType& fm) {
        const PhysicalProperties_string_op::Rock rock;
        const PhysicalProperties_string_op::Fracture_string frac;
        const PhysicalProperties_string_op::Water water;

        fm.getOrCreateMatrixScalar(rock.k_xx_tag, 1.0e-13);
        fm.getOrCreateMatrixScalar(rock.k_yy_tag, 1.0e-13);
        fm.getOrCreateMatrixScalar(rock.k_zz_tag, 1.0e-13);
        fm.getOrCreateMatrixScalar(rock.lambda_tag, 2.0);
        fm.getOrCreateMatrixScalar(rock.phi_tag, 0.2);

        fm.getOrCreateMatrixScalar(rock.rho_tag, 2600.0);
        fm.getOrCreateMatrixScalar(rock.cp_tag, 1000.0);
        fm.getOrCreateMatrixScalar(rock.c_r_tag, 0.0);

        fm.getOrCreateFractureScalar(frac.k_t_tag, 1.0e-11);
        fm.getOrCreateFractureScalar(frac.k_n_tag, 1.0e-12);
        fm.getOrCreateFractureScalar(frac.aperture_tag, 1.0e-3);
        fm.getOrCreateFractureScalar(frac.lambda_tag, 2.0);
        fm.getOrCreateFractureScalar(frac.phi_tag, 0.2);

        fm.getOrCreateFractureScalar(frac.rho_tag, 2600.0);
        fm.getOrCreateFractureScalar(frac.cp_tag, 1000.0);
        fm.getOrCreateFractureScalar(frac.c_r_tag, 0.0);

        fm.getOrCreateFractureScalar(water.k_tag, 0.6);
    }

} // namespace FIM_Engine

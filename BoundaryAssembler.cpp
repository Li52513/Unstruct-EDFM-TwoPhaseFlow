/**
 * @file BoundaryAssembler.cpp
 * @brief Boundary and leakoff assembly facade implementation
 */
#include "BoundaryAssembler.h"
#include "FVM_Ops_AD.h"
#include "Well_WellControlTypes.h"
#include "SolverContrlStrName_op.h"
#include "ADVar.hpp"
#include "AD_FluidEvaluator.h"
#include "CapRelPerm_HD_AD.h"
#include "FIM_TransientEngine/StateSync.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <initializer_list>
#include <limits>
#include <stdexcept>

namespace {
    constexpr double kBoundaryCoeffEps = 1.0e-14;
    constexpr double kEps = 1.0e-12; // denominator floor

    struct WellFieldTags {
        std::string pressure;
        std::string k_xx;
        std::string k_yy;
        std::string k_zz;
        std::string mob_density_w;
        std::string mob_density_g;
        std::string rho_w;
        std::string rho_g;
        std::string lambda_w_mob;
        std::string lambda_g_mob;
        std::string h_w;
        std::string h_g;
        std::string temperature;
        std::string saturation;
    };

    inline WellFieldTags BuildWellFieldTags() {
        WellFieldTags tags;
        const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        const auto sCfg = PhysicalProperties_string_op::SaturationEquation_String::FIM();
        const PhysicalProperties_string_op::Rock rock;
        const PhysicalProperties_string_op::Water water;
        const PhysicalProperties_string_op::CO2 gas;
        const PhysicalProperties_string_op::WellAuxiliaryFields aux;
        tags.pressure = pCfg.pressure_field;
        tags.k_xx = rock.k_xx_tag;
        tags.k_yy = rock.k_yy_tag;
        tags.k_zz = rock.k_zz_tag;
        tags.mob_density_w = aux.mob_density_w;
        tags.mob_density_g = aux.mob_density_g;
        tags.rho_w = water.rho_tag;
        tags.rho_g = gas.rho_tag;
        tags.lambda_w_mob = water.lambda_w_tag;
        tags.lambda_g_mob = gas.lambda_g_tag;
        tags.h_w = water.h_tag;
        tags.h_g = gas.h_tag;
        tags.temperature = tCfg.temperatue_field;
        tags.saturation = sCfg.saturation;
        return tags;
    }

    // Helper utilities
    inline bool isValidEqIdx(int eqIdx, const std::vector<double>& residual, const std::vector<double>& jacobianDiag) {
        return eqIdx >= 0 && eqIdx < static_cast<int>(residual.size()) && eqIdx < static_cast<int>(jacobianDiag.size());
    }

    inline bool isValidLocalIdx(int localIdx, const std::shared_ptr<volScalarField>& field) {
        return field && localIdx >= 0 && localIdx < static_cast<int>(field->data.size());
    }

    inline double safeGetFieldValue(const std::shared_ptr<volScalarField>& field, int localIdx, double default_val = 0.0) {
        if (isValidLocalIdx(localIdx, field)) return field->data[localIdx];
        return default_val;
    }

    inline int get3DMatrixOffset(const MeshManager_3D& mgr) {
        return mgr.fracture_network().getSolverIndexOffset();
    }

    inline void ResolveIndex(int input_id, WellTargetDomain domain, int n_matrix_dofs, int& solverIdx, int& localIdx) {
        if (domain == WellTargetDomain::Matrix) {
            localIdx = input_id;
            solverIdx = input_id;
        }
        else {
            if (input_id >= n_matrix_dofs) {
                solverIdx = input_id;
                localIdx = input_id - n_matrix_dofs;
            }
            else {
                localIdx = input_id;
                solverIdx = input_id + n_matrix_dofs;
            }
        }
    }

    inline void ResolveIndexWithSpace(
        int input_id,
        CompletionIdSpace id_space,
        WellTargetDomain domain,
        int n_matrix_dofs,
        int& solverIdx,
        int& localIdx)
    {
        if (id_space == CompletionIdSpace::SolverIndex) {
            solverIdx = input_id;
            localIdx = (domain == WellTargetDomain::Matrix) ? input_id : (input_id - n_matrix_dofs);
            return;
        }

        if (id_space == CompletionIdSpace::FractureLocalIndex) {
            if (domain == WellTargetDomain::Matrix) {
                solverIdx = input_id;
                localIdx = input_id;
            }
            else {
                localIdx = input_id;
                solverIdx = input_id + n_matrix_dofs;
            }
            return;
        }

        ResolveIndex(input_id, domain, n_matrix_dofs, solverIdx, localIdx);
    }

    inline void ResolveStepCompletionIndex(
        const WellScheduleStep& step,
        int n_matrix_dofs,
        int& solverIdx,
        int& localIdx)
    {
        if (step.completion_solver_index >= 0) {
            solverIdx = step.completion_solver_index;
            localIdx = (step.domain == WellTargetDomain::Matrix)
                ? solverIdx
                : (solverIdx - n_matrix_dofs);
            return;
        }

        ResolveIndexWithSpace(
            step.completion_id,
            step.completion_id_space,
            step.domain,
            n_matrix_dofs,
            solverIdx,
            localIdx);
    }

    inline double ResolveWellBhpValue(
        const WellScheduleStep& step,
        double p_cell,
        const std::unordered_map<std::string, double>* well_bhp_by_name)
    {
        if (well_bhp_by_name) {
            auto it = well_bhp_by_name->find(step.well_name);
            if (it != well_bhp_by_name->end() && std::isfinite(it->second)) {
                return it->second;
            }
        }
        return (step.control_mode == WellControlMode::BHP) ? step.target_value : p_cell;
    }

    /**
     * @brief Create a constant AD variable with zero gradient.
     * @tparam N AD dimension.
     * @param v Scalar value.
     * @return Constant AD variable.
     */
    template<int N>
    inline ADVar<N> MakeConstAD(double v) {
        ADVar<N> x;
        x.val = v;
        for (int k = 0; k < N; ++k) x.grad(k) = 0.0;
        return x;
    }

    inline FIM_Engine::FluidPropertyEvalContext BuildFluidEvalContext(
        const FluidPropertyEvalConfig& fluid_cfg)
    {
        return FIM_Engine::BuildFluidPropertyEvalContext(fluid_cfg);
    }

    struct PhaseSplit {
        double fw = 0.0;
        double fg = 0.0;
    };

    /**
     * @brief Check whether user provided explicit phase fractions.
     * @param step Well schedule step.
     * @param[out] split User-normalized split when available.
     * @return True when explicit fractions are provided.
     */
    inline bool TryResolveUserSpecifiedPhaseSplit(const WellScheduleStep& step, PhaseSplit& split) {
        const double fw_user = std::max(0.0, step.frac_w);
        const double fg_user = std::max(0.0, step.frac_g);
        const double sum_user = fw_user + fg_user;
        if (sum_user <= kEps) {
            return false;
        }
        split.fw = fw_user / sum_user;
        split.fg = fg_user / sum_user;
        return true;
    }

    // Total-mode split priority:
    // 1) explicit user fractions,
    // 2) local mobility-proportional split,
    // 3) fallback to water-only (1,0) when mobility is unavailable.
    inline PhaseSplit ResolveTotalPhaseSplit(const WellScheduleStep& step, double mob_w, double mob_g) {
        PhaseSplit split;
        if (TryResolveUserSpecifiedPhaseSplit(step, split)) {
            return split;
        }

        const double mw = std::max(0.0, mob_w);
        const double mg = std::max(0.0, mob_g);
        const double sum_mob = mw + mg;
        if (sum_mob > kEps) {
            split.fw = mw / sum_mob;
            split.fg = mg / sum_mob;
            return split;
        }

        split.fw = 1.0;
        split.fg = 0.0;
        return split;
    }

    /**
     * @brief Resolve injection composition split for Total mode.
     * @details
     * Priority: explicit user fractions -> component_mode selector -> injection_is_co2 hint.
     * Falls back to water-only when no explicit information is provided.
     * @param step Well schedule step.
     * @param has_w_eq Whether water equation is available.
     * @param has_g_eq Whether gas equation is available.
     * @return Injection phase split.
     */
    inline PhaseSplit ResolveInjectionPhaseSplit(const WellScheduleStep& step, bool has_w_eq, bool has_g_eq) {
        if (has_w_eq && !has_g_eq) return PhaseSplit{ 1.0, 0.0 };
        if (!has_w_eq && has_g_eq) return PhaseSplit{ 0.0, 1.0 };
        if (!has_w_eq && !has_g_eq) return PhaseSplit{ 0.0, 0.0 };

        PhaseSplit user_split;
        if (TryResolveUserSpecifiedPhaseSplit(step, user_split)) {
            return user_split;
        }

        if (step.component_mode == WellComponentMode::Water) return PhaseSplit{ 1.0, 0.0 };
        if (step.component_mode == WellComponentMode::Gas) return PhaseSplit{ 0.0, 1.0 };

        return step.injection_is_co2 ? PhaseSplit{ 0.0, 1.0 } : PhaseSplit{ 1.0, 0.0 };
    }

    inline PhaseSplit ResolvePhaseSplit(const WellScheduleStep& step, double mob_w, double mob_g, bool has_w_eq, bool has_g_eq) {
        // Single-equation mode: map all well mass to the available mass row.
        if (has_w_eq && !has_g_eq) return PhaseSplit{ 1.0, 0.0 };
        if (!has_w_eq && has_g_eq) return PhaseSplit{ 0.0, 1.0 };
        if (!has_w_eq && !has_g_eq) return PhaseSplit{ 0.0, 0.0 };

        if (step.component_mode == WellComponentMode::Water) return PhaseSplit{ 1.0, 0.0 };
        if (step.component_mode == WellComponentMode::Gas) return PhaseSplit{ 0.0, 1.0 };
        return ResolveTotalPhaseSplit(step, mob_w, mob_g);
    }

    /**
     * @brief Clamp saturation to constitutive-safe interval for kr/Pc evaluation.
     * @tparam N AD dimension.
     * @param sw Water saturation AD variable.
     * @param vg van Genuchten parameters.
     * @return Clamped saturation preserving interior gradients.
     */
    template<int N>
    inline ADVar<N> ClampSwForConstitutiveWell(const ADVar<N>& sw, const CapRelPerm::VGParams& vg) {
        const double eps = 1.0e-8;
        const double lower = std::max(0.0, std::min(1.0, vg.Swr + eps));
        const double upper_raw = std::max(0.0, std::min(1.0, 1.0 - vg.Sgr - eps));
        const double upper = std::max(lower, upper_raw);
        if (sw.val < lower) return MakeConstAD<N>(lower);
        if (sw.val > upper) return MakeConstAD<N>(upper);
        return sw;
    }

    /**
     * @brief Standard-condition densities used for optional std-volume rate conversion.
     */
    struct StandardConditionDensities {
        double rho_w = 1.0;
        double rho_g = 1.0;
        bool enabled = false;
    };

    /**
     * @brief Evaluate standard-condition densities for water and CO2.
     * @param step Well schedule step.
     * @return Standard-condition density bundle.
     * @throws std::invalid_argument Invalid standard-state pressure/temperature.
     * @throws std::runtime_error Failed or nonphysical density evaluation.
     */
    inline StandardConditionDensities EvaluateStandardConditionDensities(
        const WellScheduleStep& step,
        const FIM_Engine::FluidPropertyEvalContext& fluid_ctx)
    {
        StandardConditionDensities out;
        if (step.rate_target_type != WellRateTargetType::StdVolumeRate) {
            return out;
        }

        if (fluid_ctx.config.enable_two_phase_constant) {
            out.rho_w = fluid_ctx.config.water.rho;
            out.rho_g = fluid_ctx.config.gas.rho;
            if (!std::isfinite(out.rho_w) || !std::isfinite(out.rho_g) ||
                out.rho_w <= 0.0 || out.rho_g <= 0.0) {
                throw std::runtime_error("Invalid two-phase constant standard-condition densities.");
            }
            out.enabled = true;
            return out;
        }

        if (fluid_ctx.config.enable_single_phase_constant) {
            out.rho_w = fluid_ctx.config.single_phase_is_co2
                ? fluid_ctx.config.gas.rho
                : fluid_ctx.config.water.rho;
            out.rho_g = fluid_ctx.config.gas.rho;
            if (!std::isfinite(out.rho_w) || !std::isfinite(out.rho_g) ||
                out.rho_w <= 0.0 || out.rho_g <= 0.0) {
                throw std::runtime_error("Invalid single-phase constant standard-condition densities.");
            }
            out.enabled = true;
            return out;
        }

        if (!(step.std_pressure > 0.0) || !(step.std_temperature > 0.0)) {
            throw std::invalid_argument("Standard-condition conversion requires std_pressure>0 and std_temperature>0.");
        }

        ADVar<1> P_std(step.std_pressure);
        ADVar<1> T_std(step.std_temperature);
        const auto props_w_std = FIM_Engine::EvalTwoPhaseWaterFluid<1>(fluid_ctx, P_std, T_std);
        const auto props_g_std = FIM_Engine::EvalTwoPhaseGasFluid<1>(fluid_ctx, P_std, T_std);

        if (!std::isfinite(props_w_std.rho.val) || !std::isfinite(props_g_std.rho.val) ||
            props_w_std.rho.val <= 0.0 || props_g_std.rho.val <= 0.0) {
            throw std::runtime_error("Failed to evaluate positive finite standard-condition densities.");
        }

        out.rho_w = props_w_std.rho.val;
        out.rho_g = props_g_std.rho.val;
        out.enabled = true;
        return out;
    }

    /**
     * @brief Build a constant phase mass source from target split.
     * @tparam N AD dimension.
     * @param step Well schedule step.
     * @param split_const Constant split in [0,1].
     * @param std_dens Standard-condition densities (used only for std-volume target).
     * @param water_phase True for water phase, false for gas phase.
     * @return Phase mass source AD variable.
     */
    template<int N>
    inline ADVar<N> BuildMassRateBySplitConst(
        const WellScheduleStep& step,
        double split_const,
        const StandardConditionDensities& std_dens,
        bool water_phase)
    {
        const double split = std::max(0.0, split_const);
        const double rho_std = water_phase ? std_dens.rho_w : std_dens.rho_g;
        const double q_mass = (step.rate_target_type == WellRateTargetType::StdVolumeRate)
            ? (step.target_value * split * rho_std)
            : (step.target_value * split);
        return MakeConstAD<N>(q_mass);
    }

    /**
     * @brief Build an AD phase mass source from AD split.
     * @tparam N AD dimension.
     * @param step Well schedule step.
     * @param split_ad AD split variable.
     * @param std_dens Standard-condition densities (used only for std-volume target).
     * @param water_phase True for water phase, false for gas phase.
     * @return Phase mass source AD variable.
     */
    template<int N>
    inline ADVar<N> BuildMassRateBySplitAD(
        const WellScheduleStep& step,
        const ADVar<N>& split_ad,
        const StandardConditionDensities& std_dens,
        bool water_phase)
    {
        if (step.rate_target_type == WellRateTargetType::StdVolumeRate) {
            const double rho_std = water_phase ? std_dens.rho_w : std_dens.rho_g;
            return MakeConstAD<N>(step.target_value * rho_std) * split_ad;
        }
        return MakeConstAD<N>(step.target_value) * split_ad;
    }
    inline const char* ToBoundaryTypeLabel(BoundarySetting::BoundaryType type) {
        switch (type) {
        case BoundarySetting::BoundaryType::Dirichlet: return "Dirichlet";
        case BoundarySetting::BoundaryType::Neumann:   return "Neumann";
        case BoundarySetting::BoundaryType::Robin:     return "Robin";
        default:                                       return "Unknown";
        }
    }
    inline BoundaryTagStats& TouchBoundaryTagStats(
        BoundaryAssemblyStats& stats,
        int physicalTag,
        BoundarySetting::BoundaryType type)
    {
        const std::string key = std::to_string(physicalTag) + "_" + ToBoundaryTypeLabel(type);
        return stats.perTagType[key];
    }

    struct BoundaryFieldTags {
        std::string pressure;
        std::string temperature;
        std::string saturation;
        std::array<std::string, 4> pressure_aliases;
        std::array<std::string, 3> temperature_aliases;
        std::array<std::string, 4> saturation_aliases;
        std::string k_xx;
        std::string k_yy;
        std::string k_zz;
        std::string phi;
        std::string lambda_rock;
        std::string lambda_eff;
        std::string fr_k_t;
        std::string fr_phi;
        std::string fr_lambda;
    };

    inline BoundaryFieldTags BuildBoundaryFieldTags() {
        BoundaryFieldTags tags;
        const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        const auto sCfg = PhysicalProperties_string_op::SaturationEquation_String::FIM();
        const PhysicalProperties_string_op::Rock rock;
        const PhysicalProperties_string_op::Fracture_string frac;
        const PhysicalProperties_string_op::EffectiveProps eff;
        const PhysicalProperties_string_op::LegacyFieldAliasPolicy aliases;
        tags.pressure = pCfg.pressure_field;
        tags.temperature = tCfg.temperatue_field;
        tags.saturation = sCfg.saturation;
        tags.pressure_aliases = aliases.pressure_aliases;
        tags.temperature_aliases = aliases.temperature_aliases;
        tags.saturation_aliases = aliases.saturation_aliases;
        tags.k_xx = rock.k_xx_tag;
        tags.k_yy = rock.k_yy_tag;
        tags.k_zz = rock.k_zz_tag;
        tags.phi = rock.phi_tag;
        tags.lambda_rock = rock.lambda_tag;
        tags.lambda_eff = eff.lambda_eff_tag;
        tags.fr_k_t = frac.k_t_tag;
        tags.fr_phi = frac.phi_tag;
        tags.fr_lambda = frac.lambda_tag;
        return tags;
    }

    struct DomainBounds {
        double xmin = 0.0;
        double xmax = 0.0;
        double ymin = 0.0;
        double ymax = 0.0;
        double zmin = 0.0;
        double zmax = 0.0;
    };

    inline DomainBounds ComputeDomainBounds(const Mesh& mesh) {
        DomainBounds b;
        b.xmin = b.ymin = b.zmin = std::numeric_limits<double>::infinity();
        b.xmax = b.ymax = b.zmax = -std::numeric_limits<double>::infinity();
        for (const auto& kv : mesh.getNodesMap()) {
            const auto& c = kv.second.coord;
            b.xmin = std::min(b.xmin, c.m_x);
            b.xmax = std::max(b.xmax, c.m_x);
            b.ymin = std::min(b.ymin, c.m_y);
            b.ymax = std::max(b.ymax, c.m_y);
            b.zmin = std::min(b.zmin, c.m_z);
            b.zmax = std::max(b.zmax, c.m_z);
        }
        if (!std::isfinite(b.xmin)) {
            b = DomainBounds{};
        }
        return b;
    }

    inline int InferBoundaryTagFromPoint(const Vector& p, const DomainBounds& b, bool enableZTag) {
        const double Lx = std::max(1.0, std::abs(b.xmax - b.xmin));
        const double Ly = std::max(1.0, std::abs(b.ymax - b.ymin));
        const double Lz = std::max(1.0, std::abs(b.zmax - b.zmin));
        const double tolx = std::max(1.0e-8, 1.0e-7 * Lx);
        const double toly = std::max(1.0e-8, 1.0e-7 * Ly);
        const double tolz = std::max(1.0e-8, 1.0e-7 * Lz);

        if (std::abs(p.m_x - b.xmin) <= tolx) return MeshTags::LEFT;
        if (std::abs(p.m_x - b.xmax) <= tolx) return MeshTags::RIGHT;
        if (std::abs(p.m_y - b.ymin) <= toly) return MeshTags::BOTTOM;
        if (std::abs(p.m_y - b.ymax) <= toly) return MeshTags::TOP;
        if (enableZTag) {
            if (std::abs(p.m_z - b.zmin) <= tolz) return MeshTags::TAG_FRONT;
            if (std::abs(p.m_z - b.zmax) <= tolz) return MeshTags::TAG_BACK;
        }
        return -1;
    }

    inline bool IsTemperatureFieldName(const BoundaryFieldTags& tags, const std::string& fieldName, int dofOffset) {
        if (fieldName == tags.temperature) return true;
        return dofOffset == 2;
    }

    inline bool IsSaturationFieldName(const BoundaryFieldTags& tags, const std::string& fieldName, int dofOffset) {
        if (fieldName == tags.saturation) return true;
        return dofOffset == 1 && fieldName != tags.temperature;
    }

    inline int ResolvePrimaryJacobianSlot(const BoundaryFieldTags& tags, const std::string& fieldName, int dofOffset) {
        if (IsTemperatureFieldName(tags, fieldName, dofOffset)) return 2;
        if (IsSaturationFieldName(tags, fieldName, dofOffset)) return 1;
        return 0;
    }

    inline ADVar<3> MakeStateAD(double val, int gradIdx) {
        ADVar<3> x = MakeConstAD<3>(val);
        if (gradIdx >= 0 && gradIdx < 3) {
            x.grad(gradIdx) = 1.0;
        }
        return x;
    }

    inline ADVar<3> ClampSwFractureForConstitutive(const ADVar<3>& sw) {
        constexpr double kOuterSlope = 1.0e-6;
        if (sw.val < 0.0) {
            ADVar<3> out = MakeConstAD<3>(0.0);
            for (int k = 0; k < 3; ++k) out.grad(k) = kOuterSlope * sw.grad(k);
            return out;
        }
        if (sw.val > 1.0) {
            ADVar<3> out = MakeConstAD<3>(1.0);
            for (int k = 0; k < 3; ++k) out.grad(k) = kOuterSlope * sw.grad(k);
            return out;
        }
        return sw;
    }

    inline void AddBoundaryEqContribution(
        int eqIdx,
        const ADVar<3>& q,
        std::vector<double>& residualRef,
        std::vector<std::array<double, 3>>& jacobianFull)
    {
        if (eqIdx < 0 || eqIdx >= static_cast<int>(residualRef.size()) || eqIdx >= static_cast<int>(jacobianFull.size())) {
            return;
        }
        residualRef[eqIdx] += q.val;
        jacobianFull[eqIdx][0] += q.grad(0);
        jacobianFull[eqIdx][1] += q.grad(1);
        jacobianFull[eqIdx][2] += q.grad(2);
    }

    struct BoundaryLocalState {
        ADVar<3> P = MakeConstAD<3>(0.0);
        ADVar<3> Sw = MakeConstAD<3>(1.0);
        ADVar<3> T = MakeConstAD<3>(300.0);
        ADVar<3> mob_w = MakeConstAD<3>(0.0);
        ADVar<3> mob_g = MakeConstAD<3>(0.0);
        ADVar<3> h_w = MakeConstAD<3>(0.0);
        ADVar<3> h_g = MakeConstAD<3>(0.0);
        ADVar<3> lambda_eff = MakeConstAD<3>(1.0);
        ADVar<3> k_n = MakeConstAD<3>(1.0e-13);
        bool has_sw = false;
        bool valid = false;
    };

    template<typename FieldManagerType>
    inline BoundaryLocalState BuildBoundaryLocalState(
        FieldManagerType& fm,
        bool fracture_domain,
        int localIdx,
        const Vector& faceNormal,
        const FIM_Engine::FluidPropertyEvalContext& fluid_ctx,
        const CapRelPerm::VGParams& vg,
        const CapRelPerm::RelPermParams& rp,
        const BoundaryFieldTags& tags)
    {
        BoundaryLocalState s;

        auto getScalar = [&](const std::string& name) -> std::shared_ptr<volScalarField> {
            return fracture_domain ? fm.getFractureScalar(name) : fm.getMatrixScalar(name);
        };

        auto getScalarWithAliases = [&](const std::string& primary,
            const auto& aliases) -> std::shared_ptr<volScalarField> {
                auto fld = getScalar(primary);
                if (fld) return fld;
                for (const auto& alias : aliases) {
                    if (alias.empty() || primary == alias) continue;
                    fld = getScalar(alias);
                    if (fld) return fld;
                }
                return nullptr;
            };

        // Keep legacy test/workflow compatibility: support historical field aliases.
        auto fP = getScalarWithAliases(tags.pressure, tags.pressure_aliases);
        auto fT = getScalarWithAliases(tags.temperature, tags.temperature_aliases);
        auto fSw = getScalarWithAliases(tags.saturation, tags.saturation_aliases);
        if (!isValidLocalIdx(localIdx, fP)) {
            return s;
        }

        s.P = MakeStateAD(fP->data[localIdx], 0);
        s.T = MakeStateAD(safeGetFieldValue(fT, localIdx, 300.0), 2);
        s.has_sw = isValidLocalIdx(localIdx, fSw);
        s.Sw = MakeStateAD(s.has_sw ? fSw->data[localIdx] : 1.0, s.has_sw ? 1 : -1);

        ADVar<3> Sw_const = s.Sw;
        ADVar<3> Pc = MakeConstAD<3>(0.0);
        ADVar<3> krw = MakeConstAD<3>(1.0);
        ADVar<3> krg = MakeConstAD<3>(0.0);

        if (s.has_sw) {
            if (fracture_domain) {
                Sw_const = ClampSwFractureForConstitutive(Sw_const);
                krw = Sw_const;
                krg = MakeConstAD<3>(1.0) - Sw_const;
                Pc = MakeConstAD<3>(0.0);
            }
            else {
                Sw_const = ClampSwForConstitutiveWell<3>(s.Sw, vg);
                CapRelPerm::kr_Mualem_vG<3>(Sw_const, vg, rp, krw, krg);
                Pc = CapRelPerm::pc_vG<3>(Sw_const, vg);
            }
        }
        const AD_Fluid::ADFluidProperties<3> propsW = s.has_sw
            ? FIM_Engine::EvalTwoPhaseWaterFluid<3>(fluid_ctx, s.P, s.T)
            : FIM_Engine::EvalSinglePhaseFluid<3>(fluid_ctx, s.P, s.T);
        const AD_Fluid::ADFluidProperties<3> propsG = s.has_sw
            ? FIM_Engine::EvalTwoPhaseGasFluid<3>(fluid_ctx, s.P + Pc, s.T)
            : FIM_Engine::EvalTwoPhaseGasFluid<3>(fluid_ctx, s.P, s.T);

        s.mob_w = (propsW.rho * krw) / propsW.mu;
        s.mob_g = s.has_sw ? ((propsG.rho * krg) / propsG.mu) : MakeConstAD<3>(0.0);
        s.h_w = propsW.h;
        s.h_g = propsG.h;

        auto fPhi = getScalar(fracture_domain ? tags.fr_phi : tags.phi);
        auto fLamRock = getScalar(fracture_domain ? tags.fr_lambda : tags.lambda_rock);
        auto fLamEff = getScalar(tags.lambda_eff);

        const double phi_val = std::max(0.0, std::min(1.0, safeGetFieldValue(fPhi, localIdx, 0.2)));
        const double lam_r = std::max(0.0, safeGetFieldValue(fLamRock, localIdx, 2.0));
        if (isValidLocalIdx(localIdx, fLamEff)) {
            s.lambda_eff = MakeConstAD<3>(std::max(0.0, fLamEff->data[localIdx]));
        }
        else if (s.has_sw) {
            s.lambda_eff = MakeConstAD<3>((1.0 - phi_val) * lam_r) + MakeConstAD<3>(phi_val) * (Sw_const * propsW.k + (MakeConstAD<3>(1.0) - Sw_const) * propsG.k);
        }
        else {
            s.lambda_eff = MakeConstAD<3>((1.0 - phi_val) * lam_r) + MakeConstAD<3>(phi_val) * propsW.k;
        }

        double k_n_val = 1.0e-13;
        if (fracture_domain) {
            auto fk = getScalar(tags.fr_k_t);
            k_n_val = std::max(1.0e-20, safeGetFieldValue(fk, localIdx, 1.0e-12));
        }
        else {
            const double kxx = std::max(0.0, safeGetFieldValue(getScalar(tags.k_xx), localIdx, 1.0e-13));
            const double kyy = std::max(0.0, safeGetFieldValue(getScalar(tags.k_yy), localIdx, 1.0e-13));
            const double kzz = std::max(0.0, safeGetFieldValue(getScalar(tags.k_zz), localIdx, 1.0e-13));
            const double nx = faceNormal.m_x;
            const double ny = faceNormal.m_y;
            const double nz = faceNormal.m_z;
            k_n_val = nx * nx * kxx + ny * ny * kyy + nz * nz * kzz;
            if (k_n_val <= 0.0) {
                k_n_val = std::max({ kxx, kyy, kzz, 1.0e-20 });
            }
        }
        s.k_n = MakeConstAD<3>(std::max(1.0e-20, k_n_val));

        s.valid = true;
        return s;
    }

    inline ADVar<3> ComputePhaseBoundaryMassFlux(
        const BoundarySetting::BoundaryConditionManager* phaseBC,
        int physicalTag,
        const Vector& faceMidpoint,
        double area,
        double T_geom,
        const ADVar<3>& primaryVar,
        const ADVar<3>& k_n,
        const ADVar<3>& mobility)
    {
        ADVar<3> q = MakeConstAD<3>(0.0);
        if (!phaseBC || !phaseBC->HasBC(physicalTag)) {
            return q;
        }

        const auto bc = phaseBC->GetBCCoefficients(physicalTag, faceMidpoint);
        if (bc.type == BoundarySetting::BoundaryType::Dirichlet) {
            if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                return q;
            }
            const ADVar<3> raw = FVM_Ops::Op_Boundary_Dirichlet_AD<3, ADVar<3>>(T_geom, primaryVar, bc.c / bc.a);
            const ADVar<3> q_darcy = FVM_Ops::Op_Boundary_ScaleFlux_AD<3, ADVar<3>>(raw, k_n);
            return q_darcy * mobility;
        }

        if (bc.type == BoundarySetting::BoundaryType::Robin) {
            if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                return q;
            }
            const double C_L = -bc.a;
            const double far_val = bc.c / bc.a;
            const ADVar<3> raw = MakeConstAD<3>(area) * FVM_Ops::Op_Leakoff_Source_AD<3, ADVar<3>>(true, C_L, primaryVar, far_val);
            const ADVar<3> q_darcy = FVM_Ops::Op_Boundary_ScaleFlux_AD<3, ADVar<3>>(raw, k_n);
            return q_darcy * mobility;
        }

        return FVM_Ops::Op_Boundary_Neumann_AD<3, ADVar<3>>(area, bc.c);
    }

    inline bool TryGetInflowSwFraction(
        const BoundarySetting::BoundaryConditionManager* saturationBC,
        int physicalTag,
        const Vector& boundaryPoint,
        BoundaryAssemblyStats& stats,
        double& sw_bc_out)
    {
        if (!saturationBC || !saturationBC->HasBC(physicalTag)) {
            return false;
        }

        const auto bc = saturationBC->GetBCCoefficients(physicalTag, boundaryPoint);
        if (bc.type != BoundarySetting::BoundaryType::Dirichlet) {
            ++stats.satDriveIgnored;
            return false;
        }
        if (std::abs(bc.a) <= kBoundaryCoeffEps) {
            ++stats.satDriveIgnored;
            return false;
        }

        sw_bc_out = std::clamp(bc.c / bc.a, 0.0, 1.0);
        return true;
    }

    inline bool ComputeTemperatureConductionFlux(
        const BoundarySetting::BoundaryConditionManager* temperatureBC,
        int physicalTag,
        const Vector& boundaryPoint,
        double area,
        double T_geom,
        const BoundaryLocalState& state,
        ADVar<3>& q_cond)
    {
        q_cond = MakeConstAD<3>(0.0);
        if (!temperatureBC || !temperatureBC->HasBC(physicalTag)) {
            return false;
        }

        const auto bc = temperatureBC->GetBCCoefficients(physicalTag, boundaryPoint);
        if (bc.type == BoundarySetting::BoundaryType::Dirichlet) {
            if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                return false;
            }
            const ADVar<3> raw = FVM_Ops::Op_Boundary_Dirichlet_AD<3, ADVar<3>>(T_geom, state.T, bc.c / bc.a);
            q_cond = FVM_Ops::Op_Boundary_ScaleFlux_AD<3, ADVar<3>>(raw, state.lambda_eff);
            return true;
        }

        if (bc.type == BoundarySetting::BoundaryType::Robin) {
            if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                return false;
            }
            const double C_L = -bc.a;
            const double far_val = bc.c / bc.a;
            const ADVar<3> raw = MakeConstAD<3>(area) * FVM_Ops::Op_Leakoff_Source_AD<3, ADVar<3>>(true, C_L, state.T, far_val);
            q_cond = FVM_Ops::Op_Boundary_ScaleFlux_AD<3, ADVar<3>>(raw, state.lambda_eff);
            return true;
        }

        q_cond = FVM_Ops::Op_Boundary_Neumann_AD<3, ADVar<3>>(area, bc.c);
        return true;
    }

    inline bool ComputePressureDrivenTwoPhaseMassFlux(
        const BoundarySetting::BoundaryConditionManager* pressureBC,
        const BoundarySetting::BoundaryConditionManager* saturationBC,
        int physicalTag,
        const Vector& boundaryPoint,
        double area,
        double T_geom,
        const BoundaryLocalState& state,
        BoundaryAssemblyStats& stats,
        ADVar<3>& q_w,
        ADVar<3>& q_g)
    {
        q_w = MakeConstAD<3>(0.0);
        q_g = MakeConstAD<3>(0.0);
        if (!pressureBC || !pressureBC->HasBC(physicalTag)) {
            return false;
        }

        const auto bc = pressureBC->GetBCCoefficients(physicalTag, boundaryPoint);
        if (bc.type == BoundarySetting::BoundaryType::Dirichlet ||
            bc.type == BoundarySetting::BoundaryType::Robin) {
            ADVar<3> q_darcy = MakeConstAD<3>(0.0);
            if (bc.type == BoundarySetting::BoundaryType::Dirichlet) {
                if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                    return false;
                }
                const ADVar<3> raw = FVM_Ops::Op_Boundary_Dirichlet_AD<3, ADVar<3>>(T_geom, state.P, bc.c / bc.a);
                q_darcy = FVM_Ops::Op_Boundary_ScaleFlux_AD<3, ADVar<3>>(raw, state.k_n);
            }
            else {
                if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                    return false;
                }
                const double C_L = -bc.a;
                const double far_val = bc.c / bc.a;
                const ADVar<3> raw = MakeConstAD<3>(area) * FVM_Ops::Op_Leakoff_Source_AD<3, ADVar<3>>(true, C_L, state.P, far_val);
                q_darcy = FVM_Ops::Op_Boundary_ScaleFlux_AD<3, ADVar<3>>(raw, state.k_n);
            }

            q_w = q_darcy * state.mob_w;
            q_g = state.has_sw ? (q_darcy * state.mob_g) : MakeConstAD<3>(0.0);

            const ADVar<3> q_total_mass = q_w + q_g;
            if (state.has_sw && q_total_mass.val < 0.0) {
                double sw_in = 1.0;
                if (TryGetInflowSwFraction(saturationBC, physicalTag, boundaryPoint, stats, sw_in)) {
                    ++stats.satCompositionUsed;
                    const double sg_in = 1.0 - sw_in;
                    q_w = q_total_mass * sw_in;
                    q_g = q_total_mass * sg_in;
                }
            }
            return true;
        }

        ADVar<3> q_total = FVM_Ops::Op_Boundary_Neumann_AD<3, ADVar<3>>(area, bc.c);
        if (!state.has_sw) {
            q_w = q_total;
            return true;
        }

        double fw = 1.0;
        if (q_total.val < 0.0) {
            double sw_in = 1.0;
            if (TryGetInflowSwFraction(saturationBC, physicalTag, boundaryPoint, stats, sw_in)) {
                ++stats.satCompositionUsed;
                fw = sw_in;
            }
            else {
                const double mw = std::max(0.0, state.mob_w.val);
                const double mg = std::max(0.0, state.mob_g.val);
                const double msum = mw + mg;
                fw = (msum > 1.0e-30) ? (mw / msum) : 1.0;
            }
        }
        else {
            const double mw = std::max(0.0, state.mob_w.val);
            const double mg = std::max(0.0, state.mob_g.val);
            const double msum = mw + mg;
            fw = (msum > 1.0e-30) ? (mw / msum) : 1.0;
        }

        fw = std::clamp(fw, 0.0, 1.0);
        const double fg = 1.0 - fw;
        q_w = q_total * fw;
        q_g = q_total * fg;
        return true;
    }

    inline void AccumulateCoupledEqContribution(
        int eqIdx,
        int primarySlot,
        const ADVar<3>& q,
        std::vector<double>& residual,
        std::vector<std::array<double, 3>>& jacobianFull,
        BoundaryAssemblyStats& stats)
    {
        if (eqIdx < 0 || eqIdx >= static_cast<int>(residual.size()) || eqIdx >= static_cast<int>(jacobianFull.size())) {
            ++stats.invalidEqRows;
            return;
        }

        ++stats.visitedEqRows;
        const bool all_zero = (std::abs(q.val) <= 1.0e-16) &&
            (std::abs(q.grad(0)) <= 1.0e-16) &&
            (std::abs(q.grad(1)) <= 1.0e-16) &&
            (std::abs(q.grad(2)) <= 1.0e-16);
        if (all_zero) {
            ++stats.zeroEqRows;
            return;
        }

        AddBoundaryEqContribution(eqIdx, q, residual, jacobianFull);
        ++stats.nonzeroEqRows;
        stats.sumResidual += q.val;
        if (primarySlot >= 0 && primarySlot < 3) {
            stats.sumJacobianDiag += q.grad(primarySlot);
        }
    }
} // namespace end

BoundaryAssemblyStats BoundaryAssembler::Assemble_2D(
    MeshManager& mgr,
    const BoundarySetting::BoundaryConditionManager& bcMgr,
    int dofOffset,
    FieldManager_2D& fm,
    const std::string& fieldName,
    std::vector<double>& residual,
    std::vector<double>& jacobianDiag)
{
    BoundaryAssemblyStats stats;
    if (dofOffset < 0 || dofOffset >= 3) {
        std::cerr << "[BoundaryAssembler::Assemble_2D] Invalid dofOffset = " << dofOffset << " (expected 0..2)." << std::endl;
        return stats;
    }
    if (residual.size() != jacobianDiag.size()) {
        std::cerr << "[BoundaryAssembler::Assemble_2D] residual/jacobianDiag size mismatch: "
            << residual.size() << " vs " << jacobianDiag.size() << std::endl;
        return stats;
    }

    std::vector<std::array<double, 3>> jacobianFull(residual.size(), std::array<double, 3>{ 0.0, 0.0, 0.0 });
    stats = Assemble_2D_FullJac(
        mgr, bcMgr, dofOffset, fm, fieldName, residual, jacobianFull,
        nullptr, nullptr, FluidPropertyEvalConfig(), CapRelPerm::VGParams(), CapRelPerm::RelPermParams());

    const BoundaryFieldTags tags = BuildBoundaryFieldTags();
    const int jacSlot = ResolvePrimaryJacobianSlot(tags, fieldName, dofOffset);
    for (size_t i = 0; i < jacobianDiag.size(); ++i) {
        jacobianDiag[i] += jacobianFull[i][jacSlot];
    }
    return stats;
}

BoundaryAssemblyStats BoundaryAssembler::Assemble_2D_FullJac(
    MeshManager& mgr,
    const BoundarySetting::BoundaryConditionManager& bcMgr,
    int dofOffset,
    FieldManager_2D& fm,
    const std::string& fieldName,
    std::vector<double>& residual,
    std::vector<std::array<double, 3>>& jacobianFull,
    const BoundarySetting::BoundaryConditionManager* coupledPressureBC,
    const BoundarySetting::BoundaryConditionManager* coupledSaturationBC,
    const FluidPropertyEvalConfig& fluid_cfg,
    const CapRelPerm::VGParams& vg,
    const CapRelPerm::RelPermParams& rp)
{
    BoundaryAssemblyStats stats;
    const FIM_Engine::FluidPropertyEvalContext fluid_ctx = BuildFluidEvalContext(fluid_cfg);
    if (dofOffset < 0 || dofOffset >= 3) {
        std::cerr << "[BoundaryAssembler::Assemble_2D_FullJac] Invalid dofOffset = " << dofOffset << " (expected 0..2)." << std::endl;
        return stats;
    }
    if (residual.size() != jacobianFull.size()) {
        std::cerr << "[BoundaryAssembler::Assemble_2D_FullJac] residual/jacobianFull size mismatch: "
            << residual.size() << " vs " << jacobianFull.size() << std::endl;
        return stats;
    }

    const BoundaryFieldTags tags = BuildBoundaryFieldTags();
    const bool is_temp_eq = IsTemperatureFieldName(tags, fieldName, dofOffset);
    const bool is_sat_eq = IsSaturationFieldName(tags, fieldName, dofOffset);

    const auto& mesh = mgr.mesh();
    for (size_t i = 0; i < mesh.getCells().size(); ++i) {
        const auto& cell = mesh.getCells()[i];
        const int solverIdx = static_cast<int>(i);
        const int eqIdx = mgr.getEquationIndex(solverIdx, dofOffset);
        if (eqIdx < 0 || eqIdx >= static_cast<int>(residual.size())) {
            ++stats.invalidEqRows;
            continue;
        }

        for (int globalFaceID : cell.CellFaceIDs) {
            const int faceIdx = mesh.getFaceIndex(globalFaceID);
            if (faceIdx < 0) continue;
            const auto& face = mesh.getFaces()[faceIdx];
            if (!face.isBoundary() || !bcMgr.HasBC(face.physicalGroupId)) continue;

            const auto bc = bcMgr.GetBCCoefficients(face.physicalGroupId, face.midpoint);
            auto& tagStats = TouchBoundaryTagStats(stats, face.physicalGroupId, bc.type);
            ++tagStats.faces;
            ++stats.visitedEqRows;

            const BoundaryLocalState state = BuildBoundaryLocalState(
                fm, false, static_cast<int>(i), face.normal, fluid_ctx, vg, rp, tags);
            if (!state.valid) {
                ++tagStats.skipped;
                ++stats.invalidEqRows;
                continue;
            }

            const Vector d_owner = face.midpoint - cell.center;
            const double dist = std::max(std::abs(d_owner * face.normal), 1.0e-12);
            double area = std::max(face.vectorE.Mag(), 0.0);
            if (area <= 1.0e-20) area = std::max(face.length, 0.0);
            const double T_geom = area / dist;

            const ADVar<3> var_cell = is_temp_eq ? state.T : (is_sat_eq ? state.Sw : state.P);
            const ADVar<3> phys_coeff = is_temp_eq
                ? state.lambda_eff
                : (is_sat_eq ? (state.k_n * state.mob_g) : (state.k_n * state.mob_w));

            ADVar<3> raw_flux = MakeConstAD<3>(0.0);
            ADVar<3> q_flux = MakeConstAD<3>(0.0);

            if (bc.type == BoundarySetting::BoundaryType::Dirichlet) {
                if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                    ++tagStats.skipped;
                    ++stats.zeroEqRows;
                    continue;
                }
                raw_flux = FVM_Ops::Op_Boundary_Dirichlet_AD<3, ADVar<3>>(T_geom, var_cell, bc.c / bc.a);
                q_flux = FVM_Ops::Op_Boundary_ScaleFlux_AD<3, ADVar<3>>(raw_flux, phys_coeff);
            }
            else if (bc.type == BoundarySetting::BoundaryType::Neumann) {
                raw_flux = FVM_Ops::Op_Boundary_Neumann_AD<3, ADVar<3>>(area, bc.c);
                q_flux = raw_flux;
            }
            else if (bc.type == BoundarySetting::BoundaryType::Robin) {
                if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                    ++tagStats.skipped;
                    ++stats.zeroEqRows;
                    continue;
                }
                const double C_L = -bc.a;
                const double far_val = bc.c / bc.a;
                raw_flux = MakeConstAD<3>(area) * FVM_Ops::Op_Leakoff_Source_AD<3, ADVar<3>>(true, C_L, var_cell, far_val);
                // Legacy Robin/Leakoff semantics in Assemble_*_FullJac:
                // q = area * C_L * (phi - phi_far), without extra k/mob/lambda scaling.
                q_flux = raw_flux;
            }

            if (is_temp_eq) {
                const ADVar<3> q_mass_w = ComputePhaseBoundaryMassFlux(
                    coupledPressureBC, face.physicalGroupId, face.midpoint, area, T_geom, state.P, state.k_n, state.mob_w);
                ADVar<3> q_mass_g = MakeConstAD<3>(0.0);
                if (state.has_sw) {
                    q_mass_g = ComputePhaseBoundaryMassFlux(
                        coupledSaturationBC, face.physicalGroupId, face.midpoint, area, T_geom, state.Sw, state.k_n, state.mob_g);
                }
                q_flux += q_mass_w * state.h_w + q_mass_g * state.h_g;
            }

            AddBoundaryEqContribution(eqIdx, q_flux, residual, jacobianFull);

            const double diag_add = q_flux.grad(dofOffset);
            if (std::abs(q_flux.val) <= 1.0e-16 &&
                std::abs(q_flux.grad(0)) <= 1.0e-16 &&
                std::abs(q_flux.grad(1)) <= 1.0e-16 &&
                std::abs(q_flux.grad(2)) <= 1.0e-16) {
                ++tagStats.skipped;
                ++stats.zeroEqRows;
            }
            else {
                ++tagStats.applied;
                ++stats.nonzeroEqRows;
                tagStats.sumR += std::abs(q_flux.val);
                tagStats.sumDiag += std::abs(diag_add);
            }

            stats.sumResidual += q_flux.val;
            stats.sumJacobianDiag += diag_add;
            stats.matrixBCCount++;
        }
    }

    const int nMat = mgr.getMatrixDOFCount();
    const DomainBounds bounds = ComputeDomainBounds(mesh);
    const auto& frNet = mgr.fracture_network();
    const auto& orderedElems = frNet.getOrderedFractureElements();

    for (const FractureElement* elemPtr : orderedElems) {
        if (!elemPtr) continue;
        const FractureElement& elem = *elemPtr;

        const int solverIdx = elem.solverIndex;
        if (solverIdx < nMat) continue;

        const int localIdx = solverIdx - nMat;
        const int eqIdx = mgr.getEquationIndex(solverIdx, dofOffset);
        if (eqIdx < 0 || eqIdx >= static_cast<int>(residual.size()) || localIdx < 0) {
            ++stats.invalidEqRows;
            continue;
        }

        if (elem.parentFractureID < 0 || elem.parentFractureID >= static_cast<int>(frNet.fractures.size())) {
            ++stats.invalidEqRows;
            continue;
        }
        const auto& frac = frNet.fractures[elem.parentFractureID];
        if (elem.id <= 0 || elem.id >= static_cast<int>(frac.intersections.size())) {
            ++stats.invalidEqRows;
            continue;
        }

        const Vector p0 = frac.intersections[elem.id - 1].point;
        const Vector p1 = frac.intersections[elem.id].point;

        Vector tangent = p1 - p0;
        if (tangent.Mag() <= 1.0e-12) {
            tangent = frac.end - frac.start;
        }
        const double tanMag = std::max(tangent.Mag(), 1.0e-12);
        const Vector tanUnit = tangent / tanMag;

        const double dist = std::max(0.5 * std::max(elem.length, 0.0), 1.0e-12);
        const double area = std::max(elem.aperture, 1.0e-12);

        for (int endId = 0; endId < 2; ++endId) {
            const bool isStart = (endId == 0);
            const Vector boundaryPoint = isStart ? p0 : p1;
            const Vector boundaryNormal = tanUnit * (isStart ? -1.0 : 1.0);
            const bool isTerminalEndpoint =
                (isStart && elem.id == 1) ||
                (!isStart && elem.id == static_cast<int>(frac.intersections.size()) - 1);

            const int inferredTag = InferBoundaryTagFromPoint(boundaryPoint, bounds, false);
            if (inferredTag < 0) {
                if (isTerminalEndpoint) {
                    ++stats.inferredTagFail;
                }
                continue;
            }
            if (!bcMgr.HasBC(inferredTag)) continue;

            const auto bc = bcMgr.GetBCCoefficients(inferredTag, boundaryPoint);
            auto& tagStats = TouchBoundaryTagStats(stats, inferredTag, bc.type);
            ++tagStats.faces;
            ++stats.visitedEqRows;

            const BoundaryLocalState state = BuildBoundaryLocalState(
                fm, true, localIdx, boundaryNormal, fluid_ctx, vg, rp, tags);
            if (!state.valid) {
                ++tagStats.skipped;
                ++stats.invalidEqRows;
                continue;
            }

            const double T_geom = area / dist;
            const ADVar<3> var_cell = is_temp_eq ? state.T : (is_sat_eq ? state.Sw : state.P);
            const ADVar<3> phys_coeff = is_temp_eq
                ? state.lambda_eff
                : (is_sat_eq ? (state.k_n * state.mob_g) : (state.k_n * state.mob_w));

            ADVar<3> raw_flux = MakeConstAD<3>(0.0);
            ADVar<3> q_flux = MakeConstAD<3>(0.0);

            if (bc.type == BoundarySetting::BoundaryType::Dirichlet) {
                if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                    ++tagStats.skipped;
                    ++stats.zeroEqRows;
                    continue;
                }
                raw_flux = FVM_Ops::Op_Boundary_Dirichlet_AD<3, ADVar<3>>(T_geom, var_cell, bc.c / bc.a);
                q_flux = FVM_Ops::Op_Boundary_ScaleFlux_AD<3, ADVar<3>>(raw_flux, phys_coeff);
            }
            else if (bc.type == BoundarySetting::BoundaryType::Neumann) {
                raw_flux = FVM_Ops::Op_Boundary_Neumann_AD<3, ADVar<3>>(area, bc.c);
                q_flux = raw_flux;
            }
            else if (bc.type == BoundarySetting::BoundaryType::Robin) {
                if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                    ++tagStats.skipped;
                    ++stats.zeroEqRows;
                    continue;
                }
                const double C_L = -bc.a;
                const double far_val = bc.c / bc.a;
                raw_flux = MakeConstAD<3>(area) * FVM_Ops::Op_Leakoff_Source_AD<3, ADVar<3>>(true, C_L, var_cell, far_val);
                q_flux = raw_flux;
            }

            if (is_temp_eq) {
                const ADVar<3> q_mass_w = ComputePhaseBoundaryMassFlux(
                    coupledPressureBC, inferredTag, boundaryPoint, area, T_geom, state.P, state.k_n, state.mob_w);
                ADVar<3> q_mass_g = MakeConstAD<3>(0.0);
                if (state.has_sw) {
                    q_mass_g = ComputePhaseBoundaryMassFlux(
                        coupledSaturationBC, inferredTag, boundaryPoint, area, T_geom, state.Sw, state.k_n, state.mob_g);
                }
                q_flux += q_mass_w * state.h_w + q_mass_g * state.h_g;
            }

            AddBoundaryEqContribution(eqIdx, q_flux, residual, jacobianFull);
            const double diag_add = q_flux.grad(dofOffset);

            if (std::abs(q_flux.val) <= 1.0e-16 &&
                std::abs(q_flux.grad(0)) <= 1.0e-16 &&
                std::abs(q_flux.grad(1)) <= 1.0e-16 &&
                std::abs(q_flux.grad(2)) <= 1.0e-16) {
                ++tagStats.skipped;
                ++stats.zeroEqRows;
            }
            else {
                ++tagStats.applied;
                ++stats.nonzeroEqRows;
                tagStats.sumR += std::abs(q_flux.val);
                tagStats.sumDiag += std::abs(diag_add);
            }

            stats.sumResidual += q_flux.val;
            stats.sumJacobianDiag += diag_add;
            stats.fractureBCCount++;
        }
    }

    return stats;
}

BoundaryAssemblyStats BoundaryAssembler::Assemble_3D(
    MeshManager_3D& mgr,
    const BoundarySetting::BoundaryConditionManager& bcMgr,
    int dofOffset,
    FieldManager_3D& fm,
    const std::string& fieldName,
    std::vector<double>& residual,
    std::vector<double>& jacobianDiag)
{
    BoundaryAssemblyStats stats;
    if (dofOffset < 0 || dofOffset >= 3) {
        std::cerr << "[BoundaryAssembler::Assemble_3D] Invalid dofOffset = " << dofOffset << " (expected 0..2)." << std::endl;
        return stats;
    }
    if (residual.size() != jacobianDiag.size()) {
        std::cerr << "[BoundaryAssembler::Assemble_3D] residual/jacobianDiag size mismatch: "
            << residual.size() << " vs " << jacobianDiag.size() << std::endl;
        return stats;
    }

    std::vector<std::array<double, 3>> jacobianFull(residual.size(), std::array<double, 3>{ 0.0, 0.0, 0.0 });
    stats = Assemble_3D_FullJac(
        mgr, bcMgr, dofOffset, fm, fieldName, residual, jacobianFull,
        nullptr, nullptr, FluidPropertyEvalConfig(), CapRelPerm::VGParams(), CapRelPerm::RelPermParams());

    const BoundaryFieldTags tags = BuildBoundaryFieldTags();
    const int jacSlot = ResolvePrimaryJacobianSlot(tags, fieldName, dofOffset);
    for (size_t i = 0; i < jacobianDiag.size(); ++i) {
        jacobianDiag[i] += jacobianFull[i][jacSlot];
    }
    return stats;
}

BoundaryAssemblyStats BoundaryAssembler::Assemble_3D_FullJac(
    MeshManager_3D& mgr,
    const BoundarySetting::BoundaryConditionManager& bcMgr,
    int dofOffset,
    FieldManager_3D& fm,
    const std::string& fieldName,
    std::vector<double>& residual,
    std::vector<std::array<double, 3>>& jacobianFull,
    const BoundarySetting::BoundaryConditionManager* coupledPressureBC,
    const BoundarySetting::BoundaryConditionManager* coupledSaturationBC,
    const FluidPropertyEvalConfig& fluid_cfg,
    const CapRelPerm::VGParams& vg,
    const CapRelPerm::RelPermParams& rp)
{
    BoundaryAssemblyStats stats;
    const FIM_Engine::FluidPropertyEvalContext fluid_ctx = BuildFluidEvalContext(fluid_cfg);
    if (dofOffset < 0 || dofOffset >= 3) {
        std::cerr << "[BoundaryAssembler::Assemble_3D_FullJac] Invalid dofOffset = " << dofOffset << " (expected 0..2)." << std::endl;
        return stats;
    }
    if (residual.size() != jacobianFull.size()) {
        std::cerr << "[BoundaryAssembler::Assemble_3D_FullJac] residual/jacobianFull size mismatch: "
            << residual.size() << " vs " << jacobianFull.size() << std::endl;
        return stats;
    }

    const BoundaryFieldTags tags = BuildBoundaryFieldTags();
    const bool is_temp_eq = IsTemperatureFieldName(tags, fieldName, dofOffset);
    const bool is_sat_eq = IsSaturationFieldName(tags, fieldName, dofOffset);

    const auto& mesh = mgr.mesh();
    for (size_t i = 0; i < mesh.getCells().size(); ++i) {
        const auto& cell = mesh.getCells()[i];
        const int solverIdx = static_cast<int>(i);
        const int eqIdx = mgr.getEquationIndex(solverIdx, dofOffset);
        if (eqIdx < 0 || eqIdx >= static_cast<int>(residual.size())) {
            ++stats.invalidEqRows;
            continue;
        }

        for (int globalFaceID : cell.CellFaceIDs) {
            const int faceIdx = mesh.getFaceIndex(globalFaceID);
            if (faceIdx < 0) continue;
            const auto& face = mesh.getFaces()[faceIdx];
            if (!face.isBoundary() || !bcMgr.HasBC(face.physicalGroupId)) continue;

            const auto bc = bcMgr.GetBCCoefficients(face.physicalGroupId, face.midpoint);
            auto& tagStats = TouchBoundaryTagStats(stats, face.physicalGroupId, bc.type);
            ++tagStats.faces;
            ++stats.visitedEqRows;

            const BoundaryLocalState state = BuildBoundaryLocalState(
                fm, false, static_cast<int>(i), face.normal, fluid_ctx, vg, rp, tags);
            if (!state.valid) {
                ++tagStats.skipped;
                ++stats.invalidEqRows;
                continue;
            }

            const Vector d_owner = face.midpoint - cell.center;
            const double dist = std::max(std::abs(d_owner * face.normal), 1.0e-12);
            double area = std::max(face.vectorE.Mag(), 0.0);
            if (area <= 1.0e-20) area = std::max(face.length, 0.0);
            const double T_geom = area / dist;

            const ADVar<3> var_cell = is_temp_eq ? state.T : (is_sat_eq ? state.Sw : state.P);
            const ADVar<3> phys_coeff = is_temp_eq
                ? state.lambda_eff
                : (is_sat_eq ? (state.k_n * state.mob_g) : (state.k_n * state.mob_w));

            ADVar<3> raw_flux = MakeConstAD<3>(0.0);
            ADVar<3> q_flux = MakeConstAD<3>(0.0);

            if (bc.type == BoundarySetting::BoundaryType::Dirichlet) {
                if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                    ++tagStats.skipped;
                    ++stats.zeroEqRows;
                    continue;
                }
                raw_flux = FVM_Ops::Op_Boundary_Dirichlet_AD<3, ADVar<3>>(T_geom, var_cell, bc.c / bc.a);
                q_flux = FVM_Ops::Op_Boundary_ScaleFlux_AD<3, ADVar<3>>(raw_flux, phys_coeff);
            }
            else if (bc.type == BoundarySetting::BoundaryType::Neumann) {
                raw_flux = FVM_Ops::Op_Boundary_Neumann_AD<3, ADVar<3>>(area, bc.c);
                q_flux = raw_flux;
            }
            else if (bc.type == BoundarySetting::BoundaryType::Robin) {
                if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                    ++tagStats.skipped;
                    ++stats.zeroEqRows;
                    continue;
                }
                const double C_L = -bc.a;
                const double far_val = bc.c / bc.a;
                raw_flux = MakeConstAD<3>(area) * FVM_Ops::Op_Leakoff_Source_AD<3, ADVar<3>>(true, C_L, var_cell, far_val);
                q_flux = raw_flux;
            }

            if (is_temp_eq) {
                const ADVar<3> q_mass_w = ComputePhaseBoundaryMassFlux(
                    coupledPressureBC, face.physicalGroupId, face.midpoint, area, T_geom, state.P, state.k_n, state.mob_w);
                ADVar<3> q_mass_g = MakeConstAD<3>(0.0);
                if (state.has_sw) {
                    q_mass_g = ComputePhaseBoundaryMassFlux(
                        coupledSaturationBC, face.physicalGroupId, face.midpoint, area, T_geom, state.Sw, state.k_n, state.mob_g);
                }
                q_flux += q_mass_w * state.h_w + q_mass_g * state.h_g;
            }

            AddBoundaryEqContribution(eqIdx, q_flux, residual, jacobianFull);

            const double diag_add = q_flux.grad(dofOffset);
            if (std::abs(q_flux.val) <= 1.0e-16 &&
                std::abs(q_flux.grad(0)) <= 1.0e-16 &&
                std::abs(q_flux.grad(1)) <= 1.0e-16 &&
                std::abs(q_flux.grad(2)) <= 1.0e-16) {
                ++tagStats.skipped;
                ++stats.zeroEqRows;
            }
            else {
                ++tagStats.applied;
                ++stats.nonzeroEqRows;
                tagStats.sumR += std::abs(q_flux.val);
                tagStats.sumDiag += std::abs(diag_add);
            }

            stats.sumResidual += q_flux.val;
            stats.sumJacobianDiag += diag_add;
            stats.matrixBCCount++;
        }
    }

    const int nMat = get3DMatrixOffset(mgr);
    const DomainBounds bounds = ComputeDomainBounds(mesh);
    for (const auto& edge : mgr.fracture_network().getGlobalEdges()) {
        if (edge.neighborCell_solverIndex != -1) continue;
        if (edge.ownerCell_solverIndex < nMat) continue;

        const int inferredTag = InferBoundaryTagFromPoint(edge.midpoint, bounds, true);
        if (inferredTag < 0) {
            ++stats.inferredTagFail;
            continue;
        }
        if (!bcMgr.HasBC(inferredTag)) continue;

        const int solverIdx = edge.ownerCell_solverIndex;
        const int localIdx = solverIdx - nMat;
        const int eqIdx = mgr.getEquationIndex(solverIdx, dofOffset);
        if (eqIdx < 0 || eqIdx >= static_cast<int>(residual.size()) || localIdx < 0) {
            ++stats.invalidEqRows;
            continue;
        }

        const auto bc = bcMgr.GetBCCoefficients(inferredTag, edge.midpoint);
        auto& tagStats = TouchBoundaryTagStats(stats, inferredTag, bc.type);
        ++tagStats.faces;
        ++stats.visitedEqRows;

        const BoundaryLocalState state = BuildBoundaryLocalState(
            fm, true, localIdx, edge.normal, fluid_ctx, vg, rp, tags);
        if (!state.valid) {
            ++tagStats.skipped;
            ++stats.invalidEqRows;
            continue;
        }

        double dist = std::abs(edge.ownerToNeighbor * edge.normal);
        if (dist <= 1.0e-12) {
            dist = edge.ownerToNeighbor.Mag();
        }
        dist = std::max(dist, 1.0e-12);
        const double area = std::max(edge.length, 0.0);
        const double T_geom = area / dist;

        const ADVar<3> var_cell = is_temp_eq ? state.T : (is_sat_eq ? state.Sw : state.P);
        const ADVar<3> phys_coeff = is_temp_eq
            ? state.lambda_eff
            : (is_sat_eq ? (state.k_n * state.mob_g) : (state.k_n * state.mob_w));

        ADVar<3> raw_flux = MakeConstAD<3>(0.0);
        ADVar<3> q_flux = MakeConstAD<3>(0.0);

        if (bc.type == BoundarySetting::BoundaryType::Dirichlet) {
            if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                ++tagStats.skipped;
                ++stats.zeroEqRows;
                continue;
            }
            raw_flux = FVM_Ops::Op_Boundary_Dirichlet_AD<3, ADVar<3>>(T_geom, var_cell, bc.c / bc.a);
            q_flux = FVM_Ops::Op_Boundary_ScaleFlux_AD<3, ADVar<3>>(raw_flux, phys_coeff);
        }
        else if (bc.type == BoundarySetting::BoundaryType::Neumann) {
            raw_flux = FVM_Ops::Op_Boundary_Neumann_AD<3, ADVar<3>>(area, bc.c);
            q_flux = raw_flux;
        }
        else if (bc.type == BoundarySetting::BoundaryType::Robin) {
            if (std::abs(bc.a) <= kBoundaryCoeffEps) {
                ++tagStats.skipped;
                ++stats.zeroEqRows;
                continue;
            }
            const double C_L = -bc.a;
            const double far_val = bc.c / bc.a;
            raw_flux = MakeConstAD<3>(area) * FVM_Ops::Op_Leakoff_Source_AD<3, ADVar<3>>(true, C_L, var_cell, far_val);
            q_flux = raw_flux;
        }

        if (is_temp_eq) {
            const ADVar<3> q_mass_w = ComputePhaseBoundaryMassFlux(
                coupledPressureBC, inferredTag, edge.midpoint, area, T_geom, state.P, state.k_n, state.mob_w);
            ADVar<3> q_mass_g = MakeConstAD<3>(0.0);
            if (state.has_sw) {
                q_mass_g = ComputePhaseBoundaryMassFlux(
                    coupledSaturationBC, inferredTag, edge.midpoint, area, T_geom, state.Sw, state.k_n, state.mob_g);
            }
            q_flux += q_mass_w * state.h_w + q_mass_g * state.h_g;
        }

        AddBoundaryEqContribution(eqIdx, q_flux, residual, jacobianFull);
        const double diag_add = q_flux.grad(dofOffset);

        if (std::abs(q_flux.val) <= 1.0e-16 &&
            std::abs(q_flux.grad(0)) <= 1.0e-16 &&
            std::abs(q_flux.grad(1)) <= 1.0e-16 &&
            std::abs(q_flux.grad(2)) <= 1.0e-16) {
            ++tagStats.skipped;
            ++stats.zeroEqRows;
        }
        else {
            ++tagStats.applied;
            ++stats.nonzeroEqRows;
            tagStats.sumR += std::abs(q_flux.val);
            tagStats.sumDiag += std::abs(diag_add);
        }

        stats.sumResidual += q_flux.val;
        stats.sumJacobianDiag += diag_add;
        stats.fractureBCCount++;
    }

    return stats;
}

BoundaryAssemblyStats BoundaryAssembler::Assemble_2D_CoupledN3_FullJac(
    MeshManager& mgr,
    FieldManager_2D& fm,
    const BoundarySetting::BoundaryConditionManager* pressureBC,
    const BoundarySetting::BoundaryConditionManager* saturationBC,
    const BoundarySetting::BoundaryConditionManager* temperatureBC,
    int dofOffset_P,
    int dofOffset_W,
    int dofOffset_G,
    int dofOffset_E,
    std::vector<double>& residual,
    std::vector<std::array<double, 3>>& jacobianFull,
    const FluidPropertyEvalConfig& fluid_cfg,
    const CapRelPerm::VGParams& vg,
    const CapRelPerm::RelPermParams& rp)
{
    BoundaryAssemblyStats stats;
    const FIM_Engine::FluidPropertyEvalContext fluid_ctx = BuildFluidEvalContext(fluid_cfg);
    if (residual.size() != jacobianFull.size()) {
        std::cerr << "[BoundaryAssembler::Assemble_2D_CoupledN3_FullJac] residual/jacobianFull size mismatch: "
            << residual.size() << " vs " << jacobianFull.size() << std::endl;
        return stats;
    }
    if (dofOffset_P < 0 || dofOffset_P >= 3 ||
        dofOffset_W < 0 || dofOffset_W >= 3 ||
        dofOffset_G < 0 || dofOffset_G >= 3 ||
        dofOffset_E < 0 || dofOffset_E >= 3) {
        std::cerr << "[BoundaryAssembler::Assemble_2D_CoupledN3_FullJac] Invalid dof offsets." << std::endl;
        return stats;
    }

    const BoundaryFieldTags tags = BuildBoundaryFieldTags();
    const auto& mesh = mgr.mesh();

    for (size_t i = 0; i < mesh.getCells().size(); ++i) {
        const auto& cell = mesh.getCells()[i];
        const int solverIdx = static_cast<int>(i);
        for (int globalFaceID : cell.CellFaceIDs) {
            const int faceIdx = mesh.getFaceIndex(globalFaceID);
            if (faceIdx < 0) continue;
            const auto& face = mesh.getFaces()[faceIdx];
            if (!face.isBoundary()) continue;

            const int tag = face.physicalGroupId;
            const bool hasPressure = pressureBC && pressureBC->HasBC(tag);
            const bool hasTemp = temperatureBC && temperatureBC->HasBC(tag);
            if (!hasPressure && !hasTemp) continue;

            const BoundaryLocalState state = BuildBoundaryLocalState(
                fm, false, static_cast<int>(i), face.normal, fluid_ctx, vg, rp, tags);
            if (!state.valid) {
                ++stats.invalidEqRows;
                continue;
            }

            const Vector d_owner = face.midpoint - cell.center;
            const double dist = std::max(std::abs(d_owner * face.normal), 1.0e-12);
            double area = std::max(face.vectorE.Mag(), 0.0);
            if (area <= 1.0e-20) area = std::max(face.length, 0.0);
            const double T_geom = area / dist;

            ADVar<3> q_w = MakeConstAD<3>(0.0);
            ADVar<3> q_g = MakeConstAD<3>(0.0);
            ADVar<3> q_cond = MakeConstAD<3>(0.0);
            const bool massApplied = ComputePressureDrivenTwoPhaseMassFlux(
                pressureBC, saturationBC, tag, face.midpoint, area, T_geom, state, stats, q_w, q_g);
            const bool condApplied = ComputeTemperatureConductionFlux(
                temperatureBC, tag, face.midpoint, area, T_geom, state, q_cond);
            if (!massApplied && !condApplied) continue;

            ADVar<3> q_energy = q_cond;
            if (massApplied) {
                q_energy += q_w * state.h_w + q_g * state.h_g;
            }

            const int eqIdxW = mgr.getEquationIndex(solverIdx, dofOffset_W);
            const int eqIdxG = mgr.getEquationIndex(solverIdx, dofOffset_G);
            const int eqIdxE = mgr.getEquationIndex(solverIdx, dofOffset_E);
            if (massApplied) {
                AccumulateCoupledEqContribution(eqIdxW, dofOffset_P, q_w, residual, jacobianFull, stats);
                AccumulateCoupledEqContribution(eqIdxG, dofOffset_G, q_g, residual, jacobianFull, stats);
            }
            AccumulateCoupledEqContribution(eqIdxE, dofOffset_E, q_energy, residual, jacobianFull, stats);
            ++stats.matrixBCCount;
        }
    }

    const int nMat = mgr.getMatrixDOFCount();
    const DomainBounds bounds = ComputeDomainBounds(mesh);
    const auto& frNet = mgr.fracture_network();
    const auto& orderedElems = frNet.getOrderedFractureElements();
    for (const FractureElement* elemPtr : orderedElems) {
        if (!elemPtr) continue;
        const FractureElement& elem = *elemPtr;
        const int solverIdx = elem.solverIndex;
        if (solverIdx < nMat) continue;

        const int localIdx = solverIdx - nMat;
        if (localIdx < 0) {
            ++stats.invalidEqRows;
            continue;
        }

        if (elem.parentFractureID < 0 || elem.parentFractureID >= static_cast<int>(frNet.fractures.size())) {
            ++stats.invalidEqRows;
            continue;
        }
        const auto& frac = frNet.fractures[elem.parentFractureID];
        if (elem.id <= 0 || elem.id >= static_cast<int>(frac.intersections.size())) {
            ++stats.invalidEqRows;
            continue;
        }

        const Vector p0 = frac.intersections[elem.id - 1].point;
        const Vector p1 = frac.intersections[elem.id].point;
        Vector tangent = p1 - p0;
        if (tangent.Mag() <= 1.0e-12) tangent = frac.end - frac.start;
        const double tanMag = std::max(tangent.Mag(), 1.0e-12);
        const Vector tanUnit = tangent / tanMag;

        const double dist = std::max(0.5 * std::max(elem.length, 0.0), 1.0e-12);
        const double area = std::max(elem.aperture, 1.0e-12);
        const double T_geom = area / dist;

        for (int endId = 0; endId < 2; ++endId) {
            const bool isStart = (endId == 0);
            const Vector boundaryPoint = isStart ? p0 : p1;
            const Vector outwardNormal = isStart ? (tanUnit * -1.0) : tanUnit;

            const int inferredTag = InferBoundaryTagFromPoint(boundaryPoint, bounds, false);
            if (inferredTag < 0) {
                ++stats.inferredTagFail;
                continue;
            }

            const bool hasPressure = pressureBC && pressureBC->HasBC(inferredTag);
            const bool hasTemp = temperatureBC && temperatureBC->HasBC(inferredTag);
            if (!hasPressure && !hasTemp) continue;

            const BoundaryLocalState state = BuildBoundaryLocalState(
                fm, true, localIdx, outwardNormal, fluid_ctx, vg, rp, tags);
            if (!state.valid) {
                ++stats.invalidEqRows;
                continue;
            }

            ADVar<3> q_w = MakeConstAD<3>(0.0);
            ADVar<3> q_g = MakeConstAD<3>(0.0);
            ADVar<3> q_cond = MakeConstAD<3>(0.0);
            const bool massApplied = ComputePressureDrivenTwoPhaseMassFlux(
                pressureBC, saturationBC, inferredTag, boundaryPoint, area, T_geom, state, stats, q_w, q_g);
            const bool condApplied = ComputeTemperatureConductionFlux(
                temperatureBC, inferredTag, boundaryPoint, area, T_geom, state, q_cond);
            if (!massApplied && !condApplied) continue;

            ADVar<3> q_energy = q_cond;
            if (massApplied) {
                q_energy += q_w * state.h_w + q_g * state.h_g;
            }

            const int eqIdxW = mgr.getEquationIndex(solverIdx, dofOffset_W);
            const int eqIdxG = mgr.getEquationIndex(solverIdx, dofOffset_G);
            const int eqIdxE = mgr.getEquationIndex(solverIdx, dofOffset_E);
            if (massApplied) {
                AccumulateCoupledEqContribution(eqIdxW, dofOffset_P, q_w, residual, jacobianFull, stats);
                AccumulateCoupledEqContribution(eqIdxG, dofOffset_G, q_g, residual, jacobianFull, stats);
            }
            AccumulateCoupledEqContribution(eqIdxE, dofOffset_E, q_energy, residual, jacobianFull, stats);
            ++stats.fractureBCCount;
        }
    }

    return stats;
}

BoundaryAssemblyStats BoundaryAssembler::Assemble_3D_CoupledN3_FullJac(
    MeshManager_3D& mgr,
    FieldManager_3D& fm,
    const BoundarySetting::BoundaryConditionManager* pressureBC,
    const BoundarySetting::BoundaryConditionManager* saturationBC,
    const BoundarySetting::BoundaryConditionManager* temperatureBC,
    int dofOffset_P,
    int dofOffset_W,
    int dofOffset_G,
    int dofOffset_E,
    std::vector<double>& residual,
    std::vector<std::array<double, 3>>& jacobianFull,
    const FluidPropertyEvalConfig& fluid_cfg,
    const CapRelPerm::VGParams& vg,
    const CapRelPerm::RelPermParams& rp)
{
    BoundaryAssemblyStats stats;
    const FIM_Engine::FluidPropertyEvalContext fluid_ctx = BuildFluidEvalContext(fluid_cfg);
    if (residual.size() != jacobianFull.size()) {
        std::cerr << "[BoundaryAssembler::Assemble_3D_CoupledN3_FullJac] residual/jacobianFull size mismatch: "
            << residual.size() << " vs " << jacobianFull.size() << std::endl;
        return stats;
    }
    if (dofOffset_P < 0 || dofOffset_P >= 3 ||
        dofOffset_W < 0 || dofOffset_W >= 3 ||
        dofOffset_G < 0 || dofOffset_G >= 3 ||
        dofOffset_E < 0 || dofOffset_E >= 3) {
        std::cerr << "[BoundaryAssembler::Assemble_3D_CoupledN3_FullJac] Invalid dof offsets." << std::endl;
        return stats;
    }

    const BoundaryFieldTags tags = BuildBoundaryFieldTags();
    const auto& mesh = mgr.mesh();
    for (size_t i = 0; i < mesh.getCells().size(); ++i) {
        const auto& cell = mesh.getCells()[i];
        const int solverIdx = static_cast<int>(i);
        for (int globalFaceID : cell.CellFaceIDs) {
            const int faceIdx = mesh.getFaceIndex(globalFaceID);
            if (faceIdx < 0) continue;
            const auto& face = mesh.getFaces()[faceIdx];
            if (!face.isBoundary()) continue;

            const int tag = face.physicalGroupId;
            const bool hasPressure = pressureBC && pressureBC->HasBC(tag);
            const bool hasTemp = temperatureBC && temperatureBC->HasBC(tag);
            if (!hasPressure && !hasTemp) continue;

            const BoundaryLocalState state = BuildBoundaryLocalState(
                fm, false, static_cast<int>(i), face.normal, fluid_ctx, vg, rp, tags);
            if (!state.valid) {
                ++stats.invalidEqRows;
                continue;
            }

            const Vector d_owner = face.midpoint - cell.center;
            const double dist = std::max(std::abs(d_owner * face.normal), 1.0e-12);
            double area = std::max(face.vectorE.Mag(), 0.0);
            if (area <= 1.0e-20) area = std::max(face.length, 0.0);
            const double T_geom = area / dist;

            ADVar<3> q_w = MakeConstAD<3>(0.0);
            ADVar<3> q_g = MakeConstAD<3>(0.0);
            ADVar<3> q_cond = MakeConstAD<3>(0.0);
            const bool massApplied = ComputePressureDrivenTwoPhaseMassFlux(
                pressureBC, saturationBC, tag, face.midpoint, area, T_geom, state, stats, q_w, q_g);
            const bool condApplied = ComputeTemperatureConductionFlux(
                temperatureBC, tag, face.midpoint, area, T_geom, state, q_cond);
            if (!massApplied && !condApplied) continue;

            ADVar<3> q_energy = q_cond;
            if (massApplied) {
                q_energy += q_w * state.h_w + q_g * state.h_g;
            }

            const int eqIdxW = mgr.getEquationIndex(solverIdx, dofOffset_W);
            const int eqIdxG = mgr.getEquationIndex(solverIdx, dofOffset_G);
            const int eqIdxE = mgr.getEquationIndex(solverIdx, dofOffset_E);
            if (massApplied) {
                AccumulateCoupledEqContribution(eqIdxW, dofOffset_P, q_w, residual, jacobianFull, stats);
                AccumulateCoupledEqContribution(eqIdxG, dofOffset_G, q_g, residual, jacobianFull, stats);
            }
            AccumulateCoupledEqContribution(eqIdxE, dofOffset_E, q_energy, residual, jacobianFull, stats);
            ++stats.matrixBCCount;
        }
    }

    const int nMat = get3DMatrixOffset(mgr);
    const DomainBounds bounds = ComputeDomainBounds(mesh);
    for (const auto& edge : mgr.fracture_network().getGlobalEdges()) {
        if (edge.neighborCell_solverIndex != -1) continue;
        if (edge.ownerCell_solverIndex < nMat) continue;

        const int inferredTag = InferBoundaryTagFromPoint(edge.midpoint, bounds, true);
        if (inferredTag < 0) {
            ++stats.inferredTagFail;
            continue;
        }

        const bool hasPressure = pressureBC && pressureBC->HasBC(inferredTag);
        const bool hasTemp = temperatureBC && temperatureBC->HasBC(inferredTag);
        if (!hasPressure && !hasTemp) continue;

        const int solverIdx = edge.ownerCell_solverIndex;
        const int localIdx = solverIdx - nMat;
        if (localIdx < 0) {
            ++stats.invalidEqRows;
            continue;
        }

        const BoundaryLocalState state = BuildBoundaryLocalState(
            fm, true, localIdx, edge.normal, fluid_ctx, vg, rp, tags);
        if (!state.valid) {
            ++stats.invalidEqRows;
            continue;
        }

        double dist = std::abs(edge.ownerToNeighbor * edge.normal);
        if (dist <= 1.0e-12) dist = edge.ownerToNeighbor.Mag();
        dist = std::max(dist, 1.0e-12);
        const double area = std::max(edge.length, 0.0);
        const double T_geom = area / dist;

        ADVar<3> q_w = MakeConstAD<3>(0.0);
        ADVar<3> q_g = MakeConstAD<3>(0.0);
        ADVar<3> q_cond = MakeConstAD<3>(0.0);
        const bool massApplied = ComputePressureDrivenTwoPhaseMassFlux(
            pressureBC, saturationBC, inferredTag, edge.midpoint, area, T_geom, state, stats, q_w, q_g);
        const bool condApplied = ComputeTemperatureConductionFlux(
            temperatureBC, inferredTag, edge.midpoint, area, T_geom, state, q_cond);
        if (!massApplied && !condApplied) continue;

        ADVar<3> q_energy = q_cond;
        if (massApplied) {
            q_energy += q_w * state.h_w + q_g * state.h_g;
        }

        const int eqIdxW = mgr.getEquationIndex(solverIdx, dofOffset_W);
        const int eqIdxG = mgr.getEquationIndex(solverIdx, dofOffset_G);
        const int eqIdxE = mgr.getEquationIndex(solverIdx, dofOffset_E);
        if (massApplied) {
            AccumulateCoupledEqContribution(eqIdxW, dofOffset_P, q_w, residual, jacobianFull, stats);
            AccumulateCoupledEqContribution(eqIdxG, dofOffset_G, q_g, residual, jacobianFull, stats);
        }
        AccumulateCoupledEqContribution(eqIdxE, dofOffset_E, q_energy, residual, jacobianFull, stats);
        ++stats.fractureBCCount;
    }

    return stats;
}

BoundaryAssemblyStats BoundaryAssembler::Assemble_Wells_2D(
    MeshManager& mgr, FieldManager_2D& fm, const std::vector<WellScheduleStep>& active_steps,
    int dofOffset_P, int dofOffset_W, int dofOffset_G, int dofOffset_E,
    std::vector<double>& residual, std::vector<double>& jacobianDiag)
{
    BoundaryAssemblyStats stats;
    if (residual.size() != jacobianDiag.size() || dofOffset_P < 0 || dofOffset_P >= 3) {
        std::cerr << "[Assemble_Wells_2D] Error: Invalid Vector Size or P-Offset.\n";
        return stats;
    }
    if ((dofOffset_W >= 3) || (dofOffset_G >= 3) || (dofOffset_E >= 3)) {
        std::cerr << "[Assemble_Wells_2D] Error: One or more equation offsets exceed ADVar<3> bounds.\n";
        return stats;
    }

    const auto& mesh = mgr.mesh();
    int nMat = mgr.getMatrixDOFCount();
    const WellFieldTags tags = BuildWellFieldTags();
    const FIM_Engine::FluidPropertyEvalContext fluid_ctx = BuildFluidEvalContext(FluidPropertyEvalConfig());

    for (const auto& step : active_steps) {
        int solverIdx = -1, localIdx = -1;
        ResolveStepCompletionIndex(step, nMat, solverIdx, localIdx);

        std::shared_ptr<volScalarField> pField_P, pField_Kxx, pField_Kyy;
        std::shared_ptr<volScalarField> pField_MobDenW, pField_MobDenG, pField_RhoW, pField_RhoG, pField_LamW, pField_LamG;
        std::shared_ptr<volScalarField> pField_Hw, pField_Hg;

        if (step.domain == WellTargetDomain::Matrix) {
            if (localIdx < 0 || localIdx >= static_cast<int>(mesh.getCells().size())) {
                std::cerr << "[WellSkip][2D][" << step.well_name << "] Invalid Matrix completion_id=" << localIdx << "\n";
                continue;
            }
            pField_P = fm.getMatrixScalar(tags.pressure);
            pField_Kxx = fm.getMatrixScalar(tags.k_xx);
            pField_Kyy = fm.getMatrixScalar(tags.k_yy);
            pField_MobDenW = fm.getMatrixScalar(tags.mob_density_w);
            pField_MobDenG = fm.getMatrixScalar(tags.mob_density_g);
            pField_RhoW = fm.getMatrixScalar(tags.rho_w);
            pField_RhoG = fm.getMatrixScalar(tags.rho_g);
            pField_LamW = fm.getMatrixScalar(tags.lambda_w_mob);
            pField_LamG = fm.getMatrixScalar(tags.lambda_g_mob);
            pField_Hw = fm.getMatrixScalar(tags.h_w);
            pField_Hg = fm.getMatrixScalar(tags.h_g);
        }
        else {
            pField_P = fm.getFractureScalar(tags.pressure);
            pField_MobDenW = fm.getFractureScalar(tags.mob_density_w);
            pField_MobDenG = fm.getFractureScalar(tags.mob_density_g);
            pField_RhoW = fm.getFractureScalar(tags.rho_w);
            pField_RhoG = fm.getFractureScalar(tags.rho_g);
            pField_LamW = fm.getFractureScalar(tags.lambda_w_mob);
            pField_LamG = fm.getFractureScalar(tags.lambda_g_mob);
            pField_Hw = fm.getFractureScalar(tags.h_w);
            pField_Hg = fm.getFractureScalar(tags.h_g);
        }

        if (!isValidLocalIdx(localIdx, pField_P)) {
            std::cerr << "[WellSkip][2D][" << step.well_name << "] Missing pressure field.\n";
            continue;
        }

        double mob_den_w = 0.0, mob_den_g = 0.0;
        if (pField_MobDenW) {
            mob_den_w = safeGetFieldValue(pField_MobDenW, localIdx, 0.0);
        }
        else if (pField_RhoW && pField_LamW) {
            mob_den_w = safeGetFieldValue(pField_RhoW, localIdx) * safeGetFieldValue(pField_LamW, localIdx);
        }
        else {
            std::cerr << "[WellSkip][2D][" << step.well_name << "] Missing water mobility data.\n";
            continue;
        }

        if (pField_MobDenG) {
            mob_den_g = safeGetFieldValue(pField_MobDenG, localIdx, 0.0);
        }
        else if (pField_RhoG && pField_LamG) {
            mob_den_g = safeGetFieldValue(pField_RhoG, localIdx) * safeGetFieldValue(pField_LamG, localIdx);
        }
        else {
            mob_den_g = 0.0;
        }

        double WI = step.wi_override;
        if (WI < 0.0) {
            if (step.domain == WellTargetDomain::Fracture) {
                throw std::runtime_error(std::string("[Assemble_Wells_2D] Fracture well requires explicit wi_override: ") + step.well_name);
            }
            double V = mesh.getCells()[localIdx].volume;
            double kxx = safeGetFieldValue(pField_Kxx, localIdx, 1e-13);
            double kyy = safeGetFieldValue(pField_Kyy, localIdx, 1e-13);
            double k_eff = std::sqrt(kxx * kyy);
            double r_eq = std::sqrt(V / 3.141592653589793);

            if (V <= kEps || step.rw <= 0.0 || r_eq <= step.rw || k_eff <= 0.0) {
                std::cerr << "[WellSkip][2D][" << step.well_name << "] Invalid geom for auto-WI (V=" << V << ", r_eq=" << r_eq << ").\n";
                continue;
            }
            double denom = std::log(r_eq / step.rw) + step.skin;
            if (denom <= kEps) {
                std::cerr << "[WellSkip][2D][" << step.well_name << "] WI denominator <= eps.\n";
                continue;
            }
            WI = (2.0 * 3.141592653589793 * k_eff * 1.0) / denom;
        }

        ADVar<3> P_cell = MakeConstAD<3>(pField_P->data[localIdx]);
        P_cell.grad(dofOffset_P) = 1.0;

        ADVar<3> q_mass_w = MakeConstAD<3>(0.0), q_mass_g = MakeConstAD<3>(0.0);
        const bool has_w_eq = (dofOffset_W >= 0);
        const bool has_g_eq = (dofOffset_G >= 0);
        const bool total_mode = has_w_eq && has_g_eq && (step.component_mode == WellComponentMode::Total);

        if (step.control_mode == WellControlMode::BHP) {
            const ADVar<3> dP = P_cell - MakeConstAD<3>(step.target_value);
            if (total_mode) {
                if (dP.val < 0.0) {
                    const PhaseSplit inj_split = ResolveInjectionPhaseSplit(step, has_w_eq, has_g_eq);
                    const ADVar<3> q_total = MakeConstAD<3>(WI * (mob_den_w + mob_den_g)) * dP;
                    q_mass_w = MakeConstAD<3>(inj_split.fw) * q_total;
                    q_mass_g = MakeConstAD<3>(inj_split.fg) * q_total;
                }
                else {
                    q_mass_w = FVM_Ops::Op_Well_BHP_Source_AD<3, ADVar<3>>(WI * mob_den_w, P_cell, step.target_value);
                    q_mass_g = FVM_Ops::Op_Well_BHP_Source_AD<3, ADVar<3>>(WI * mob_den_g, P_cell, step.target_value);
                }
            }
            else {
                const PhaseSplit split = ResolvePhaseSplit(step, mob_den_w, mob_den_g, has_w_eq, has_g_eq);
                if (split.fw > 0.0) {
                    q_mass_w = FVM_Ops::Op_Well_BHP_Source_AD<3, ADVar<3>>(WI * split.fw * mob_den_w, P_cell, step.target_value);
                }
                if (split.fg > 0.0) {
                    q_mass_g = FVM_Ops::Op_Well_BHP_Source_AD<3, ADVar<3>>(WI * split.fg * mob_den_g, P_cell, step.target_value);
                }
            }
        }
        else if (step.control_mode == WellControlMode::Rate) {
            const StandardConditionDensities std_dens = EvaluateStandardConditionDensities(step, fluid_ctx);
            if (total_mode) {
                PhaseSplit user_split;
                const bool has_user_split = TryResolveUserSpecifiedPhaseSplit(step, user_split);
                if (step.target_value < 0.0) {
                    const PhaseSplit inj_split = ResolveInjectionPhaseSplit(step, has_w_eq, has_g_eq);
                    q_mass_w = BuildMassRateBySplitConst<3>(step, inj_split.fw, std_dens, true);
                    q_mass_g = BuildMassRateBySplitConst<3>(step, inj_split.fg, std_dens, false);
                }
                else if (has_user_split) {
                    q_mass_w = BuildMassRateBySplitConst<3>(step, user_split.fw, std_dens, true);
                    q_mass_g = BuildMassRateBySplitConst<3>(step, user_split.fg, std_dens, false);
                }
                else {
                    const PhaseSplit split = ResolveTotalPhaseSplit(step, mob_den_w, mob_den_g);
                    q_mass_w = BuildMassRateBySplitConst<3>(step, split.fw, std_dens, true);
                    q_mass_g = BuildMassRateBySplitConst<3>(step, split.fg, std_dens, false);
                }
            }
            else {
                const PhaseSplit split = ResolvePhaseSplit(step, mob_den_w, mob_den_g, has_w_eq, has_g_eq);
                if (split.fw > 0.0) {
                    q_mass_w = BuildMassRateBySplitConst<3>(step, split.fw, std_dens, true);
                }
                if (split.fg > 0.0) {
                    q_mass_g = BuildMassRateBySplitConst<3>(step, split.fg, std_dens, false);
                }
            }
        }

        if (dofOffset_W >= 0) {
            int eqIdx_W = mgr.getEquationIndex(solverIdx, dofOffset_W);
            if (isValidEqIdx(eqIdx_W, residual, jacobianDiag)) {
                residual[eqIdx_W] += q_mass_w.val;
                jacobianDiag[eqIdx_W] += q_mass_w.grad(dofOffset_P);
                stats.sumResidual += q_mass_w.val;
                stats.sumJacobianDiag += q_mass_w.grad(dofOffset_P);
            }
        }
        if (dofOffset_G >= 0) {
            int eqIdx_G = mgr.getEquationIndex(solverIdx, dofOffset_G);
            if (isValidEqIdx(eqIdx_G, residual, jacobianDiag)) {
                residual[eqIdx_G] += q_mass_g.val;
                jacobianDiag[eqIdx_G] += q_mass_g.grad(dofOffset_P);
                stats.sumResidual += q_mass_g.val;
                stats.sumJacobianDiag += q_mass_g.grad(dofOffset_P);
            }
        }
        if (dofOffset_E >= 0 && isValidLocalIdx(localIdx, pField_Hw)) {
            double hw = pField_Hw->data[localIdx];
            double hg = safeGetFieldValue(pField_Hg, localIdx, 0.0);
            if (step.injection_temperature > 0.0) {
                const double p_eval = (step.control_mode == WellControlMode::BHP) ? step.target_value : pField_P->data[localIdx];
                if (q_mass_w.val < 0.0) {
                    ADVar<1> P_inj(p_eval), T_inj(step.injection_temperature);
                    hw = FIM_Engine::EvalSinglePhaseFluid<1>(fluid_ctx, P_inj, T_inj).h.val;
                }
                if (q_mass_g.val < 0.0) {
                    ADVar<1> P_inj(p_eval), T_inj(step.injection_temperature);
                    hg = FIM_Engine::EvalTwoPhaseGasFluid<1>(fluid_ctx, P_inj, T_inj).h.val;
                }
            }
            ADVar<3> q_energy = FVM_Ops::Op_Well_Energy_Source_AD<3, ADVar<3>>(q_mass_w, hw, q_mass_g, hg);
            int eqIdx_E = mgr.getEquationIndex(solverIdx, dofOffset_E);
            if (isValidEqIdx(eqIdx_E, residual, jacobianDiag)) {
                residual[eqIdx_E] += q_energy.val;
                jacobianDiag[eqIdx_E] += q_energy.grad(dofOffset_P);
                stats.sumResidual += q_energy.val;
                stats.sumJacobianDiag += q_energy.grad(dofOffset_P);
            }
        }

        if (step.domain == WellTargetDomain::Matrix) stats.matrixBCCount++;
        else stats.fractureBCCount++;
    }
    return stats;
}
BoundaryAssemblyStats BoundaryAssembler::Assemble_Wells_3D(
    MeshManager_3D& mgr, FieldManager_3D& fm, const std::vector<WellScheduleStep>& active_steps,
    int dofOffset_P, int dofOffset_W, int dofOffset_G, int dofOffset_E,
    std::vector<double>& residual, std::vector<double>& jacobianDiag)
{
    BoundaryAssemblyStats stats;
    if (residual.size() != jacobianDiag.size() || dofOffset_P < 0 || dofOffset_P >= 3) {
        std::cerr << "[Assemble_Wells_3D] Error: Invalid Vector Size or P-Offset.\n";
        return stats;
    }
    if ((dofOffset_W >= 3) || (dofOffset_G >= 3) || (dofOffset_E >= 3)) {
        std::cerr << "[Assemble_Wells_3D] Error: One or more equation offsets exceed ADVar<3> bounds.\n";
        return stats;
    }

    const auto& mesh = mgr.mesh();
    int nMat = get3DMatrixOffset(mgr);
    const WellFieldTags tags = BuildWellFieldTags();
    const FIM_Engine::FluidPropertyEvalContext fluid_ctx = BuildFluidEvalContext(FluidPropertyEvalConfig());

    for (const auto& step : active_steps) {
        int solverIdx = -1, localIdx = -1;
        ResolveStepCompletionIndex(step, nMat, solverIdx, localIdx);

        std::shared_ptr<volScalarField> pField_P, pField_Kxx, pField_Kyy, pField_Kzz;
        std::shared_ptr<volScalarField> pField_MobDenW, pField_MobDenG, pField_RhoW, pField_RhoG, pField_LamW, pField_LamG;
        std::shared_ptr<volScalarField> pField_Hw, pField_Hg;

        if (step.domain == WellTargetDomain::Matrix) {
            if (localIdx < 0 || localIdx >= static_cast<int>(mesh.getCells().size())) {
                std::cerr << "[WellSkip][3D][" << step.well_name << "] Invalid Matrix completion_id.\n";
                continue;
            }
            pField_P = fm.getMatrixScalar(tags.pressure);
            pField_Kxx = fm.getMatrixScalar(tags.k_xx);
            pField_Kyy = fm.getMatrixScalar(tags.k_yy);
            pField_Kzz = fm.getMatrixScalar(tags.k_zz);
            pField_MobDenW = fm.getMatrixScalar(tags.mob_density_w);
            pField_MobDenG = fm.getMatrixScalar(tags.mob_density_g);
            pField_RhoW = fm.getMatrixScalar(tags.rho_w);
            pField_RhoG = fm.getMatrixScalar(tags.rho_g);
            pField_LamW = fm.getMatrixScalar(tags.lambda_w_mob);
            pField_LamG = fm.getMatrixScalar(tags.lambda_g_mob);
            pField_Hw = fm.getMatrixScalar(tags.h_w);
            pField_Hg = fm.getMatrixScalar(tags.h_g);
        }
        else {
            pField_P = fm.getFractureScalar(tags.pressure);
            pField_MobDenW = fm.getFractureScalar(tags.mob_density_w);
            pField_MobDenG = fm.getFractureScalar(tags.mob_density_g);
            pField_RhoW = fm.getFractureScalar(tags.rho_w);
            pField_RhoG = fm.getFractureScalar(tags.rho_g);
            pField_LamW = fm.getFractureScalar(tags.lambda_w_mob);
            pField_LamG = fm.getFractureScalar(tags.lambda_g_mob);
            pField_Hw = fm.getFractureScalar(tags.h_w);
            pField_Hg = fm.getFractureScalar(tags.h_g);
        }

        if (!isValidLocalIdx(localIdx, pField_P)) {
            std::cerr << "[WellSkip][3D][" << step.well_name << "] Missing pressure field.\n";
            continue;
        }

        double mob_den_w = 0.0, mob_den_g = 0.0;
        if (pField_MobDenW) {
            mob_den_w = safeGetFieldValue(pField_MobDenW, localIdx, 0.0);
        }
        else if (pField_RhoW && pField_LamW) {
            mob_den_w = safeGetFieldValue(pField_RhoW, localIdx) * safeGetFieldValue(pField_LamW, localIdx);
        }
        else {
            std::cerr << "[WellSkip][3D][" << step.well_name << "] Missing mobility data.\n";
            continue;
        }

        if (pField_MobDenG) {
            mob_den_g = safeGetFieldValue(pField_MobDenG, localIdx, 0.0);
        }
        else if (pField_RhoG && pField_LamG) {
            mob_den_g = safeGetFieldValue(pField_RhoG, localIdx) * safeGetFieldValue(pField_LamG, localIdx);
        }

        double WI = step.wi_override;
        if (WI < 0.0) {
            if (step.domain == WellTargetDomain::Fracture) {
                throw std::runtime_error(std::string("[Assemble_Wells_3D] Fracture well requires explicit wi_override: ") + step.well_name);
            }
            const auto& cell = mesh.getCells()[localIdx];
            double V = cell.volume;

            double L = step.L_override;
            WellAxis axis = step.well_axis;

            if (L < 0.0 && axis == WellAxis::None) {
                std::cerr << "[Assemble_Wells_3D] Warning: well_axis=None with L_override<0. Defaulting to Z-axis.\n";
                axis = WellAxis::Z;
            }

            if (L < 0.0) {
                double min_ax = 1e20, max_ax = -1e20;
                for (int nodeID : cell.CellNodeIDs) {
                    const auto& pt = mesh.getNodesMap().at(nodeID).coord;
                    double val = (axis == WellAxis::X) ? pt.m_x : (axis == WellAxis::Y) ? pt.m_y : pt.m_z;
                    if (val < min_ax) min_ax = val;
                    if (val > max_ax) max_ax = val;
                }
                L = max_ax - min_ax;
            }

            if (V <= kEps || step.rw <= 0.0 || L <= kEps) {
                std::cerr << "[WellSkip][3D][" << step.well_name << "] Invalid geometry (V/L<=eps).\n";
                continue;
            }

            double r_eq = std::sqrt(V / (3.141592653589793 * L));
            double kxx = safeGetFieldValue(pField_Kxx, localIdx, 1e-13);
            double kyy = safeGetFieldValue(pField_Kyy, localIdx, 1e-13);
            double kzz = safeGetFieldValue(pField_Kzz, localIdx, 1e-13);

            double k_plane = 1e-13;
            if (axis == WellAxis::X) k_plane = std::sqrt(kyy * kzz);
            else if (axis == WellAxis::Y) k_plane = std::sqrt(kxx * kzz);
            else k_plane = std::sqrt(kxx * kyy);

            if (r_eq <= step.rw || k_plane <= 0.0) {
                std::cerr << "[WellSkip][3D][" << step.well_name << "] Auto-WI degenerated (r_eq <= rw or k<=0).\n";
                continue;
            }

            double denom = std::log(r_eq / step.rw) + step.skin;
            if (denom <= kEps) {
                std::cerr << "[WellSkip][3D][" << step.well_name << "] WI denominator <= eps.\n";
                continue;
            }
            WI = (2.0 * 3.141592653589793 * k_plane * L) / denom;
        }

        ADVar<3> P_cell = MakeConstAD<3>(pField_P->data[localIdx]);
        P_cell.grad(dofOffset_P) = 1.0;

        ADVar<3> q_mass_w = MakeConstAD<3>(0.0), q_mass_g = MakeConstAD<3>(0.0);
        const bool has_w_eq = (dofOffset_W >= 0);
        const bool has_g_eq = (dofOffset_G >= 0);
        const bool total_mode = has_w_eq && has_g_eq && (step.component_mode == WellComponentMode::Total);

        if (step.control_mode == WellControlMode::BHP) {
            const ADVar<3> dP = P_cell - MakeConstAD<3>(step.target_value);
            if (total_mode) {
                if (dP.val < 0.0) {
                    const PhaseSplit inj_split = ResolveInjectionPhaseSplit(step, has_w_eq, has_g_eq);
                    const ADVar<3> q_total = MakeConstAD<3>(WI * (mob_den_w + mob_den_g)) * dP;
                    q_mass_w = MakeConstAD<3>(inj_split.fw) * q_total;
                    q_mass_g = MakeConstAD<3>(inj_split.fg) * q_total;
                }
                else {
                    q_mass_w = FVM_Ops::Op_Well_BHP_Source_AD<3, ADVar<3>>(WI * mob_den_w, P_cell, step.target_value);
                    q_mass_g = FVM_Ops::Op_Well_BHP_Source_AD<3, ADVar<3>>(WI * mob_den_g, P_cell, step.target_value);
                }
            }
            else {
                const PhaseSplit split = ResolvePhaseSplit(step, mob_den_w, mob_den_g, has_w_eq, has_g_eq);
                if (split.fw > 0.0) {
                    q_mass_w = FVM_Ops::Op_Well_BHP_Source_AD<3, ADVar<3>>(WI * split.fw * mob_den_w, P_cell, step.target_value);
                }
                if (split.fg > 0.0) {
                    q_mass_g = FVM_Ops::Op_Well_BHP_Source_AD<3, ADVar<3>>(WI * split.fg * mob_den_g, P_cell, step.target_value);
                }
            }
        }
        else if (step.control_mode == WellControlMode::Rate) {
            const StandardConditionDensities std_dens = EvaluateStandardConditionDensities(step, fluid_ctx);
            if (total_mode) {
                PhaseSplit user_split;
                const bool has_user_split = TryResolveUserSpecifiedPhaseSplit(step, user_split);
                if (step.target_value < 0.0) {
                    const PhaseSplit inj_split = ResolveInjectionPhaseSplit(step, has_w_eq, has_g_eq);
                    q_mass_w = BuildMassRateBySplitConst<3>(step, inj_split.fw, std_dens, true);
                    q_mass_g = BuildMassRateBySplitConst<3>(step, inj_split.fg, std_dens, false);
                }
                else if (has_user_split) {
                    q_mass_w = BuildMassRateBySplitConst<3>(step, user_split.fw, std_dens, true);
                    q_mass_g = BuildMassRateBySplitConst<3>(step, user_split.fg, std_dens, false);
                }
                else {
                    const PhaseSplit split = ResolveTotalPhaseSplit(step, mob_den_w, mob_den_g);
                    q_mass_w = BuildMassRateBySplitConst<3>(step, split.fw, std_dens, true);
                    q_mass_g = BuildMassRateBySplitConst<3>(step, split.fg, std_dens, false);
                }
            }
            else {
                const PhaseSplit split = ResolvePhaseSplit(step, mob_den_w, mob_den_g, has_w_eq, has_g_eq);
                if (split.fw > 0.0) {
                    q_mass_w = BuildMassRateBySplitConst<3>(step, split.fw, std_dens, true);
                }
                if (split.fg > 0.0) {
                    q_mass_g = BuildMassRateBySplitConst<3>(step, split.fg, std_dens, false);
                }
            }
        }

        if (dofOffset_W >= 0) {
            int eqIdx_W = mgr.getEquationIndex(solverIdx, dofOffset_W);
            if (isValidEqIdx(eqIdx_W, residual, jacobianDiag)) {
                residual[eqIdx_W] += q_mass_w.val;
                jacobianDiag[eqIdx_W] += q_mass_w.grad(dofOffset_P);
                stats.sumResidual += q_mass_w.val;
                stats.sumJacobianDiag += q_mass_w.grad(dofOffset_P);
            }
        }
        if (dofOffset_G >= 0) {
            int eqIdx_G = mgr.getEquationIndex(solverIdx, dofOffset_G);
            if (isValidEqIdx(eqIdx_G, residual, jacobianDiag)) {
                residual[eqIdx_G] += q_mass_g.val;
                jacobianDiag[eqIdx_G] += q_mass_g.grad(dofOffset_P);
                stats.sumResidual += q_mass_g.val;
                stats.sumJacobianDiag += q_mass_g.grad(dofOffset_P);
            }
        }
        if (dofOffset_E >= 0 && isValidLocalIdx(localIdx, pField_Hw)) {
            double hw = pField_Hw->data[localIdx];
            double hg = safeGetFieldValue(pField_Hg, localIdx, 0.0);
            if (step.injection_temperature > 0.0) {
                const double p_eval = (step.control_mode == WellControlMode::BHP) ? step.target_value : pField_P->data[localIdx];
                if (q_mass_w.val < 0.0) {
                    ADVar<1> P_inj(p_eval), T_inj(step.injection_temperature);
                    hw = FIM_Engine::EvalSinglePhaseFluid<1>(fluid_ctx, P_inj, T_inj).h.val;
                }
                if (q_mass_g.val < 0.0) {
                    ADVar<1> P_inj(p_eval), T_inj(step.injection_temperature);
                    hg = FIM_Engine::EvalTwoPhaseGasFluid<1>(fluid_ctx, P_inj, T_inj).h.val;
                }
            }
            ADVar<3> q_energy = FVM_Ops::Op_Well_Energy_Source_AD<3, ADVar<3>>(q_mass_w, hw, q_mass_g, hg);
            int eqIdx_E = mgr.getEquationIndex(solverIdx, dofOffset_E);
            if (isValidEqIdx(eqIdx_E, residual, jacobianDiag)) {
                residual[eqIdx_E] += q_energy.val;
                jacobianDiag[eqIdx_E] += q_energy.grad(dofOffset_P);
                stats.sumResidual += q_energy.val;
                stats.sumJacobianDiag += q_energy.grad(dofOffset_P);
            }
        }

        if (step.domain == WellTargetDomain::Matrix) stats.matrixBCCount++;
        else stats.fractureBCCount++;
    }
    return stats;
}
namespace {
    inline ADVar<3> MakeConstAD3(double v) {
        ADVar<3> x;
        x.val = v;
        for (int k = 0; k < 3; ++k) x.grad(k) = 0.0;
        return x;
    }

    inline void AddWellEqContribution(
        int eqIdx,
        const ADVar<3>& q,
        std::vector<std::array<double, 3>>& jacobianFull,
        std::vector<double>& residualRef)
    {
        if (eqIdx < 0 || eqIdx >= static_cast<int>(residualRef.size()) || eqIdx >= static_cast<int>(jacobianFull.size())) {
            return;
        }
        residualRef[eqIdx] += q.val;
        jacobianFull[eqIdx][0] += q.grad(0);
        jacobianFull[eqIdx][1] += q.grad(1);
        jacobianFull[eqIdx][2] += q.grad(2);
    }

    inline void FillWellEqLinearization(
        WellEquationLinearization& out,
        int eq_idx,
        int eq_dof,
        const ADVar<3>& q,
        double dq_dPbh)
    {
        out.valid = (eq_idx >= 0);
        out.eq_index = eq_idx;
        out.eq_dof = eq_dof;
        out.q = q.val;
        out.dq_dcell = { { q.grad(0), q.grad(1), q.grad(2) } };
        out.dq_dPbh = dq_dPbh;
    }
}

BoundaryAssemblyStats BoundaryAssembler::Assemble_Wells_2D_FullJac(
    MeshManager& mgr, FieldManager_2D& fm, const std::vector<WellScheduleStep>& active_steps,
    int dofOffset_P, int dofOffset_W, int dofOffset_G, int dofOffset_E,
    std::vector<double>& residual, std::vector<std::array<double, 3>>& jacobianFull,
    const FluidPropertyEvalConfig& fluid_cfg, const CapRelPerm::VGParams& vg, const CapRelPerm::RelPermParams& rp,
    std::vector<WellCompletionLinearization>* completionLinearizations,
    const std::unordered_map<std::string, double>* wellBhpByName)
{
    BoundaryAssemblyStats stats;
    if (residual.size() != jacobianFull.size() || dofOffset_P < 0 || dofOffset_P >= 3) {
        std::cerr << "[Assemble_Wells_2D_FullJac] Error: Invalid vector size or P-offset.\n";
        return stats;
    }
    if ((dofOffset_W >= 3) || (dofOffset_G >= 3) || (dofOffset_E >= 3)) {
        std::cerr << "[Assemble_Wells_2D_FullJac] Error: One or more equation offsets exceed ADVar<3> bounds.\n";
        return stats;
    }

    const auto& mesh = mgr.mesh();
    const int nMat = mgr.getMatrixDOFCount();
    const WellFieldTags tags = BuildWellFieldTags();
    const FIM_Engine::FluidPropertyEvalContext fluid_ctx = BuildFluidEvalContext(fluid_cfg);
    if (completionLinearizations) completionLinearizations->clear();

    for (size_t step_idx = 0; step_idx < active_steps.size(); ++step_idx) {
        const auto& step = active_steps[step_idx];
        int solverIdx = -1, localIdx = -1;
        ResolveStepCompletionIndex(step, nMat, solverIdx, localIdx);

        std::shared_ptr<volScalarField> pField_P, pField_T, pField_Sw;
        std::shared_ptr<volScalarField> pField_Kxx, pField_Kyy;

        if (step.domain == WellTargetDomain::Matrix) {
            if (localIdx < 0 || localIdx >= static_cast<int>(mesh.getCells().size())) {
                throw std::runtime_error("[Assemble_Wells_2D_FullJac] Matrix completion index out of range.");
            }
            pField_P = fm.getMatrixScalar(tags.pressure);
            pField_T = fm.getMatrixScalar(tags.temperature);
            pField_Sw = fm.getMatrixScalar(tags.saturation);
            pField_Kxx = fm.getMatrixScalar(tags.k_xx);
            pField_Kyy = fm.getMatrixScalar(tags.k_yy);
        }
        else {
            pField_P = fm.getFractureScalar(tags.pressure);
            pField_T = fm.getFractureScalar(tags.temperature);
            pField_Sw = fm.getFractureScalar(tags.saturation);
        }

        if (!isValidLocalIdx(localIdx, pField_P)) {
            throw std::runtime_error("[Assemble_Wells_2D_FullJac] Completion points to invalid pressure field index.");
        }

        double WI = step.wi_override;
        if (WI < 0.0) {
            if (step.domain == WellTargetDomain::Fracture) {
                throw std::runtime_error(std::string("[Assemble_Wells_2D_FullJac] Fracture well requires explicit wi_override: ") + step.well_name);
            }

            const double V = mesh.getCells()[localIdx].volume;
            const double kxx = safeGetFieldValue(pField_Kxx, localIdx, 1e-13);
            const double kyy = safeGetFieldValue(pField_Kyy, localIdx, 1e-13);
            const double k_eff = std::sqrt(kxx * kyy);
            const double r_eq = std::sqrt(V / 3.141592653589793);
            if (V <= kEps || step.rw <= 0.0 || r_eq <= step.rw || k_eff <= 0.0) {
                throw std::runtime_error("[Assemble_Wells_2D_FullJac] Invalid geometry/properties for automatic WI evaluation.");
            }

            const double denom = std::log(r_eq / step.rw) + step.skin;
            if (denom <= kEps) {
                throw std::runtime_error("[Assemble_Wells_2D_FullJac] Invalid WI denominator for automatic WI evaluation.");
            }
            WI = (2.0 * 3.141592653589793 * k_eff * 1.0) / denom;
        }

        ADVar<3> P_cell = MakeConstAD3(pField_P->data[localIdx]);
        P_cell.grad(0) = 1.0;

        ADVar<3> T_cell = MakeConstAD3(safeGetFieldValue(pField_T, localIdx, 300.0));
        T_cell.grad(2) = 1.0;

        const bool has_sw = isValidLocalIdx(localIdx, pField_Sw);
        ADVar<3> Sw_cell = MakeConstAD3(has_sw ? pField_Sw->data[localIdx] : 1.0);
        if (has_sw) Sw_cell.grad(1) = 1.0;

        ADVar<3> Sw_const = Sw_cell;
        ADVar<3> Pc_cell = MakeConstAD3(0.0);
        ADVar<3> krw = MakeConstAD3(1.0);
        ADVar<3> krg = MakeConstAD3(0.0);
        if (has_sw) {
            if (step.domain == WellTargetDomain::Fracture) {
                if (Sw_const.val < 0.0) Sw_const = MakeConstAD3(0.0);
                else if (Sw_const.val > 1.0) Sw_const = MakeConstAD3(1.0);
                krw = Sw_const;
                krg = MakeConstAD3(1.0) - Sw_const;
            }
            else {
                Sw_const = ClampSwForConstitutiveWell<3>(Sw_cell, vg);
                Pc_cell = CapRelPerm::pc_vG<3>(Sw_const, vg);
                CapRelPerm::kr_Mualem_vG<3>(Sw_const, vg, rp, krw, krg);
            }
        }
        else if (fluid_ctx.config.single_phase_is_co2) {
            krw = MakeConstAD3(0.0);
            krg = MakeConstAD3(1.0);
        }

        AD_Fluid::ADFluidProperties<3> propsW{};
        AD_Fluid::ADFluidProperties<3> propsG{};
        if (has_sw) {
            propsW = FIM_Engine::EvalTwoPhaseWaterFluid<3>(fluid_ctx, P_cell, T_cell);
            propsG = FIM_Engine::EvalTwoPhaseGasFluid<3>(fluid_ctx, P_cell + Pc_cell, T_cell);
        }
        else {
            propsW = FIM_Engine::EvalSinglePhaseFluid<3>(fluid_ctx, P_cell, T_cell);
            propsG = fluid_ctx.config.single_phase_is_co2
                ? propsW
                : FIM_Engine::EvalTwoPhaseGasFluid<3>(fluid_ctx, P_cell, T_cell);
        }

        ADVar<3> mob_w = (propsW.rho * krw) / propsW.mu;
        ADVar<3> mob_g = (propsG.rho * krg) / propsG.mu;

        const bool has_w_eq = (dofOffset_W >= 0);
        const bool has_g_eq = (dofOffset_G >= 0);
        const bool total_mode = has_w_eq && has_g_eq && (step.component_mode == WellComponentMode::Total);

        const double pbh_eval = ResolveWellBhpValue(step, P_cell.val, wellBhpByName);
        ADVar<3> dP = P_cell - MakeConstAD3(pbh_eval);

        ADVar<3> q_mass_w = MakeConstAD3(0.0);
        ADVar<3> q_mass_g = MakeConstAD3(0.0);
        double dq_mass_w_dPbh = 0.0;
        double dq_mass_g_dPbh = 0.0;

        if (total_mode) {
            if (dP.val < 0.0) {
                const PhaseSplit inj_split = ResolveInjectionPhaseSplit(step, has_w_eq, has_g_eq);
                ADVar<3> q_total = MakeConstAD3(WI) * (mob_w + mob_g) * dP;
                q_mass_w = MakeConstAD3(inj_split.fw) * q_total;
                q_mass_g = MakeConstAD3(inj_split.fg) * q_total;
                const double dq_total_dPbh = -WI * (mob_w.val + mob_g.val);
                dq_mass_w_dPbh = inj_split.fw * dq_total_dPbh;
                dq_mass_g_dPbh = inj_split.fg * dq_total_dPbh;
            }
            else {
                q_mass_w = MakeConstAD3(WI) * mob_w * dP;
                q_mass_g = MakeConstAD3(WI) * mob_g * dP;
                dq_mass_w_dPbh = -WI * mob_w.val;
                dq_mass_g_dPbh = -WI * mob_g.val;
            }
        }
        else {
            const PhaseSplit split = ResolvePhaseSplit(step, mob_w.val, mob_g.val, has_w_eq, has_g_eq);
            if (split.fw > 0.0) {
                q_mass_w = MakeConstAD3(WI * split.fw) * mob_w * dP;
                dq_mass_w_dPbh = -WI * split.fw * mob_w.val;
            }
            if (split.fg > 0.0) {
                q_mass_g = MakeConstAD3(WI * split.fg) * mob_g * dP;
                dq_mass_g_dPbh = -WI * split.fg * mob_g.val;
            }
        }

        ADVar<3> h_w_well = propsW.h;
        ADVar<3> h_g_well = propsG.h;
        if (step.injection_temperature > 0.0) {
            ADVar<3> T_inj = MakeConstAD3(step.injection_temperature);
            ADVar<3> P_inj_const = MakeConstAD3(pbh_eval);
            if (q_mass_w.val < 0.0) {
                h_w_well = has_sw
                    ? FIM_Engine::EvalTwoPhaseWaterFluid<3>(fluid_ctx, P_inj_const, T_inj).h
                    : FIM_Engine::EvalSinglePhaseFluid<3>(fluid_ctx, P_inj_const, T_inj).h;
            }
            if (q_mass_g.val < 0.0) {
                h_g_well = (has_sw || !fluid_ctx.config.single_phase_is_co2)
                    ? FIM_Engine::EvalTwoPhaseGasFluid<3>(fluid_ctx, P_inj_const, T_inj).h
                    : FIM_Engine::EvalSinglePhaseFluid<3>(fluid_ctx, P_inj_const, T_inj).h;
            }
        }

        ADVar<3> q_total = q_mass_w + q_mass_g;
        ADVar<3> q_energy = (q_mass_w * h_w_well) + (q_mass_g * h_g_well);
        const double dq_total_dPbh = dq_mass_w_dPbh + dq_mass_g_dPbh;
        const double dq_energy_dPbh = dq_mass_w_dPbh * h_w_well.val + dq_mass_g_dPbh * h_g_well.val;

        const int eqIdx_W = (dofOffset_W >= 0) ? mgr.getEquationIndex(solverIdx, dofOffset_W) : -1;
        const int eqIdx_G = (dofOffset_G >= 0) ? mgr.getEquationIndex(solverIdx, dofOffset_G) : -1;
        const int eqIdx_E = (dofOffset_E >= 0) ? mgr.getEquationIndex(solverIdx, dofOffset_E) : -1;

        if (dofOffset_W >= 0) AddWellEqContribution(eqIdx_W, q_mass_w, jacobianFull, residual);
        if (dofOffset_G >= 0) AddWellEqContribution(eqIdx_G, q_mass_g, jacobianFull, residual);
        if (dofOffset_E >= 0) AddWellEqContribution(eqIdx_E, q_energy, jacobianFull, residual);

        if (completionLinearizations) {
            WellCompletionLinearization rec;
            rec.step_index = step_idx;
            rec.well_name = step.well_name;
            rec.completion_solver_index = (step.completion_solver_index >= 0) ? step.completion_solver_index : solverIdx;
            rec.domain = step.domain;
            rec.pbh_eval = pbh_eval;
            rec.q_total = q_total.val;
            rec.dqtotal_dcell = { { q_total.grad(0), q_total.grad(1), q_total.grad(2) } };
            rec.dqtotal_dPbh = dq_total_dPbh;
            FillWellEqLinearization(rec.water, eqIdx_W, dofOffset_W, q_mass_w, dq_mass_w_dPbh);
            FillWellEqLinearization(rec.gas, eqIdx_G, dofOffset_G, q_mass_g, dq_mass_g_dPbh);
            FillWellEqLinearization(rec.energy, eqIdx_E, dofOffset_E, q_energy, dq_energy_dPbh);
            completionLinearizations->push_back(rec);
        }

        if (step.domain == WellTargetDomain::Matrix) stats.matrixBCCount++;
        else stats.fractureBCCount++;
    }

    return stats;
}
BoundaryAssemblyStats BoundaryAssembler::Assemble_Wells_3D_FullJac(
    MeshManager_3D& mgr, FieldManager_3D& fm, const std::vector<WellScheduleStep>& active_steps,
    int dofOffset_P, int dofOffset_W, int dofOffset_G, int dofOffset_E,
    std::vector<double>& residual, std::vector<std::array<double, 3>>& jacobianFull,
    const FluidPropertyEvalConfig& fluid_cfg, const CapRelPerm::VGParams& vg, const CapRelPerm::RelPermParams& rp,
    std::vector<WellCompletionLinearization>* completionLinearizations,
    const std::unordered_map<std::string, double>* wellBhpByName)
{
    BoundaryAssemblyStats stats;
    if (residual.size() != jacobianFull.size() || dofOffset_P < 0 || dofOffset_P >= 3) {
        std::cerr << "[Assemble_Wells_3D_FullJac] Error: Invalid vector size or P-offset.\n";
        return stats;
    }
    if ((dofOffset_W >= 3) || (dofOffset_G >= 3) || (dofOffset_E >= 3)) {
        std::cerr << "[Assemble_Wells_3D_FullJac] Error: One or more equation offsets exceed ADVar<3> bounds.\n";
        return stats;
    }

    const auto& mesh = mgr.mesh();
    const int nMat = get3DMatrixOffset(mgr);
    const WellFieldTags tags = BuildWellFieldTags();
    const FIM_Engine::FluidPropertyEvalContext fluid_ctx = BuildFluidEvalContext(fluid_cfg);
    if (completionLinearizations) completionLinearizations->clear();

    for (size_t step_idx = 0; step_idx < active_steps.size(); ++step_idx) {
        const auto& step = active_steps[step_idx];
        int solverIdx = -1, localIdx = -1;
        ResolveStepCompletionIndex(step, nMat, solverIdx, localIdx);

        std::shared_ptr<volScalarField> pField_P, pField_T, pField_Sw;
        std::shared_ptr<volScalarField> pField_Kxx, pField_Kyy, pField_Kzz;

        if (step.domain == WellTargetDomain::Matrix) {
            if (localIdx < 0 || localIdx >= static_cast<int>(mesh.getCells().size())) {
                throw std::runtime_error("[Assemble_Wells_3D_FullJac] Matrix completion index out of range.");
            }
            pField_P = fm.getMatrixScalar(tags.pressure);
            pField_T = fm.getMatrixScalar(tags.temperature);
            pField_Sw = fm.getMatrixScalar(tags.saturation);
            pField_Kxx = fm.getMatrixScalar(tags.k_xx);
            pField_Kyy = fm.getMatrixScalar(tags.k_yy);
            pField_Kzz = fm.getMatrixScalar(tags.k_zz);
        }
        else {
            pField_P = fm.getFractureScalar(tags.pressure);
            pField_T = fm.getFractureScalar(tags.temperature);
            pField_Sw = fm.getFractureScalar(tags.saturation);
        }

        if (!isValidLocalIdx(localIdx, pField_P)) {
            throw std::runtime_error("[Assemble_Wells_3D_FullJac] Completion points to invalid pressure field index.");
        }

        double WI = step.wi_override;
        if (WI < 0.0) {
            if (step.domain == WellTargetDomain::Fracture) {
                throw std::runtime_error(std::string("[Assemble_Wells_3D_FullJac] Fracture well requires explicit wi_override: ") + step.well_name);
            }

            const auto& cell = mesh.getCells()[localIdx];
            const double V = cell.volume;

            double L = step.L_override;
            WellAxis axis = step.well_axis;
            if (L < 0.0 && axis == WellAxis::None) axis = WellAxis::Z;

            if (L < 0.0) {
                double min_ax = 1e20, max_ax = -1e20;
                for (int nodeID : cell.CellNodeIDs) {
                    const auto& pt = mesh.getNodesMap().at(nodeID).coord;
                    const double val = (axis == WellAxis::X) ? pt.m_x : (axis == WellAxis::Y) ? pt.m_y : pt.m_z;
                    if (val < min_ax) min_ax = val;
                    if (val > max_ax) max_ax = val;
                }
                L = max_ax - min_ax;
            }

            if (V <= kEps || step.rw <= 0.0 || L <= kEps) {
                throw std::runtime_error("[Assemble_Wells_3D_FullJac] Invalid geometry/properties for automatic WI evaluation.");
            }

            const double r_eq = std::sqrt(V / (3.141592653589793 * L));
            const double kxx = safeGetFieldValue(pField_Kxx, localIdx, 1e-13);
            const double kyy = safeGetFieldValue(pField_Kyy, localIdx, 1e-13);
            const double kzz = safeGetFieldValue(pField_Kzz, localIdx, 1e-13);

            double k_plane = 1e-13;
            if (axis == WellAxis::X) k_plane = std::sqrt(kyy * kzz);
            else if (axis == WellAxis::Y) k_plane = std::sqrt(kxx * kzz);
            else k_plane = std::sqrt(kxx * kyy);

            if (r_eq <= step.rw || k_plane <= 0.0) {
                throw std::runtime_error("[Assemble_Wells_3D_FullJac] Degenerated WI geometry in automatic WI evaluation.");
            }

            const double denom = std::log(r_eq / step.rw) + step.skin;
            if (denom <= kEps) {
                throw std::runtime_error("[Assemble_Wells_3D_FullJac] Invalid WI denominator for automatic WI evaluation.");
            }
            WI = (2.0 * 3.141592653589793 * k_plane * L) / denom;
        }

        ADVar<3> P_cell = MakeConstAD3(pField_P->data[localIdx]);
        P_cell.grad(0) = 1.0;

        ADVar<3> T_cell = MakeConstAD3(safeGetFieldValue(pField_T, localIdx, 300.0));
        T_cell.grad(2) = 1.0;

        const bool has_sw = isValidLocalIdx(localIdx, pField_Sw);
        ADVar<3> Sw_cell = MakeConstAD3(has_sw ? pField_Sw->data[localIdx] : 1.0);
        if (has_sw) Sw_cell.grad(1) = 1.0;

        ADVar<3> Sw_const = Sw_cell;
        ADVar<3> Pc_cell = MakeConstAD3(0.0);
        ADVar<3> krw = MakeConstAD3(1.0);
        ADVar<3> krg = MakeConstAD3(0.0);
        if (has_sw) {
            if (step.domain == WellTargetDomain::Fracture) {
                if (Sw_const.val < 0.0) Sw_const = MakeConstAD3(0.0);
                else if (Sw_const.val > 1.0) Sw_const = MakeConstAD3(1.0);
                krw = Sw_const;
                krg = MakeConstAD3(1.0) - Sw_const;
            }
            else {
                Sw_const = ClampSwForConstitutiveWell<3>(Sw_cell, vg);
                Pc_cell = CapRelPerm::pc_vG<3>(Sw_const, vg);
                CapRelPerm::kr_Mualem_vG<3>(Sw_const, vg, rp, krw, krg);
            }
        }
        else if (fluid_ctx.config.single_phase_is_co2) {
            krw = MakeConstAD3(0.0);
            krg = MakeConstAD3(1.0);
        }

        AD_Fluid::ADFluidProperties<3> propsW{};
        AD_Fluid::ADFluidProperties<3> propsG{};
        if (has_sw) {
            propsW = FIM_Engine::EvalTwoPhaseWaterFluid<3>(fluid_ctx, P_cell, T_cell);
            propsG = FIM_Engine::EvalTwoPhaseGasFluid<3>(fluid_ctx, P_cell + Pc_cell, T_cell);
        }
        else {
            propsW = FIM_Engine::EvalSinglePhaseFluid<3>(fluid_ctx, P_cell, T_cell);
            propsG = fluid_ctx.config.single_phase_is_co2
                ? propsW
                : FIM_Engine::EvalTwoPhaseGasFluid<3>(fluid_ctx, P_cell, T_cell);
        }

        ADVar<3> mob_w = (propsW.rho * krw) / propsW.mu;
        ADVar<3> mob_g = (propsG.rho * krg) / propsG.mu;

        const bool has_w_eq = (dofOffset_W >= 0);
        const bool has_g_eq = (dofOffset_G >= 0);
        const bool total_mode = has_w_eq && has_g_eq && (step.component_mode == WellComponentMode::Total);

        const double pbh_eval = ResolveWellBhpValue(step, P_cell.val, wellBhpByName);
        ADVar<3> dP = P_cell - MakeConstAD3(pbh_eval);

        ADVar<3> q_mass_w = MakeConstAD3(0.0);
        ADVar<3> q_mass_g = MakeConstAD3(0.0);
        double dq_mass_w_dPbh = 0.0;
        double dq_mass_g_dPbh = 0.0;

        if (total_mode) {
            if (dP.val < 0.0) {
                const PhaseSplit inj_split = ResolveInjectionPhaseSplit(step, has_w_eq, has_g_eq);
                ADVar<3> q_total = MakeConstAD3(WI) * (mob_w + mob_g) * dP;
                q_mass_w = MakeConstAD3(inj_split.fw) * q_total;
                q_mass_g = MakeConstAD3(inj_split.fg) * q_total;
                const double dq_total_dPbh = -WI * (mob_w.val + mob_g.val);
                dq_mass_w_dPbh = inj_split.fw * dq_total_dPbh;
                dq_mass_g_dPbh = inj_split.fg * dq_total_dPbh;
            }
            else {
                q_mass_w = MakeConstAD3(WI) * mob_w * dP;
                q_mass_g = MakeConstAD3(WI) * mob_g * dP;
                dq_mass_w_dPbh = -WI * mob_w.val;
                dq_mass_g_dPbh = -WI * mob_g.val;
            }
        }
        else {
            const PhaseSplit split = ResolvePhaseSplit(step, mob_w.val, mob_g.val, has_w_eq, has_g_eq);
            if (split.fw > 0.0) {
                q_mass_w = MakeConstAD3(WI * split.fw) * mob_w * dP;
                dq_mass_w_dPbh = -WI * split.fw * mob_w.val;
            }
            if (split.fg > 0.0) {
                q_mass_g = MakeConstAD3(WI * split.fg) * mob_g * dP;
                dq_mass_g_dPbh = -WI * split.fg * mob_g.val;
            }
        }

        ADVar<3> h_w_well = propsW.h;
        ADVar<3> h_g_well = propsG.h;
        if (step.injection_temperature > 0.0) {
            ADVar<3> T_inj = MakeConstAD3(step.injection_temperature);
            ADVar<3> P_inj_const = MakeConstAD3(pbh_eval);
            if (q_mass_w.val < 0.0) {
                h_w_well = has_sw
                    ? FIM_Engine::EvalTwoPhaseWaterFluid<3>(fluid_ctx, P_inj_const, T_inj).h
                    : FIM_Engine::EvalSinglePhaseFluid<3>(fluid_ctx, P_inj_const, T_inj).h;
            }
            if (q_mass_g.val < 0.0) {
                h_g_well = (has_sw || !fluid_ctx.config.single_phase_is_co2)
                    ? FIM_Engine::EvalTwoPhaseGasFluid<3>(fluid_ctx, P_inj_const, T_inj).h
                    : FIM_Engine::EvalSinglePhaseFluid<3>(fluid_ctx, P_inj_const, T_inj).h;
            }
        }

        ADVar<3> q_total = q_mass_w + q_mass_g;
        ADVar<3> q_energy = (q_mass_w * h_w_well) + (q_mass_g * h_g_well);
        const double dq_total_dPbh = dq_mass_w_dPbh + dq_mass_g_dPbh;
        const double dq_energy_dPbh = dq_mass_w_dPbh * h_w_well.val + dq_mass_g_dPbh * h_g_well.val;

        const int eqIdx_W = (dofOffset_W >= 0) ? mgr.getEquationIndex(solverIdx, dofOffset_W) : -1;
        const int eqIdx_G = (dofOffset_G >= 0) ? mgr.getEquationIndex(solverIdx, dofOffset_G) : -1;
        const int eqIdx_E = (dofOffset_E >= 0) ? mgr.getEquationIndex(solverIdx, dofOffset_E) : -1;

        if (dofOffset_W >= 0) AddWellEqContribution(eqIdx_W, q_mass_w, jacobianFull, residual);
        if (dofOffset_G >= 0) AddWellEqContribution(eqIdx_G, q_mass_g, jacobianFull, residual);
        if (dofOffset_E >= 0) AddWellEqContribution(eqIdx_E, q_energy, jacobianFull, residual);

        if (completionLinearizations) {
            WellCompletionLinearization rec;
            rec.step_index = step_idx;
            rec.well_name = step.well_name;
            rec.completion_solver_index = (step.completion_solver_index >= 0) ? step.completion_solver_index : solverIdx;
            rec.domain = step.domain;
            rec.pbh_eval = pbh_eval;
            rec.q_total = q_total.val;
            rec.dqtotal_dcell = { { q_total.grad(0), q_total.grad(1), q_total.grad(2) } };
            rec.dqtotal_dPbh = dq_total_dPbh;
            FillWellEqLinearization(rec.water, eqIdx_W, dofOffset_W, q_mass_w, dq_mass_w_dPbh);
            FillWellEqLinearization(rec.gas, eqIdx_G, dofOffset_G, q_mass_g, dq_mass_g_dPbh);
            FillWellEqLinearization(rec.energy, eqIdx_E, dofOffset_E, q_energy, dq_energy_dPbh);
            completionLinearizations->push_back(rec);
        }

        if (step.domain == WellTargetDomain::Matrix) stats.matrixBCCount++;
        else stats.fractureBCCount++;
    }

    return stats;
}

#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "BCAdapter.h"
#include "DiffusionCentral.h"
#include "ConvectionUpwind_Flux.h"
#include "Solver_TimeLoopUtils.h"
#include "IMPES_TimeTermAssemblerandSolver.h"
#include "Solver_AssemblerCOO.h"
#include "FVM_WellCoupling_TwoPhase.h"
#include "IMPES_CommonUtils.h"
#include "TwoPhaseWells_StrictRate.h"

namespace IMPES
{

    struct PressureAssemblyConfig
    {
        std::string operator_tag = "p_impes";
        std::string pressure_field = "p_w";
        std::string pressure_old_field = "p_w_old";
        std::string pressure_prev_field = "p_w_prev";
        std::string phi_field = "phi_r";
        std::string lambda_total_field = "lambda_t";
        std::string rho_total_field = "rho_t";
        std::string rho_total_old_field = "rho_t_old";
        std::string drho_dp_field = "drho_t_dp";
        std::string drho_water_field = "drho_w_dp";
        std::string drho_gas_field = "drho_g_dp";
        std::string temperature_eval_field = "T";
        std::string saturation_field = "s_w";
        std::string rho_w_field = "rho_w";
        std::string rho_g_field = "rho_g";
        std::string kxx_field = "kxx";
        std::string kyy_field = "kyy";
        std::string kzz_field = "kzz";
        std::string total_mass_flux_name = "mf_total";
        std::string total_vol_flux_name = "Qf_total";
        std::string total_velocity_name = "ufn_total";
        double rock_compressibility = 0.0;
        Vector gravity = { 0.0, 0.0, 0.0 };
        bool enable_buoyancy = false;
        bool auto_update_total_density = true;
        int gradient_smoothing = 0;
        bool use_constant_phase_density = true; // true: rely on existing density fields (from ppm) and skip table lookups
    };

    struct PressureAssemblyResult
    {
        SparseSystemCOO system;
        std::vector<int> cell_lid;
        std::shared_ptr<faceScalarField> total_mass_flux;
        std::shared_ptr<faceScalarField> total_vol_flux;
        std::shared_ptr<faceScalarField> total_face_velocity;
    };

    namespace detail
    {
        inline bool ensureTotalDensityField(
            MeshManager& mgr,
            FieldRegistry& reg,
            const PressureAssemblyConfig& cfg)
        {
            if (!cfg.auto_update_total_density || cfg.rho_total_field.empty()) return true;

            auto sw = reg.get<volScalarField>(cfg.saturation_field);
            auto rho_w = reg.get<volScalarField>(cfg.rho_w_field);
            auto rho_g = reg.get<volScalarField>(cfg.rho_g_field);
            if (!sw || !rho_w || !rho_g) return reg.has(cfg.rho_total_field);

            const size_t n = mgr.mesh().getCells().size();
            if (!cfg.rho_total_old_field.empty())
            {
                auto existing = reg.get<volScalarField>(cfg.rho_total_field);
                if (existing)
                {
                    auto rho_old = reg.getOrCreate<volScalarField>(cfg.rho_total_old_field, n, 0.0);
                    rho_old->data = existing->data;
                }
            }

            auto rho_t = reg.getOrCreate<volScalarField>(cfg.rho_total_field, n, 0.0);
            std::shared_ptr<volScalarField> drho_total;
            auto drho_w = cfg.drho_water_field.empty() ? nullptr : reg.get<volScalarField>(cfg.drho_water_field);
            auto drho_g = cfg.drho_gas_field.empty() ? nullptr : reg.get<volScalarField>(cfg.drho_gas_field);
            if (!cfg.drho_dp_field.empty() && drho_w && drho_g)
            {
                drho_total = reg.getOrCreate<volScalarField>(cfg.drho_dp_field, n, 0.0);
                std::fill(drho_total->data.begin(), drho_total->data.end(), 0.0);
            }

            const auto& cells = mgr.mesh().getCells();
            const auto& id2idx = mgr.mesh().getCellId2Index();
            for (const auto& c : cells)
            {
                if (c.id < 0) continue;
                const size_t i = id2idx.at(c.id);
                const double sw_i = clampValue((*sw)[i], 0.0, 1.0);
                const double rho_w_i = std::max((*rho_w)[i], 0.0);
                const double rho_g_i = std::max((*rho_g)[i], 0.0);
                (*rho_t)[i] = sw_i * rho_w_i + (1.0 - sw_i) * rho_g_i;

                if (drho_total && drho_w && drho_g)
                {
                    const double dw = (*drho_w)[i];
                    const double dg = (*drho_g)[i];
                    (*drho_total)[i] = sw_i * dw + (1.0 - sw_i) * dg;
                }
            }
            return true;
        }

        inline void ensureFieldExists(
            MeshManager& mgr,
            FieldRegistry& reg,
            const std::string& name,
            double init_value = 0.0)
        {
            if (name.empty() || reg.has(name)) return;
            const size_t n = mgr.mesh().getCells().size();
            reg.getOrCreate<volScalarField>(name, n, init_value);
        }

        inline bool buildTotalFluxFields(
            MeshManager& mgr,
            FieldRegistry& reg,
            FaceFieldRegistry& freg,
            const PressureBCAdapter& pbc,
            const PressureAssemblyConfig& cfg,
            const OperatorFieldNames& nm)
        {
            return FVM::Convection::buildFlux_Darcy_Mass(
                mgr, reg, freg,
                nm.a_f_diff, nm.s_f_diff,
                cfg.pressure_field,
                cfg.rho_total_field,
                cfg.total_mass_flux_name,
                cfg.total_vol_flux_name,
                cfg.total_velocity_name,
                &pbc,
                cfg.enable_buoyancy);
        }

        inline bool updatePhaseDensityAndDerivative_const(
            MeshManager& /*mgr*/,
            FieldRegistry& reg,
            const PressureAssemblyConfig& /*cfg*/,
            const std::string& phase,
            const std::string& rho_field,
            const std::string& drho_field)
        {
            if (rho_field.empty()) return true;
            auto rho = reg.get<volScalarField>(rho_field);
            if (!rho)
            {
                std::cerr << "[IMPES][Pressure] constant-property mode missing density field '"
                          << rho_field << "' for phase " << phase << ".\n";
                return false;
            }
            if (!drho_field.empty())
            {
                auto drho = reg.getOrCreate<volScalarField>(drho_field, rho->data.size(), 0.0);
                std::fill(drho->data.begin(), drho->data.end(), 0.0);
            }
            return true;
        }

        inline bool updatePhaseDensityAndDerivative_variable(
            MeshManager& mgr,
            FieldRegistry& reg,
            const PressureAssemblyConfig& cfg,
            const std::string& phase,
            const std::string& rho_field,
            const std::string& drho_field)
        {
            if (rho_field.empty() || drho_field.empty()) return true;
            return computeRhoAndDrhoDpAt(
                mgr, reg,
                cfg.pressure_field,
                cfg.temperature_eval_field,
                phase,
                rho_field,
                drho_field);
        }

        inline bool updatePhaseDensityAndDerivative(
            MeshManager& mgr,
            FieldRegistry& reg,
            const PressureAssemblyConfig& cfg,
            const std::string& phase,
            const std::string& rho_field,
            const std::string& drho_field)
        {
            if (cfg.use_constant_phase_density)
            {
                return updatePhaseDensityAndDerivative_const(mgr, reg, cfg, phase, rho_field, drho_field);
            }
            return updatePhaseDensityAndDerivative_variable(mgr, reg, cfg, phase, rho_field, drho_field);
        }
    } // namespace detail

    inline bool assemblePressureSystem(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& pbc,
        const std::vector<WellDOF_TwoPhase>& wells,
        double dt,
        const PressureAssemblyConfig& cfg,
        PressureAssemblyResult& result)
    {
        if (!(dt > 0.0))
        {
            std::cerr << "[IMPES][Pressure] invalid time step.\n";
            return false;
        }

        if (!detail::updatePhaseDensityAndDerivative(mgr, reg, cfg, "water", cfg.rho_w_field, cfg.drho_water_field))
        {
            std::cerr << "[IMPES][Pressure] failed to update water density/derivative.\n";
            return false;
        }
        if (!detail::updatePhaseDensityAndDerivative(mgr, reg, cfg, "co2", cfg.rho_g_field, cfg.drho_gas_field))
        {
            std::cerr << "[IMPES][Pressure] failed to update CO2 density/derivative.\n";
            return false;
        }

        if (!detail::ensureTotalDensityField(mgr, reg, cfg))
        {
            std::cerr << "[IMPES][Pressure] failed to build total density field.\n";
            return false;
        }

        detail::ensureFieldExists(mgr, reg, cfg.rho_total_old_field);
        detail::ensureFieldExists(mgr, reg, cfg.drho_dp_field);

        if (!reg.has(cfg.lambda_total_field))
        {
            std::cerr << "[IMPES][Pressure] missing total mobility field '" << cfg.lambda_total_field << "'.\n";
            return false;
        }

        const OperatorFieldNames nm = makeNames(cfg.operator_tag);

        if (!TimeTerm_IMPES_Pressure(
            mgr, reg, dt,
            cfg.phi_field,
            cfg.pressure_old_field,
            cfg.pressure_field,
            cfg.rho_total_old_field,
            cfg.rho_total_field,
            cfg.drho_dp_field,
            cfg.rock_compressibility,
            nm.a_time,
            nm.b_time))
        {
            std::cerr << "[IMPES][Pressure] IMPES time-term assembly failed.\n";
            return false;
        }

        std::vector<std::string> mobility_tokens;
        mobility_tokens.emplace_back("kxx:" + cfg.kxx_field);
        mobility_tokens.emplace_back("kyy:" + cfg.kyy_field);
        mobility_tokens.emplace_back("kzz:" + cfg.kzz_field);
        mobility_tokens.emplace_back(cfg.lambda_total_field);

        if (!FVM::Diffusion::build_FaceCoeffs_Central(
            mgr, reg, freg,
            nm.a_f_diff, nm.s_f_diff,
            cfg.pressure_field,
            mobility_tokens,
            cfg.rho_total_field,
            FVM::Diffusion::RhoFaceMethod::Linear,
            cfg.gravity,
            pbc,
            cfg.enable_buoyancy,
            cfg.gradient_smoothing))
        {
            std::cerr << "[IMPES][Pressure] diffusion operator build failed.\n";
            return false;
        }

        SparseSystemCOO sys;
        if (!assemble_COO(mgr, reg, freg, "ddt+diffusion", nm, &sys))
        {
            std::cerr << "[IMPES][Pressure] assemble_COO failed.\n";
            return false;
        }

        int N = 0;
        result.cell_lid = buildUnknownMap(mgr.mesh(), N);
        if (static_cast<size_t>(N) != result.cell_lid.size())
        {
            std::cerr << "[IMPES][Pressure] inconsistent unknown map.\n";
            return false;
        }

        int desired_n = sys.n;
        for (const auto& well : wells)
        {
            desired_n = std::max(desired_n, well.lid + 1);
        }
        sys.expandTo(desired_n);

        for (const auto& well : wells)
        {
            FVM::TwoPhaseWellCoupling::couple_well_to_pressure_equation(
                sys, mgr, reg, well, result.cell_lid);
        }

        if (!cfg.total_mass_flux_name.empty())
        {
            if (!detail::buildTotalFluxFields(mgr, reg, freg, pbc, cfg, nm))
            {
                std::cerr << "[IMPES][Pressure] failed to compute total mass flux.\n";
                return false;
            }
            result.total_mass_flux = freg.get<faceScalarField>(cfg.total_mass_flux_name);
            result.total_vol_flux = freg.get<faceScalarField>(cfg.total_vol_flux_name);
            result.total_face_velocity = freg.get<faceScalarField>(cfg.total_velocity_name);
        }
        else
        {
            result.total_mass_flux.reset();
            result.total_vol_flux.reset();
            result.total_face_velocity.reset();
        }

        result.system = std::move(sys);
        return true;
    }
} // namespace IMPES


namespace IMPES_revised
{
    /**
     * @brief Configuration parameters for the revised IMPES two-phase pressure assembly.
     */
    struct PressureAssemblyConfig
    {
        // ======== 输入字段：时间项需要 ========
        /// rock porosity φ_r
        std::string phi_field = "phi_r";
        /// rock compressibility field c_r(p,T)
        std::string c_r_field = "c_r";
        /// total compressibility c_t = c_r + (Sw*rho_w*c_w + Sg*rho_g*c_g)/rho_t
        std::string c_t_field = "c_t";
        /// total density at current eval step rho_t
        std::string rho_total_field = "rho_t";
        /// total density at previous time level rho_t^n
        std::string rho_total_old_field = "rho_t_old";
        /// previous time level water pressure p_w^n
        std::string pressure_old_field = "p_w_old";
        /// current eval water pressure p_w (unknown of this step)
        std::string pressure_field = "p_w";

        // 可选：单相旧接口里用的 drho/dp，保留做诊断或输出
        std::string drho_dp_field = "drho_t_dp";

        // ======== 其它两相物性字段（目前时间项不直接用，但其它模块可能会用到） ========
        std::string c_w_field = "c_w";
        std::string c_g_field = "c_g";
        std::string rho_w_field = "rho_w";
        std::string rho_g_field = "rho_g";
        std::string saturation_field = "s_w";

        // ======== 扩散算子相关字段 ========
        std::string operator_tag = "p_impes_rev";
        std::string pressure_prev_field = "p_w_prev";
        std::string capillary_pressure_field = "Pc";
        std::string lambda_w_field = "lambda_w";
        std::string lambda_g_field = "lambda_g";
        std::string kxx_field = "kxx";
        std::string kyy_field = "kyy";
        std::string kzz_field = "kzz";
        VGParams      vg_params;
        RelPermParams relperm_params;

        // ======== 导出/派生场 ========
        std::string lambda_mass_field = "lambda_mass_total";     ///< λ_m = λ_w ρ_w + λ_g ρ_g
        std::string lambda_capillary_field = "lambda_capillary_mass"; ///< λ_cap = λ_g ρ_g
        std::string lambda_gravity_field = "lambda_gravity_mass";   ///< λ_gra = λ_w ρ_w² + λ_g ρ_g²
        std::string gravity_dummy_field = "gravity_dummy_scalar";
        std::string unity_scalar_field = "unity_scalar";
        std::string rho_eff_gravity_field = "rho_eff_gravity";
        std::string rho_coeff_field = "rho_coeff_mass";
        std::string rho_buoy_field = "rho_buoy_ratio";
        double      rho_buoy_epsilon = 1e-30;

        std::string total_mass_flux_name = "mf_total";
        std::string total_vol_flux_name = "Qf_total";
        std::string total_velocity_name = "ufn_total";
        std::string capillary_correction_flux_name = "mf_capillary_corr";
        std::string gravity_correction_flux_name = "mf_gravity_corr";

        Vector gravity = { 0.0, -9.8, 0.0 };
        bool   enable_buoyancy = true;
        int    gradient_smoothing = 0;
    };

    using PressureAssemblyResult = IMPES::PressureAssemblyResult;
    using TwoPhasePressureResult = PressureAssemblyResult;

    namespace detailed
    {
        inline std::shared_ptr<volScalarField> ensureField(
            FieldRegistry& reg,
            const std::string& name,
            size_t n,
            double init_value = 0.0)
        {
            if (name.empty()) return {};
            return reg.getOrCreate<volScalarField>(name, n, init_value);
        }

        inline bool gatherFieldToVector(
            Mesh& mesh,
            const FieldRegistry& reg,
            const std::string& field_name,
            const std::vector<int>& lid_of_cell,
            int nUnknowns,
            std::vector<double>& vec)
        {
            auto fld = reg.get<volScalarField>(field_name);
            if (!fld)
            {
                std::cerr << "[IMPES_revised][Pressure] missing field '" << field_name << "'.\n";
                return false;
            }
            vec.assign(nUnknowns, 0.0);
            const auto& cells = mesh.getCells();
            const auto& id2idx = mesh.getCellId2Index();
            for (const auto& c : cells)
            {
                if (c.id < 0) continue;
                const size_t idx = id2idx.at(c.id);
                const int lid = lid_of_cell[idx];
                if (lid < 0) continue;
                vec[lid] = (*fld)[idx];
            }
            return true;
        }

        inline void applyCOO(const SparseSystemCOO& sys, const std::vector<double>& x, std::vector<double>& y)
        {
            y.assign(sys.n, 0.0);
            for (const auto& t : sys.A)
            {
                if (t.c >= static_cast<int>(x.size())) continue;
                y[t.r] += t.v * x[t.c];
            }
            for (int i = 0; i < sys.n && i < static_cast<int>(sys.b.size()); ++i)
            {
                y[i] -= sys.b[i];
            }
        }

        inline bool computeDerivedMobilityFields(
            MeshManager& mgr,
            FieldRegistry& reg,
            const PressureAssemblyConfig& cfg)
        {
            auto lambda_w = reg.get<volScalarField>(cfg.lambda_w_field);
            auto lambda_g = reg.get<volScalarField>(cfg.lambda_g_field);
            auto rho_w = reg.get<volScalarField>(cfg.rho_w_field);
            auto rho_g = reg.get<volScalarField>(cfg.rho_g_field);
            if (!lambda_w || !lambda_g || !rho_w || !rho_g)
            {
                std::cerr << "[IMPES_revised][Pressure] missing lambda/rho fields.\n";
                return false;
            }

            const size_t n = mgr.mesh().getCells().size();
            auto lambda_mass = ensureField(reg, cfg.lambda_mass_field, n, 0.0);
            auto lambda_cap = ensureField(reg, cfg.lambda_capillary_field, n, 0.0);
            auto lambda_grav = ensureField(reg, cfg.lambda_gravity_field, n, 0.0);
            auto gravity_dummy = ensureField(reg, cfg.gravity_dummy_field, n, 0.0);
            auto unity = ensureField(reg, cfg.unity_scalar_field, n, 1.0);
            auto rho_eff = ensureField(reg, cfg.rho_eff_gravity_field, n, 0.0);
            auto rho_coeff = ensureField(reg, cfg.rho_coeff_field, n, 0.0);
            auto rho_buoy = ensureField(reg, cfg.rho_buoy_field, n, 0.0);
            if (!lambda_mass || !lambda_cap || !lambda_grav || !gravity_dummy || !unity || !rho_eff || !rho_coeff || !rho_buoy)
            {
                std::cerr << "[IMPES_revised][Pressure] failed to allocate derived fields.\n";
                return false;
            }

            const auto& cells = mgr.mesh().getCells();
            const auto& id2idx = mgr.mesh().getCellId2Index();
            const double tiny = std::max(1e-30, cfg.rho_buoy_epsilon);
            for (const auto& c : cells)
            {
                if (c.id < 0) continue;
                const size_t i = id2idx.at(c.id);
                const double lw = std::max((*lambda_w)[i], 0.0);
                const double lg = std::max((*lambda_g)[i], 0.0);
                const double rw = std::max((*rho_w)[i], 0.0);
                const double rg = std::max((*rho_g)[i], 0.0);

                const double lam_mass = lw * rw + lg * rg;
                const double lam_cap = lg * rg;
                const double lam_grav = lw * rw * rw + lg * rg * rg;
                const double rhoBuoy = (lam_mass > tiny) ? (lam_grav / std::max(lam_mass, tiny)) : 0.0;

                (*lambda_mass)[i] = lam_mass;
                (*lambda_cap)[i] = lam_cap;
                (*lambda_grav)[i] = lam_grav;
                (*gravity_dummy)[i] = 0.0;
                (*unity)[i] = 1.0;
                (*rho_eff)[i] = rhoBuoy;
                (*rho_coeff)[i] = lam_mass;
                (*rho_buoy)[i] = rhoBuoy;
            }
            return true;
        }

        inline bool buildAndApplyOperator(
            MeshManager& mgr,
            FieldRegistry& reg,
            FaceFieldRegistry& freg,
            const PressureBCAdapter& pbc,
            const OperatorFieldNames& nm,
            const std::vector<std::string>& mobility_tokens,
            const std::string& rho_coeff_field,
            const std::string& rho_buoy_field,
            const std::string& x_field,
            const PressureAssemblyConfig& cfg,
            bool enable_buoyancy,
            const std::vector<double>& eval_vec,
            std::vector<double>& result_vec)
        {
            if (!FVM::Diffusion::TwoPhaseDiffusionTemplate::build_FaceCoeffs_Central(
                mgr, reg, freg,
                nm.a_f_diff, nm.s_f_diff,
                x_field,
                mobility_tokens,
                rho_coeff_field,
                rho_buoy_field,
                FVM::Diffusion::RhoFaceMethod::Linear,
                cfg.gravity,
                pbc,
                enable_buoyancy,
                cfg.gradient_smoothing))
            {
                std::cerr << "[IMPES_revised][Pressure] helper diffusion build failed.\n";
                return false;
            }

            SparseSystemCOO helper;
            if (!assemble_COO(mgr, reg, freg, "diffusion", nm, &helper))
            {
                std::cerr << "[IMPES_revised][Pressure] helper assemble failed.\n";
                return false;
            }
            applyCOO(helper, eval_vec, result_vec);
            return true;
        }

        inline bool computeFaceFluxFromOperator(
            MeshManager& mgr,
            FieldRegistry& reg,
            FaceFieldRegistry& freg,
            const OperatorFieldNames& nm,
            const std::string& x_field,
            const std::string& flux_name)
        {
            if (flux_name.empty()) return true;

            auto aF = freg.get<faceScalarField>(nm.a_f_diff.c_str());
            auto sF = freg.get<faceScalarField>(nm.s_f_diff.c_str());
            if (!aF || !sF)
            {
                std::cerr << "[IMPES_revised][Pressure] missing operator face fields for flux export.\n";
                return false;
            }

            auto flux = freg.getOrCreate<faceScalarField>(flux_name.c_str(), aF->data.size(), 0.0);
            std::shared_ptr<volScalarField> xFld_sp;
            const volScalarField* xFld = nullptr;
            if (!x_field.empty())
            {
                xFld_sp = reg.get<volScalarField>(x_field.c_str());
                if (!xFld_sp)
                {
                    std::cerr << "[IMPES_revised][Pressure] missing field '" << x_field
                              << "' while exporting flux '" << flux_name << "'.\n";
                    return false;
                }
                xFld = xFld_sp.get();
            }

            const Mesh& mesh = mgr.mesh();
            const auto& faces = mesh.getFaces();
            const auto& id2idx = mesh.getCellId2Index();

            auto cellValue = [&](int cellId)->double
                {
                    if (!xFld || cellId < 0) return 0.0;
                    auto it = id2idx.find(cellId);
                    if (it == id2idx.end()) return 0.0;
                    const size_t idx = static_cast<size_t>(it->second);
                    if (idx >= xFld->data.size()) return 0.0;
                    return (*xFld)[idx];
                };

            for (const auto& F : faces)
            {
                const int iF = F.id - 1;
                const double xP = cellValue(F.ownerCell);
                double faceFlux = 0.0;
                if (!F.isBoundary() && F.neighborCell >= 0)
                {
                    const double xN = cellValue(F.neighborCell);
                    faceFlux = (*aF)[iF] * (xP - xN) - (*sF)[iF];
                }
                else
                {
                    faceFlux = (*aF)[iF] * xP - (*sF)[iF];
                }
                (*flux)[iF] = faceFlux;
            }
            return true;
        }
    } // namespace detailed

    inline bool assemblePressureTwoPhase(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const PressureBCAdapter& pbc,
        const std::vector<WellDOF_TwoPhase>& wells,
        double dt,
        const PressureAssemblyConfig& cfg,
        PressureAssemblyResult& result)
    {
        if (dt <= 0.0)
        {
            std::cerr << "[IMPES_revised][Pressure] invalid dt.\n";
            return false;
        }

        if (!detailed::computeDerivedMobilityFields(mgr, reg, cfg))
        {
            return false;
        }

        const OperatorFieldNames nm = makeNames(cfg.operator_tag);

       if (!TimeTerm_IMPES_Pressure(
            mgr, reg, dt,
            cfg.phi_field,            // φ_r
            cfg.pressure_old_field,   // p_w^n
            cfg.pressure_field,       // p_w (eval)
            cfg.rho_total_old_field,  // ρ_t^n
            cfg.rho_total_field,      // ρ_t
            cfg.c_t_field,            // c_t
            cfg.c_r_field,            // c_r
            nm.a_time,
            nm.b_time))
        {
            std::cerr << "[IMPES_revised][Pressure] time-term assembly failed.\n";
            return false;
        }

        std::vector<std::string> mobility_tokens{
            "kxx:" + cfg.kxx_field,
            "kyy:" + cfg.kyy_field,
            "kzz:" + cfg.kzz_field
        };

        if (!FVM::Diffusion::TwoPhaseDiffusionTemplate::build_FaceCoeffs_Central(
            mgr, reg, freg,
            nm.a_f_diff, nm.s_f_diff,
            cfg.pressure_field,
            mobility_tokens,
            cfg.rho_coeff_field,
            "",
            FVM::Diffusion::RhoFaceMethod::Linear,
            cfg.gravity,
            pbc,
            false,
            cfg.gradient_smoothing))
        {
            std::cerr << "[IMPES_revised][Pressure] main diffusion build failed.\n";
            return false;
        }

        SparseSystemCOO sys;
        if (!assemble_COO(mgr, reg, freg, "ddt+diffusion", nm, &sys))
        {
            std::cerr << "[IMPES_revised][Pressure] assemble_COO failed.\n";
            return false;
        }

        int N = 0;
        result.cell_lid = buildUnknownMap(mgr.mesh(), N);
        if (sys.n < N) sys.expandTo(N);

        std::vector<double> eval_vec;
        if (!detailed::gatherFieldToVector(mgr.mesh(), reg, cfg.capillary_pressure_field, result.cell_lid, N, eval_vec))
        {
            std::cerr << "[IMPES_revised][Pressure] failed to gather Pc.\n";
            return false;
        }

        OperatorFieldNames nmCap = makeNames(cfg.operator_tag + "_cap");
        std::vector<double> capContribution;
        if (!detailed::buildAndApplyOperator(
            mgr, reg, freg, pbc,
            nmCap,
            mobility_tokens,
            cfg.lambda_capillary_field,
            "",
            cfg.capillary_pressure_field,
            cfg,
            false,
            eval_vec,
            capContribution))
        {
            return false;
        }
        for (int r = 0; r < N; ++r)
        {
            sys.addb(r, -capContribution[r]);
        }
        if (!cfg.capillary_correction_flux_name.empty())
        {
            if (!detailed::computeFaceFluxFromOperator(
                mgr, reg, freg,
                nmCap,
                cfg.capillary_pressure_field,
                cfg.capillary_correction_flux_name))
            {
                return false;
            }
        }

        OperatorFieldNames nmGrav = makeNames(cfg.operator_tag + "_grav");
        std::vector<double> zeroVec(N, 0.0), gravContribution;
        if (!detailed::buildAndApplyOperator(
            mgr, reg, freg, pbc,
            nmGrav,
            mobility_tokens,
            cfg.rho_coeff_field,
            cfg.rho_buoy_field,
            cfg.gravity_dummy_field,
            cfg,
            cfg.enable_buoyancy,
            zeroVec,
            gravContribution))
        {
            return false;
        }
        for (int r = 0; r < N; ++r)
        {
            sys.addb(r, -gravContribution[r]);
        }
        if (!cfg.gravity_correction_flux_name.empty())
        {
            if (!detailed::computeFaceFluxFromOperator(
                mgr, reg, freg,
                nmGrav,
                cfg.gravity_dummy_field,
                cfg.gravity_correction_flux_name))
            {
                return false;
            }
        }

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
                cfg.vg_params,        // 或 cfg.vg / cfg.vg_matrix，看你 InitConfig 怎么命名
                cfg.relperm_params    // 同上，对应 RelPermParams
            );
        }

        if (!cfg.total_mass_flux_name.empty())
        {
            if (!FVM::Convection::buildFlux_Darcy_Mass(
                mgr, reg, freg,
                nm.a_f_diff, nm.s_f_diff,
                cfg.pressure_field,
                cfg.rho_eff_gravity_field,
                cfg.total_mass_flux_name,
                cfg.total_vol_flux_name,
                cfg.total_velocity_name,
                &pbc,
                false))
            {
                std::cerr << "[IMPES_revised][Pressure] failed to compute total flux.\n";
                return false;
            }
            result.total_mass_flux = freg.get<faceScalarField>(cfg.total_mass_flux_name);
            result.total_vol_flux = freg.get<faceScalarField>(cfg.total_vol_flux_name);
            result.total_face_velocity = freg.get<faceScalarField>(cfg.total_velocity_name);
        }
        else
        {
            result.total_mass_flux.reset();
            result.total_vol_flux.reset();
            result.total_face_velocity.reset();
        }

        result.system = std::move(sys);
        return true;
    }
} // namespace IMPES_revised





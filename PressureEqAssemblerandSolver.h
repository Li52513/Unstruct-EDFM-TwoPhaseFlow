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
#include "TwoPhaseWells_StrictRate.h"
#include "TimeTermAssemblerandSolver.h"
#include "SolverContrlStrName.h"

namespace IMPES_Iteration
{
    struct PressureAssemblyConfig
    {
        PressureEquation_String             P_Eq_str;
        std::string operator_tag =          P_Eq_str.operator_tag;                   // pressure operator tag for nm
        std::string pressure_field =        P_Eq_str.pressure_field;                 // current eval pressure field
        std::string pressure_old_field =    P_Eq_str.pressure_old_field;
        std::string pressure_prev_field =   P_Eq_str.pressure_prev_field;
        std::string pressure_g =            P_Eq_str.pressure_g;                     //current CO2 pressure
        std::string Pc_field =              P_Eq_str.Pc_field;
        Vector gravity = { 0.0, 0.0, 0.0 };
        bool enable_buoyancy = false;
        int gradient_smoothing = 0;
        TwoPhase_VG_Parameters              VG_Parameter;

        //扩散项离散系数临时储存名称
        std::string rho_coeff_field =       P_Eq_str.rho_coeff_field;                  ///< ρ_coeff = λ_w ρ_w + λ_g ρ_g
        std::string rho_capillary_field =   P_Eq_str.rho_capillary_field;               ///< ρ_cap = λ_g ρ_g
        std::string rho_gravity_field =     P_Eq_str.rho_gravity_field;                ///< ρ_gra = (λ_w ρ_w² + λ_g ρ_g²)/(λ_w ρ_w + λ_g ρ_g)
		std::string rho_mix_field =         P_Eq_str.rho_mix_field; 			       ///< ρ_mix = ρ_w s_w + ρ_g s_g
        std::string lambda_gravity_field =  P_Eq_str.lambda_gravity_field;             ///< λ_gra = λ_w ρ_w + λ_g ρ_g
        std::string gravity_dummy_field =   P_Eq_str.gravity_dummy_field;              /// 用于占位的重力场
    };

    struct PressureAssemblyResult
    {
        SparseSystemCOO system;
        std::vector<int> cell_lid;
    };


    //工具函数
    namespace detailed
    {
        /**
        * brief : 确保场存在且尺寸正确
        **/
        inline std::shared_ptr<volScalarField> ensureField(
            FieldRegistry& reg,
            const std::string& name,
            size_t n,
            double init_value = 0.0)
        {
            if (name.empty()) return {};
            return reg.getOrCreate<volScalarField>(name, n, init_value);
        }

        /**
        * brief: 从场中收集数据到向量
        **/
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

        /**
        * brief: 计算向量y=Ax-b
        **/
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

        /**
        * brief: 计算扩散项的通量系数
        **/
        inline bool computeDerivedMobilityFields
        (
            MeshManager& mgr,
            FieldRegistry& reg,
            const PressureAssemblyConfig& cfg
        )
        {
            //输入场
            auto lambda_g = reg.get<volScalarField>(TwoPhase::CO2().lambda_g_tag);
            auto lambda_w = reg.get<volScalarField>(TwoPhase::Water().lambda_w_tag);
            auto lambda_mass = reg.get<volScalarField>(TwoPhase::Auxiliaryparameters().lambda_mass_tag);
            auto rho_g = reg.get<volScalarField>(TwoPhase::CO2().rho_tag);
            auto rho_w = reg.get<volScalarField>(TwoPhase::Water().rho_tag);
            auto kxx = reg.get<volScalarField>("kxx");
            auto kyy = reg.get<volScalarField>("kyy");
            auto kzz = reg.get<volScalarField>("kzz");
            auto s_w = reg.get<volScalarField>(SaturationEquation_String().saturation);

            if (!lambda_g || !lambda_w || !rho_g || !rho_w
                || !kxx || !kyy || !kzz || !lambda_mass || !s_w)
            {
                std::cerr << "[IMPES_revised][Pressure] missing lambda/rho/perm/s_w fields.\n";
                return false;
            }

            //计算得到的场
            const size_t n = mgr.mesh().getCells().size();
            auto rho_coeff = ensureField(reg, cfg.rho_coeff_field, n, 0.0);
            auto rho_capillary = ensureField(reg, cfg.rho_capillary_field, n, 0.0);
            auto rho_grav = ensureField(reg, cfg.rho_gravity_field, n, 0.0);
			auto rho_mix = ensureField(reg, cfg.rho_mix_field, n, 0.0);
            auto lambda_gravity = ensureField(reg, cfg.lambda_gravity_field, n, 0.0);
            auto gravity_dummy = ensureField(reg, cfg.gravity_dummy_field, n, 0.0);
            if (!rho_coeff || !rho_capillary || !rho_grav || !rho_mix || !gravity_dummy || !lambda_gravity)
            {
                std::cerr << "[IMPES_revised][Pressure] failed to allocate derived fields.\n";
                return false;
            }

            ///1）第一项，主压力扩散项 \nabla\cdot\left\{\left[-\mathbf{k}\left(\lambda_{w}\rho_{w}+\lambda_{g}\rho_{g}\right)\right]\nabla p_{w}\right\}
            TwoPhase::calculateteTotalMobilityField(mgr, reg); //计算\lambda_{mass}=\left(\lambda_{w}\rho_{w}+\lambda_{g}\rho_{g}\right)，为(*lambda_mass)[i]赋值并更新

            ///2）第二项，毛细压力扩散项 \nabla\cdot\left\{\left[-\mathbf{k}\lambda_{g}\rho_{g}\right]\nabla P_{c}\right\}
            ///3）第三项，重力项 \nabla\cdot\left\{\left[-\mathbf{k}\left(\lambda_{w}\rho_{w}^{2}+\lambda_{g}\rho_{g}^{2}\right)\right]\mathbf{g}\right\}
            const auto& cells = mgr.mesh().getCells();
            const auto& id2idx = mgr.mesh().getCellId2Index();
            const double tiny = 1e-30;
            for (const auto& c : cells)
            {
                if (c.id < 0) continue;
                const size_t i = id2idx.at(c.id);
                const double lw = std::max((*lambda_w)[i], 0.0);
                const double lg = std::max((*lambda_g)[i], 0.0);
                const double rw = std::max((*rho_w)[i], 0.0);
                const double rg = std::max((*rho_g)[i], 0.0);
                double sw = (*s_w)[i];
                sw = std::max(0.0, std::min(1.0, sw));
                const double sg = 1.0 - sw;
                

                const double rhoCoeff = std::max((*lambda_mass)[i], 0.0);
                const double rhoCapillary = lg * rg;

                const double rhoGravity = lw * rw * rw + lg * rg * rg;
                const double rhoBuoy = (rhoCoeff > tiny) ? (rhoGravity / std::max(rhoCoeff, tiny)) : 0.0;
                const double lambda_Buoy = std::max((*lambda_mass)[i], 0.0);

                (*rho_mix)[i] = rw * sw + rg * sg;
                (*rho_coeff)[i] = rhoCoeff;
                (*rho_capillary)[i] = rhoCapillary;
                (*rho_grav)[i] = rhoBuoy;
                (*gravity_dummy)[i] = 0.0;
                (*lambda_gravity)[i] = lambda_Buoy;
            }
            return true;
        }

        /// 一个 Null BCProvider，永远不返回 ABC
        struct NullPressureBCAdapter {
            bool getABC(int, double&, double&, double&) const { return false; }
        };

        inline bool buildAndApplyOperator
        (
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
            std::vector<double>& result_vec
        )
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

	}//namespace detailed

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

        ///1）计算扩散项的通量系数
        if (!detailed::computeDerivedMobilityFields(mgr, reg, cfg))
        {
            return false;
        }

		//—— 组装压力方程矩阵 ——//
		///2）构造算子字段名称
        const OperatorFieldNames nm = makeNames(cfg.operator_tag);
		///3）组装矩阵
		
        ///—— 组装时间项 ——//
        
        if (!TimeTerm_IMPES_Pressure(mgr, reg, dt, cfg.pressure_old_field, cfg.pressure_field, nm.a_time, nm.b_time))
        {
			std::cerr << "[IMPES_revised][Pressure] IMPES time-term assembly failed.\n";
			return false;
        }

		///—— 组装主扩散项 ——//
        std::vector<std::string> mobility_tokens = {"kxx:kxx", "kyy:kyy", "kzz:kzz"};

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
			std::cerr << "[IMPES_revised][Pressure] diffusion operator build failed.\n";
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
        if (!detailed::gatherFieldToVector(mgr.mesh(), reg, cfg.Pc_field, result.cell_lid, N, eval_vec))
        {
            std::cerr << "[IMPES_revised][Pressure] failed to gather Pc.\n";
            return false;
        }
        PressureBC::Registry nullReg;
        PressureBCAdapter   nullPbc{ nullReg };

        ///—— 组装毛细压力扩散项 ——//
        OperatorFieldNames nmCap = makeNames(cfg.operator_tag + "_cap");
		std::vector<double> capContribution;

        if (!detailed::buildAndApplyOperator(
            mgr, reg, freg, nullPbc,
            nmCap,
            mobility_tokens,
            cfg.rho_capillary_field,
            "",
            cfg.Pc_field,
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
		///—— 组装重力项 ——//
        OperatorFieldNames nmGrav = makeNames(cfg.operator_tag + "_grav");
        std::vector<double> zeroVec1(N, 0.0), gravContribution;
        if (!detailed::buildAndApplyOperator(
            mgr, reg, freg, nullPbc,
            nmGrav,
            mobility_tokens,
            cfg.rho_coeff_field,
            cfg.rho_gravity_field,
            cfg.gravity_dummy_field,
            cfg,
            cfg.enable_buoyancy,
            zeroVec1,
            gravContribution))
        {
            return false;
        }
        for (int r = 0; r < N; ++r)
        {
            sys.addb(r, -gravContribution[r]);
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
                cfg.VG_Parameter.vg_params,        // 或 cfg.vg / cfg.vg_matrix，看你 InitConfig 怎么命名
                cfg.VG_Parameter.relperm_params    // 同上，对应 RelPermParams
            );
        }
        result.system = std::move(sys);
        return true;
    }	    
}
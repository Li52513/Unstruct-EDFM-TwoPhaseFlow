#pragma once
#include <string>
#include "CapRelPerm.h"

namespace IMPES_Iteration
{
    struct PressureEquation_String
    {
        std::string operator_tag = "p_w_IMPES"; // pressure operator tag for nm
        std::string pressure_field = "p_w";    // current eval pressure field
        std::string pressure_old_field = "p_w_old";
        std::string pressure_prev_field = "p_w_prev";
        std::string pressure_g = "p_g";   //current CO2 pressure
        std::string Pc_field = "Pc";

        //扩散项离散系数临时储存名称
        std::string rho_coeff_field = "rho_coeff_mass";              ///< ρ_coeff = λ_w ρ_w + λ_g ρ_g
        std::string rho_capillary_field = "rho_capillary_mass";      ///< ρ_cap = λ_g ρ_g
        std::string rho_gravity_field = "rho_gravity_mass";          ///< ρ_gra = (λ_w ρ_w² + λ_g ρ_g²)/(λ_w ρ_w + λ_g ρ_g)
        std::string rho_mix_field = "rho_mix"; 			             ///< ρ_mix = ρ_w s_w + ρ_g s_g
        std::string lambda_gravity_field = "lambda_gravity_mass";   ///< λ_gra = λ_w ρ_w + λ_g ρ_g
        std::string gravity_dummy_field = "gravity_dummy_scalar";   /// 用于占位的重力场
    };

    struct FaceMassRate_String
    {
        ///储存通量名称的场
        std::string total_mass_flux_name = "mf_total";
        std::string capillary_correction_flux_name = "mf_capillary_corr";
        std::string gravity_correction_flux_name = "mf_gravity_corr";
        std::string total_vol_flux_name = "Qf_total";
        std::string total_velocity_name = "ufn_total";
    };

    struct SaturationEquation_String
    {
        std::string saturation =        "s_w";              ///当前时间层的水相饱和度（正在求解的）
        std::string saturation_old =    "s_w_old";          ///上一时间步的水相饱和度 （已知解，用于时间项）
        std::string saturation_prev =   "s_w_prev";         ///外迭代 / RK2 stage 备用的饱和度拷贝（可选）  
        std::string water_source_field = "";                /// 水相源汇项，单位 [kg/s]，例如井源的水相质量源 若为空字符串则视为无显式体源项
    };

    struct FluxSplitConfig_String
    {
        std::string water_mass_flux = "mf_w";          // Output: water-phase mass flux
        std::string gas_mass_flux = "mf_g";            // Output: gas-phase (CO2) mass flux
        std::string fractional_flow_face = "fw_face";  // Optional: water fractional flow on faces (can be empty)   
    };

    struct TwoPhase_VG_Parameters
    {
        VGParams	                        vg_params;
        RelPermParams                       relperm_params;
    };

    /// 饱和度时间步控制方案
    enum class SatTimeControlScheme
    {
        SimpleCFL,   ///< 现有方案: CFL ≈ dt |F_div| / (φ V ρ), dt_suggest 来自 ΔS_max 条件
        RedondoLike  ///< Redondo 风格: 用 λ_w, λ_g, dPc/dSw, f_w', Σ|q_tot| 构造 dt_CFL
    };

}


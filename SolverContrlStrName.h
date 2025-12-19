#pragma once
#include <string>
#include "CapRelPerm.h"


namespace PhysicalProperties_string
{
    struct Water
    {
        // Basic properties
        std::string rho_tag = "rho_w";
        std::string mu_tag = "mu_w";
        std::string k_tag = "k_w";
        std::string cp_tag = "cp_w";
        std::string drho_w_dp_tag = "drho_wdp";
        std::string c_w_tag = "c_w";
        // old timen layer
        std::string rho_old_tag = "rho_w_old";
        std::string mu_old_tag = "mu_w_old";
        std::string k_old_tag = "k_w_old";
        std::string cp_old_tag = "cp_w_old";
        std::string drho_w_dp_old_tag = "drho_wdp_old";
        std::string c_w_old_tag = "c_w_old";
        // Two-phase properties
        std::string k_rw_tag = "k_rw";				//relative permeability of water
        std::string lambda_w_tag = "lambda_w";		//mobility of water
        std::string k_rw_old_tag = "k_rw_old";				//relative permeability of water
        std::string lambda_w_old_tag = "lambda_w_old";		//mobility of water
    };

    struct CO2
    {
        // Basic properties
        std::string rho_tag = "rho_g";
        std::string mu_tag = "mu_g";
        std::string k_tag = "k_g";
        std::string cp_tag = "cp_g";
        std::string drho_g_dp_tag = "drho_gdp";
        std::string c_g_tag = "c_g";
        // old timen layer
        std::string rho_old_tag = "rho_g_old";
        std::string mu_old_tag = "mu_g_old";
        std::string k_old_tag = "k_g_old";
        std::string cp_old_tag = "cp_g_old";
        std::string drho_g_dp_old_tag = "drho_gdp_old";
        std::string c_g_old_tag = "c_g_old";
        // Two-phase properties
        std::string k_rg_tag = "k_rg";				        //relative permeability 
        std::string lambda_g_tag = "lambda_g";		        //mobility
		std::string k_rg_old_tag = "k_rg_old";				//relative permeability
		std::string lambda_g_old_tag = "lambda_g_old";		//mobility
    };

    struct Rock
    {
        // Bsic properties
        std::string phi_tag = "phi_r";	
        std::string phi_old_tag = "phi_r_old";      //porosity (old time layer)//porosity
        std::string rho_tag = "rho_r";				//rock density
        std::string cp_tag = "cp_r";				//rock specific heat capacity
        std::string lambda_tag = "lambda_r";				//rock thermal conductivity
        std::string k_xx_tag = "kxx";
        std::string k_yy_tag = "kyy";
        std::string k_zz_tag = "kzz";
        std::string c_r_tag = "c_r";				//rock compressibility
    };

    struct SinglePhase_case
    {
        std::string C_eff_tag = "C_eff";
        std::string lambda_eff_tag = "lambda_eff";
        std::string C_eff_old_tag = "C_eff_old";
        std::string lambda_eff_old_tag = "lambda_eff_old";
    };

    struct TwoPhase_case
    {
        std::string Pc_tag = "Pc";					//capillary pressure
        std::string dPc_dSw_tag = "dPc_dSw";		//derivative of capillary pressure to water saturation
        std::string lambda_mass_tag = "lambda_mass";// equation : lambda_mass = lambda_w * rho_w + lambda_g * rho_g
    };
}

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

namespace FC_P_IMPES_I
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
}


namespace SinglePhase 
{
    struct PhysicalParameters_String
    {
        //Fluid in matrix
        std::string rho_fluid_field = PhysicalProperties_string::CO2().rho_tag;
        std::string rho_fluid_old_field = PhysicalProperties_string::CO2().rho_old_tag;
        std::string drho_dp_tag = PhysicalProperties_string::CO2().drho_g_dp_tag;
        std::string cp_fluid_field = PhysicalProperties_string::CO2().cp_tag;
        std::string C_eff_field = PhysicalProperties_string::SinglePhase_case().C_eff_tag;
        std::string C_eff_old_field = PhysicalProperties_string::SinglePhase_case().C_eff_old_tag;

        //Matrix Rock
        std::string phi_tag = PhysicalProperties_string::Rock().phi_tag;				//porosity
        std::string c_r_tag = PhysicalProperties_string::Rock().c_r_tag;				//rock compressibility
    };
    
    struct PressureEquation_String
    {
        std::string operator_tag = "p_g_singlePhase";   // pressure operator tag for nm
        std::string pressure_field = "p_g";             // current eval pressure field
        std::string pressure_old_field = "p_g_old";
        std::string pressure_prev_field = "p_g_prev";
    };

    struct TemperatureEquation_String
    {
        std::string operator_tag = "T_singlePhase";   // Temperature operator tag for nm
        std::string temperatue_field = "T";
        std::string temperatue_old_field = "T_old";
        std::string temperatue_prev_field = "T_prev";
    };

    struct FaceMassRate_String
    {
        ///储存通量名称的场
        std::string total_mass_flux_name = "mf_total";
        std::string total_vol_flux_name = "Qf_total";
        std::string total_velocity_name = "ufn_total";
    };
}






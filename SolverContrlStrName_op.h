#pragma once

#include <string>
#include <vector>

// 引入 VG 模型及相渗参数定义
#include "CapRelPerm_HD.h"

/**
 * @namespace PhysicalProperties_string
 * @brief 物理属性与求解器控制字符串定义的统一管理命名空间 (Refactored Version)
 * @details
 * 采用 "合并同类项" 策略：
 * 1. 移除了 SinglePhase, IMPES_Iteration 等子命名空间。
 * 2. 将分散的 PressureEquation_String 等合并为通用结构体。
 * 3. 通过静态方法提供不同求解器的预设配置。
 */
namespace PhysicalProperties_string_op
{
    // =========================================================
    // Part 1: 物质固有属性 (Material Intrinsic Properties)
    // =========================================================

    struct CO2
    {
        std::string rho_tag = "rho_g";
        std::string mu_tag = "mu_g";
        std::string k_tag = "lambda_g";   // 导热系数
        std::string cp_tag = "cp_g";
        std::string h_tag = "h_g";
        std::string drho_dp_tag = "drho_g_dp";
        std::string c_g_tag = "c_g";

        // Old time layer
        std::string rho_old_tag = "rho_g_old";
        std::string mu_old_tag = "mu_g_old";
        std::string k_old_tag = "lambda_g_old";
        std::string cp_old_tag = "cp_g_old";
        std::string drho_dp_old_tag = "drho_g_dp_old";
        std::string c_g_old_tag = "c_g_old";
        std::string h_g_old_tag = "h_g_old";

        // Two-phase properties
        std::string k_rg_tag = "k_rg";
        std::string lambda_g_tag = "lambda_g_mob"; // mobility (避免与导热系数混淆)
        std::string k_rg_old_tag = "k_rg_old";
        std::string lambda_g_old_tag = "lambda_g_mob_old";
    };

    struct Water
    {
        std::string rho_tag = "rho_w";
        std::string mu_tag = "mu_w";
        std::string k_tag = "lambda_w";   // 导热系数
        std::string cp_tag = "cp_w";
        std::string h_tag = "h_w";
        std::string drho_dp_tag = "drho_w_dp";
        std::string c_w_tag = "c_w";

        // Old time layer
        std::string rho_old_tag = "rho_w_old";
        std::string mu_old_tag = "mu_w_old";
        std::string k_old_tag = "lambda_w_old";
        std::string cp_old_tag = "cp_w_old";
        std::string drho_dp_old_tag = "drho_w_dp_old";
        std::string c_w_old_tag = "c_w_old";

        // Two-phase properties
        std::string k_rw_tag = "k_rw";
        std::string lambda_w_tag = "lambda_w_mob";
        std::string k_rw_old_tag = "k_rw_old";
        std::string lambda_w_old_tag = "lambda_w_mob_old";
    };

    struct Rock
    {
        std::string k_xx_tag = "K_xx";
        std::string k_yy_tag = "K_yy";
        std::string k_zz_tag = "K_zz";
        std::string phi_tag = "phi";
        std::string phi_old_tag = "phi_old";
        std::string rho_tag = "rho_r";
        std::string cp_tag = "cp_r";
        std::string lambda_tag = "lambda_r";
        std::string c_r_tag = "c_r";
    };

    // [Renamed & Merged] 裂缝物性
    struct Fracture_string
    {
        std::string k_n_tag = "fr_k_n";
        std::string k_t_tag = "fr_k_t";
        std::string aperture_tag = "fr_aperture";
        std::string phi_tag = "fr_phi";
        std::string rho_tag = "fr_rho";
        std::string cp_tag = "fr_cp";
        std::string lambda_tag = "fr_lambda";
        std::string c_r_tag = "fr_c_r";
    };

    // [Merged] 有效/混合物性 (原 SinglePhase_case)
    struct EffectiveProps
    {
        std::string C_eff_tag = "C_eff";      // 有效热容
        std::string C_eff_old_tag = "C_eff_old";
        std::string lambda_eff_tag = "lambda_eff"; // 有效导热系数
        std::string lambda_eff_old_tag = "lambda_eff_old";
    };

    // [Merged] 两相辅助变量 (原 TwoPhase_case)
    struct TwoPhaseAux
    {
        std::string Pc_tag = "Pc";
        std::string dPc_dSw_tag = "dPc_dSw";
        std::string lambda_mass_tag = "lambda_mass"; // Total mass mobility
    };


    // =========================================================
    // Part 2: 求解器控制变量 (Merged Structs)
    // =========================================================

    /**
     * @struct PressureEquation_String
     * @brief [合并同类项] 通用压力方程控制字符串
     * @details 包含 SinglePhase 和 IMPES 所需的所有压力相关字段
     */
    struct PressureEquation_String
    {
        // --- Common Fields ---
        std::string operator_tag = "p_op";
        std::string pressure_field = "p";      // Current Pressure (查表用主变量)
        std::string pressure_old_field = "p_old";
        std::string pressure_prev_field = "p_prev";

        // --- IMPES / TwoPhase Specific ---
        std::string pressure_g = "p_g";    // Gas phase pressure
        std::string pressure_g_old = "p_g_old";
        std::string pressure_g_prev = "p_g_prev";
        std::string Pc_field = "Pc";

        // IMPES Coefficients 扩散项离散系数临时储存名称
        std::string rho_coeff_field = "rho_coeff_mass";
        std::string rho_capillary_field = "rho_capillary_mass";
        std::string rho_gravity_field = "rho_gravity_mass";
        std::string rho_mix_field = "rho_mix";
        std::string lambda_gravity_field = "lambda_gravity_mass";
        std::string gravity_dummy_field = "gravity_dummy_scalar";

        // --- Factory Methods ---

        static PressureEquation_String SinglePhase() {
            PressureEquation_String s;
            s.operator_tag = "p_g_singlePhase";
            s.pressure_field = "p_g";
            s.pressure_old_field = "p_g_old";
            s.pressure_prev_field = "p_g_prev";
            return s;
        }

        static PressureEquation_String IMPES() {
            PressureEquation_String s;
            s.operator_tag = "p_w_IMPES";
            s.pressure_field = "p_w"; // 或者 "p_global"，视您的主变量而定
            s.pressure_old_field = "p_w_old";
            s.pressure_prev_field = "p_w_prev";
            // ... IMPES specific fields use defaults or can be customized here
            return s;
        }
    };

    /**
     * @struct TemperatureEquation_String
     * @brief [合并同类项] 通用温度方程控制字符串
     */
    struct TemperatureEquation_String
    {
        std::string operator_tag = "T_op";
        std::string temperatue_field = "T";      // 查表用温度变量
        std::string temperatue_old_field = "T_old";
        std::string temperatue_prev_field = "T_prev";

        static TemperatureEquation_String SinglePhase() {
            TemperatureEquation_String s;
            s.operator_tag = "T_singlePhase";
            return s;
        }

        static TemperatureEquation_String IMPES() {
            TemperatureEquation_String s;
            s.operator_tag = "T_IMPES";
            s.temperatue_field = "T_res"; // IMPES 下可能叫 T_res
            s.temperatue_old_field = "T_res_old";
            return s;
        }
    };

    /**
     * @struct FaceMassRate_String
     * @brief [合并同类项] 通用通量字符串
     */
    struct FaceMassRate_String
    {
        std::string total_mass_flux_name = "mf_total";
        std::string total_vol_flux_name = "Qf_total";
        std::string total_velocity_name = "ufn_total";

        // IMPES Specific
        std::string capillary_correction_flux_name = "mf_capillary_corr";
        std::string gravity_correction_flux_name = "mf_gravity_corr";

        static FaceMassRate_String SinglePhase() {
            return FaceMassRate_String();
        }

        static FaceMassRate_String IMPES() {
            return FaceMassRate_String(); // 默认值即涵盖
        }
    };

    /**
     * @struct SaturationEquation_String
     * @brief 饱和度方程 (仅 IMPES 使用)
     */
    struct SaturationEquation_String
    {
        std::string saturation = "s_w";
        std::string saturation_old = "s_w_old";
        std::string saturation_prev = "s_w_prev";
        std::string water_source_field = "src_w";
    };

    /**
     * @struct FluxSplitConfig_String
     * @brief 分流量配置
     */
    struct FluxSplitConfig_String
    {
        std::string water_mass_flux = "mf_w";
        std::string gas_mass_flux = "mf_g";
        std::string fractional_flow_face = "fw_face";
    };

    struct TwoPhase_VG_Parameters
    {
        CapRelPerm::VGParams      vg_params;
        CapRelPerm::RelPermParams relperm_params;
    };

    enum class SatTimeControlScheme
    {
        SimpleCFL,
        RedondoLike
    };
}
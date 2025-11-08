#pragma once
#include <string>
#include <algorithm>
#include <cmath>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "FieldAcessForDiscre.h"
#include "Solver_AssemblerCOO.h"
#include "BCAdapter.h"
#include "TemperatureBCAdapter.h"
namespace FVM {
	namespace Convection {

        /**
        * @brief 路线 A1：达西驱动（质量通量形态，复用已有 a_f/s_f）
        *
        * 需求：FaceField 中已存在 a_f 与 s_f（含密度的“质量通量系数/源”），
        *      与装配约定 b += s_f 一致：mf = a_f*Δp − s_f。
        *
         * 输出：创建/写入面场 mf, Qf, ufn（三者均覆盖写）。
        *
        * @param a_name  质量形 a_f（如 "a_f_Diff" 或 "a_f_Diff_p_w"）
        * @param s_name  质量形 s_f（如 "s_f_Diff" 或 "s_f_Diff_p_w"）
        * @param p_name  压力场名（如 "p", "p_w", "p_g"）
        * @param rho_name 密度体场名（用于 Qf = mf/ρ_face 与边界处理）
        * @param mf_name 输出质量通量面场名（默认 "mf"）
        * @param Qf_name 输出体积通量面场名（默认 "Qf"）
        * @param ufn_name 输出法向速度面场名（默认 "ufn"）
        */
        bool buildFlux_Darcy_Mass(
            MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
            const std::string& a_name = "a_f_Diff",
            const std::string& s_name = "s_f_Diff",
            const std::string& p_name = "p",
            const std::string& rho_name = "rho",
            const std::string& mf_name = "mf",
            const std::string& Qf_name = "Qf",
            const std::string& ufn_name = "ufn",
            const PressureBCAdapter* bc = nullptr,
            bool clampDirichletBackflow = false);


        /**
        * @brief 路线 A2：达西驱动（体积通量形态，先 Qf 后乘 ρ）
        *
        * 需求：FaceField 中已存在“体积形” a_vol/s_vol（不含密度），
        *      Qf = a_vol*Δp + s_vol；mf = ρ_up * Qf（内部上风，边界取 owner）。
        *
        * 输出：创建/写入面场 mf, Qf, ufn。
         *
        * @param a_vol   体积形 a_f（如 "a_f_Diff_vol"）
        * @param s_vol   体积形 s_f（如 "s_f_Diff_vol"）
        * @param p_name  压力场名
        * @param rho_name 密度体场名（用于上风密度与边界）
        * @param mf_name 输出质量通量面场名
        * @param Qf_name 输出体积通量面场名
        * @param ufn_name 输出法向速度面场名
        */
        bool buildFlux_Darcy_Vol
        (
            MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
            const std::string& a_vol = "a_f_Diff_vol",
            const std::string& s_vol = "s_f_Diff_vol",
            const std::string& p_name = "p",
            const std::string& rho_name = "rho",
            const std::string& mf_name = "mf",
            const std::string& Qf_name = "Qf",
            const std::string& ufn_name = "ufn",
            const PressureBCAdapter* bc = nullptr,
            const TemperatureBCAdapter* tbc = nullptr)
            ;

        /**
        * @brief 路线 B：外给速度（常量或体场 U），直接构造通量
        *
        * 逻辑：U_f 线性插值到面 → u_n = U_f·n → Qf = u_n * |A| → mf = ρ_up * Qf。
        * 若 U_name 为空，使用常量 U_const。
        *
        * 输出：创建/写入面场 mf, Qf, ufn。
        *
        * @param U_name   体矢量速度场名（可留空）
        * @param U_const  常量速度（当 U_name 为空或缺失时使用）
        * @param rho_name 密度体场名（内部上风，边界取 owner）
        * @param mf_name  输出质量通量面场名
        * @param Qf_name  输出体积通量面场名
        * @param ufn_name 输出法向速度面场名
        */
        bool buildFlux_FromVelocity
        (
            MeshManager& mgr, const FieldRegistry& reg, FaceFieldRegistry& freg,
            const std::string& U_name = "",
            Vector U_const = { 1, 0, 0 },
            const std::string& rho_name = "rho",
            const std::string& mf_name = "mf",
            const std::string& Qf_name = "Qf",
            const std::string& ufn_name = "ufn",
            const PressureBCAdapter * bc = nullptr);



    }

}

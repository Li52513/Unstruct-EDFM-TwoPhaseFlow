#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "Solver_AssemblerCOO.h"
#include "CapRelPerm.h"   // VGParams, RelPermParams, kr_Mualem_vG
#include "WellConfig_TwoPhase.h"    // WellDOF_TwoPhase (mode, role, target, Tin, etc.)

namespace FVM {
    namespace TwoPhaseWellsStrict {

        /// 小常数，避免除零
        static constexpr double kTiny = 1e-30;

        /**
         * @brief 注入混合物在井底条件下的流性参数（质量形式）。
         *
         * - krw, krg: 由 VG+Mualem 模型基于井底给定饱和度 s_w_bh 计算
         * - lambda_w, lambda_g: 体积动度（不含渗透率；渗透率包含在 WI 中）
         * - lambda_mass_w, lambda_mass_g: 质量动度 = λ * ρ
         * - lambda_mass: 总质量动度 = λ_mass_w + λ_mass_g
         * - fw_mass: 注入混合物的水相质量分率
         */
        struct InjectionMixtureProps
        {
            double krw = 0.0;
            double krg = 0.0;
            double lambda_w = 0.0;
            double lambda_g = 0.0;
            double lambda_mass_w = 0.0;
            double lambda_mass_g = 0.0;
            double lambda_mass = 0.0;
            double fw_mass = 0.0;
        };

        /**
         * @brief 根据井底给定状态计算注入混合物的质量动度和质量分率。
         *
         * 使用统一的 VG + Mualem 相对渗透率模型，保证井源与体方程在物性上的一致性。
         *
         * @param well      两相井 DOF（包含 s_w_bh, mu_w_inj, mu_g_inj, rho_w_inj, rho_g_inj 等）
         * @param vg_params van Genuchten 参数
         * @param rp_params 相对渗透率参数（Mualem）
         */
        inline InjectionMixtureProps computeInjectionMixtureProps(
            const WellDOF_TwoPhase& well,
            const VGParams& vg_params,
            const RelPermParams& rp_params)
        {
            InjectionMixtureProps props;

            // 1) VG + Mualem 模型：基于井底给定饱和度 s_w_bh 计算 krw, krg
            kr_Mualem_vG(well.s_w_bh, vg_params, rp_params,
                props.krw, props.krg);

            // 2) 体积动度（不含渗透率；渗透率在 WI 中）
            props.lambda_w = props.krw / std::max(well.mu_w_inj, 1e-12);
            props.lambda_g = props.krg / std::max(well.mu_g_inj, 1e-12);

            // 3) 质量动度：λ_mass = λ * ρ
            props.lambda_mass_w = props.lambda_w * well.rho_w_inj;
            props.lambda_mass_g = props.lambda_g * well.rho_g_inj;
            props.lambda_mass = props.lambda_mass_w + props.lambda_mass_g;

            if (props.lambda_mass > kTiny) {
                props.fw_mass = props.lambda_mass_w / props.lambda_mass;
            }
            else {
                props.fw_mass = 0.0;
            }

            return props;
        }

        /**
         * @brief 将两相井耦合到压力方程，使用“质量动度”实现严格的 Rate 约束。
         *
         * 对 Injector + Mode::Rate：
         *   - 使用注入混合物质量动度 λ_mass,inj 作为 PI 的动度部分
         *   - 构造离散方程：Σ WI_i λ_mass,inj (p_bh - p_i) = well.target (kg/s)
         *   - 其中 well.target 明确解释为“总注入质量流量” (kg/s)
         *
         * 对 Producer 或 Pressure 模式：
         *   - 使用当前储层场 λ_w, λ_g, ρ_w, ρ_g 构造局部质量动度进行耦合
         *
         * @param sys          压力方程稀疏系统 (A, b)
         * @param mgr          网格管理器
         * @param reg          场变量注册表（提供 lambda_w, lambda_g, rho_w, rho_g, mask, WI）
         * @param well         当前井的 DOF 配置（包含 mode, role, target, lid 等）
         * @param lid_of_cell  cell.id -> 全局自由度编号的映射（长度=单元数）
         * @param vg_params    VG 参数（用于注入井 kr 计算）
         * @param rp_params    相对渗透率参数
         */
        inline void couple_well_to_pressure_equation_strict_rate(
            SparseSystemCOO& sys,
            MeshManager& mgr,
            const FieldRegistry& reg,
            const WellDOF_TwoPhase& well,
            const std::vector<int>& lid_of_cell,
            const VGParams& vg_params,
            const RelPermParams& rp_params)
        {
            Mesh& mesh = mgr.mesh();
            const int Nc = static_cast<int>(mesh.getCells().size());
            const int well_lid = well.lid;

            if (well_lid < 0) {
                std::cerr << "[TwoPhaseWellsStrict] Invalid well LID for "
                    << well.name << ", skipping.\n";
                return;
            }

            // 0) 井自己的 DOF 行：Pressure 模式直接钉住 P_bh
            if (well.mode == WellDOF_TwoPhase::Mode::Pressure) {
                sys.addA(well_lid, well_lid, 1.0);
                sys.addb(well_lid, well.target);   ///< P_bh = target (Pa)
            }

            // 1) 必要字段：完井掩码 & 几何井指数 WI（不含动度）
            auto mask_field = reg.get<volScalarField>(well.mask_field.c_str());
            auto WI_field = reg.get<volScalarField>(well.PI_field_w.c_str()); // 存 WI

            if (!mask_field || !WI_field) {
                std::cerr << "[TwoPhaseWellsStrict] Missing mask or WI field for well "
                    << well.name << "\n";
                return;
            }

            // 2) 储层相动度与密度场（用于 Producer / Pressure 模式）
            auto lambda_w_field = reg.get<volScalarField>("lambda_w");
            auto lambda_g_field = reg.get<volScalarField>("lambda_g");
            auto rho_w_field = reg.get<volScalarField>("rho_w");
            auto rho_g_field = reg.get<volScalarField>("rho_g");

            if (!lambda_w_field || !lambda_g_field || !rho_w_field || !rho_g_field) {
                std::cerr << "[TwoPhaseWellsStrict] Missing fields "
                    "'lambda_w', 'lambda_g', 'rho_w', 'rho_g'.\n";
                return;
            }

            const bool isInjector = (well.role == WellDOF_TwoPhase::Role::Injector);

            // 3) 对注入井预计算混合物质量动度（位置无关）
            InjectionMixtureProps injProps;
            if (isInjector) {
                injProps = computeInjectionMixtureProps(well, vg_params, rp_params);
                if (injProps.lambda_mass <= kTiny) {
                    std::cerr << "[TwoPhaseWellsStrict] Injection mixture has zero "
                        "mass mobility for well " << well.name << ", skip.\n";
                    return;
                }
            }

            // 4) 遍历完井单元，装配 Cell 行和 Well 行
            double rate_diag_sum = 0.0; ///< Rate 模式井行的对角项 Σ PI_i

            const size_t nCells = std::min(
                static_cast<size_t>(Nc),
                std::min(mask_field->data.size(), WI_field->data.size())
            );

            for (size_t cidx = 0; cidx < nCells; ++cidx)
            {
                if ((*mask_field)[cidx] <= 1e-12) continue;  // 未完井单元跳过

                const double WI_i = (*WI_field)[cidx];
                if (WI_i <= 0.0) continue;

                // 4.1 确定该单元的质量动度 λ_mass_i
                double lambda_mass_i = 0.0;

                if (isInjector) {
                    // 注入井：使用井底混合物质量动度（不随位置变）
                    lambda_mass_i = injProps.lambda_mass;
                }
                else {
                    // 生产井：使用储层相动度+密度构造质量动度
                    const double lambda_w_res = (*lambda_w_field)[cidx];
                    const double lambda_g_res = (*lambda_g_field)[cidx];
                    const double rho_w_res = (*rho_w_field)[cidx];
                    const double rho_g_res = (*rho_g_field)[cidx];

                    const double lambda_mass_w_res = lambda_w_res * rho_w_res;
                    const double lambda_mass_g_res = lambda_g_res * rho_g_res;
                    lambda_mass_i = lambda_mass_w_res + lambda_mass_g_res;
                }

                if (lambda_mass_i <= kTiny) continue;

                const double PI_mass_i = WI_i * lambda_mass_i;
                const int cell_lid = lid_of_cell[cidx];
                if (cell_lid < 0) continue;

                // 4.2 Cell 行：+PI_mass_i * p_i - PI_mass_i * p_bh
                sys.addA(cell_lid, cell_lid, PI_mass_i);

                if (well.mode == WellDOF_TwoPhase::Mode::Pressure) {
                    // P_bh 已知，-PI * P_bh -> RHS
                    sys.addb(cell_lid, PI_mass_i * well.target);
                }
                else {
                    // Rate 模式：p_bh 是未知 DOF (well_lid)
                    sys.addA(cell_lid, well_lid, -PI_mass_i);

                    // 井行方程 Σ PI_i (p_bh - p_i) = target_mass
                    // 其中 Σ(PI_i * p_bh) 对应井行对角；Σ(-PI_i * p_i) 是井行的非对角
                    rate_diag_sum += PI_mass_i;
                    sys.addA(well_lid, cell_lid, -PI_mass_i);
                }
            }

            // 5) Rate 模式井行对角和 RHS
            if (well.mode == WellDOF_TwoPhase::Mode::Rate) {
                if (rate_diag_sum <= kTiny) {
                    std::cerr << "[TwoPhaseWellsStrict] Rate-mode well " << well.name
                        << " has zero total PI, skip.\n";
                    return;
                }
                sys.addA(well_lid, well_lid, rate_diag_sum);
                sys.addb(well_lid, well.target);   ///< target 解释为 Σ(dotM_i) = target (kg/s)
            }
        }

        //======================================================================//
        // 2. 统一计算完井单元上的相质量流率（供饱和度/温度方程复用）
        //======================================================================//

        /**
         * @brief 基于已经求解出的 p_bh 和单元压力 p_cell，计算完井单元上的水/气相质量流率。
         *
         * 对 Injector：
         *   - 使用注入混合物质量动度 λ_mass,inj 和质量分率 f_w^mass
         *   - dotM_tot_i = WI_i λ_mass,inj (p_bh - p_cell), 若 <0 则截断为 0（只允许流入）
         *   - dotM_w_i   =  f_w^mass dotM_tot_i
         *   - dotM_g_i   = (1-f_w^mass) dotM_tot_i
         *
         * 对 Producer：
         *   - 使用储层质量动度 λ_mass,res = λ_w ρ_w + λ_g ρ_g
         *   - 分相质量分率 f_w^mass,res = λ_w ρ_w / λ_mass,res
         *   - dotM_tot_i = WI_i λ_mass,res (p_bh - p_cell)，对生产井为负（流出单元）
         *   - dotM_w_i   = f_w^mass,res dotM_tot_i
         *   - dotM_g_i   = (1-f_w^mass,res) dotM_tot_i
         *
         * @param mgr                  网格管理器
         * @param reg                  场变量注册表（提供 p_field, λ_w, λ_g, ρ_w, ρ_g）
         * @param well                 井 DOF（包含 p_bh）
         * @param vg_params            VG 参数
         * @param rp_params            相对渗透率参数
         * @param pressure_field_name  单元压力场名称（通常为 p_w）
         * @param[out] Mw_cell         单元水相质量流率（kg/s），流入单元为正
         * @param[out] Mg_cell         单元气相质量流率（kg/s），流入单元为正
         */
        inline bool compute_well_phase_mass_rates_strict(
            MeshManager& mgr,
            const FieldRegistry& reg,
            const WellDOF_TwoPhase& well,
            const VGParams& vg_params,
            const RelPermParams& rp_params,
            const std::string& pressure_field_name,
            std::vector<double>& Mw_cell,
            std::vector<double>& Mg_cell)
        {
            Mesh& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const int Nc = static_cast<int>(cells.size());

            auto p_field = reg.get<volScalarField>(pressure_field_name.c_str());
            if (!p_field) {
                std::cerr << "[TwoPhaseWellsStrict] Missing pressure field '"
                    << pressure_field_name << "'.\n";
                return false;
            }

            auto mask_field = reg.get<volScalarField>(well.mask_field.c_str());
            auto WI_field = reg.get<volScalarField>(well.PI_field_w.c_str());
            if (!mask_field || !WI_field) {
                std::cerr << "[TwoPhaseWellsStrict] Missing mask/WI field for well "
                    << well.name << ".\n";
                return false;
            }

            auto lambda_w_field = reg.get<volScalarField>("lambda_w");
            auto lambda_g_field = reg.get<volScalarField>("lambda_g");
            auto rho_w_field = reg.get<volScalarField>("rho_w");
            auto rho_g_field = reg.get<volScalarField>("rho_g");

            if (!lambda_w_field || !lambda_g_field || !rho_w_field || !rho_g_field) {
                std::cerr << "[TwoPhaseWellsStrict] Missing lambda/rho fields.\n";
                return false;
            }

            Mw_cell.assign(Nc, 0.0);
            Mg_cell.assign(Nc, 0.0);

            const bool isInjector = (well.role == WellDOF_TwoPhase::Role::Injector);
            InjectionMixtureProps injProps;
            if (isInjector) {
                injProps = computeInjectionMixtureProps(well, vg_params, rp_params);
                if (injProps.lambda_mass <= kTiny) {
                    std::cerr << "[TwoPhaseWellsStrict] Injection mixture has zero "
                        << "mass mobility for well " << well.name << ".\n";
                    return false;
                }
            }

            const double p_bh = well.p_bh;  ///< 假定外层已用压力解更新 well.p_bh

            const auto& id2idx = mesh.getCellId2Index();

            // 用于总质量流量 check（可选）
            double Mtot_well = 0.0;

            for (int ic = 0; ic < Nc; ++ic)
            {
                const auto& c = cells[ic];
                if (c.id < 0) continue;
                const size_t i = id2idx.at(c.id);

                if (i >= mask_field->data.size() ||
                    i >= WI_field->data.size()) continue;

                if ((*mask_field)[i] <= 1e-12) continue; // 非完井

                const double WI_i = (*WI_field)[i];
                if (WI_i <= 0.0) continue;

                const double p_cell = (*p_field)[i];

                // 构造局部质量动度 λ_mass_i
                double lambda_mass_i = 0.0;
                double fw_mass_i = 0.0;

                if (isInjector) {
                    lambda_mass_i = injProps.lambda_mass;
                    fw_mass_i = injProps.fw_mass;
                }
                else {
                    const double lambda_w_res = (*lambda_w_field)[i];
                    const double lambda_g_res = (*lambda_g_field)[i];
                    const double rho_w_res = (*rho_w_field)[i];
                    const double rho_g_res = (*rho_g_field)[i];

                    const double lambda_mass_w_res = lambda_w_res * rho_w_res;
                    const double lambda_mass_g_res = lambda_g_res * rho_g_res;
                    lambda_mass_i = lambda_mass_w_res + lambda_mass_g_res;

                    if (lambda_mass_i > kTiny) {
                        fw_mass_i = lambda_mass_w_res / lambda_mass_i;
                    }
                    else {
                        fw_mass_i = 0.0;
                    }
                }

                if (lambda_mass_i <= kTiny) continue;

                double dotM_tot_i = WI_i * lambda_mass_i * (p_bh - p_cell);

                if (isInjector) {
                    // 注入井：只允许流入地层（流入为正）
                    if (dotM_tot_i < 0.0) dotM_tot_i = 0.0;
                }
                else {
                    // 生产井：dotM_tot_i 应为负（流出单元）
                    if (dotM_tot_i > 0.0) dotM_tot_i = 0.0;
                }

                const double dotM_w_i = fw_mass_i * dotM_tot_i;
                const double dotM_g_i = dotM_tot_i - dotM_w_i;

                // 约定：流入单元为正，流出单元为负
                Mw_cell[ic] += dotM_w_i;
                Mg_cell[ic] += dotM_g_i;

                Mtot_well += dotM_tot_i;
            }

            // 可选：输出 check，特别在 Injector + Rate 模式下，Mtot_well 应接近 well.target
            std::cout << "[TwoPhaseWellsStrict] Well '" << well.name
                << "' total mass rate (computed) = " << Mtot_well
                << " kg/s, target = " << well.target << " kg/s\n";

            return true;
        }

        //======================================================================//
        // 3. 饱和度方程的井源项：水相质量源（显式更新用）
        //======================================================================//

        /**
         * @brief 为所有井构建水相质量源数组 wellSources，用于显式饱和度更新。
         *
         * 约定：
         *   - wellSources[i] 单位为 kg/s，表示“流入单元 i 的水相质量流率”
         *   - 注入井：Mw_cell > 0，生产井：Mw_cell < 0
         *
         * @param mgr          网格管理器
         * @param reg          场注册表
         * @param wells        所有两相井 DOF
         * @param vg_params    VG 参数
         * @param rp_params    相对渗透率参数
         * @param pressure_field_name 压力场名称（如 "p_w"）
         * @param[out] wellSources    单元水相质量源（长度=单元数）
         */
        inline bool build_saturation_well_sources_strict(
            MeshManager& mgr,
            const FieldRegistry& reg,
            const std::vector<WellDOF_TwoPhase>& wells,
            const VGParams& vg_params,
            const RelPermParams& rp_params,
            const std::string& pressure_field_name,
            std::vector<double>& wellSources)
        {
            Mesh& mesh = mgr.mesh();
            const int Nc = static_cast<int>(mesh.getCells().size());
            wellSources.assign(Nc, 0.0);

            std::vector<double> Mw_cell, Mg_cell;

            for (const auto& well : wells)
            {
                if (!compute_well_phase_mass_rates_strict(
                    mgr, reg, well,
                    vg_params, rp_params,
                    pressure_field_name,
                    Mw_cell, Mg_cell))
                {
                    return false;
                }

                // 将每口井的 Mw_cell 叠加到全场 wellSources 中
                for (int ic = 0; ic < Nc; ++ic) {
                    wellSources[ic] += Mw_cell[ic];
                }
            }
            return true;
        }

        //======================================================================//
        // 4. 温度方程井源项：使用同一批质量流量 + 给定注入温度
        //======================================================================//

        /**
         * @brief 将井与温度方程耦合，使用“严格 Rate + Tin” 的质量/能量源。
         *
         * 注入井：
         *   - 使用 compute_well_phase_mass_rates_strict 得到 Mw_cell, Mg_cell (kg/s)
         *   - 对每个完井单元 i，添加热源：Q_h,i = Mw_i cp_w_inj Tin + Mg_i cp_g_inj Tin
         *     （单位 W）
         *
         * 生产井（可选）：
         *   - 使用 Mw_cell, Mg_cell (<0) 构造抽能项：
         *       diag += -(Mw_i cp_w_cell + Mg_i cp_g_cell)
         *     即：M_tot cp T 项线性化到对角，不额外加 RHS
         *
         * @param sysT                 温度方程线性系统
         * @param mgr                  网格管理器
         * @param reg                  场注册表（提供 cp_w, cp_g）
         * @param wells                井 DOF
         * @param vg_params            VG 参数
         * @param rp_params            相对渗透率参数
         * @param pressure_field_name  压力场名称（用于质量流率计算）
         * @param cp_w_field_name      单元水相 cp 字段名（如 "cp_w"）
         * @param cp_g_field_name      单元气相 cp 字段名（如 "cp_g"）
         * @param lid_of_cell          cell.id -> 温度系统自由度编号
         * @param thickness            2D 模型的“厚度” [m]，3D 可设为 1.0
         */
        inline bool couple_wells_to_temperature_equation_strict(
            SparseSystemCOO& sysT,
            MeshManager& mgr,
            const FieldRegistry& reg,
            const std::vector<WellDOF_TwoPhase>& wells,
            const VGParams& vg_params,
            const RelPermParams& rp_params,
            const std::string& pressure_field_name,
            const std::string& cp_w_field_name,
            const std::string& cp_g_field_name,
            const std::vector<int>& lid_of_cell
            )
        {
            Mesh& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const int Nc = static_cast<int>(cells.size());

            auto cp_w_field = reg.get<volScalarField>(cp_w_field_name.c_str());
            auto cp_g_field = reg.get<volScalarField>(cp_g_field_name.c_str());

            if (!cp_w_field || !cp_g_field) {
                std::cerr << "[TwoPhaseWellsStrict] Missing cp fields '"
                    << cp_w_field_name << "', '" << cp_g_field_name << "'.\n";
                return false;
            }

            std::vector<double> Mw_cell(Nc, 0.0), Mg_cell(Nc, 0.0);

            // 逐井耦合（注入/生产分别处理）
            for (const auto& well : wells)
            {
                if (!compute_well_phase_mass_rates_strict(
                    mgr, reg, well,
                    vg_params, rp_params,
                    pressure_field_name,
                    Mw_cell, Mg_cell))
                {
                    return false;
                }

                const bool isInjector = (well.role == WellDOF_TwoPhase::Role::Injector);
                const auto& id2idx = mesh.getCellId2Index();

                for (int ic = 0; ic < Nc; ++ic)
                {
                    const auto& c = cells[ic];
                    if (c.id < 0) continue;
                    const size_t i = id2idx.at(c.id);

                    const int row = lid_of_cell[ic];
                    if (row < 0) continue;

                    const double Mw_i = Mw_cell[ic];  // kg/s，流入单元为正
                    const double Mg_i = Mg_cell[ic];

                    if (std::abs(Mw_i) < 1e-20 && std::abs(Mg_i) < 1e-20)
                        continue;

                    if (isInjector) {
                        // 注入井：固定温度 Tin 的热源项，Q_h = m_dot * cp_inj * Tin
                        const double Q_h =
                            Mw_i * well.cp_w_inj * well.Tin  +
                            Mg_i * well.cp_g_inj * well.Tin;

                        sysT.addb(row, Q_h);
                    }
                    else {
                        // 生产井：抽能，写成对角线性项：-(Mw cp_w + Mg cp_g) * T
                        const double cp_w_i = (*cp_w_field)[i];
                        const double cp_g_i = (*cp_g_field)[i];

                        // 注意 Mw_i/Mg_i 对生产井为负数 => diag_contrib > 0
                        const double diag_contrib =
                            -(Mw_i * cp_w_i + Mg_i * cp_g_i);

                        if (std::abs(diag_contrib) > 0.0) {
                            sysT.addA(row, row, diag_contrib);
                        }
                    }
                }
            }

            return true;
        }

    } // namespace TwoPhaseWellsStrict
} // namespace FVM

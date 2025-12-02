#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "CapRelPerm.h"

#include "SolverContrlStrName.h"
//#include "FluxSplitterandSolver.h"

namespace IMPES_Iteration
{
    /**
     * \brief 配置两相饱和度输运所需的场名称与 VG/相对渗透率参数。
     *
     *  该配置不负责更新物性，仅约定：
     *  - 使用哪个单元标量场作为水相饱和度 S_w；
     *  - 历史时间层 S_w^n 存在哪里；
     *  - （可选）上一外迭代或 RK2 中间层 S_w^{prev} 存在哪里；
     *  - 水相质量通量面场 mf_w 的名称；
     *  - 孔隙度与水相密度的场名称，用于时间项与显式更新；
     *  - VG/相对渗透率模型参数（供外部更新物性时使用）。
     */
    struct SaturationTransportConfig
    {
        SaturationEquation_String           S_Eq_str;
        std::string saturation =            S_Eq_str.saturation;                 ///当前时间层的水相饱和度（正在求解的）
        std::string saturation_old =        S_Eq_str.saturation_old;         ///上一时间步的水相饱和度 （已知解，用于时间项）
        std::string saturation_prev =       S_Eq_str.saturation_prev;       ///外迭代 / RK2 stage 备用的饱和度拷贝（可选）  

        std::string water_mass_flux =       S_Eq_str.water_mass_flux;          // Output: water-phase mass flux
        std::string water_source_field;         /// （可选）水相源汇项，单位 [kg/s]，例如井源的水相质量源 若为空字符串则视为无显式体源项

        TwoPhase_VG_Parameters              VG_Parameter;

        /// 时间步控制开关：SimpleCFL / RedondoLike
        SatTimeControlScheme time_control_scheme = SatTimeControlScheme::SimpleCFL;

        /// Redondo 风格 & ΔS_max 条件使用的参数（SimpleCFL 也会用 dS_max）
        double CFL_safety = 0.9;  ///< C_CFL 安全系数, dt_CFL = CFL_safety * ...
        double dS_max = 0.2;  ///< 单步允许的最大 |ΔS_w|, 用于 ΔS_max 条件
    };

    /**
     * \brief 单个时间步内饱和度推进的诊断数据。
     *
     * max_CFL / max_dS 可用于时间步自适应控制；
     * suggested_dt 可由当前步给出下一个时间步的建议上限。
     */
    struct SaturationStepStats
    {
        double max_CFL = 0.0;   ///< 对流 CFL 最大值（基于水相质量通量）
        double max_dS = 0.0;   ///< 单元内 |ΔS_w| 最大值
        double suggested_dt = 1e100; ///< 由 CFL 与 ΔS_max 约束给出的 dt 建议上限
    };

    namespace detail
    {
        /**
         * \brief clamp S_w into effective range [Swr, 1-Sgr].
         */
        inline double clampSw(double sw, const VGParams& vg)
        {
            const double smin = vg.Swr;
            const double smax = 1.0 - vg.Sgr;
            return std::max(smin, std::min(smax, sw));
        }

        /**
         * \brief 计算给定 S_w 下的质量型分流函数 f_w(S_w)。
         *
         * 使用 Mualem-van Genuchten 相对渗透率模型：
         *   - krw, krg 由 kr_Mualem_vG(sw, vg, rp, krw, krg) 给出；
         *   - λ_w = krw / μ_w, λ_g = krg / μ_g；
         *   - f_w = (λ_w ρ_w) / (λ_w ρ_w + λ_g ρ_g).
         *
         * \param sw     当前水相饱和度（未夹取）
         * \param vg     Van Genuchten 参数
         * \param rp     相对渗透率参数
         * \param mu_w   水相黏度
         * \param mu_g   气相黏度
         * \param rho_w  水相密度
         * \param rho_g  气相密度
         * \return       质量型水相分流函数 f_w \in [0,1]
         */
        inline double fractionalFlow_w_mass(
            double sw,
            const VGParams& vg,
            const RelPermParams& rp,
            double mu_w,
            double mu_g,
            double rho_w,
            double rho_g)
        {
            const double sw_eff = clampSw(sw, vg);
            double krw = 0.0, krg = 0.0;
            kr_Mualem_vG(sw_eff, vg, rp, krw, krg);

            mu_w = std::max(mu_w, 1e-20);
            mu_g = std::max(mu_g, 1e-20);
            rho_w = std::max(rho_w, 0.0);
            rho_g = std::max(rho_g, 0.0);

            const double lamW = krw / mu_w;
            const double lamG = krg / mu_g;

            const double num = lamW * rho_w;
            const double den = std::max(num + lamG * rho_g, 1e-30);

            double fw = num / den;
            if (fw < 0.0) fw = 0.0;
            if (fw > 1.0) fw = 1.0;
            return fw;
        }

        /**
         * \brief 计算 f_w(S_w) 及其导数 df_w/dS_w（有限差分近似）。
         *
         * 使用中心差分：
         *   df/dS ≈ (f(S+δS) - f(S-δS)) / (S+δS - S-δS)
         * 若接近残余/上限饱和度导致 S+δS 与 S-δS 太接近，则退化为前/后向差分，
         * 极端情况下导数视为 0（对 CFL 约束是安全的）。
         */
        inline void fw_and_dfw_mass(
            double sw,
            const VGParams& vg,
            const RelPermParams& rp,
            double mu_w,
            double mu_g,
            double rho_w,
            double rho_g,
            double& fw,
            double& dfw_dSw)
        {
            const double sw_eff = clampSw(sw, vg);
            const double smin = vg.Swr;
            const double smax = 1.0 - vg.Sgr;
            const double epsS = 1e-4;

            double sL = std::max(smin, sw_eff - epsS);
            double sR = std::min(smax, sw_eff + epsS);

            // 当前点的 f_w
            fw = fractionalFlow_w_mass(sw_eff, vg, rp, mu_w, mu_g, rho_w, rho_g);

            const double tiny = 1e-12;
            dfw_dSw = 0.0;

            if (sR - sL < tiny)
            {
                // 几乎到端点，导数退为 0，对稳定性是保守的
                return;
            }

            const double fwL = fractionalFlow_w_mass(sL, vg, rp, mu_w, mu_g, rho_w, rho_g);
            const double fwR = fractionalFlow_w_mass(sR, vg, rp, mu_w, mu_g, rho_w, rho_g);

            dfw_dSw = (fwR - fwL) / (sR - sL);
        }
    } // namespace detail

    /**
     * \brief 显式 Euler 形式推进两相水相饱和度（单一步长）。
     *
     * 离散形式（质量型）：
     * \f[
     *   \phi_i V_i \rho_{w,i} \frac{S_{w,i}^{n+1} - S_{w,i}^n}{\Delta t}
     *   + \sum_{f \in \partial i} F_{w,f}^{(n)}
     *   = Q_{w,i}^{(n)},
     * \f]
     * 其中 \f$F_{w,f}\f$ 为面水相质量通量，采用 owner->neighbor 为正方向。
     * 于是显式更新为：
     * \f[
     *   S_{w,i}^{n+1}
     *   = S_{w,i}^n
     *   + \frac{\Delta t}{\phi_i V_i \rho_{w,i}}
     *     \big( Q_{w,i}^{(n)} - \sum_{f} F_{w,f}^{(n)} \big).
     * \f]
     *
     * @param mgr     网格管理器（提供单元/体积等）
     * @param reg     单元场注册表（S_w、phi_r、rho_w、源项等）
     * @param freg    面场注册表（mf_w）
     * @param cfg     饱和度输运配置（场名等）
     * @param dt      时间步长
     * @param stats   输出本时间步的 CFL / ΔS_{max} / 建议 dt
     * @return true 表示推进成功；false 表示关键字段缺失或几何异常
     */
    inline bool advanceSaturationExplicit_Euler(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const SaturationTransportConfig& cfg,
        double dt,
        SaturationStepStats& stats)
    {
        stats = SaturationStepStats{};

        // ---- 0) 基本检查 ---- //
        if (dt <= 0.0) {
            std::cerr << "[Saturation] advanceSaturationExplicit_Euler: dt <= 0.\n";
            return false;
        }

        auto s_w = reg.get<volScalarField>(cfg.saturation);
        auto s_w_old = reg.get<volScalarField>(cfg.saturation_old);
        if (!s_w || !s_w_old) {
            std::cerr << "[Saturation] missing saturation fields '"
                << cfg.saturation << "' or '" << cfg.saturation_old << "'.\n";
            return false;
        }

        auto phi = reg.get<volScalarField>(TwoPhase::Rock().phi_tag);
        if (!phi) {
            std::cerr << " missing porosity field  <<.\n";
            return false;
        }

        auto rho_w = reg.get<volScalarField>(TwoPhase::Water().rho_tag);
        if (!rho_w)
        {
            std::cerr << " missing density of water field  <<.\n";
            return false;
        }

        auto mf_w = freg.get<faceScalarField>(cfg.water_mass_flux);
        if (!mf_w)
        {
            std::cerr << " missing face fieldof water mass flux  <<.\n";
            return false;
        }

        ///耦合井源模型 提供接口
        std::shared_ptr<volScalarField> q_w;
        if (!cfg.water_source_field.empty()) {
            q_w = reg.get<volScalarField>(cfg.water_source_field.c_str());
            if (!q_w) {
                std::cerr << "[Saturation] requested water source field '"
                    << cfg.water_source_field << "' not found.\n";
                return false;
            }
        }

        Mesh& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& faces = mesh.getFaces();
        const auto& id2idx = mesh.getCellId2Index();

        const size_t nCells = cells.size();
        const size_t nFaces = faces.size();
        if (s_w->data.size() != nCells || s_w_old->data.size() != nCells) {
            std::cerr << "[Saturation] saturation field size mismatch with mesh cells.\n";
            return false;
        }
        if (phi->data.size() != nCells || rho_w->data.size() != nCells) {
            std::cerr << "[Saturation] porosity/rho_w field size mismatch with mesh cells.\n";
            return false;
        } 

        // ---- 1) 构造每个单元的水相质量通量散度 F_div[i] ---- //
        // F_div[i] = sum_{faces} sigma(i,f) * Fw_f
        // 其中 sigma(i,f) = +1 (owner) / -1 (neighbor)

        volScalarField F_div("div_mf_w", nCells, 0.0);

        std::vector<double> sumAbsQtot; // Redondo 方案才会用到
        sumAbsQtot.assign(nCells, 0.0);

        // 如果要用 Redondo, 提前拿 total_vol_flux face 场
        PressureEquation_String P_names;
        std::shared_ptr<faceScalarField> Qf_total;
        if (cfg.time_control_scheme == SatTimeControlScheme::RedondoLike) 
        {
            Qf_total = freg.get<faceScalarField>(P_names.total_vol_flux_name);
            if (!Qf_total) {
                std::cerr << "[Saturation] RedondoLike requires total_vol_flux face field '"
                    << P_names.total_vol_flux_name << "'.\n";
                return false;
            }
        }

        for (const auto& F : faces)
        {
            const int iF = F.id - 1;
            if (iF < 0 || static_cast<size_t>(iF) >= nFaces) continue;

            const double Fw_face = (*mf_w)[iF];

            const int ownerId = F.ownerCell;
            const int neighId = F.neighborCell;

            if (ownerId >= 0)
            {
                const size_t iOwner = id2idx.at(ownerId);
                F_div[iOwner] += Fw_face; // owner -> neighbor 视为外流为正
                if (Qf_total) sumAbsQtot[iOwner] += std::fabs((*Qf_total)[iF]);
            }
            if (neighId >= 0) 
            {
                const size_t iNeigh = id2idx.at(neighId);
                F_div[iNeigh] -= Fw_face; // 对 neighbor 而言，这是流入，取负号
                if (Qf_total) sumAbsQtot[iNeigh] += std::fabs((*Qf_total)[iF]);
            }
        }

        // ---- 2) 显式更新 S_w，并统计 CFL 与 ΔS_max ---- //
        const double CFL_safety = cfg.CFL_safety;
        const double dS_target = cfg.dS_max;
        const double tiny = 1e-30;

        double max_CFL = 0.0;
        double max_dS = 0.0;
        double min_dt_dS = 1e100;
        double min_dt_CFL = 1e100; // SimpleCFL 下不会用到

        // Redondo 模式下才取这些字段
        std::shared_ptr<volScalarField> lambda_w, lambda_g, rho_g, mu_w, mu_g, dPc_dSw;
        if (cfg.time_control_scheme == SatTimeControlScheme::RedondoLike)
        {
            lambda_w = reg.get<volScalarField>(TwoPhase::Water().lambda_w_tag);
            lambda_g = reg.get<volScalarField>(TwoPhase::CO2().lambda_g_tag);
            rho_g = reg.get<volScalarField>(TwoPhase::CO2().rho_tag);
            mu_w = reg.get<volScalarField>(TwoPhase::Water().mu_tag);
            mu_g = reg.get<volScalarField>(TwoPhase::CO2().mu_tag);
            dPc_dSw = reg.get<volScalarField>(TwoPhase::Auxiliaryparameters().dPc_dSw_tag);
            if (!lambda_w || !lambda_g || !rho_g || !mu_w || !mu_g || !dPc_dSw)
            {
                std::cerr << "[Saturation] RedondoLike requires lambda/mu/rho_g/dPc_dSw fields.\n";
                return false;
            }
        }

        const auto& vg = cfg.VG_Parameter.vg_params;
        const auto& rp = cfg.VG_Parameter.relperm_params;

        for (size_t ic = 0; ic < nCells; ++ic) {
            const auto& cell = cells[ic];
            if (cell.id < 0) continue;

            const double V = cell.volume;
            const double phi_i = std::max((*phi)[ic], 1e-12);
            const double rho_i = std::max((*rho_w)[ic], 1e-12);
            const double Sw_n = (*s_w)[ic];
            const double F_div_i = F_div[ic];
            const double Qw_i = q_w ? (*q_w)[ic] : 0.0;

            // ---- 2.1 显式更新 S_w (与旧版完全一样) ---- //
            const double coef = dt / (phi_i * V * rho_i);
            const double dS = coef * (Qw_i - F_div_i);

            double Sw_new = Sw_n + dS;
            Sw_new = std::max(0.0, std::min(1.0, Sw_new));
            Sw_new = clamp(Sw_new, cfg.VG_Parameter.vg_params.Swr, 1 - cfg.VG_Parameter.vg_params.Sgr);
            (*s_w)[ic] = Sw_new;

            const double dS_abs = std::fabs(Sw_new - Sw_n);
            if (dS_abs > max_dS) max_dS = dS_abs;

            // ---- 2.2 ΔS_max 条件 → dt_ΔS,i (两种方案共用) ---- //
            const double denom_loc = std::fabs(Qw_i - F_div_i);
            if (denom_loc > tiny) {
                const double dt_loc = dS_target * phi_i * V * rho_i / denom_loc;
                if (dt_loc < min_dt_dS) min_dt_dS = dt_loc;
            }

            // ---- 2.3 不同时间步控制方案 ---- //
            if (cfg.time_control_scheme == SatTimeControlScheme::SimpleCFL) {
                // 旧方案: CFL_i ≈ dt |F_div| / (φ V ρ), dt_suggest 只看 ΔS_max
                const double CFL_i = std::fabs(F_div_i) * dt / (phi_i * V * rho_i);
                if (CFL_i > max_CFL) max_CFL = CFL_i;
            }
            else { // RedondoLike
                const double lamW_i = std::max((*lambda_w)[ic], 0.0);
                const double lamG_i = std::max((*lambda_g)[ic], 0.0);
                const double lamSum = lamW_i + lamG_i;
                double Psi_i = 0.0;
                if (lamSum > tiny) Psi_i = (lamW_i * lamG_i) / lamSum;

                const double pc_prime = std::fabs((*dPc_dSw)[ic]);
                const double sum_qtot = sumAbsQtot[ic];

                const double rho_g_i = std::max((*rho_g)[ic], 1e-12);
                const double mu_w_i = std::max((*mu_w)[ic], 1e-20);
                const double mu_g_i = std::max((*mu_g)[ic], 1e-20);

                double fw_i = 0.0, dfw_dSw_i = 0.0;
                detail::fw_and_dfw_mass(
                    Sw_n, vg, rp,
                    mu_w_i, mu_g_i,
                    rho_i, rho_g_i,
                    fw_i, dfw_dSw_i);

                const double cap_term = 2.0 * Psi_i * pc_prime;
                const double adv_term = 4.0 * std::fabs(dfw_dSw_i) * sum_qtot;
                const double denom_CFL = cap_term + adv_term;

                if (denom_CFL > tiny) {
                    const double dt_CFL_i = CFL_safety * phi_i * V / denom_CFL;
                    if (dt_CFL_i < min_dt_CFL) min_dt_CFL = dt_CFL_i;

                    const double CFL_i = dt / dt_CFL_i; // Redondo 风格 CFL
                    if (CFL_i > max_CFL) max_CFL = CFL_i;
                }
            }
        }

        // ---- 3) 整体 dt 建议 ---- //
        double dt_suggest = 1e100;
        if (cfg.time_control_scheme == SatTimeControlScheme::SimpleCFL) {
            dt_suggest = min_dt_dS; // 完全保持原行为
        }
        else {
            if (min_dt_CFL < dt_suggest) dt_suggest = min_dt_CFL;
            if (min_dt_dS < dt_suggest) dt_suggest = std::min(dt_suggest, min_dt_dS);
        }

        stats.max_CFL = max_CFL;
        stats.max_dS = max_dS;
        stats.suggested_dt = dt_suggest;
        return true;

    }

    // 预留：RK2 版本接口（后续实现）
   /**
    * \brief 两阶段 RK2 形式推进水相饱和度（接口预留，待实现）。
    *
    * 典型过程：
    *  1. 以 S_w^n 更新物性，解压力 → 得到通量，计算 dS/dt^{(1)}；
    *  2. 用 S_w^{(1)} = S_w^n + dt * dS/dt^{(1)} 更新物性，再解压力 → dS/dt^{(2)}；
    *  3. S_w^{n+1} = S_w^n + dt/2 * (dS/dt^{(1)} + dS/dt^{(2)}).
    *
    * 具体实现将结合现有 assemblePressureTwoPhase 与 splitTwoPhaseMassFlux。
    */
    inline bool advanceSaturation_RK2_placeholder()
    {
        std::cerr << "[Saturation] advanceSaturation_RK2: not implemented yet.\n";
        return false;
    }
}
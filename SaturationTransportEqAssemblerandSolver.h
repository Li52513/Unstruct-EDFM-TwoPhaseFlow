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
#include "FluxSplitterandSolver.h"

namespace IMPES_Iteration
{
    /// 饱和度时间积分格式
    enum class SatTimeIntegrationScheme
    {
        ExplicitEuler,  ///< 一阶显式 Euler
        HeunRK2         ///< 二阶 Heun / RK2
    };

    /// Explicit saturation update: choose which phase mass balance drives the saturation update.
    /// - Water: conserves water-phase mass (when no clamping occurs).
    /// - Gas:   conserves gas-phase (CO2) mass (when no clamping occurs).
    enum class SaturationPrimaryPhase
    {
        Water,
        Gas
    };

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
        std::string water_source_field;                                      /// （可选）水相源汇项，单位 [kg/s]，例如井源的水相质量源 若为空字符串则视为无显式体源项

        std::string gas_source_field;                                        /// (optional) gas/CO2 mass source [kg/s] (e.g. well source)

        TwoPhase_VG_Parameters              VG_Parameter;

        /// 时间积分格式（新加）
        SatTimeIntegrationScheme time_integration_scheme = SatTimeIntegrationScheme::HeunRK2;
        SaturationPrimaryPhase primary_phase = SaturationPrimaryPhase::Water;
        /// 时间步控制开关：SimpleCFL / RedondoLike
        SatTimeControlScheme time_control_scheme = SatTimeControlScheme::RedondoLike;

        /// Redondo 风格 & ΔS_max 条件使用的参数（SimpleCFL 也会用 dS_max）
        double CFL_safety = 0.8;  ///< C_CFL 安全系数, dt_CFL = CFL_safety * ...
        double dS_max = 0.1;  ///< 单步允许的最大 |ΔS_w|, 用于 ΔS_max 条件

        bool useTraditionMethod = false;
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
        int    n_clamped_cells = 0; ///< 本步发生饱和度夹取(clamp)的单元数（夹取会破坏严格守恒）
        double max_abs_clamp_deltaS = 0.0; ///< |S_clamped - S_raw| 的最大值
    };

    namespace detail
    {
        /**
         * \brief clamp S_w into effective range [Swr, 1-Sgr].
         */
        inline bool vgBoundsValid(const VGParams& vg, double eps = 1e-12)
        {
            return (vg.Swr >= 0.0) && (vg.Sgr >= 0.0) && (vg.Swr + vg.Sgr <= 1.0 - 10.0 * eps);
        }

        inline double clampSw(double sw, const VGParams& vg, double eps = 1e-12)
        {
            if (!vgBoundsValid(vg, eps)) {
                // 参数不合法：退回到 [0,1]，同时建议你在日志里报错
                return std::max(0.0, std::min(1.0, sw));
            }
            const double lo = vg.Swr + eps;
            const double hi = 1.0 - vg.Sgr - eps;
            return std::max(lo, std::min(hi, sw));
        }

        inline double clampSg(double sg, const VGParams& vg, double eps = 1e-12)
        {
            if (!vgBoundsValid(vg, eps)) {
                return std::max(0.0, std::min(1.0, sg));
            }
            const double lo = vg.Sgr + eps;
            const double hi = 1.0 - vg.Swr - eps;
            return std::max(lo, std::min(hi, sg));
        }
    } // namespace detail

    /**
     * \brief 显式 Euler 形式推进两相水相饱和度（单一步长），
     *        只依赖 splitTwoPhaseMassFlux 给出的水相质量通量 mf_w。
     *
     * 离散形式（质量型）：
     * \f[
     *   \phi_i V_i \rho_{w,i} \frac{S_{w,i}^{n+1} - S_{w,i}^n}{\Delta t}
     *   + \sum_{f \in \partial i} F_{w,f}^{(n)}
     *   = Q_{w,i}^{(n)},
     * \f]
     * 其中 \f$F_{w,f}\f$ 为面水相质量通量（owner→neighbor 为正方向）。
     * 显式更新为：
     * \f[
     *   S_{w,i}^{n+1}
     *   = S_{w,i}^n
     *   + \frac{\Delta t}{\phi_i V_i \rho_{w,i}}
     *     \big( Q_{w,i}^{(n)} - \sum_{f} F_{w,f}^{(n)} \big).
     * \f]
     *
     * 时间步控制：
     *  - 记录本步 max_dS = max_i |S_{w,i}^{n+1} - S_{w,i}^n|;
     *  - 记录简单 CFL：
     *      CFL_i = |F_div_i| * dt / (phi_i * V_i * rho_{w,i});
     *  - 由两类约束给出建议 dt 上限：
     *      dt_dS,i   = dS_max * phi_i * V_i * rho_{w,i} / |Qw_i - F_div_i|;
     *      dt_CFL,i  = CFL_safety * phi_i * V_i * rho_{w,i} / |F_div_i|;
     *    最终 suggested_dt = min_i( dt_dS,i, dt_CFL,i ).
     *
     * @param mgr     网格管理器（提供单元/体积等）
     * @param reg     单元场注册表（S_w、phi_r、rho_w、源项等）
     * @param freg    面场注册表（mf_w）
     * @param cfg     饱和度输运配置（场名、VG 参数、dS_max、CFL_safety）
     * @param FluxCfg 通量拆分配置（提供 water_mass_flux 名称）
     * @param dt      时间步长
     * @param stats   输出本时间步的 CFL / ΔS_max / 建议 dt
     * @return true  表示推进成功；false 表示关键字段缺失或几何异常
     */
    inline bool advanceSaturationExplicit_Euler(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const SaturationTransportConfig& cfg,
        const FluxSplitConfig& FluxCfg,
        double dt,
        SaturationStepStats& stats)
    {
        stats = SaturationStepStats{};

        // ---- 0) 基本检查 ---- //
        if (dt <= 0.0) {
            std::cerr << "[Saturation] advanceSaturationExplicit_Euler: dt <= 0.\n";
            return false;
        }

        // 当前时间层 S_w 与上一时间步 S_w_old（只做尺寸检查与回滚用）
        auto s_w = reg.get<volScalarField>(cfg.saturation);
        auto s_w_old = reg.get<volScalarField>(cfg.saturation_old);
        if (!s_w || !s_w_old) {
            std::cerr << "[Saturation] missing saturation fields '"
                << cfg.saturation << "' or '" << cfg.saturation_old << "'.\n";
            return false;
        }

        // 孔隙度与水相密度（时间项）
        auto phi = reg.get<volScalarField>(PhysicalProperties_string::Rock().phi_tag);
        auto phi_old = reg.get<volScalarField>(PhysicalProperties_string::Rock().phi_old_tag);
        if (!phi || !phi_old) {
            std::cerr << "[Saturation] missing porosity field (current/old).\n";
            return false;
        }

        const bool advance_water = (cfg.primary_phase == SaturationPrimaryPhase::Water);
        const bool advance_gas = (cfg.primary_phase == SaturationPrimaryPhase::Gas);

        // density for the selected primary phase
        std::shared_ptr<volScalarField> rho = nullptr;
        std::shared_ptr<volScalarField> rho_old = nullptr;
        if (advance_water) {
            rho = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_tag);
            rho_old = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_old_tag);
        }
        else if (advance_gas) {
            rho = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_tag);
            rho_old = reg.get<volScalarField>(PhysicalProperties_string::CO2().rho_old_tag);
        }
        else {
            std::cerr << "[Saturation] unknown cfg.primary_phase.\n";
            return false;
        }
        if (!rho || !rho_old) {
            std::cerr << "[Saturation] missing density field (current/old) for primary phase.\n";
            return false;
        }

        // 水相质量通量面场（由 splitTwoPhaseMassFlux 提供）
        std::shared_ptr<faceScalarField> mf;
        if (advance_water) mf = freg.get<faceScalarField>(FluxCfg.water_mass_flux);
        if (advance_gas)   mf = freg.get<faceScalarField>(FluxCfg.gas_mass_flux);
        if (!mf) {
            std::cerr << "[Saturation] missing face field of primary mass flux.\n";
            return false;
        }

        // 可选：单元体水相质量源（井源等），单位 [kg/s]
        std::shared_ptr<volScalarField> q;
        if (advance_water && !cfg.water_source_field.empty()) {
            q = reg.get<volScalarField>(cfg.water_source_field.c_str());
            if (!q) {
                std::cerr << "[Saturation] requested water source field '"
                    << cfg.water_source_field << "' not found.\n";
                return false;
            }
        }
        if (advance_gas && !cfg.gas_source_field.empty()) {
            q = reg.get<volScalarField>(cfg.gas_source_field.c_str());
            if (!q) {
                std::cerr << "[Saturation] requested gas source field '"
                    << cfg.gas_source_field << "' not found.\n";
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
        if (phi->data.size() != nCells || rho->data.size() != nCells || rho_old->data.size() != nCells) {
            std::cerr << "[Saturation] porosity/rho/rho_old field size mismatch with mesh cells.\n";
            return false;
        }
        if (mf->data.size() != nFaces) {
            std::cerr << "[Saturation] primary mass flux face field size mismatch with mesh faces.\n";
            return false;
        }
        if (q && q->data.size() != nCells) {
            std::cerr << "[Saturation] primary source field size mismatch with mesh cells.\n";
            return false;
        }

        // ---- 1) 构造每个单元的相质量通量散度 F_div[i] ---- //
        // F_div[i] = sum_{faces} sigma(i,f) * Fw_f
        // 其中 sigma(i,f) = +1 (owner) / -1 (neighbor)
        volScalarField F_div("div_mf_primary", nCells, 0.0);

        for (const auto& F : faces)
        {
            const int iF = F.id - 1;
            if (iF < 0 || static_cast<size_t>(iF) >= nFaces) continue;

            const double F_face = (*mf)[iF];
            const int ownerId = F.ownerCell;
            const int neighId = F.neighborCell;

            if (ownerId >= 0)
            {
                const size_t iOwner = id2idx.at(ownerId);
                // owner -> neighbor 为正：对 owner 是“流出”贡献
                F_div[iOwner] += F_face;
            }
            if (neighId >= 0)
            {
                const size_t iNeigh = id2idx.at(neighId);
                // 对 neighbor 而言 same Fw_face 是“流入”，取负号
                F_div[iNeigh] -= F_face;
            }
        }

        // ---- 2) 显式更新 S_w，并统计 CFL 与 ΔS_max、建议 dt ---- //
        const double CFL_safety = cfg.CFL_safety;
        const double dS_target = cfg.dS_max;
        const double tiny = 1e-30;
        const auto& vg = cfg.VG_Parameter.vg_params;

        double max_CFL = 0.0;
        double max_dS = 0.0;
        double min_dt_dS = 1e100;
        double min_dt_CFL = 1e100;

        for (size_t ic = 0; ic < nCells; ++ic)
        {
            const auto& cell = cells[ic];
            if (cell.id < 0) continue;

            const double V = cell.volume;
            if (V <= tiny) {
                std::cerr << "[Saturation] cell volume too small at cell id="
                    << cell.id << "\n";
                continue;
            }

            const double phi_np1 = std::max((*phi)[ic], 1e-12);
            const double phi_n = std::max((*phi_old)[ic], 1e-12);
            const double rho_np1 = std::max((*rho)[ic], 1e-12);
            const double rho_n = std::max((*rho_old)[ic], 1e-12);
            const double Sw_n = (*s_w)[ic];

            const double F_div_i = F_div[ic];
            const double Q_i = q ? (*q)[ic] : 0.0;

            // 2.1 显式更新（质量守恒形式）：
            //   M_w^{n+1} = M_w^n + dt*(Qw - div(Fw)),
            //   S_w^{n+1} = M_w^{n+1} / (phi * V * rho_w^{n+1})
            //
            // 这样在 rho_w 随压力变化时，也能保证（不发生 clamp 时）离散水相质量严格守恒。
            double Sw_raw = Sw_n;
            double Sw_new = Sw_n;
            if (cfg.useTraditionMethod)
            {
                // 传统显示更新（不严格守恒）
				const double coef = dt / (phi_np1 * V * rho_np1);
				const double dS = coef * (Q_i - F_div_i);
				Sw_raw = Sw_n + dS;
                Sw_new = detail::clampSw(Sw_raw, vg);
            }
            else
            {
				// 新方法：基于相质量守恒的显示更新
                if (advance_water)
                {
                    const double Mw_n = phi_n * V * rho_n * Sw_n;
                    const double Mw_np1 = Mw_n + dt * (Q_i - F_div_i);

                    const double denom = std::max(phi_np1 * V * rho_np1, tiny);
                    Sw_raw = Mw_np1 / denom;
                    Sw_new = detail::clampSw(Sw_raw, vg);
                }
                else
                {
                    const double Sg_n = 1.0 - Sw_n;
                    const double Mg_n = phi_n * V * rho_n * Sg_n;
                    const double Mg_np1 = Mg_n + dt * (Q_i - F_div_i);

                    const double denom = std::max(phi_np1 * V * rho_np1, tiny);
                    const double Sg_raw = Mg_np1 / denom;
                    const double Sg_new = detail::clampSg(Sg_raw, vg);

                    Sw_raw = 1.0 - Sg_raw;
                    Sw_new = 1.0 - Sg_new;

                }
            }
            const double clamp_delta = std::abs(Sw_new - Sw_raw);
            if (clamp_delta > 1e-4) {
                ++stats.n_clamped_cells;
                if (clamp_delta > stats.max_abs_clamp_deltaS) stats.max_abs_clamp_deltaS = clamp_delta;
            }

            (*s_w)[ic] = Sw_new;

            // 2.2 统计 |ΔS| 与 ΔS_max 约束对应的 dt
            const double dS_abs = std::fabs(Sw_new - Sw_n);
            if (dS_abs > max_dS) max_dS = dS_abs;

            const double denom_loc = std::fabs(Q_i - F_div_i);
            if (denom_loc > tiny) {
                const double dt_loc = dS_target * phi_np1 * V * rho_np1 / denom_loc;
                if (dt_loc < min_dt_dS) min_dt_dS = dt_loc;
            }

            // 2.3 简单 CFL：基于对流散度 F_div
            const double CFL_i = std::fabs(F_div_i) * dt / (phi_np1 * V * rho_np1);
            if (CFL_i > max_CFL) max_CFL = CFL_i;

            const double adv_denom = std::fabs(F_div_i);
            if (adv_denom > tiny) {
                const double dt_CFL_i = CFL_safety * phi_np1 * V * rho_np1 / adv_denom;
                if (dt_CFL_i < min_dt_CFL) min_dt_CFL = dt_CFL_i;
            }
        }

        // ---- 3) 根据两类约束给出整体 dt 建议 ---- //
        double dt_suggest = 1e100;
        if (min_dt_dS < dt_suggest) dt_suggest = min_dt_dS;
        if (min_dt_CFL < dt_suggest) dt_suggest = std::min(dt_suggest, min_dt_CFL);

        // 如果没有任何约束起作用，就退化为当前 dt
        if (!(dt_suggest > 0.0 && dt_suggest < 1e90)) {
            dt_suggest = dt;
        }

        stats.max_CFL = max_CFL;
        stats.max_dS = max_dS;
        stats.suggested_dt = dt_suggest;

        // 最后一刀：确保全场在有效区间（防止某些分支漏 clamp）
        for (size_t ic = 0; ic < nCells; ++ic) {
            (*s_w)[ic] = detail::clampSw((*s_w)[ic], cfg.VG_Parameter.vg_params);
        }

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
    /**
        * \brief Heun / RK2 形式推进水相饱和度（单一步长），
        *        需要在两个 stage 之间重建水相质量通量 mf_w。
        *
        * 形式上：
        *   dS/dt = R(S),  R(S) = (Qw - div Fw) / (phi V rho_w)
        *   k1    = R(S^n)
        *   S*    = S^n + dt * k1
        *   k2    = R(S*)
        *   S^{n+1} = S^n + dt/2 * (k1 + k2)
        *
        * 使用约定：
        *  - 调用本函数前，调用者应基于当前 S_w 构建好 mf_w（FluxCfg.water_mass_flux）。
        *  - rebuildFlux(stageId) 会在 stage=1 时被调用一次，此时 cfg.saturation
        *    中已经保存了预测场 S*；rebuildFlux 负责：
        *      * 更新相对渗透率/动度等物性；
        *      * 拆分两相通量，更新 water_mass_flux 面场；
        *      * 如需，可同步更新水相体源项 q_w（例如井源 Qw_well）。
        *
        * @tparam FluxRebuilder  可调用对象，签名为 bool(int stageId)
        */
   
    inline bool advanceSaturationHeun_RK2(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const SaturationTransportConfig& cfg,
        const FluxSplitConfig& FluxCfg,
        double dt,
        SaturationStepStats& stats)
    {
        stats = SaturationStepStats{};

        if (dt <= 0.0) {
            std::cerr << "[Saturation] advanceSaturationHeun_RK2: dt <= 0.\n";
            return false;
        }

        // --- 0) 基本场访问与尺寸检查 --- //
        auto s_w = reg.get<volScalarField>(cfg.saturation);
        auto s_w_old = reg.get<volScalarField>(cfg.saturation_old);   // 仅作尺寸检查/回滚用
        auto s_w_prev = reg.get<volScalarField>(cfg.saturation_prev);  // 可选：存储 S^n

        if (!s_w || !s_w_old) {
            std::cerr << "[Saturation] advanceSaturationHeun_RK2: missing saturation fields '"
                << cfg.saturation << "' or '" << cfg.saturation_old << "'.\n";
            return false;
        }

        auto phi = reg.get<volScalarField>(PhysicalProperties_string::Rock().phi_tag);
        auto rho_w = reg.get<volScalarField>(PhysicalProperties_string::Water().rho_tag);
        if (!phi || !rho_w) {
            std::cerr << "[Saturation] advanceSaturationHeun_RK2: missing porosity or rho_w.\n";
            return false;
        }

        auto mf_w = freg.get<faceScalarField>(FluxCfg.water_mass_flux);
        if (!mf_w) {
            std::cerr << "[Saturation] advanceSaturationHeun_RK2: missing water mass flux field '"
                << FluxCfg.water_mass_flux << "'.\n";
            return false;
        }

        std::shared_ptr<volScalarField> q_w;
        if (!cfg.water_source_field.empty()) {
            q_w = reg.get<volScalarField>(cfg.water_source_field.c_str());
            if (!q_w) {
                std::cerr << "[Saturation] advanceSaturationHeun_RK2: requested water source field '"
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
            std::cerr << "[Saturation] advanceSaturationHeun_RK2: saturation size mismatch.\n";
            return false;
        }
        if (phi->data.size() != nCells || rho_w->data.size() != nCells) {
            std::cerr << "[Saturation] advanceSaturationHeun_RK2: porosity/rho_w size mismatch.\n";
            return false;
        }
        if (mf_w->data.size() != nFaces) {
            std::cerr << "[Saturation] advanceSaturationHeun_RK2: mf_w size mismatch.\n";
            return false;
        }
        if (q_w && q_w->data.size() != nCells) {
            std::cerr << "[Saturation] advanceSaturationHeun_RK2: q_w size mismatch.\n";
            return false;
        }

        // --- 0.1) 备份当前 S^n（到局部数组 + 可选的 saturation_prev） --- //
        std::vector<double> S_old(nCells, 0.0);
        for (size_t ic = 0; ic < nCells; ++ic) {
            S_old[ic] = (*s_w)[ic];
        }
        if (s_w_prev) {
            s_w_prev->data = S_old; // 方便外部需要访问 S^n 时使用
        }

        // --- 通用小工具：从当前 mf_w 构建 F_div --- //
        auto build_div_from_mf = [&](faceScalarField& mf, volScalarField& F_div) {
            std::fill(F_div.data.begin(), F_div.data.end(), 0.0);
            for (const auto& F : faces) {
                const int iF = F.id - 1;
                if (iF < 0 || static_cast<size_t>(iF) >= nFaces) continue;

                const double Fw_face = mf[iF];
                const int ownerId = F.ownerCell;
                const int neighId = F.neighborCell;

                if (ownerId >= 0) {
                    const size_t iOwner = id2idx.at(ownerId);
                    F_div[iOwner] += Fw_face;  // owner->neighbor 为正，视为 owner 的“流出”
                }
                if (neighId >= 0) {
                    const size_t iNeigh = id2idx.at(neighId);
                    F_div[iNeigh] -= Fw_face;  // 对 neighbor 则为“流入”
                }
            }
            };

        const double CFL_safety = cfg.CFL_safety;
        const double dS_target = cfg.dS_max;
        const double tiny = 1e-30;
        // --- Stage 1: 基于 S^n 和 mf_w^n 计算 k1 并得到预测 S* --- //
        volScalarField F_div_stage1("div_mf_w_stage1", nCells, 0.0);
        build_div_from_mf(*mf_w, F_div_stage1);

        std::vector<double> k1(nCells, 0.0);
        const auto& vg = cfg.VG_Parameter.vg_params;


        double max_CFL = 0.0;
        double max_dS_stage1 = 0.0;
        double min_dt_dS = 1e100;
        double min_dt_CFL = 1e100;

        for (size_t ic = 0; ic < nCells; ++ic)
        {
            const auto& cell = cells[ic];
            if (cell.id < 0) continue;

            const double V = cell.volume;
            if (V <= tiny) {
                std::cerr << "[Saturation] advanceSaturationHeun_RK2: cell volume too small at id="
                    << cell.id << "\n";
                continue;
            }

            const double phi_i = std::max((*phi)[ic], 1e-12);
            const double rho_i = std::max((*rho_w)[ic], 1e-12);
            const double Sw_n = S_old[ic];

            const double F_div_i = F_div_stage1[ic];
            const double Qw_i = q_w ? (*q_w)[ic] : 0.0;

            const double mass_coeff = phi_i * V * rho_i;
            const double RHS1 = (Qw_i - F_div_i);
            const double rate1 = RHS1 / std::max(mass_coeff, tiny); // dS/dt
            k1[ic] = rate1;

            // 预测 S*
            double S_star = Sw_n + dt * rate1;
            S_star = detail::clampSw(S_star, vg);

            (*s_w)[ic] = S_star;

            // ΔS_stage1 仅用于诊断，不是最终 ΔS
            const double dS_stage1 = std::fabs(S_star - Sw_n);
            if (dS_stage1 > max_dS_stage1) max_dS_stage1 = dS_stage1;

            // dt_dS 约束：基于 |RHS1|
            const double denom_loc = std::fabs(RHS1);
            if (denom_loc > tiny) {
                const double dt_loc = dS_target * mass_coeff / denom_loc;
                if (dt_loc < min_dt_dS) min_dt_dS = dt_loc;
            }

            // CFL 约束：基于 |F_div|
            const double CFL_i = std::fabs(F_div_i) * dt / std::max(mass_coeff, tiny);
            if (CFL_i > max_CFL) max_CFL = CFL_i;

            const double adv_denom = std::fabs(F_div_i);
            if (adv_denom > tiny) {
                const double dt_CFL_i = CFL_safety * mass_coeff / adv_denom;
                if (dt_CFL_i < min_dt_CFL) min_dt_CFL = dt_CFL_i;
            }
        }
        mf_w = freg.get<faceScalarField>(FluxCfg.water_mass_flux);
        if (!mf_w || mf_w->data.size() != nFaces) {
            std::cerr << "[Saturation] advanceSaturationHeun_RK2: mf_w after rebuild invalid.\n";
            return false;
        }

        // --- Stage 2: 基于 S* 和 mf_w^* 计算 k2，并合成 S^{n+1} --- //
        volScalarField F_div_stage2("div_mf_w_stage2", nCells, 0.0);
        build_div_from_mf(*mf_w, F_div_stage2);

        std::vector<double> k2(nCells, 0.0);

        double max_dS_final = 0.0;  // 最终 |S^{n+1} - S^n|
        double min_dt_dS2 = 1e100;
        double min_dt_CFL2 = 1e100;
        double max_CFL2 = 0.0;

        for (size_t ic = 0; ic < nCells; ++ic)
        {
            const auto& cell = cells[ic];
            if (cell.id < 0) continue;

            const double V = cell.volume;
            if (V <= tiny) continue;

            const double phi_i = std::max((*phi)[ic], 1e-12);
            const double rho_i = std::max((*rho_w)[ic], 1e-12);
            const double Sw_n = S_old[ic];

            const double F_div_i2 = F_div_stage2[ic];
            const double Qw_i = q_w ? (*q_w)[ic] : 0.0;

            const double mass_coeff = phi_i * V * rho_i;
            const double RHS2 = (Qw_i - F_div_i2);
            const double rate2 = RHS2 / std::max(mass_coeff, tiny);
            k2[ic] = rate2;

            // Heun 合成：S^{n+1} = S^n + dt/2 * (k1 + k2)
            double S_new = Sw_n + 0.5 * dt * (k1[ic] + k2[ic]);
            S_new = std::max(0.0, std::min(1.0, S_new));

            (*s_w)[ic] = S_new;

            const double dS_final = std::fabs(S_new - Sw_n);
            if (dS_final > max_dS_final) max_dS_final = dS_final;

            const double denom_loc = std::max(std::fabs(RHS2), tiny);
            const double dt_loc = dS_target * mass_coeff / denom_loc;
            if (dt_loc < min_dt_dS2) min_dt_dS2 = dt_loc;

            const double CFL_i2 = std::fabs(F_div_i2) * dt / std::max(mass_coeff, tiny);
            if (CFL_i2 > max_CFL2) max_CFL2 = CFL_i2;

            const double adv_denom2 = std::fabs(F_div_i2);
            if (adv_denom2 > tiny) {
                const double dt_CFL_i2 = CFL_safety * mass_coeff / adv_denom2;
                if (dt_CFL_i2 < min_dt_CFL2) min_dt_CFL2 = dt_CFL_i2;
            }
        }

        // --- 3) 汇总 CFL / ΔS / dt 建议 --- //
        stats.max_CFL = std::max(max_CFL, max_CFL2);
        stats.max_dS = max_dS_final;

        double dt_suggest = 1e100;
        if (min_dt_dS < dt_suggest) dt_suggest = min_dt_dS;
        if (min_dt_dS2 < dt_suggest) dt_suggest = min_dt_dS2;
        if (min_dt_CFL < dt_suggest) dt_suggest = min_dt_CFL;
        if (min_dt_CFL2 < dt_suggest) dt_suggest = min_dt_CFL2;

        if (!(dt_suggest > 0.0 && dt_suggest < 1e90)) {
            dt_suggest = dt; // 回退：无有效约束时保持当前 dt
        }
        stats.suggested_dt = dt_suggest;

        // 最后一刀：确保全场在有效区间（防止某些分支漏 clamp）
        for (size_t ic = 0; ic < nCells; ++ic) {
            (*s_w)[ic] = detail::clampSw((*s_w)[ic], cfg.VG_Parameter.vg_params);
        }

        return true;
    }
}

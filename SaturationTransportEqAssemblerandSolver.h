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
            }
            if (neighId >= 0) 
            {
                const size_t iNeigh = id2idx.at(neighId);
                F_div[iNeigh] -= Fw_face; // 对 neighbor 而言，这是流入，取负号
            }
        }

        // ---- 2) 显式更新 S_w，并统计 CFL 与 ΔS_max ---- //
        const double CFL_target = 0.5;     // 缺省 CFL 阈值，可在外部配置
        const double dS_target = 0.1;     // 缺省 ΔS_max 阈值（用于建议时间步）

        double dt_suggest = 1e100;
        double max_CFL = 0.0;
        double max_dS = 0.0;

        for (size_t ic = 0; ic < nCells; ++ic)
        {
            const auto& cell = cells[ic];
            if (cell.id < 0) continue;
            const double V = cell.volume;
            if (V <= 0.0)
            {
                std::cerr << "[Saturation] cell volume <= 0, cell.id = "
                    << cell.id << "\n";
                return false;
            }
            const double phi_i = std::max((*phi)[ic], 1e-12);
            const double rho_i = std::max((*rho_w)[ic], 1e-12);
            const double Sw_old = (*s_w_old)[ic];
            double Sw_n = (*s_w)[ic];

            const double F_div_i = F_div[ic];
            const double Qw_i = q_w ? (*q_w)[ic] : 0.0;

            // 显式更新公式：
           // Sw^{n+1} = Sw^n + dt / (phi V rho) * (Qw - F_div)
            const double coef = dt / (phi_i * V * rho_i);
            const double dS = coef * (Qw_i - F_div_i);

            double Sw_new = Sw_n + dS;

            // 简单夹取：可结合 VG 残余饱和度参数进一步改进
            Sw_new = std::max(0.0, std::min(1.0, Sw_new));

            (*s_w)[ic] = Sw_new;

            const double dS_abs = std::fabs(Sw_new - Sw_n);
            if (dS_abs > max_dS) max_dS = dS_abs;

            // 简单 CFL 估计：CFL_i ≈ dt * |F_div_i| / (phi V rho)
            const double CFL_i = std::fabs(F_div_i) * dt / (phi_i * V * rho_i);
            if (CFL_i > max_CFL) max_CFL = CFL_i;

            // 根据 ΔS_target 给出局部 dt 建议：dt_i ≈ dS_target * phi V rho / |Qw - F_div|
            const double denom_loc = std::fabs(Qw_i - F_div_i);
            if (denom_loc > 0.0) {
                const double dt_loc = dS_target * phi_i * V * rho_i / denom_loc;
                if (dt_loc < dt_suggest) dt_suggest = dt_loc;
            }
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
#pragma once
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FVM_WellDOF_TwoPhase.h"
#include "CapRelPerm.h"            // VGParams, RelPermParams
#include "TwoPhaseWells_StrictRate.h" // PeacemanTwoPhase + TwoPhaseWellsStrict

namespace IMPES_Iteration
{
    /**
     * @brief 基于 Peaceman 严格 Rate 两相井，构建水相质量源场 Qw_well [kg/s]。
     *
     * 约定：
     *  - Qw_well[i] > 0 : 向单元 i 注入水相
     *  - Qw_well[i] < 0 : 从单元 i 采出水相
     *
     * @param mgr        网格管理器
     * @param reg        场注册表（提供 p_w, lambda_w, lambda_g, rho_w, rho_g 等）
     * @param wells      两相井 DOF 列表（p_bh 已由压力解更新）
     * @param vg_params  VG 参数
     * @param rp_params  相对渗透率参数
     * @param pw_field   单元压力场名（通常为 "p_w"）
     * @param Qw_name    输出水相源场名，例如 "Qw_well"
     */
    inline bool buildWaterSourceFieldFromWells(
        MeshManager& mgr,
        FieldRegistry& reg,
        const std::vector<WellDOF_TwoPhase>& wells,
        const VGParams& vg_params,
        const RelPermParams& rp_params,
        const std::string& pw_field,
        const std::string& Qw_name)
    {
        Mesh& mesh = mgr.mesh();
        const int Nc = static_cast<int>(mesh.getCells().size());

        // 1) 调用井模块，生成 cell 级水相质量源数组
        std::vector<double> wellSources;
        if (!FVM::TwoPhaseWellsStrict::build_saturation_well_sources_strict(
            mgr, reg, wells,
            vg_params, rp_params,
            pw_field,
            wellSources))
        {
            std::cerr << "[IMPES][Well] build_saturation_well_sources_strict failed.\n";
            return false;
        }

        if ((int)wellSources.size() != Nc) {
            std::cerr << "[IMPES][Well] wellSources size mismatch.\n";
            return false;
        }

        // 2) 写入 FieldRegistry 的体源场
        auto Qw = reg.getOrCreate<volScalarField>(Qw_name.c_str(), Nc, 0.0);
        if (!Qw) {
            std::cerr << "[IMPES][Well] failed to create Qw_well field '"
                << Qw_name << "'.\n";
            return false;
        }
        Qw->data.assign(Nc, 0.0);
        for (int ic = 0; ic < Nc; ++ic) {
            (*Qw)[ic] = wellSources[ic]; // kg/s, 注入正，采出负
        }

        return true;
    }

    /// \brief 调试用：打印单口生产井当前的 Peaceman 质量指数和预期流量。
    ///
    /// \param mgr   网格管理器
    /// \param reg   场注册表（提供 p_w, lambda_w, lambda_g, rho_w, rho_g, mask, WI）
    /// \param well  目标井（通常为生产井，mode=Pressure）
    /// \param pw_field_name  压力场名称（如 "p_w"）
    ///
    /// 调用时机：压力方程收敛、well.p_bh 已更新之后。
inline void debugPrintPeacemanPIAndRate(
    MeshManager& mgr,
    const FieldRegistry& reg,
    const WellDOF_TwoPhase& well,
    const std::string& pw_field_name)
{
    Mesh& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();
    const int Nc = static_cast<int>(cells.size());

    auto pF   = reg.get<volScalarField>(pw_field_name.c_str());
    auto mask = reg.get<volScalarField>(well.mask_field.c_str());
    auto WI   = reg.get<volScalarField>(well.PI_field_w.c_str());

    auto lambda_w = reg.get<volScalarField>(TwoPhase::Water().lambda_w_tag);
    auto lambda_g = reg.get<volScalarField>(TwoPhase::CO2().lambda_g_tag);
    auto rho_w    = reg.get<volScalarField>(TwoPhase::Water().rho_tag);
    auto rho_g    = reg.get<volScalarField>(TwoPhase::CO2().rho_tag);

    if (!pF || !mask || !WI || !lambda_w || !lambda_g || !rho_w || !rho_g) {
        std::cerr << "[DebugPI] missing fields for well " << well.name << "\n";
        return;
    }

    const double p_bh = well.p_bh; // Pressure 模式下应等于 well.target

    double sum_PI_mass = 0.0;
    double sum_dotM    = 0.0;
    double sum_deltap  = 0.0;
    int    nPerf       = 0;

    for (int ic = 0; ic < Nc; ++ic) {
        const auto& c = cells[ic];
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);

        if (i >= mask->data.size() || i >= WI->data.size()) continue;
        if ((*mask)[i] <= 1e-12) continue;

        const double WI_i = (*WI)[i];
        if (WI_i <= 0.0) continue;

        const double p_cell = (*pF)[i];

        const double lambda_w_res = (*lambda_w)[i];
        const double lambda_g_res = (*lambda_g)[i];
        const double rho_w_res    = (*rho_w)[i];
        const double rho_g_res    = (*rho_g)[i];

        const double lambda_mass_i =
            lambda_w_res * rho_w_res + lambda_g_res * rho_g_res;

        if (lambda_mass_i <= 1e-30) continue;

        const double PI_mass_i = WI_i * lambda_mass_i;
        const double dp        = p_bh - p_cell;
        const double dotM_i    = PI_mass_i * dp;

        sum_PI_mass += PI_mass_i;
        sum_dotM    += dotM_i;
        sum_deltap  += dp;
        ++nPerf;
    }

    if (nPerf == 0) {
        std::cout << "[DebugPI] well " << well.name
                  << " has no perforated cells.\n";
        return;
    }

    const double avg_dp = sum_deltap / static_cast<double>(nPerf);

    std::cout << "[DebugPI] Well '" << well.name << "' (mode=Pressure)\n"
              << "  p_bh       = " << p_bh << " Pa\n"
              << "  nPerf      = " << nPerf << "\n"
              << "  Sum PI_mass= " << sum_PI_mass << " (kg/(s·Pa))\n"
              << "  Avg Δp     = " << avg_dp << " Pa\n"
              << "  Est Q_mass = " << sum_dotM << " kg/s\n";
}

} // namespace IMPES_Iteration

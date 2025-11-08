#pragma once
#include <string>
#include <cmath>
#include <iostream>
#include <numbers>
#include <memory> 
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "VolField.h"
#include "FVM_SourceTerm_InjectionMask.h"  // InjectionCircle + 已构建好的 "inj_mask"

namespace FVM {
    namespace SourceTerm {

        // —— 小工具：π（跨平台而不依赖 M_PI 宏）——
        inline double pi_value() {
#if __cpp_lib_math_constants >= 201907L
            return std::numbers::pi_v<double>;
#else
            return 3.141592653589793238462643383279502884;
#endif
        }

        // ========== 选项与枚举 ==========
        struct PeacemanPIOptions {
            double thickness = 1.0;    // 2D 单位厚度
            double skin = 0.0;    // 皮肤因子 s（可为负；数值将兜底）
            double beta_re = 0.21;   // re = beta_re * sqrt(A_cell)
            double min_log_arg = 1e-12;  // ln 入参兜底
            bool   clamp_nonpos = true;  // 负/NaN/异常导纳钳 0
            bool   verbose;  // 打印日志
        };

        enum class FluidPropertyMode { Constant, Field };

        // 从体场估计常量的配置
        struct PIConstantEstimateOpts {
            // 字段名（按你工程中的约定）
            std::string k_iso_name = "k";    // 若存在等向渗透率场，优先
            std::string kxx_name = "kxx";  // 各向异性主值
            std::string kyy_name = "kyy";
            std::string rho_name = "rho_g";
            std::string mu_name = "mu_g";

            bool   volume_weight = true;   // 2D: 用 cell.volume(面积) 做权重；3D: 用体积
            bool   use_mask = false;  // 若 true，仅使用 mask>0 区域加权
            std::string mask_name = "inj_mask";

            double eps_k = 1e-30;  // 渗透率/对数/黏度兜底
            double k_default = 1e-14;
            double rho_default = 800.0;
            double mu_default = 1.48e-5;

            bool   verbose;
        };

        // ========== 接口 1：从体场估计 k_const / rho_const / mu_const ==========
        inline bool estimate_PI_constants_from_fields(
            ::MeshManager& mgr,
            ::FieldRegistry& reg,
            double& k_const,
            double& rho_const,
            double& mu_const,
            const PIConstantEstimateOpts& opt = PIConstantEstimateOpts{}
        ) {
            ::Mesh& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const auto& id2i = mesh.getCellId2Index();

            // 取字段指针（可能为 nullptr）
            auto kIsoF = reg.get<::volScalarField>(opt.k_iso_name);
            auto kxxF = reg.get<::volScalarField>(opt.kxx_name);
            auto kyyF = reg.get<::volScalarField>(opt.kyy_name);
            auto rhoF = reg.get<::volScalarField>(opt.rho_name);
            auto muF = reg.get<::volScalarField>(opt.mu_name);

            const bool hasIso = (bool)kIsoF;
            const bool hasKxy = (bool)kxxF && (bool)kyyF;
            const bool hasRho = (bool)rhoF;
            const bool hasMu = (bool)muF;

            std::shared_ptr<::volScalarField> maskF;
            if (opt.use_mask) maskF = reg.get<::volScalarField>(opt.mask_name);

            double Wsum = 0.0, sum_k = 0.0, sum_rho = 0.0, sum_mu = 0.0;
            size_t cnt = 0;

            for (const auto& c : cells) {
                if (c.id < 0) continue;                 // 跳过 ghost/无效
                const size_t i = id2i.at(c.id);

                // 可选的 mask 过滤
                if (maskF && (*maskF)[i] <= 0.0) continue;

                // 权重：2D->面积；3D->体积
                double w = opt.volume_weight ? std::max(c.volume, 0.0) : 1.0;
                if (w <= 0.0) continue;

                // k_eff：优先使用 k（等向）；否则 sqrt(kxx*kyy)
                double keff_i = opt.k_default;
                if (hasIso) {
                    keff_i = std::max((*kIsoF)[i], opt.eps_k);
                }
                else if (hasKxy) {
                    const double kx = std::max((*kxxF)[i], opt.eps_k);
                    const double ky = std::max((*kyyF)[i], opt.eps_k);
                    keff_i = std::sqrt(kx * ky);
                }

                const double rho_i = hasRho ? (*rhoF)[i] : opt.rho_default;
                const double mu_i = hasMu ? (*muF)[i] : opt.mu_default;

                sum_k += w * keff_i;
                sum_rho += w * rho_i;
                sum_mu += w * mu_i;
                Wsum += w;
                ++cnt;
            }

            if (Wsum <= 0.0 || cnt == 0) {
                k_const = opt.k_default;
                rho_const = opt.rho_default;
                mu_const = opt.mu_default;
                if (opt.verbose) {
                    std::cerr << "[PI-const] fields missing/empty; fallback to defaults: "
                        << "k=" << k_const << ", rho=" << rho_const << ", mu=" << mu_const << "\n";
                }
                return true;
            }

            k_const = sum_k / Wsum;
            rho_const = sum_rho / Wsum;
            mu_const = sum_mu / Wsum;

            if (opt.verbose) {
                std::cout << "[PI-const] k_const=" << k_const
                    << " (m^2), rho_const=" << rho_const
                    << " (kg/m^3), mu_const=" << mu_const
                    << " (Pa·s)  weight=" << (opt.volume_weight ? "volume" : "uniform")
                    << (opt.use_mask ? " [masked]" : "") << "\n";
            }
            return true;
        }

        // ========== 接口 2：构建 Peaceman-like 质量导纳 PI ==========
        inline bool build_peaceman_PI
        (
            ::MeshManager& mgr,
            ::FieldRegistry& reg,
            const wellCircle& well,           // 圆井几何（用于 r_w；re 用单元面积）
            const std::string& inj_mask_name,      // "inj_mask"（Step1 已生成）
            const std::string& out_PI_name,        // "PI"
            double k_const, double rho_const, double mu_const,   // 常物性（或下方传场）
            FluidPropertyMode fp_mode = FluidPropertyMode::Constant,
            const std::shared_ptr<::volScalarField>  rho_field = nullptr,
            const std::shared_ptr<::volScalarField>  mu_field = nullptr,
            const PeacemanPIOptions& opt = PeacemanPIOptions{}
        ) {
            ::Mesh& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const size_t nC = cells.size();

            // 读取 mask；若无则按 0 处理
            auto maskPtr = reg.get<::volScalarField>(inj_mask_name);
            if (!maskPtr) {
                std::cerr << "[PeacemanPI] mask '" << inj_mask_name << "' not found; PI=0.\n";
                auto PIPtr = reg.getOrCreate<::volScalarField>(out_PI_name, nC, 0.0);
                PIPtr->data.assign(nC, 0.0);
                return true;
            }
            const auto& mask = maskPtr->data;

            // 输出 PI
            auto PIPtr = reg.getOrCreate<::volScalarField>(out_PI_name, nC, 0.0);
            auto& PI = PIPtr->data;
            PI.assign(nC, 0.0);

            const double r_w = std::max(well.r, 1e-12);
            const double two_pi_h = 2.0 * pi_value() * opt.thickness;
            constexpr double eps_mu = 1e-30;
            constexpr double eps_denom = 1e-12;

            size_t n_set = 0;
            for (size_t i = 0; i < nC; ++i) {
                if (mask[i] <= 0.0) { PI[i] = 0.0; continue; }

                const auto& c = cells[i];
                const double A = std::max(c.volume, 0.0);   // 2D: 面积；3D: 体积
                if (!(A > 0.0)) { PI[i] = 0.0; continue; }

                const double re = std::max(opt.beta_re * std::sqrt(A), r_w * (1.0 + 1e-9));

                // 局部流体性质（可选场模式）
                double rho = rho_const, mu = mu_const;
                if (fp_mode == FluidPropertyMode::Field) {
                    if (rho_field) rho = rho_field->data[i];
                    if (mu_field) mu = mu_field->data[i];
                }
                const double mu_safe = std::max(mu, eps_mu);
                const double logTerm = std::log(std::max(re / r_w, opt.min_log_arg));
                const double denom_raw = logTerm + opt.skin;
                const double denom_safe = (denom_raw >= eps_denom ? denom_raw : eps_denom);

                double pi_i = (rho * k_const * two_pi_h) / (mu_safe * denom_safe);

                // 面积/体积占比缩放（边界切分的粗略考虑：mask 由 Step1 决定）
                pi_i *= mask[i];

                if (opt.clamp_nonpos && !(pi_i > 0.0)) pi_i = 0.0;
                PI[i] = pi_i;
                if (pi_i > 0.0) ++n_set;
            }

            if (opt.verbose) {
                std::cout << "[PeacemanPI] '" << out_PI_name << "' built: " << n_set
                    << " active / " << nC << " cells. (h=" << opt.thickness
                    << ", skin=" << opt.skin << ", beta_re=" << opt.beta_re << ")\n";
            }
            return true;
        }

    } // namespace SourceTerm
} // namespace FVM

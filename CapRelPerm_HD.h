/**
 * @file CapRelPerm.h
 * @brief 毛管力(Capillary Pressure)与相对渗透率(Relative Permeability)模型定义
 * @details 实现了基于 Van Genuchten (VG) 模型的本构关系。
 * 为了保证数值稳定性 (Numerical Stability)，在残余饱和度附近引入了平滑扩展 (Linearization)，
 * 防止导数爆炸 (Derivative Explosion)，符合工业级模拟器标准 (e.g., Eclipse, CMG)。
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>

namespace CapRelPerm{

    // =========================================================
    // 全局常量定义
    // =========================================================

    /// @brief 极小值，用于防止除零错误
    constexpr double kTiny = 1e-12;

    /// @brief 最大毛管力限制 (Pa)，防止物理量溢出
    constexpr double kPcMax = 1e7;

    /// @brief VG 模型平滑处理的有效饱和度阈值
    /// @details 当 Se < kSe_Eps 时，采用线性延伸策略代替指数公式
    constexpr double kSe_Eps = 1.0e-3;

    // =========================================================
    // 参数结构体定义
    // =========================================================

    /**
     * @struct VGParams
     * @brief Van Genuchten 模型参数
     * @details 描述毛管力与饱和度的关系: Se = (1 + (alpha * Pc)^n)^(-m)
     */
    struct VGParams
    {
        double alpha = 1.0 / 5e4;   ///< 进气值倒数 [1/Pa], 约为 1/Pc_entry
        double n = 2.0;             ///< 孔径分布指数, 必须 > 1.0
        double Swr = 0.25;          ///< 束缚水饱和度 (Irreducible Water Saturation)
        double Sgr = 0.05;          ///< 残余气饱和度 (Residual Gas Saturation)

        /// @brief 计算 Mualem 参数 m = 1 - 1/n
        double m() const { return 1.0 - 1.0 / n; }
    };

    /**
     * @struct RelPermParams
     * @brief 相对渗透率参数 (Mualem 模型)
     */
    struct RelPermParams
    {
        double L = 0.5;             ///< 孔隙连通性因子 (通常为 0.5)
    };

    // =========================================================
    // 辅助函数
    // =========================================================

    /**
     * @brief 检查 VG 参数是否物理有效
     * @param vg VG参数结构体
     * @return true 有效, false 无效
     */
    inline bool vg_params_valid(const VGParams& vg)
    {
        if (vg.alpha <= 0.0) return false;
        if (vg.n <= 1.0) return false;
        if (vg.Swr < 0.0 || vg.Sgr < 0.0) return false;
        if (vg.Swr + vg.Sgr >= 1.0 - kTiny) return false;
        return true;
    }

    /**
     * @brief 计算有效饱和度 (Effective Saturation)
     * @details Se = (Sw - Swr) / (1 - Swr - Sgr)
     * @param Sw 当前水相饱和度
     * @param vg VG参数
     * @return double 有效饱和度 Se (未截断)
     */
    inline double calculate_Se(double Sw, const VGParams& vg)
    {
        double denom = 1.0 - vg.Swr - vg.Sgr;
        if (std::abs(denom) < kTiny) return (Sw >= vg.Swr) ? 1.0 : 0.0; // 防御性保护
        return (Sw - vg.Swr) / denom;
    }

    // =========================================================
    // 核心物理模型实现 (含平滑处理)
    // =========================================================

    /**
     * @brief 计算毛管力 Pc(Sw) - 工业级鲁棒版本
     * @details
     * 1. 当 Se >= kSe_Eps 时，使用标准 VG 公式: Pc = (1/alpha) * (Se^(-1/m) - 1)^(1/n)
     * 2. 当 Se < kSe_Eps 时，使用切线线性延伸 (Linear Extension)，防止 Pc -> inf。
     * 3. 结果限制在 [0, kPcMax] 范围内。
     * * @param Sw 水相饱和度
     * @param vg VG参数
     * @return double 毛管力 (Pa)
     */
    inline double pc_vG(double Sw, const VGParams& vg)
    {
        double Se = calculate_Se(Sw, vg);

        // 1. 物理上限截断 (完全润湿状态)
        if (Se >= 1.0) return 0.0;

        // 2. 正常 VG 计算区间 (Se >= Threshold)
        if (Se >= kSe_Eps) {
            double m = vg.m();
            double val = std::pow(Se, -1.0 / m) - 1.0;
            // 防止 val 因数值误差微小于0
            val = std::max(val, 0.0);
            double pc = (1.0 / vg.alpha) * std::pow(val, 1.0 / vg.n);
            return std::min(pc, kPcMax);
        }

        // 3. 线性化扩展区间 (Smoothing Region)
        // Pc(Se) ~ Pc(Eps) + Slope * (Se - Eps)
        // 目的：将无穷大的奇点替换为有限斜率的直线
        else {
            double Se_star = kSe_Eps;
            double m = vg.m();

            // 计算参考点 Pc(Eps)
            double val_star = std::pow(Se_star, -1.0 / m) - 1.0;
            double pc_star = (1.0 / vg.alpha) * std::pow(val_star, 1.0 / vg.n);

            // 计算参考点斜率 dPc/dSe | Eps
            // dPc/dSe = (1/alpha) * (1/n) * (val_star)^(1/n - 1) * (-1/m) * Se^(-1/m - 1)
            double term1 = (1.0 / (vg.alpha * vg.n));
            double term2 = std::pow(val_star, (1.0 / vg.n) - 1.0);
            double term3 = (-1.0 / m) * std::pow(Se_star, (-1.0 / m) - 1.0);
            double slope = term1 * term2 * term3; // 斜率应为负大值

            // 线性延伸
            double pc_linear = pc_star + slope * (Se - Se_star);

            // 安全截断
            if (pc_linear < 0.0) return 0.0;
            if (pc_linear > kPcMax) return kPcMax;

            return pc_linear;
        }
    }

    /**
     * @brief 计算毛管力对饱和度的导数 dPc/dSw - 工业级鲁棒版本
     * @details 对应 pc_vG 的逻辑，确保导数连续且有限。
     * * @param Sw 水相饱和度
     * @param vg VG参数
     * @return double dPc/dSw (Pa)
     */
    inline double d_pc_vG(double Sw, const VGParams& vg)
    {
        double Se = calculate_Se(Sw, vg);
        double dSe_dSw = 1.0 / (1.0 - vg.Swr - vg.Sgr);

        // 1. 饱和区导数为 0
        if (Se >= 1.0) return 0.0;

        // 2. 正常区间导数
        if (Se >= kSe_Eps) {
            double m = vg.m();
            double se_pow = std::pow(Se, -1.0 / m); // Se^(-1/m)
            double val = se_pow - 1.0;

            if (val <= kTiny) return 0.0; // 避免除零

            // dPc/dSe 公式推导:
            // Pc = (1/a) * val^(1/n)
            // dPc/dSe = (1/a) * (1/n) * val^(1/n - 1) * d(val)/dSe
            // d(val)/dSe = (-1/m) * Se^(-1/m - 1)

            double termA = (1.0 / (vg.alpha * vg.n)) * std::pow(val, (1.0 / vg.n) - 1.0);
            double termB = (-1.0 / m) * std::pow(Se, (-1.0 / m) - 1.0);

            double dPc_dSe = termA * termB;
            return dPc_dSe * dSe_dSw;
        }

        // 3. 线性化区间导数 (常数斜率)
        else {
            // 计算 kSe_Eps 处的斜率 (同 pc_vG 中的逻辑)
            double Se_star = kSe_Eps;
            double m = vg.m();
            double val_star = std::pow(Se_star, -1.0 / m) - 1.0;

            // 检查是否已经触顶 kPcMax
            // 如果 pc_vG 返回了 kPcMax，且此时 Se 进一步减小导致 Pc 更大，则物理上导数应截断?
            // 工业处理：如果线性延伸值 > kPcMax，通常认为进入了刚性边界，导数可设为0或保持斜率。
            // 为了数值收敛，保持斜率通常更好，能指示牛顿法回到有效区间。
            // 但如果 pc_vG 已经被 clamp 住了，导数如果非零，可能导致不一致。
            // 此处策略：计算线性斜率，如果 pc_linear > kPcMax，则视为截断区，导数为0。

            double term1 = (1.0 / (vg.alpha * vg.n));
            double term2 = std::pow(val_star, (1.0 / vg.n) - 1.0);
            double term3 = (-1.0 / m) * std::pow(Se_star, (-1.0 / m) - 1.0);
            double slope_Se = term1 * term2 * term3;

            // 检查值是否超限
            double pc_star = (1.0 / vg.alpha) * std::pow(val_star, 1.0 / vg.n);
            double pc_linear = pc_star + slope_Se * (Se - Se_star);

            if (pc_linear > kPcMax || pc_linear < 0.0) {
                return 0.0;
            }

            return slope_Se * dSe_dSw;
        }
    }

    /**
     * @brief 计算相对渗透率 (Mualem-vG Model)
     * @details
     * krw = sqrt(Se) * (1 - (1 - Se^(1/m))^m)^2
     * krg = sqrt(1 - Se) * (1 - Se^(1/m))^(2m)
     * 增加了对 Se < 0 和 Se > 1 的安全处理。
     * * @param Sw 水相饱和度
     * @param vg VG参数
     * @param rp 相渗参数(L)
     * @param[out] krw 水相相对渗透率
     * @param[out] krg 气相相对渗透率
     */
    inline void kr_Mualem_vG(double Sw, const VGParams& vg, const RelPermParams& rp, double& krw, double& krg)
    {
        double Se = calculate_Se(Sw, vg);

        // 边界截断
        if (Se >= 1.0) {
            krw = 1.0;
            krg = 0.0;
            return;
        }
        if (Se <= 0.0) {
            krw = 0.0;
            krg = 1.0;
            return;
        }

        double m = vg.m();

        // Krw 计算
        // Term1 = Se^0.5 (假设 L=0.5)
        // Term2 = (1 - (1 - Se^(1/m))^m)^2
        double Se_inv_m = std::pow(Se, 1.0 / m); // Se^(1/m)
        double one_minus_Se_inv_m = 1.0 - Se_inv_m;
        // 保护幂运算底数非负
        if (one_minus_Se_inv_m < 0.0) one_minus_Se_inv_m = 0.0;

        double term_m = std::pow(one_minus_Se_inv_m, m);
        double bracket_w = 1.0 - term_m;
        krw = std::pow(Se, rp.L) * (bracket_w * bracket_w);

        // Krg 计算
        // Krg = (1 - Se)^L * (1 - Se^(1/m))^(2m)
        // 注意：这里的公式依据常见的 Mualem CO2-Brine 形式，可能有变体。
        // 标准 Mualem 气相: Krg = sqrt(1-Se) * (1 - Se^(1/m))^(2m)
        // 此处假设 rp.L 也适用于气相的指数
        double bracket_g = 1.0 - Se_inv_m;
        if (bracket_g < 0.0) bracket_g = 0.0;

        krg = std::pow(1.0 - Se, rp.L) * std::pow(bracket_g, 2.0 * m);
    }

    /**
     * @brief 计算相对渗透率导数 (dKr/dSw)
     * @details 对 Mualem-vG 模型进行解析求导。
     * * @param Sw 水相饱和度
     * @param vg VG参数
     * @param rp 相渗参数
     * @param[out] dkrw_dSw 水相相渗导数
     * @param[out] dkrg_dSw 气相相渗导数
     */
    inline void d_kr_Mualem_vG(double Sw, const VGParams& vg, const RelPermParams& rp, double& dkrw_dSw, double& dkrg_dSw)
    {
        double Se = calculate_Se(Sw, vg);
        double dSe_dSw = 1.0 / (1.0 - vg.Swr - vg.Sgr);

        // 边界区导数为0 (或可考虑单侧导数，此处简化为0)
        if (Se >= 1.0 - kTiny || Se <= kTiny) {
            dkrw_dSw = 0.0;
            dkrg_dSw = 0.0;
            return;
        }

        double m = vg.m();
        double inv_m = 1.0 / m;

        // --- dKrw/dSe 推导 ---
        // krw = Se^L * B^2, where B = 1 - (1 - Se^inv_m)^m
        // Let X = 1 - Se^inv_m
        // B = 1 - X^m
        // dKrw/dSe = L * Se^(L-1) * B^2 + Se^L * 2 * B * dB/dSe

        double Se_pow_L = std::pow(Se, rp.L);
        double Se_inv_m = std::pow(Se, inv_m);
        double X = 1.0 - Se_inv_m;
        if (X < 0) X = 0;

        double X_pow_m = std::pow(X, m);
        double B = 1.0 - X_pow_m;

        // dB/dSe = -m * X^(m-1) * dX/dSe
        // dX/dSe = -inv_m * Se^(inv_m - 1)
        double dX_dSe = -inv_m * std::pow(Se, inv_m - 1.0);
        double dB_dSe = -m * std::pow(X, m - 1.0) * dX_dSe;

        double dKrw_dSe = rp.L * std::pow(Se, rp.L - 1.0) * B * B +
            Se_pow_L * 2.0 * B * dB_dSe;

        dkrw_dSw = dKrw_dSe * dSe_dSw;

        // --- dKrg/dSe 推导 ---
        // krg = (1-Se)^L * (1 - Se^inv_m)^(2m)
        // Let Y = 1 - Se^inv_m = X
        // krg = (1-Se)^L * Y^(2m)
        // dKrg/dSe = -L*(1-Se)^(L-1) * Y^(2m) + (1-Se)^L * 2m * Y^(2m-1) * dY/dSe

        double term1 = std::pow(1.0 - Se, rp.L);
        double term2 = std::pow(X, 2.0 * m);

        double dTerm1_dSe = -rp.L * std::pow(1.0 - Se, rp.L - 1.0);
        // dY/dSe is same as dX/dSe calculated above
        double dTerm2_dSe = 2.0 * m * std::pow(X, 2.0 * m - 1.0) * dX_dSe;

        double dKrg_dSe = dTerm1_dSe * term2 + term1 * dTerm2_dSe;

        dkrg_dSw = dKrg_dSe * dSe_dSw;
    }

} // namespace CapRelPerm
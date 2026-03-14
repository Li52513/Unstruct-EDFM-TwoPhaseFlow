#pragma once

#include "MeshManager.h"
#include "3D_MeshManager.h"
#include "2D_FieldManager.h"
#include "3D_FieldManager.h"
#include "FIM_StateMap.h"
#include "FIM_BlockSparseMatrix.h"
#include "FIM_GlobalAssembler.h"
#include "FIM_ConnectionManager.h"
#include "FIM_TopologyBuilder2D.h"
#include "FIM_TopologyBuilder3D.h"
#include "TransmissibilitySolver_2D.h"
#include "TransmissibilitySolver_3D.h"
#include "BoundaryAssembler.h"
#include "AD_FluidEvaluator.h"
#include "CapRelPerm_HD_AD.h"
#include "FVM_Ops_AD.h"
#include "2D_PostProcess.h"
#include "3D_PostProcess.h"
#include "SolverContrlStrName_op.h"
#include "FIM_TransientSupport.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iomanip>
#include <type_traits>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#define MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define MKDIR(path) mkdir(path, 0777)
#endif

    namespace FIM_Engine {

    /**
     * @brief V3 诊断级别枚举
     */
    enum class DiagLevel { Off, Summary, Hotspot, Forensic };

    /**
     * @struct EqContrib
     * @brief 方程贡献拆解容器，独立记录各个物理过程对特定方程残差和 Jacobian 对角线的贡献
     */
    struct EqContrib {
        double R_acc = 0.0, R_flux = 0.0, R_well = 0.0, R_bc = 0.0;
        double D_acc = 0.0, D_flux = 0.0, D_well = 0.0, D_bc = 0.0;

        /** @brief 重置当前方程的贡献记录 */
        void reset() {
            R_acc = R_flux = R_well = R_bc = 0.0;
            D_acc = D_flux = D_well = D_bc = 0.0;
        }

        /** @brief 获取当前方程的总残差值 */
        double R_total() const { return R_acc + R_flux + R_well + R_bc; }

        /** @brief 获取当前方程的 Jacobian 主对角线总值 */
        double D_total() const { return D_acc + D_flux + D_well + D_bc; }
    };

    // =====================================================================
    // 求解器与初始条件模块 (核心数据结构)
    // =====================================================================

    enum class SolverRoute { FIM, IMPES };
    enum class LinearSolverType { SparseLU, BiCGSTAB };

    /** @brief 瞬态求解初始条件 */
    struct InitialConditions {
        double P_init = 2.0e5;   ///< 初始地层压力 (Pa)
        double T_init = 300.0;   ///< 初始地层温度 (K)
        double Sw_init = 0.2;    ///< 初始含水饱和度 (-)
    };

    /** @brief 瞬态全隐式求解器控制参数 */
    struct TransientSolverParams {
        int max_steps = 50;                     ///< 模拟最大推进步数

        // 时间步长控制
        double dt_init = 1.0;                   ///< 初始时间步长（建议0.05~1.0s，应对瞬态冲击效应）
        double dt_min = 1e-4;                   ///< 最小允许时间步长 (s)
        double dt_max = 86400.0;                ///< 最大允许时间步长 (s)

        // 非线性控制
        int max_newton_iter = 8;                ///< Day6 计划要求每步最大 8~15 次非线性迭代
        double abs_res_tol = 1e-6;              ///< 绝对残差收敛容差（建议不用做最终判定）

        // 停滞(Stagnation)与回退控制
        double stagnation_growth_tol = 0.995;   ///< 停滞判定容忍度
        double stagnation_abs_res_tol = 1.0e6;  ///< 停滞时允许的最大绝对残差 (单相建议 1e4, 两相建议 1e6)
        double stagnation_min_drop = 2.0e-3;    ///< 停滞时要求的最小相对残差下降量
        double best_iter_growth_trigger = 1.5;  ///< Trigger best-iterate takeover when current residual grows beyond this factor vs in-step best.
        int best_iter_guard_min_iter = 3;       ///< Minimum Newton iterations before enabling best-iterate takeover.
        bool enable_best_iter_guard = false;    ///< Disable best-iterate takeover by default for strict acceptance.
        bool enable_dt_floor_best = false;      ///< Disable dt-floor best-state acceptance by default.
        bool enable_stagnation_accept = false;  ///< Disable stagnation-based acceptance by default.
        bool enable_dt_floor_hold = false;      ///< Emergency-only fallback at dt_min; keep false for strict acceptance.

        // 双收敛准则: 相对残差 + 相对更新量
        double rel_res_tol = 1.0e-3;            ///< Relative residual threshold (max_res / res_iter1).
        double rel_update_tol = 1.0e-6;         ///< Relative state-update infinity norm threshold.

        // Armijo 回溯线搜索
        bool enable_armijo_line_search = true;  ///< Enable residual-tested backtracking line search.
        int armijo_max_backtracks = 8;          ///< Max backtracking count per Newton step.
        double armijo_beta = 0.5;               ///< Backtracking shrink factor in (0,1).
        double armijo_c1 = 1.0e-4;              ///< Armijo sufficient decrease coefficient.

        // 行缩放（左缩放）: 以 |diag_acc| 与 |diag(A)| 的特征量为基准
        bool enable_row_scaling = true;         ///< Enable row scaling before linear solve.
        double row_scale_floor = 1.0;           ///< Characteristic lower bound in scaling denominator.
        double row_scale_min = 1.0e-12;         ///< Min allowed row scale factor.
        double row_scale_max = 1.0e+12;         ///< Max allowed row scale factor.

        // 状态截断与保护(Limiter & Damping)
        double max_dP = 1.0e4;                  ///< [防发散] 牛顿步最大允许压力变化量(Pa)
        double max_dT = 2.0;                    ///< [防发散] 牛顿步最大允许温度变化量(K)
        double max_dSw = 0.05;                  ///< 牛顿步最大允许饱和度变化量(-)
        double min_alpha = 1.0e-8;              ///< line-search/damping lower bound; avoid forced large update at stiff states
        bool clamp_state_to_eos_bounds = false; ///< hard-clip P/T into EOS table bounds; disabled by default to avoid edge pinning

        // 线性求解器配置
        LinearSolverType lin_solver = LinearSolverType::SparseLU; ///< [默认SparseLU] 应对高度非对称病态系统的求解器
        double bicgstab_droptol = 1e-2;                           ///< 使用 BiCGSTAB 时使用的 ILUT 截断容差
        double well_source_sign = 1.0;                            ///< Well sign for outflow-positive well operators in residual: R += Q_out.

        // ==========================================================
        // V3 发散定位诊断系统专属配置 (默认开启 Summary 模式)
        // ==========================================================
        DiagLevel diag_level = DiagLevel::Summary;   ///< 诊断级别
        int diag_print_every_iter = 1;               ///< Summary 常驻输出的迭代步频次
        double diag_blowup_factor = 5.0;             ///< 残差爆炸触发阈值 (res_new > res_old * factor)
        int diag_hot_repeat_iters = 3;               ///< 连续同一 Hotspot 卡死触发深挖的次数
        double diag_hot_res_change_tol = 1e-2;       ///< Hotspot 残差变化停滞判定容忍度
        int diag_clamp_trigger = 20;                 ///< Limiter 触发风暴警告的阈值次数
        int diag_max_hot_conn = 5;                   ///< 触发深挖时输出的 Top-K 拓扑连接数
        int diag_max_clamp_dump = 10;                ///< 事故快照中最多倾倒的 Limiter 事件明细数
        double diag_flux_spike_factor = 10.0;        ///< 全局通量尖峰倍数报警阈值
        double diag_eos_near_bound_ratio = 0.02;     ///< 系统中越界网格达到此比例触发 EOS 全局警告
        bool diag_incident_once_per_step = true;     ///< 防止刷屏，同一时间步最多生成一次事故快照 (JSON)
    };

    // =====================================================================
    // SFINAE 模板推导：用于在编译期安全判断求解器是否为迭代法 (是否有 iterations() 接口)
    // =====================================================================
    template <typename T, typename = void>
    struct is_iterative_solver : std::false_type {};

    template <typename T>
    struct is_iterative_solver<T, std::void_t<decltype(std::declval<T>().iterations())>> : std::true_type {};

    // =====================================================================
    // V3 补丁：3D 拓扑与交互对启动预检查 (Pre-check)
    // =====================================================================
    template<typename MeshMgrType>
    inline void Run3DDiagnosticPrecheck(MeshMgrType& mgr, const std::vector<Connection>& conns, const TransientSolverParams& params) {
        if constexpr (std::is_same_v<MeshMgrType, MeshManager_3D>) {
            if (params.diag_level == DiagLevel::Off) return;

            std::cout << "\n=========================================================\n"
                << "[PRE3D-DIAG] V3 Diagnostic Pre-check Starting...\n";

            // 1. 校验 InteractionPairs (基岩-裂缝交互对)
            int invalid_idx = 0, invalid_area = 0, invalid_dist = 0;
            const auto& pairs = mgr.getInteractionPairs();
            for (const auto& p : pairs) {
                if (p.matrixSolverIndex < 0 || p.fracCellSolverIndex < 0) invalid_idx++;
                if (p.intersectionArea <= 1e-12) invalid_area++;
                if (p.distMatrixToFracPlane <= 1e-12) invalid_dist++;
            }
            std::cout << "  [PRE3D-PAIR] Total Pairs: " << pairs.size() << "\n"
                << "               Invalid Index: " << invalid_idx << "\n"
                << "               Area <= eps  : " << invalid_area << "\n"
                << "               Dist <= eps  : " << invalid_dist << "\n";

            // 2. 校验全局拓扑连接 (Connections)
            int mm = 0, mf = 0, ff = 0, fi = 0;
            int neg_t_flow = 0, zero_t_flow = 0;
            int neg_t_heat = 0, zero_t_heat = 0;

            for (const auto& c : conns) {
                if (c.type == ConnectionType::Matrix_Matrix) mm++;
                else if (c.type == ConnectionType::Matrix_Fracture) mf++;
                else if (c.type == ConnectionType::Fracture_Fracture) ff++;
                else if (c.type == ConnectionType::Fracture_Internal) fi++;

                if (c.T_Flow < 0.0) neg_t_flow++;
                if (std::abs(c.T_Flow) < 1e-16) zero_t_flow++;

                if (c.T_Heat < 0.0) neg_t_heat++;
                if (std::abs(c.T_Heat) < 1e-16) zero_t_heat++;
            }
            std::cout << "  [PRE3D-CONN] Connections Count -> MM: " << mm << ", MF(NNC): " << mf
                << ", FF: " << ff << ", FI: " << fi << "\n"
                << "               T_Flow Anomaly  -> Negative: " << neg_t_flow << ", Zero: " << zero_t_flow << "\n"
                << "               T_Heat Anomaly  -> Negative: " << neg_t_heat << ", Zero: " << zero_t_heat << "\n"
                << "=========================================================\n\n";
        }
    }


    /**
     * @brief Day6 可选外部模块注入: 自定义初始化器与边界条件设置
     * @details
     * - property_initializer: 在矩阵系统建立之前执行，用于函数式注入基岩与裂缝物性
     * - pressure_bc / saturation_bc / temperature_bc: 对应场变量的边界条件设置
     */
    inline const char* ConnectionTypeLabel(ConnectionType type) {
        switch (type) {
        case ConnectionType::Matrix_Matrix: return "MM";
        case ConnectionType::Matrix_Fracture: return "MF";
        case ConnectionType::Fracture_Fracture: return "FF";
        case ConnectionType::Fracture_Internal: return "FI";
        default: return "Unknown";
        }
    }

    template <typename MeshMgrType, typename FieldMgrType>
    struct TransientOptionalModules {
        std::function<void(MeshMgrType&, FieldMgrType&)> property_initializer;
        const BoundarySetting::BoundaryConditionManager* pressure_bc = nullptr;
        const BoundarySetting::BoundaryConditionManager* saturation_bc = nullptr;
        const BoundarySetting::BoundaryConditionManager* temperature_bc = nullptr;
    };

    // =====================================================================
    // Traits 辅助方法提取器
    // =====================================================================

    inline int MatrixBlockCount(const MeshManager& mgr) { return mgr.getMatrixDOFCount(); }
    inline int MatrixBlockCount(const MeshManager_3D& mgr) { return mgr.fracture_network().getSolverIndexOffset(); }

    inline void MakePath(const std::string& caseName) {
        MKDIR("Test"); MKDIR("Test/Transient"); MKDIR("Test/Transient/Day6"); MKDIR(("Test/Transient/Day6/" + caseName).c_str());
    }

    template <typename FieldMgrType, typename MeshMgrType, int N>
    inline void SyncStateToFieldManager(const FIM_StateMap<N>& state, FieldMgrType& fm, const MeshMgrType& mgr) {
        int nMat = MatrixBlockCount(mgr);
        int nTotal = mgr.getTotalDOFCount();

        const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        const auto satCfg = PhysicalProperties_string_op::SaturationEquation_String::FIM();
        const PhysicalProperties_string_op::Water wCfg;
        const PhysicalProperties_string_op::CO2 gCfg;

        auto f_pw = fm.getOrCreateMatrixScalar(pCfg.pressure_field, 0.0);
        auto f_T = fm.getOrCreateMatrixScalar(tCfg.temperatue_field, 0.0);
        auto f_rhow = fm.getOrCreateMatrixScalar(wCfg.rho_tag, 0.0);
        auto f_hw = fm.getOrCreateMatrixScalar(wCfg.h_tag, 0.0);
        auto f_lamw_mob = fm.getOrCreateMatrixScalar(wCfg.lambda_w_tag, 0.0);

        auto frac_pw = fm.getOrCreateFractureScalar(pCfg.pressure_field, 0.0);
        auto frac_T = fm.getOrCreateFractureScalar(tCfg.temperatue_field, 0.0);
        auto frac_rhow = fm.getOrCreateFractureScalar(wCfg.rho_tag, 0.0);
        auto frac_hw = fm.getOrCreateFractureScalar(wCfg.h_tag, 0.0);
        auto frac_lamw_mob = fm.getOrCreateFractureScalar(wCfg.lambda_w_tag, 0.0);

        auto f_P_viz = fm.getOrCreateMatrixScalar("P", 0.0);
        auto frac_P_viz = fm.getOrCreateFractureScalar("P", 0.0);

        std::shared_ptr<volScalarField> f_sw, f_Sw_viz, f_rhog, f_hg, f_lamg_mob;
        std::shared_ptr<volScalarField> frac_sw, frac_Sw_viz, frac_rhog, frac_hg, frac_lamg_mob;

        if constexpr (N == 3) {
            f_sw = fm.getOrCreateMatrixScalar(satCfg.saturation, 0.0);
            f_Sw_viz = fm.getOrCreateMatrixScalar("S_w", 0.0);
            f_rhog = fm.getOrCreateMatrixScalar(gCfg.rho_tag, 0.0);
            f_hg = fm.getOrCreateMatrixScalar(gCfg.h_tag, 0.0);
            f_lamg_mob = fm.getOrCreateMatrixScalar(gCfg.lambda_g_tag, 0.0);

            frac_sw = fm.getOrCreateFractureScalar(satCfg.saturation, 0.0);
            frac_Sw_viz = fm.getOrCreateFractureScalar("S_w", 0.0);
            frac_rhog = fm.getOrCreateFractureScalar(gCfg.rho_tag, 0.0);
            frac_hg = fm.getOrCreateFractureScalar(gCfg.h_tag, 0.0);
            frac_lamg_mob = fm.getOrCreateFractureScalar(gCfg.lambda_g_tag, 0.0);
        }

        for (int i = 0; i < nTotal; ++i) {
            double p = state.P[i];
            double t = state.T[i];
            ADVar<N> P_ad(p), T_ad(t);
            auto propsW = AD_Fluid::Evaluator::evaluateWater<N>(P_ad, T_ad);

            double sw = (N == 3) ? state.Sw[i] : 1.0;
            double rho_w = propsW.rho.val;
            double mu_w = propsW.mu.val;
            double krw = 1.0, krg = 0.0;

            if constexpr (N == 3) {
                CapRelPerm::VGParams vg; vg.alpha = 1e-5; vg.n = 2.0; vg.Swr = 0.0; vg.Sgr = 0.0;
                CapRelPerm::RelPermParams rp; rp.L = 0.5;
                ADVar<N> krw_ad, krg_ad;
                CapRelPerm::kr_Mualem_vG<N>(ADVar<N>(sw), vg, rp, krw_ad, krg_ad);
                krw = krw_ad.val; krg = krg_ad.val;
            }

            double lambda_w_mob = krw / std::max(mu_w, 1e-18);

            if (i < nMat) {
                (*f_pw)[i] = p; (*f_T)[i] = t; (*f_rhow)[i] = rho_w; (*f_hw)[i] = propsW.h.val; (*f_lamw_mob)[i] = lambda_w_mob; (*f_P_viz)[i] = p;
                if constexpr (N == 3) {
                    auto propsG = AD_Fluid::Evaluator::evaluateCO2<N>(P_ad, T_ad);
                    (*f_sw)[i] = sw; (*f_Sw_viz)[i] = sw; (*f_rhog)[i] = propsG.rho.val; (*f_hg)[i] = propsG.h.val; (*f_lamg_mob)[i] = krg / std::max(propsG.mu.val, 1e-18);
                }
            }
            else {
                int fi = i - nMat;
                (*frac_pw)[fi] = p; (*frac_T)[fi] = t; (*frac_rhow)[fi] = rho_w; (*frac_hw)[fi] = propsW.h.val; (*frac_lamw_mob)[fi] = lambda_w_mob; (*frac_P_viz)[fi] = p;
                if constexpr (N == 3) {
                    auto propsG = AD_Fluid::Evaluator::evaluateCO2<N>(P_ad, T_ad);
                    (*frac_sw)[fi] = sw; (*frac_Sw_viz)[fi] = sw; (*frac_rhog)[fi] = propsG.rho.val; (*frac_hg)[fi] = propsG.h.val; (*frac_lamg_mob)[fi] = krg / std::max(propsG.mu.val, 1e-18);
                }
            }
        }
    }

    template <typename FieldMgrType>
    inline void InjectStaticProperties(FieldMgrType& fm) {
        const PhysicalProperties_string_op::Rock rock;
        const PhysicalProperties_string_op::Fracture_string frac;
        const PhysicalProperties_string_op::Water water;

        fm.getOrCreateMatrixScalar(rock.k_xx_tag, 1.0e-13);
        fm.getOrCreateMatrixScalar(rock.k_yy_tag, 1.0e-13);
        fm.getOrCreateMatrixScalar(rock.k_zz_tag, 1.0e-13);
        fm.getOrCreateMatrixScalar(rock.lambda_tag, 2.0);
        fm.getOrCreateMatrixScalar(rock.phi_tag, 0.2);

        fm.getOrCreateFractureScalar(frac.k_t_tag, 1.0e-11);
        fm.getOrCreateFractureScalar(frac.k_n_tag, 1.0e-12);
        fm.getOrCreateFractureScalar(frac.aperture_tag, 1.0e-3);
        fm.getOrCreateFractureScalar(frac.lambda_tag, 2.0);
        fm.getOrCreateFractureScalar(frac.phi_tag, 0.2);

        fm.getOrCreateFractureScalar(water.k_tag, 0.6);
    }

    /**
     * @brief Generic transient FIM driver used by transient scenarios. (Fully Parametrized)
     */
    template <int N, typename MeshMgrType, typename FieldMgrType>
    inline void RunGenericFIMTransient(
        const std::string& caseName,
        MeshMgrType& mgr,
        FieldMgrType& fm,
        const InitialConditions& ic,
        const std::vector<WellScheduleStep>& wells,
        const TransientSolverParams& params,
        SolverRoute route = SolverRoute::FIM,
        const TransientOptionalModules<MeshMgrType, FieldMgrType>& modules = TransientOptionalModules<MeshMgrType, FieldMgrType>()) {

        if (route == SolverRoute::IMPES) {
            throw std::runtime_error("[TODO] IMPES explicit route is reserved but currently bypassed.");
        }

        std::cout << "\n========== Starting Transient Scenario: " << caseName << " ==========\n";
        MakePath(caseName);
        InjectStaticProperties(fm);
        if (modules.property_initializer) {
            modules.property_initializer(mgr, fm);
            std::cout << "[Init] External property module injected.\n";
        }

        const int totalBlocks = mgr.getTotalDOFCount();
        const int totalEq = mgr.getTotalEquationDOFs();
        FIM_StateMap<N> state;
        state.InitSizes(totalBlocks);

        const auto pEqCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const auto tEqCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        const auto sEqCfg = PhysicalProperties_string_op::SaturationEquation_String::FIM();
        const int pressureDof = 0;
        const int saturationDof = (N == 3) ? 1 : -1;
        const int temperatureDof = (N == 3) ? 2 : 1;

        // 初始条件注入
        for (int i = 0; i < totalBlocks; ++i) {
            state.P[i] = ic.P_init;
            state.T[i] = ic.T_init;
            if constexpr (N == 3) state.Sw[i] = ic.Sw_init;
        }

        std::vector<double> vols(totalBlocks, 1.0);
        const int nMat = MatrixBlockCount(mgr);
        for (size_t i = 0; i < mgr.mesh().getCells().size(); ++i) vols[i] = mgr.mesh().getCells()[i].volume;
        for (size_t i = 0; i < mgr.fracture_network().getOrderedFractureElements().size(); ++i) {
            auto* elem = mgr.fracture_network().getOrderedFractureElements()[i];
            if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                vols[nMat + static_cast<int>(i)] = elem->length * std::max(elem->aperture, 1e-6);
            }
            else {
                const auto* frac2d = mgr.findFractureByID(mgr.fracture_network(), elem->parentFractureID);
                const double ap = frac2d ? frac2d->aperture : 1e-3;
                vols[nMat + static_cast<int>(i)] = std::max(elem->area, 1e-12) * std::max(ap, 1e-8);
            }
        }

        std::vector<Vector> blockCenters(totalBlocks, Vector(0.0, 0.0, 0.0));
        for (size_t i = 0; i < mgr.mesh().getCells().size(); ++i) blockCenters[i] = mgr.mesh().getCells()[i].center;
        if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
            const auto& orderedFrac = mgr.fracture_network().getOrderedFractureElements();
            const auto& macros = mgr.fracture_network().fractures;
            for (size_t i = 0; i < orderedFrac.size(); ++i) {
                const auto* elem = orderedFrac[i];
                Vector c(0.0, 0.0, 0.0);
                if (elem && elem->parentFractureID >= 0 && elem->parentFractureID < static_cast<int>(macros.size())) {
                    const auto& frac = macros[elem->parentFractureID];
                    double u = std::max(0.0, std::min(1.0, 0.5 * (elem->param0 + elem->param1)));
                    c = frac.start + (frac.end - frac.start) * u;
                }
                blockCenters[nMat + static_cast<int>(i)] = c;
            }
        }
        else {
            const auto& orderedFrac = mgr.fracture_network().getOrderedFractureElements();
            for (size_t i = 0; i < orderedFrac.size(); ++i) {
                blockCenters[nMat + static_cast<int>(i)] = orderedFrac[i] ? orderedFrac[i]->centroid : Vector(0.0, 0.0, 0.0);
            }
        }
        const Vector gravityVec(0.0, 0.0, -9.81);

        FIM_ConnectionManager connMgr;
        if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
            TransmissibilitySolver_2D::Calculate_Transmissibility_Matrix(mgr, fm);
            TransmissibilitySolver_2D::Calculate_Transmissibility_FractureInternal(mgr, fm);
            TransmissibilitySolver_2D::Calculate_Transmissibility_NNC(mgr, fm);
            TransmissibilitySolver_2D::Calculate_Transmissibility_FF(mgr, fm);
            FIM_TopologyBuilder2D::LoadAllConnections(connMgr, mgr, fm);
        }
        else {
            TransmissibilitySolver_3D::Calculate_Transmissibility_Matrix(mgr, fm);
            TransmissibilitySolver_3D::Calculate_Transmissibility_FractureInternal(mgr, fm);
            TransmissibilitySolver_3D::Calculate_Transmissibility_NNC(mgr, fm);
            TransmissibilitySolver_3D::Calculate_Transmissibility_FF(mgr, fm);
            FIM_TopologyBuilder3D::LoadAllConnections(connMgr, mgr, fm);
        }
        connMgr.FinalizeAndAggregate();

        FIM_BlockSparseMatrix<N> global_mat(totalBlocks);
        for (const auto& conn : connMgr.GetConnections()) {
            global_mat.AddOffDiagBlock(conn.nodeI, conn.nodeJ, Eigen::Matrix<double, N, N>::Zero());
            global_mat.AddOffDiagBlock(conn.nodeJ, conn.nodeI, Eigen::Matrix<double, N, N>::Zero());
        }
        global_mat.FreezePattern();

        double t = 0.0;
        double dt = params.dt_init;
        int total_rollbacks = 0;
        int total_limiters = 0;
        int completed_steps = 0;
        int vtk_export_count = 0;
        double final_residual = std::numeric_limits<double>::quiet_NaN();
        // State guard rails:
        // 1) always keep conservative physical lower bounds;
        // 2) optionally enforce EOS-table hard clipping (off by default).
        double p_floor = 1.0e4;
        double p_ceil = std::numeric_limits<double>::max();
        double t_floor = 273.15;
        double t_ceil = std::numeric_limits<double>::max();
        if (params.clamp_state_to_eos_bounds) {
            const auto& wt_table = WaterPropertyTable::instance();
            p_floor = std::max(p_floor, wt_table.minPressure() + 1.0);
            p_ceil = std::min(p_ceil, wt_table.maxPressure() - 1.0);
            t_floor = std::max(t_floor, wt_table.minTemperature() + 1.0e-6);
            t_ceil = std::min(t_ceil, wt_table.maxTemperature() - 1.0e-6);
            if constexpr (N == 3) {
                const auto& gt_table = CO2PropertyTable::instance();
                p_floor = std::max(p_floor, gt_table.minPressure() + 1.0);
                p_ceil = std::min(p_ceil, gt_table.maxPressure() - 1.0);
                t_floor = std::max(t_floor, gt_table.minTemperature() + 1.0e-6);
                t_ceil = std::min(t_ceil, gt_table.maxTemperature() - 1.0e-6);
            }
            if (!(p_floor < p_ceil)) { p_floor = 1.0e4; p_ceil = std::numeric_limits<double>::max(); }
            if (!(t_floor < t_ceil)) { t_floor = 273.15; t_ceil = std::numeric_limits<double>::max(); }
        }

        for (int step = 1; step <= params.max_steps; ++step) {
            bool converged = false;
            FIM_StateMap<N> old_state = state;
            FIM_StateMap<N> best_state = state;
            std::string fail_reason;
            std::string converge_mode = "none";
            int iter_used = 0;
            double res_iter1 = -1.0;
            double step_final_residual = std::numeric_limits<double>::quiet_NaN();
            double best_res = std::numeric_limits<double>::infinity();
            int best_iter = -1;

            // 【新增】：V3 诊断 - 事故快照防刷屏锁
            bool incident_dumped_this_step = false;
            int prev_hot_idx = -1;
            double prev_hot_res = -1.0;
            int hot_repeat_count = 0;
            double prev_max_mass_flux = 0.0;
            double prev_max_heat_flux = 0.0;

            // 【新增】：V3 诊断 - 在 Step 1 触发 3D 预检查
            if (step == 1) {
                Run3DDiagnosticPrecheck(mgr, connMgr.GetConnections(), params);
            }

            // Dual convergence uses previous accepted relative update.
            double last_rel_update = std::numeric_limits<double>::infinity();

            // Residual probe used by Armijo line search.
            auto compute_residual_inf_for_state = [&](const FIM_StateMap<N>& eval_state) -> double {
                FIM_BlockSparseMatrix<N> probe_mat(totalBlocks);
                probe_mat.SetZero();

                // 1) accumulation
                for (int bi = 0; bi < totalBlocks; ++bi) {
                    const double phi = 0.2, c_pr = 1000.0, rho_r = 2600.0;
                    ADVar<N> P(eval_state.P[bi]); P.grad(0) = 1.0;
                    ADVar<N> T(eval_state.T[bi]); T.grad((N == 2) ? 1 : 2) = 1.0;
                    auto pW = AD_Fluid::Evaluator::evaluateWater<N>(P, T);
                    ADVar<N> P_old(old_state.P[bi]), T_old(old_state.T[bi]);
                    auto pW_old = AD_Fluid::Evaluator::evaluateWater<N>(P_old, T_old);

                    std::vector<ADVar<N>> acc_eqs(N);
                    if constexpr (N == 2) {
                        ADVar<N> m_w = pW.rho * phi, m_w_old = pW_old.rho * phi;
                        acc_eqs[0] = (m_w - m_w_old) * (vols[bi] / dt);
                        ADVar<N> e_w = m_w * (pW.h - P / pW.rho) + ADVar<N>((1.0 - phi) * rho_r * c_pr) * T;
                        ADVar<N> e_w_old = m_w_old * (pW_old.h - P_old / pW_old.rho) + ADVar<N>((1.0 - phi) * rho_r * c_pr) * T_old;
                        acc_eqs[1] = (e_w - e_w_old) * (vols[bi] / dt);
                    }
                    else {
                        ADVar<N> Sw(eval_state.Sw[bi]); Sw.grad(1) = 1.0;
                        ADVar<N> Sg = ADVar<N>(1.0) - Sw;
                        ADVar<N> Sw_old(old_state.Sw[bi]), Sg_old = ADVar<N>(1.0) - Sw_old;
                        auto pG = AD_Fluid::Evaluator::evaluateCO2<N>(P, T);
                        auto pG_old = AD_Fluid::Evaluator::evaluateCO2<N>(P_old, T_old);
                        acc_eqs[0] = (pW.rho * phi * Sw - pW_old.rho * phi * Sw_old) * (vols[bi] / dt);
                        acc_eqs[1] = (pG.rho * phi * Sg - pG_old.rho * phi * Sg_old) * (vols[bi] / dt);
                        ADVar<N> e_fluid = pW.rho * Sw * (pW.h - P / pW.rho) + pG.rho * Sg * (pG.h - P / pG.rho);
                        ADVar<N> e_fluid_old = pW_old.rho * Sw_old * (pW_old.h - P_old / pW_old.rho) + pG_old.rho * Sg_old * (pG_old.h - P_old / pG_old.rho);
                        ADVar<N> e_rock = ADVar<N>((1.0 - phi) * rho_r * c_pr) * T;
                        ADVar<N> e_rock_old = ADVar<N>((1.0 - phi) * rho_r * c_pr) * T_old;
                        acc_eqs[2] = ((e_fluid * phi + e_rock) - (e_fluid_old * phi + e_rock_old)) * (vols[bi] / dt);
                    }
                    FIM_GlobalAssembler<N, ADVar<N>>::AssembleAccumulation(bi, acc_eqs, probe_mat);
                }

                // 2) flux
                for (const auto& conn : connMgr.GetConnections()) {
                    int i = conn.nodeI, j = conn.nodeJ;
                    auto evalFlux = [&](bool wrt_i) -> std::vector<ADVar<N>> {
                        std::vector<ADVar<N>> F(N);
                        ADVar<N> P_i(eval_state.P[i]), T_i(eval_state.T[i]), P_j(eval_state.P[j]), T_j(eval_state.T[j]);
                        const Vector& x_i = blockCenters[i]; const Vector& x_j = blockCenters[j];

                        if constexpr (N == 2) {
                            if (wrt_i) { P_i.grad(0) = 1.0; T_i.grad(1) = 1.0; }
                            else { P_j.grad(0) = 1.0; T_j.grad(1) = 1.0; }
                            auto pW_i = AD_Fluid::Evaluator::evaluateWater<N>(P_i, T_i);
                            auto pW_j = AD_Fluid::Evaluator::evaluateWater<N>(P_j, T_j);
                            ADVar<N> rho_avg_w = ADVar<N>(0.5) * (pW_i.rho + pW_j.rho);
                            ADVar<N> dPhi = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, rho_avg_w, x_i, x_j, gravityVec);
                            ADVar<N> mob_i = pW_i.rho / pW_i.mu, mob_j = pW_j.rho / pW_j.mu;
                            ADVar<N> up_mob = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi, mob_i, mob_j);
                            F[0] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mob, dPhi);
                            ADVar<N> up_h = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi, pW_i.h, pW_j.h);
                            F[1] = FVM_Ops::Compute_Heat_Flux<N, ADVar<N>>(conn.T_Heat, T_i, T_j, F[0], up_h);
                        }
                        else {
                            ADVar<N> Sw_i(eval_state.Sw[i]), Sw_j(eval_state.Sw[j]);
                            if (wrt_i) { P_i.grad(0) = 1.0; Sw_i.grad(1) = 1.0; T_i.grad(2) = 1.0; }
                            else { P_j.grad(0) = 1.0; Sw_j.grad(1) = 1.0; T_j.grad(2) = 1.0; }
                            auto pW_i = AD_Fluid::Evaluator::evaluateWater<N>(P_i, T_i), pW_j = AD_Fluid::Evaluator::evaluateWater<N>(P_j, T_j);
                            auto pG_i = AD_Fluid::Evaluator::evaluateCO2<N>(P_i, T_i), pG_j = AD_Fluid::Evaluator::evaluateCO2<N>(P_j, T_j);
                            CapRelPerm::VGParams vg; vg.alpha = 1e-5; vg.n = 2.0; vg.Swr = 0.0; vg.Sgr = 0.0;
                            CapRelPerm::RelPermParams rp; rp.L = 0.5;
                            ADVar<N> krw_i, krg_i, krw_j, krg_j;
                            CapRelPerm::kr_Mualem_vG<N>(Sw_i, vg, rp, krw_i, krg_i);
                            CapRelPerm::kr_Mualem_vG<N>(Sw_j, vg, rp, krw_j, krg_j);
                            ADVar<N> rho_avg_w = ADVar<N>(0.5) * (pW_i.rho + pW_j.rho), rho_avg_g = ADVar<N>(0.5) * (pG_i.rho + pG_j.rho);
                            ADVar<N> Pc_i = CapRelPerm::pc_vG<N>(Sw_i, vg), Pc_j = CapRelPerm::pc_vG<N>(Sw_j, vg);
                            ADVar<N> dPhi_w = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, rho_avg_w, x_i, x_j, gravityVec);
                            ADVar<N> dPhi_g = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, Pc_i, Pc_j, rho_avg_g, x_i, x_j, gravityVec);
                            ADVar<N> mobW_i = krw_i * pW_i.rho / pW_i.mu, mobW_j = krw_j * pW_j.rho / pW_j.mu;
                            ADVar<N> mobG_i = krg_i * pG_i.rho / pG_i.mu, mobG_j = krg_j * pG_j.rho / pG_j.mu;
                            ADVar<N> up_mobW = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_w, mobW_i, mobW_j);
                            ADVar<N> up_mobG = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_g, mobG_i, mobG_j);
                            F[0] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mobW, dPhi_w);
                            F[1] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mobG, dPhi_g);
                            ADVar<N> up_h_w = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_w, pW_i.h, pW_j.h);
                            ADVar<N> up_h_g = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_g, pG_i.h, pG_j.h);
                            F[2] = FVM_Ops::Compute_Heat_Flux<N, ADVar<N>>(conn.T_Heat, T_i, T_j, F[0], F[1], up_h_w, up_h_g);
                        }
                        return F;
                        };
                    auto f_wrt_i = evalFlux(true);
                    auto f_wrt_j = evalFlux(false);
                    FIM_GlobalAssembler<N, ADVar<N>>::AssembleFlux(i, j, f_wrt_i, f_wrt_j, probe_mat);
                }

                // 3) well source
                SyncStateToFieldManager(eval_state, fm, mgr);
                std::vector<double> w_res(totalEq, 0.0);
                std::vector<std::array<double, 3>> w_jac3(totalEq, std::array<double, 3>{ 0.0, 0.0, 0.0 });
                if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                    BoundaryAssembler::Assemble_Wells_2D_FullJac(mgr, fm, wells, 0, 0, (N == 3) ? 1 : -1, (N == 3) ? 2 : 1, w_res, w_jac3);
                }
                else {
                    BoundaryAssembler::Assemble_Wells_3D_FullJac(mgr, fm, wells, 0, 0, (N == 3) ? 1 : -1, (N == 3) ? 2 : 1, w_res, w_jac3);
                }
                const double kWellSourceSignProbe = params.well_source_sign;
                for (int bi = 0; bi < totalBlocks; ++bi) {
                    for (int eq = 0; eq < N; ++eq) {
                        int g_eq = mgr.getEquationIndex(bi, eq);
                        if (g_eq < 0 || g_eq >= totalEq) continue;
                        const double rWell = kWellSourceSignProbe * w_res[g_eq];
                        if (std::abs(rWell) <= 1e-16) continue;
                        probe_mat.AddResidual(bi, eq, rWell);
                    }
                }

                // 4) boundary source
                auto assembleBoundaryFieldProbe = [&](const BoundarySetting::BoundaryConditionManager* bcMgr, int dofOffset, const std::string& fieldName) {
                    if (!bcMgr || dofOffset < 0) return;
                    std::vector<double> bc_res(totalEq, 0.0);
                    std::vector<double> bc_diag(totalEq, 0.0);
                    BoundaryAssemblyStats bc_stats;
                    if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                        bc_stats = BoundaryAssembler::Assemble_2D(mgr, *bcMgr, dofOffset, fm, fieldName, bc_res, bc_diag);
                    }
                    else {
                        bc_stats = BoundaryAssembler::Assemble_3D(mgr, *bcMgr, dofOffset, fm, fieldName, bc_res, bc_diag);
                    }
                    (void)bc_stats;
                    for (int bi = 0; bi < totalBlocks; ++bi) {
                        const int eqIdx = mgr.getEquationIndex(bi, dofOffset);
                        if (eqIdx < 0 || eqIdx >= totalEq) continue;
                        const double r_bc = bc_res[eqIdx];
                        if (std::abs(r_bc) <= 1e-16) continue;
                        probe_mat.AddResidual(bi, dofOffset, r_bc);
                    }
                    };
                assembleBoundaryFieldProbe(modules.pressure_bc, pressureDof, pEqCfg.pressure_field);
                if constexpr (N == 3) {
                    assembleBoundaryFieldProbe(modules.saturation_bc, saturationDof, sEqCfg.saturation);
                }
                assembleBoundaryFieldProbe(modules.temperature_bc, temperatureDof, tEqCfg.temperatue_field);

                const Eigen::VectorXd b_probe = probe_mat.ExportEigenResidual();
                int max_probe_idx = 0;
                const double max_probe = b_probe.cwiseAbs().maxCoeff(&max_probe_idx);
                (void)max_probe_idx;
                if (!params.enable_row_scaling) {
                    return max_probe;
                }

                const Eigen::SparseMatrix<double> A_probe = probe_mat.ExportEigenSparseMatrix();
                double max_probe_scaled = 0.0;
                for (int r = 0; r < b_probe.size(); ++r) {
                    const double diag_abs = std::abs(A_probe.coeff(r, r));
                    const double denom = std::max(diag_abs, params.row_scale_floor);
                    const double scaled = std::abs(b_probe[r]) / std::max(denom, 1.0e-30);
                    if (scaled > max_probe_scaled) {
                        max_probe_scaled = scaled;
                    }
                }
                return max_probe_scaled;
                };

            for (int iter = 0; iter < params.max_newton_iter; ++iter) {
                iter_used++;
                global_mat.SetZero();

                // 【新增】：V3 诊断 - 方程贡献缓存区
                // 分配大小为所有自由度总数 (总网格数 * 每个网格的方程数 N)
                std::vector<EqContrib> eq_contribs(totalEq);

                // [V3] Diagnostics counters for summary and incident triggers
                int eos_fallback_water = 0;
                int eos_fallback_co2 = 0;
                int eos_near_bound_count = 0;
                int eos_total_samples = 0;

                double max_mass_flux = 0.0;
                double max_heat_flux = 0.0;
                int hot_mass_i = -1, hot_mass_j = -1;
                int hot_heat_i = -1, hot_heat_j = -1;
                ConnectionType hot_mass_type = ConnectionType::Matrix_Matrix;
                ConnectionType hot_heat_type = ConnectionType::Matrix_Matrix;

                // 组装积累项
                for (int bi = 0; bi < totalBlocks; ++bi) {
                    const double phi = 0.2, c_pr = 1000.0, rho_r = 2600.0;
                    ADVar<N> P(state.P[bi]); P.grad(0) = 1.0;
                    ADVar<N> T(state.T[bi]); T.grad((N == 2) ? 1 : 2) = 1.0;
                    auto pW = AD_Fluid::Evaluator::evaluateWater<N>(P, T);
                    ADVar<N> P_old(old_state.P[bi]), T_old(old_state.T[bi]);
                    auto pW_old = AD_Fluid::Evaluator::evaluateWater<N>(P_old, T_old);
                    if (params.diag_level != DiagLevel::Off) {
                        ++eos_total_samples;
                        if (pW.isFallback) ++eos_fallback_water;
                        if (pW.near_bound) ++eos_near_bound_count;
                    }

                    std::vector<ADVar<N>> acc_eqs(N);
                    if constexpr (N == 2) {
                        ADVar<N> m_w = pW.rho * phi, m_w_old = pW_old.rho * phi;
                        acc_eqs[0] = (m_w - m_w_old) * (vols[bi] / dt);
                        ADVar<N> e_w = m_w * (pW.h - P / pW.rho) + ADVar<N>((1.0 - phi) * rho_r * c_pr) * T;
                        ADVar<N> e_w_old = m_w_old * (pW_old.h - P_old / pW_old.rho) + ADVar<N>((1.0 - phi) * rho_r * c_pr) * T_old;
                        acc_eqs[1] = (e_w - e_w_old) * (vols[bi] / dt);
                    }
                    else {
                        ADVar<N> Sw(state.Sw[bi]); Sw.grad(1) = 1.0;
                        ADVar<N> Sg = ADVar<N>(1.0) - Sw;
                        ADVar<N> Sw_old(old_state.Sw[bi]), Sg_old = ADVar<N>(1.0) - Sw_old;

                        auto pG = AD_Fluid::Evaluator::evaluateCO2<N>(P, T);
                        auto pG_old = AD_Fluid::Evaluator::evaluateCO2<N>(P_old, T_old);
                        if (params.diag_level != DiagLevel::Off) {
                            ++eos_total_samples;
                            if (pG.isFallback) ++eos_fallback_co2;
                            if (pG.near_bound) ++eos_near_bound_count;
                        }

                        acc_eqs[0] = (pW.rho * phi * Sw - pW_old.rho * phi * Sw_old) * (vols[bi] / dt);
                        acc_eqs[1] = (pG.rho * phi * Sg - pG_old.rho * phi * Sg_old) * (vols[bi] / dt);

                        ADVar<N> e_fluid = pW.rho * Sw * (pW.h - P / pW.rho) + pG.rho * Sg * (pG.h - P / pG.rho);
                        ADVar<N> e_fluid_old = pW_old.rho * Sw_old * (pW_old.h - P_old / pW_old.rho) + pG_old.rho * Sg_old * (pG_old.h - P_old / pG_old.rho);
                        ADVar<N> e_rock = ADVar<N>((1.0 - phi) * rho_r * c_pr) * T;
                        ADVar<N> e_rock_old = ADVar<N>((1.0 - phi) * rho_r * c_pr) * T_old;
                        acc_eqs[2] = ((e_fluid * phi + e_rock) - (e_fluid_old * phi + e_rock_old)) * (vols[bi] / dt);
                    }
                    FIM_GlobalAssembler<N, ADVar<N>>::AssembleAccumulation(bi, acc_eqs, global_mat);
                    // Keep accumulation contributions always available for row scaling / diagnostics.
                    for (int eq = 0; eq < N; ++eq) {
                        int g_eq = mgr.getEquationIndex(bi, eq);
                        if (g_eq >= 0 && g_eq < eq_contribs.size()) {
                            eq_contribs[g_eq].R_acc += acc_eqs[eq].val;
                            eq_contribs[g_eq].D_acc += acc_eqs[eq].grad[eq];
                        }
                    }
                }

                // 组装通量项
                for (const auto& conn : connMgr.GetConnections()) {
                    int i = conn.nodeI, j = conn.nodeJ;
                    auto evalFlux = [&](bool wrt_i) -> std::vector<ADVar<N>> {
                        std::vector<ADVar<N>> F(N);
                        ADVar<N> P_i(state.P[i]), T_i(state.T[i]), P_j(state.P[j]), T_j(state.T[j]);
                        const Vector& x_i = blockCenters[i]; const Vector& x_j = blockCenters[j];

                        if constexpr (N == 2) {
                            if (wrt_i) { P_i.grad(0) = 1.0; T_i.grad(1) = 1.0; }
                            else { P_j.grad(0) = 1.0; T_j.grad(1) = 1.0; }
                            auto pW_i = AD_Fluid::Evaluator::evaluateWater<N>(P_i, T_i);
                            auto pW_j = AD_Fluid::Evaluator::evaluateWater<N>(P_j, T_j);
                            ADVar<N> rho_avg_w = ADVar<N>(0.5) * (pW_i.rho + pW_j.rho);
                            ADVar<N> dPhi = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, rho_avg_w, x_i, x_j, gravityVec);
                            ADVar<N> mob_i = pW_i.rho / pW_i.mu, mob_j = pW_j.rho / pW_j.mu;
                            ADVar<N> up_mob = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi, mob_i, mob_j);
                            F[0] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mob, dPhi);
                            ADVar<N> up_h = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi, pW_i.h, pW_j.h);
                            F[1] = FVM_Ops::Compute_Heat_Flux<N, ADVar<N>>(conn.T_Heat, T_i, T_j, F[0], up_h);
                        }
                        else {
                            ADVar<N> Sw_i(state.Sw[i]), Sw_j(state.Sw[j]);
                            if (wrt_i) { P_i.grad(0) = 1.0; Sw_i.grad(1) = 1.0; T_i.grad(2) = 1.0; }
                            else { P_j.grad(0) = 1.0; Sw_j.grad(1) = 1.0; T_j.grad(2) = 1.0; }
                            auto pW_i = AD_Fluid::Evaluator::evaluateWater<N>(P_i, T_i), pW_j = AD_Fluid::Evaluator::evaluateWater<N>(P_j, T_j);
                            auto pG_i = AD_Fluid::Evaluator::evaluateCO2<N>(P_i, T_i), pG_j = AD_Fluid::Evaluator::evaluateCO2<N>(P_j, T_j);

                            CapRelPerm::VGParams vg; vg.alpha = 1e-5; vg.n = 2.0; vg.Swr = 0.0; vg.Sgr = 0.0;
                            CapRelPerm::RelPermParams rp; rp.L = 0.5;
                            ADVar<N> krw_i, krg_i, krw_j, krg_j;
                            CapRelPerm::kr_Mualem_vG<N>(Sw_i, vg, rp, krw_i, krg_i);
                            CapRelPerm::kr_Mualem_vG<N>(Sw_j, vg, rp, krw_j, krg_j);

                            ADVar<N> rho_avg_w = ADVar<N>(0.5) * (pW_i.rho + pW_j.rho), rho_avg_g = ADVar<N>(0.5) * (pG_i.rho + pG_j.rho);
                            ADVar<N> Pc_i = CapRelPerm::pc_vG<N>(Sw_i, vg), Pc_j = CapRelPerm::pc_vG<N>(Sw_j, vg);
                            ADVar<N> dPhi_w = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, rho_avg_w, x_i, x_j, gravityVec);
                            ADVar<N> dPhi_g = FVM_Ops::Compute_Potential_Diff<N, ADVar<N>, Vector>(P_i, P_j, Pc_i, Pc_j, rho_avg_g, x_i, x_j, gravityVec);

                            ADVar<N> mobW_i = krw_i * pW_i.rho / pW_i.mu, mobW_j = krw_j * pW_j.rho / pW_j.mu;
                            ADVar<N> mobG_i = krg_i * pG_i.rho / pG_i.mu, mobG_j = krg_j * pG_j.rho / pG_j.mu;
                            ADVar<N> up_mobW = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_w, mobW_i, mobW_j);
                            ADVar<N> up_mobG = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_g, mobG_i, mobG_j);

                            F[0] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mobW, dPhi_w);
                            F[1] = FVM_Ops::Compute_Mass_Flux<N, ADVar<N>>(conn.T_Flow, up_mobG, dPhi_g);
                            ADVar<N> up_h_w = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_w, pW_i.h, pW_j.h);
                            ADVar<N> up_h_g = FVM_Ops::Op_Upwind_AD<N, ADVar<N>>(dPhi_g, pG_i.h, pG_j.h);
                            F[2] = FVM_Ops::Compute_Heat_Flux<N, ADVar<N>>(conn.T_Heat, T_i, T_j, F[0], F[1], up_h_w, up_h_g);
                        }
                        return F;
                        };
                    auto f_wrt_i = evalFlux(true);
                    auto f_wrt_j = evalFlux(false);
                    FIM_GlobalAssembler<N, ADVar<N>>::AssembleFlux(i, j, f_wrt_i, f_wrt_j, global_mat);

                    const double abs_mass_flux = [&]() {
                        if constexpr (N == 2) return std::abs(f_wrt_i[0].val);
                        else return std::max(std::abs(f_wrt_i[0].val), std::abs(f_wrt_i[1].val));
                        }();
                    const double abs_heat_flux = [&]() {
                        if constexpr (N == 2) return std::abs(f_wrt_i[1].val);
                        else return std::abs(f_wrt_i[2].val);
                        }();

                    if (abs_mass_flux > max_mass_flux) {
                        max_mass_flux = abs_mass_flux;
                        hot_mass_i = i;
                        hot_mass_j = j;
                        hot_mass_type = conn.type;
                    }
                    if (abs_heat_flux > max_heat_flux) {
                        max_heat_flux = abs_heat_flux;
                        hot_heat_i = i;
                        hot_heat_j = j;
                        hot_heat_type = conn.type;
                    }

                    // 【新增】：V3 诊断 - 记录 Flux 物理贡献
                    if (params.diag_level != DiagLevel::Off) {
                        for (int eq = 0; eq < N; ++eq) {
                            int g_eq_i = mgr.getEquationIndex(i, eq);
                            int g_eq_j = mgr.getEquationIndex(j, eq);
                            if (g_eq_i >= 0 && g_eq_i < eq_contribs.size()) {
                                eq_contribs[g_eq_i].R_flux += f_wrt_i[eq].val;
                                eq_contribs[g_eq_i].D_flux += f_wrt_i[eq].grad[eq];
                            }
                            if (g_eq_j >= 0 && g_eq_j < eq_contribs.size()) {
                                eq_contribs[g_eq_j].R_flux -= f_wrt_j[eq].val;      // 节点 j 接收 -F
                                eq_contribs[g_eq_j].D_flux -= f_wrt_j[eq].grad[eq]; // 对角线导数同理取负
                            }
                        }
                    }
                }

                // 组装井源项（全 Jacobian: d/dP, d/dSw, d/dT）
                SyncStateToFieldManager(state, fm, mgr);
                std::vector<double> w_res(totalEq, 0.0);
                std::vector<std::array<double, 3>> w_jac3(totalEq, std::array<double, 3>{ 0.0, 0.0, 0.0 });

                if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                    BoundaryAssembler::Assemble_Wells_2D_FullJac(mgr, fm, wells, 0, 0, (N == 3) ? 1 : -1, (N == 3) ? 2 : 1, w_res, w_jac3);
                }
                else {
                    BoundaryAssembler::Assemble_Wells_3D_FullJac(mgr, fm, wells, 0, 0, (N == 3) ? 1 : -1, (N == 3) ? 2 : 1, w_res, w_jac3);
                }

                double max_abs_well_dsw = 0.0;
                double max_abs_well_dt = 0.0;

                // BoundaryAssembler returns well terms in outflow-positive convention.
                // Residual here is assembled as Acc + Flux + Q_out = 0 by default.
                // Keep sign configurable for compatibility with custom conventions.
                const double kWellSourceSign = params.well_source_sign;

                for (int bi = 0; bi < totalBlocks; ++bi) {
                    for (int eq = 0; eq < N; ++eq) {
                        int g_eq = mgr.getEquationIndex(bi, eq);
                        if (g_eq < 0 || g_eq >= totalEq) continue;

                        const double dRdP = kWellSourceSign * w_jac3[g_eq][0];
                        const double dRdSw = kWellSourceSign * w_jac3[g_eq][1];
                        const double dRdT = kWellSourceSign * w_jac3[g_eq][2];
                        const double rWell = kWellSourceSign * w_res[g_eq];

                        max_abs_well_dsw = std::max(max_abs_well_dsw, std::abs(dRdSw));
                        max_abs_well_dt = std::max(max_abs_well_dt, std::abs(dRdT));

                        if (std::abs(rWell) <= 1e-16 &&
                            std::abs(dRdP) <= 1e-16 &&
                            std::abs(dRdSw) <= 1e-16 &&
                            std::abs(dRdT) <= 1e-16) {
                            continue;
                        }

                        global_mat.AddResidual(bi, eq, rWell);
                        global_mat.AddDiagJacobian(bi, eq, 0, dRdP);

                        if constexpr (N == 2) {
                            global_mat.AddDiagJacobian(bi, eq, 1, dRdT);
                        }
                        else {
                            global_mat.AddDiagJacobian(bi, eq, 1, dRdSw);
                            global_mat.AddDiagJacobian(bi, eq, 2, dRdT);
                        }

                        // 【新增】：V3 诊断 - 记录 Well 物理贡献
                        if (params.diag_level != DiagLevel::Off && g_eq >= 0 && g_eq < eq_contribs.size()) {
                            eq_contribs[g_eq].R_well += rWell;
                            if (eq == 0) eq_contribs[g_eq].D_well += dRdP;
                            else if constexpr (N == 2) { if (eq == 1) eq_contribs[g_eq].D_well += dRdT; }
                            else if constexpr (N == 3) {
                                if (eq == 1) eq_contribs[g_eq].D_well += dRdSw;
                                else if (eq == 2) eq_contribs[g_eq].D_well += dRdT;
                            }
                        }
                    }
                }

                if constexpr (N == 3) {
                    std::cout << "    [WellJac] max|dR/dSw|=" << std::scientific << max_abs_well_dsw
                        << " max|dR/dT|=" << max_abs_well_dt << "\n";
                }
                else {
                    std::cout << "    [WellJac] max|dR/dT|=" << std::scientific << max_abs_well_dt << "\n";
                }

                // 组装外部边界条件模块 (Dirichlet / Neumann / Robin)
                auto assembleBoundaryField = [&](const BoundarySetting::BoundaryConditionManager* bcMgr,
                    int dofOffset,
                    const std::string& fieldName,
                    const char* fieldLabel) {
                        if (!bcMgr || dofOffset < 0) return;

                        std::vector<double> bc_res(totalEq, 0.0);
                        std::vector<double> bc_diag(totalEq, 0.0);
                        BoundaryAssemblyStats bc_stats;

                        if constexpr (std::is_same_v<MeshMgrType, MeshManager>) {
                            bc_stats = BoundaryAssembler::Assemble_2D(mgr, *bcMgr, dofOffset, fm, fieldName, bc_res, bc_diag);
                        }
                        else {
                            bc_stats = BoundaryAssembler::Assemble_3D(mgr, *bcMgr, dofOffset, fm, fieldName, bc_res, bc_diag);
                        }

                        int appliedEq = 0;
                        for (int bi = 0; bi < totalBlocks; ++bi) {
                            const int eqIdx = mgr.getEquationIndex(bi, dofOffset);
                            if (eqIdx < 0 || eqIdx >= totalEq) continue;

                            const double r_bc = bc_res[eqIdx];
                            const double d_bc = bc_diag[eqIdx];
                            if (std::abs(r_bc) <= 1e-16 && std::abs(d_bc) <= 1e-16) continue;

                            global_mat.AddResidual(bi, dofOffset, r_bc);
                            global_mat.AddDiagJacobian(bi, dofOffset, dofOffset, d_bc);
                            ++appliedEq;

                            // 【新增】：V3 诊断 - 记录 BC 物理贡献
                            if (params.diag_level != DiagLevel::Off && eqIdx >= 0 && eqIdx < eq_contribs.size()) {
                                eq_contribs[eqIdx].R_bc += r_bc;
                                eq_contribs[eqIdx].D_bc += d_bc;
                            }
                        }
                        if (params.diag_level != DiagLevel::Off) {
                            std::cout << "    [BC-SUM] field=" << fieldLabel
                                << " faces=" << bc_stats.matrixBCCount + bc_stats.fractureBCCount
                                << " applied_eq=" << appliedEq
                                << " (visited=" << bc_stats.visitedEqRows
                                << ", nonzero=" << bc_stats.nonzeroEqRows
                                << ", zero_row=" << bc_stats.zeroEqRows
                                << ", invalid=" << bc_stats.invalidEqRows << ")\n";
                            for (const auto& kv : bc_stats.perTagType) {
                                const auto& v = kv.second;
                                std::cout << "      [BC-TAG] field=" << fieldLabel
                                    << " key=" << kv.first
                                    << " faces=" << v.faces
                                    << " applied=" << v.applied
                                    << " skipped=" << v.skipped
                                    << " |sumR|=" << v.sumR
                                    << " |sumDiag|=" << v.sumDiag << "\n";
                            }
                        }
                    };

                assembleBoundaryField(modules.pressure_bc, pressureDof, pEqCfg.pressure_field, "P");
                if constexpr (N == 3) {
                    assembleBoundaryField(modules.saturation_bc, saturationDof, sEqCfg.saturation, "Sw");
                }
                assembleBoundaryField(modules.temperature_bc, temperatureDof, tEqCfg.temperatue_field, "T");
                auto A = global_mat.ExportEigenSparseMatrix();
                auto b = global_mat.ExportEigenResidual();

                // 定位最大残差对应的方程与网格，用于定点深挖发散源
                int max_idx = 0;
                double max_res = b.cwiseAbs().maxCoeff(&max_idx);
                double max_res_scaled = max_res;
                if (params.enable_row_scaling) {
                    max_res_scaled = 0.0;
                    for (int r = 0; r < totalEq; ++r) {
                        const double diag_acc_abs = (r >= 0 && r < static_cast<int>(eq_contribs.size())) ? std::abs(eq_contribs[r].D_acc) : 0.0;
                        const double diag_abs = std::abs(A.coeff(r, r));
                        const double denom = std::max({ diag_acc_abs, diag_abs, params.row_scale_floor });
                        const double scaled = std::abs(b[r]) / std::max(denom, 1.0e-30);
                        if (scaled > max_res_scaled) {
                            max_res_scaled = scaled;
                        }
                    }
                }
                const double conv_res = params.enable_row_scaling ? max_res_scaled : max_res;

                if (iter == 0) res_iter1 = conv_res;
                step_final_residual = conv_res;
                if (conv_res < best_res) {
                    best_res = conv_res;
                    best_state = state;
                    best_iter = iter_used;
                }

                std::cout << "    [NL] step=" << step << " iter=" << iter_used
                    << " dt=" << std::scientific << dt << " res_inf=" << max_res
                    << " res_scaled=" << max_res_scaled
                    << " (at DOF=" << max_idx << ")\n";

                if (!std::isfinite(conv_res)) { fail_reason = "residual_nan_inf"; break; }
                if (conv_res < params.abs_res_tol) {
                    converged = true; converge_mode = "abs_res";
                    std::cout << "    [NL-CONVERGED] step=" << step << " mode=" << converge_mode << " res_scaled=" << conv_res << "\n";
                    break;
                }
                const double rel_res = (res_iter1 > 0.0) ? (conv_res / std::max(res_iter1, 1e-30)) : 1.0;
                if (iter_used > 1 && rel_res <= params.rel_res_tol && last_rel_update <= params.rel_update_tol) {
                    converged = true; converge_mode = "rel_res_update";
                    std::cout << "    [NL-CONVERGED] step=" << step << " mode=" << converge_mode
                        << " rel_res=" << rel_res << " rel_update=" << last_rel_update << "\n";
                    break;
                }
                // General best-iterate guard:
                // accept the in-step best state when later Newton iterations
                // clearly move away from it after an already acceptable drop.
                if (params.enable_best_iter_guard && iter_used >= params.best_iter_guard_min_iter && best_iter > 1 && res_iter1 > 0.0) {
                    const double grow_from_best = conv_res / std::max(best_res, 1e-30);
                    const double best_growth = best_res / std::max(res_iter1, 1e-30);
                    const double best_drop = 1.0 - best_growth;
                    const bool best_quality_ok = (best_res <= params.stagnation_abs_res_tol) &&
                        (best_res < res_iter1) &&
                        (best_drop >= params.stagnation_min_drop);
                    if (best_quality_ok && grow_from_best > params.best_iter_growth_trigger) {
                        state = best_state;
                        step_final_residual = best_res;
                        converged = true;
                        converge_mode = "best_iter_guard";
                        std::cout << "    [NL-CONVERGED] step=" << step
                            << " mode=" << converge_mode
                            << " best_iter=" << best_iter
                            << " best_res=" << best_res
                            << " grow_from_best=" << grow_from_best << "\n";
                        break;
                    }
                }
                // At dt floor, recover best iterate in-step when later Newton iterations
                // start blowing up after an already acceptable decrease.
                const bool at_dt_floor = (dt <= params.dt_min * (1.0 + 1.0e-12));
                if (params.enable_dt_floor_best && at_dt_floor && iter_used >= 3 && best_iter > 1) {
                    const double grow_from_best = conv_res / std::max(best_res, 1e-30);
                    const double best_growth = (res_iter1 > 0.0) ? (best_res / std::max(res_iter1, 1e-30)) : 1.0;
                    const double best_drop = 1.0 - best_growth;
                    const bool best_quality_ok = (best_res <= params.stagnation_abs_res_tol) &&
                        (best_res < res_iter1) &&
                        (best_drop >= params.stagnation_min_drop);
                    const bool accept_best = (grow_from_best > 2.0) && best_quality_ok;
                    if (accept_best) {
                        state = best_state;
                        step_final_residual = best_res;
                        converged = true;
                        converge_mode = "dt_floor_best";
                        std::cout << "    [NL-CONVERGED] step=" << step
                            << " mode=" << converge_mode
                            << " best_iter=" << best_iter
                            << " best_res=" << best_res
                            << " grow_from_best=" << grow_from_best << "\n";
                        break;
                    }
                }
                // Guarded dt-floor hold:
                // when residual at iter-1 is already within stagnation tolerance but Newton
                // corrections keep increasing residual, accept iter-1 state and move on.
                // This avoids repeated rollback loops at dt_min while still rejecting high residuals.
                if (params.enable_dt_floor_hold && at_dt_floor && iter_used >= 2 && best_iter == 1) {
                    const double grow_from_best = conv_res / std::max(best_res, 1e-30);
                    const bool best_quality_ok = (best_res <= params.stagnation_abs_res_tol);
                    const bool accept_hold = best_quality_ok && (grow_from_best > params.best_iter_growth_trigger);
                    if (accept_hold) {
                        state = best_state;
                        step_final_residual = best_res;
                        converged = true;
                        converge_mode = "dt_floor_hold";
                        std::cout << "    [NL-CONVERGED] step=" << step
                            << " mode=" << converge_mode
                            << " best_iter=" << best_iter
                            << " best_res=" << best_res
                            << " grow_from_best=" << grow_from_best << "\n";
                        break;
                    }
                }
                if (iter == params.max_newton_iter - 1) {
                    double growth = (res_iter1 > 0.0) ? (conv_res / std::max(res_iter1, 1e-30)) : 1.0;
                    const double rel_drop = 1.0 - growth;
                    const bool stagnation_ok = (growth <= params.stagnation_growth_tol) &&
                        (conv_res <= params.stagnation_abs_res_tol) &&
                        (conv_res < res_iter1) &&
                        (rel_drop >= params.stagnation_min_drop);
                    if (params.enable_stagnation_accept && stagnation_ok) {
                        converged = true; converge_mode = "stagnation";
                        std::cout << "    [NL-CONVERGED] step=" << step << " mode=" << converge_mode << " growth=" << growth << " rel_drop=" << rel_drop << " res_scaled=" << conv_res << "\n";
                    }
                    else {
                        fail_reason = "nonlinear_diverged";
                        std::cout << "    [NL-DIVERGED] step=" << step << " growth=" << growth << " rel_drop=" << rel_drop << " res_scaled=" << conv_res << "\n";
                    }
                    break;
                }

                // [修改] 参数驱动的鲁棒线性求解器设置
                Eigen::VectorXd dx;
                bool compute_ok = false;
                bool solve_ok = false;
                std::string solver_log;
                Eigen::SparseMatrix<double> A_solve = A;
                Eigen::VectorXd b_solve = b;
                bool row_scaling_applied = false;

                if (params.enable_row_scaling) {
                    Eigen::VectorXd row_scale(totalEq);
                    for (int r = 0; r < totalEq; ++r) {
                        const double diag_acc_abs = (r >= 0 && r < static_cast<int>(eq_contribs.size())) ? std::abs(eq_contribs[r].D_acc) : 0.0;
                        const double diag_abs = std::abs(A.coeff(r, r));
                        const double denom = std::max({ diag_acc_abs, diag_abs, params.row_scale_floor });
                        double s = 1.0 / std::max(denom, 1e-30);
                        s = std::max(params.row_scale_min, std::min(params.row_scale_max, s));
                        row_scale[r] = s;
                    }
                    Eigen::DiagonalMatrix<double, Eigen::Dynamic> S(row_scale);
                    A_solve = S * A;
                    b_solve = S * b;
                    row_scaling_applied = true;
                }

                if (params.lin_solver == LinearSolverType::SparseLU) {
                    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
                    solver.compute(A_solve);
                    compute_ok = (solver.info() == Eigen::Success);
                    if (compute_ok) {
                        dx = solver.solve(-b_solve);
                        solve_ok = (solver.info() == Eigen::Success);
                    }
                    solver_log = "solver=SparseLU compute_ok=" + std::string(compute_ok ? "true" : "false") +
                        " solve_ok=" + std::string(solve_ok ? "true" : "false") +
                        " nnzA=" + std::to_string(A_solve.nonZeros()) +
                        " scaled=" + std::string(row_scaling_applied ? "true" : "false") +
                        " info=" + std::to_string(solver.info());
                }
                else {
                    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
                    solver.preconditioner().setDroptol(params.bicgstab_droptol);
                    solver.compute(A_solve);
                    compute_ok = (solver.info() == Eigen::Success);
                    if (compute_ok) {
                        dx = solver.solve(-b_solve);
                        solve_ok = (solver.info() == Eigen::Success);
                    }
                    const std::string iters = compute_ok ? std::to_string(solver.iterations()) : "NA";
                    const std::string error = compute_ok ? std::to_string(solver.error()) : "NA";
                    solver_log = "solver=BiCGSTAB compute_ok=" + std::string(compute_ok ? "true" : "false") +
                        " solve_ok=" + std::string(solve_ok ? "true" : "false") +
                        " iters=" + iters +
                        " error=" + error +
                        " scaled=" + std::string(row_scaling_applied ? "true" : "false") +
                        " info=" + std::to_string(solver.info());
                }

                if (params.diag_level != DiagLevel::Off) {
                    const int printEvery = std::max(1, params.diag_print_every_iter);
                    const bool printSummary = ((iter % printEvery) == 0) || !solve_ok;
                    const double near_ratio = (eos_total_samples > 0) ? (static_cast<double>(eos_near_bound_count) / static_cast<double>(eos_total_samples)) : 0.0;
                    const bool eos_warn = (eos_fallback_water > 0) || (eos_fallback_co2 > 0) || (near_ratio >= params.diag_eos_near_bound_ratio);

                    const double mass_flux_spike = (prev_max_mass_flux > 1e-30) ? (max_mass_flux / prev_max_mass_flux) : 1.0;
                    const double heat_flux_spike = (prev_max_heat_flux > 1e-30) ? (max_heat_flux / prev_max_heat_flux) : 1.0;
                    const bool flux_spike = (mass_flux_spike > params.diag_flux_spike_factor) || (heat_flux_spike > params.diag_flux_spike_factor);

                    if (printSummary) {
                        std::cout << "    [SOLVER] " << solver_log << "\n";
                        std::cout << "    [FLUX-GLOBAL] max_mass_flux=" << max_mass_flux
                            << " (" << ConnectionTypeLabel(hot_mass_type) << ", i=" << hot_mass_i << ", j=" << hot_mass_j << ")"
                            << " max_heat_flux=" << max_heat_flux
                            << " (" << ConnectionTypeLabel(hot_heat_type) << ", i=" << hot_heat_i << ", j=" << hot_heat_j << ")\n";
                        if (eos_warn) {
                            std::cout << "    [EOS-WARN] water_fallback=" << eos_fallback_water
                                << " co2_fallback=" << eos_fallback_co2
                                << " near_bound=" << eos_near_bound_count << "/" << eos_total_samples
                                << " ratio=" << near_ratio << "\n";
                        }
                    }

                    bool trigger_hotspot = !solve_ok;
                    const double growth_factor = (res_iter1 > 0.0) ? (conv_res / std::max(res_iter1, 1e-30)) : 1.0;
                    if (iter_used > 1 && growth_factor > params.diag_blowup_factor) {
                        trigger_hotspot = true;
                    }
                    if (flux_spike) {
                        trigger_hotspot = true;
                    }

                    if (prev_hot_idx == max_idx) {
                        const double rel_change = std::abs(max_res - prev_hot_res) / std::max(std::abs(prev_hot_res), 1e-30);
                        if (rel_change <= params.diag_hot_res_change_tol) {
                            ++hot_repeat_count;
                        }
                        else {
                            hot_repeat_count = 1;
                        }
                    }
                    else {
                        hot_repeat_count = 1;
                    }
                    if (hot_repeat_count >= params.diag_hot_repeat_iters) {
                        trigger_hotspot = true;
                    }

                    if (trigger_hotspot) {
                        const int hot_block = max_idx / N;
                        const int hot_var = max_idx % N;
                        const std::string domain_type = (hot_block < MatrixBlockCount(mgr)) ? "Matrix" : "Fracture";
                        std::cout << "    >> [HOT] trigger=1 eq=" << max_idx
                            << " (block=" << hot_block << ", domain=" << domain_type << ", var=" << hot_var << ")\n";

                        const bool hot_block_valid = (hot_block >= 0 && hot_block < totalBlocks);
                        if (!hot_block_valid) {
                            std::cout << "    >> [HOT] skipped: invalid block index " << hot_block << "\n";
                        }

                        if (hot_block_valid && max_idx >= 0 && max_idx < static_cast<int>(eq_contribs.size())) {
                            std::cout << "    >> [HOT-R&J] R_total=" << eq_contribs[max_idx].R_total()
                                << " = acc(" << eq_contribs[max_idx].R_acc << ") + flux(" << eq_contribs[max_idx].R_flux
                                << ") + well(" << eq_contribs[max_idx].R_well << ") + bc(" << eq_contribs[max_idx].R_bc << ")\n";
                            std::cout << "                 D_total=" << eq_contribs[max_idx].D_total()
                                << " = acc(" << eq_contribs[max_idx].D_acc << ") + flux(" << eq_contribs[max_idx].D_flux
                                << ") + well(" << eq_contribs[max_idx].D_well << ") + bc(" << eq_contribs[max_idx].D_bc << ")\n";
                        }
                        else {
                            std::cout << "    >> [HOT-R&J] skipped: invalid eq index " << max_idx << "\n";
                        }

                        if (hot_block_valid && (!incident_dumped_this_step || !params.diag_incident_once_per_step)) {
                            nlohmann::json snap;
                            snap["DiagnosticLevel"] = "Forensic_V3";
                            snap["TimeControl"]["step"] = step;
                            snap["TimeControl"]["iter"] = iter_used;
                            snap["TimeControl"]["dt"] = dt;

                            snap["LinearSolver"] = solver_log;
                            snap["Hotspot"]["block_id"] = hot_block;
                            snap["Hotspot"]["domain"] = domain_type;
                            snap["Hotspot"]["P"] = state.P[hot_block];
                            snap["Hotspot"]["T"] = state.T[hot_block];
                            if constexpr (N == 3) {
                                snap["Hotspot"]["Sw"] = state.Sw[hot_block];
                            }
                            if (solve_ok && max_idx >= 0 && max_idx < dx.size()) {
                                snap["Hotspot"]["dx_raw"] = dx[max_idx];
                            }

                            if constexpr (N == 3) {
                                try {
                                    const double sw_val = state.Sw[hot_block];
                                    ADVar<3> Sw_ad = ADVar<3>(sw_val);
                                    Sw_ad.grad[0] = 0.0;
                                    Sw_ad.grad[1] = 1.0;
                                    Sw_ad.grad[2] = 0.0;

                                    CapRelPerm::VGParams vg;
                                    vg.alpha = 1e-5; vg.n = 2.0; vg.Swr = 0.0; vg.Sgr = 0.0;
                                    CapRelPerm::RelPermParams rp;
                                    rp.L = 0.5;

                                    ADVar<3> krw, krg;
                                    CapRelPerm::kr_Mualem_vG<3>(Sw_ad, vg, rp, krw, krg);
                                    ADVar<3> pc = CapRelPerm::pc_vG<3>(Sw_ad, vg);

                                    nlohmann::json const_snap;
                                    const_snap["Sw"] = Sw_ad.val;
                                    const_snap["krw"] = krw.val;
                                    const_snap["krg"] = krg.val;
                                    const_snap["Pc"] = pc.val;
                                    const_snap["dkrw_dSw"] = krw.grad[1];
                                    const_snap["dkrg_dSw"] = krg.grad[1];
                                    const_snap["dPc_dSw"] = pc.grad[1];
                                    if (std::isnan(krw.val) || std::isnan(pc.grad[1]) || std::isinf(pc.grad[1])) {
                                        const_snap["flag"] = "invalid_constitutive";
                                    }
                                    snap["IncidentConst"] = const_snap;
                                    std::cout << "    >> [INCIDENT-CONST] Sw=" << const_snap["Sw"]
                                        << " krw=" << const_snap["krw"]
                                        << " krg=" << const_snap["krg"]
                                        << " Pc=" << const_snap["Pc"] << "\n";
                                }
                                catch (...) {
                                    snap["IncidentConst"]["flag"] = "invalid_constitutive";
                                }
                            }

                            const std::string snap_path = "Test/Transient/Day6/" + caseName + "/fail_snapshot_step" +
                                std::to_string(step) + "_iter" + std::to_string(iter_used) + ".json";
                            std::ofstream ofs(snap_path);
                            if (ofs.is_open()) {
                                ofs << std::setw(4) << snap << std::endl;
                                incident_dumped_this_step = true;
                                std::cout << "    >> [INCIDENT] snapshot=" << snap_path << "\n";
                            }
                            else {
                                std::cout << "    >> [INCIDENT] failed to write snapshot file\n";
                            }
                        }
                    }

                    prev_hot_idx = max_idx;
                    prev_hot_res = max_res;
                    prev_max_mass_flux = max_mass_flux;
                    prev_max_heat_flux = max_heat_flux;
                }

                if (!solve_ok) { fail_reason = "linear_solve_fail"; break; }
                if (!dx.allFinite()) { fail_reason = "dx_nan_inf"; break; }

                // [修改] 应用基于时间步长感知的状态截断与保护(Damping)
                bool state_valid = true;
                double alpha = 1.0;
                // Time-step-aware damping:
                // accumulation scales with 1/dt, so update caps should shrink at least linearly with dt.
                const double dt_ref = std::max(params.dt_init, params.dt_min);
                const double dt_eff = std::max(dt, params.dt_min);
                const double damp_scale = std::max(1.0e-12, dt_eff / dt_ref);
                // Keep tiny but non-zero floors for stiff late-time rollback loops.
                const double max_dP_eff = std::max(1.0e-3, params.max_dP * damp_scale);
                const double max_dT_eff = std::max(1.0e-5, params.max_dT * damp_scale);
                const double max_dSw_eff = std::max(1.0e-6, params.max_dSw * damp_scale);
                for (int bi = 0; bi < totalBlocks; ++bi) {
                    int eqP = mgr.getEquationIndex(bi, 0), eqT = mgr.getEquationIndex(bi, (N == 3) ? 2 : 1);
                    if (eqP < 0 || eqP >= dx.size() || eqT < 0 || eqT >= dx.size()) { state_valid = false; fail_reason = "invalid_eq_index"; break; }
                    alpha = std::min(alpha, max_dP_eff / (std::abs(dx[eqP]) + 1e-14));
                    alpha = std::min(alpha, max_dT_eff / (std::abs(dx[eqT]) + 1e-14));
                    if constexpr (N == 3) {
                        int eqSw = mgr.getEquationIndex(bi, 1);
                        if (eqSw < 0 || eqSw >= dx.size()) { state_valid = false; fail_reason = "invalid_eq_index"; break; }
                        alpha = std::min(alpha, max_dSw_eff / (std::abs(dx[eqSw]) + 1e-14));
                    }
                }
                if (!state_valid) break;
                alpha = std::min(1.0, alpha);
                if (!std::isfinite(alpha) || alpha < params.min_alpha) {
                    fail_reason = "alpha_too_small";
                    break;
                }

                const auto apply_trial_update = [&](const FIM_StateMap<N>& base_state,
                    double alpha_try,
                    FIM_StateMap<N>& out_state,
                    int& limiter_added_local,
                    double& rel_update_inf) -> bool {
                        out_state = base_state;
                        limiter_added_local = 0;
                        rel_update_inf = 0.0;

                        for (int bi = 0; bi < totalBlocks; ++bi) {
                            const int eqP = mgr.getEquationIndex(bi, 0);
                            const int eqT = mgr.getEquationIndex(bi, (N == 3) ? 2 : 1);
                            if (eqP < 0 || eqP >= dx.size() || eqT < 0 || eqT >= dx.size()) return false;

                            const double dP = alpha_try * dx[eqP];
                            const double p_ref = std::max(std::abs(base_state.P[bi]), 1.0e5);
                            rel_update_inf = std::max(rel_update_inf, std::abs(dP) / p_ref);
                            out_state.P[bi] += dP;
                            if (!std::isfinite(out_state.P[bi])) return false;
                            if (params.clamp_state_to_eos_bounds) {
                                if (out_state.P[bi] < p_floor) { out_state.P[bi] = p_floor; limiter_added_local++; }
                                if (out_state.P[bi] > p_ceil) { out_state.P[bi] = p_ceil; limiter_added_local++; }
                            }

                            const double dT = alpha_try * dx[eqT];
                            const double t_ref = std::max(std::abs(base_state.T[bi]), 300.0);
                            rel_update_inf = std::max(rel_update_inf, std::abs(dT) / t_ref);

                            if constexpr (N == 3) {
                                const int eqSw = mgr.getEquationIndex(bi, 1);
                                if (eqSw < 0 || eqSw >= dx.size()) return false;
                                const double dSw = alpha_try * dx[eqSw];
                                rel_update_inf = std::max(rel_update_inf, std::abs(dSw));
                                out_state.Sw[bi] += dSw;
                                if (!std::isfinite(out_state.Sw[bi])) return false;
                                if (out_state.Sw[bi] < 0.0) { out_state.Sw[bi] = 0.0; limiter_added_local++; }
                                if (out_state.Sw[bi] > 1.0) { out_state.Sw[bi] = 1.0; limiter_added_local++; }
                            }

                            out_state.T[bi] += dT;
                            if (!std::isfinite(out_state.T[bi])) return false;
                            if (params.clamp_state_to_eos_bounds) {
                                if (out_state.T[bi] < t_floor) { out_state.T[bi] = t_floor; limiter_added_local++; }
                                if (out_state.T[bi] > t_ceil) { out_state.T[bi] = t_ceil; limiter_added_local++; }
                            }
                        }
                        return true;
                    };

                const FIM_StateMap<N> state_before_update = state;
                double accepted_alpha = alpha;
                int accepted_limiter_added = 0;
                double accepted_rel_update = std::numeric_limits<double>::infinity();
                bool update_accepted = false;

                if (params.enable_armijo_line_search) {
                    const int max_bt = std::max(1, params.armijo_max_backtracks);
                    const double bt_beta = std::max(0.1, std::min(0.95, params.armijo_beta));
                    const double c1 = std::max(1.0e-8, std::min(1.0e-2, params.armijo_c1));

                    double alpha_try = alpha;
                    for (int bt = 0; bt < max_bt; ++bt) {
                        FIM_StateMap<N> trial_state = state_before_update;
                        int trial_limiter = 0;
                        double trial_rel_update = std::numeric_limits<double>::infinity();
                        if (!apply_trial_update(state_before_update, alpha_try, trial_state, trial_limiter, trial_rel_update)) {
                            alpha_try *= bt_beta;
                            if (alpha_try < params.min_alpha) break;
                            continue;
                        }

                        const double trial_res = compute_residual_inf_for_state(trial_state);
                        const double armijo_rhs = (1.0 - c1 * alpha_try) * conv_res;
                        if (std::isfinite(trial_res) && trial_res <= armijo_rhs) {
                            state = std::move(trial_state);
                            accepted_alpha = alpha_try;
                            accepted_limiter_added = trial_limiter;
                            accepted_rel_update = trial_rel_update;
                            update_accepted = true;
                            break;
                        }

                        alpha_try *= bt_beta;
                        if (alpha_try < params.min_alpha) break;
                    }

                    if (!update_accepted) {
                        state = state_before_update;
                        fail_reason = "line_search_fail";
                        break;
                    }
                }
                else {
                    FIM_StateMap<N> trial_state = state_before_update;
                    int trial_limiter = 0;
                    double trial_rel_update = std::numeric_limits<double>::infinity();
                    if (!apply_trial_update(state_before_update, alpha, trial_state, trial_limiter, trial_rel_update)) {
                        state = state_before_update;
                        fail_reason = "state_nan_inf";
                        break;
                    }
                    state = std::move(trial_state);
                    accepted_alpha = alpha;
                    accepted_limiter_added = trial_limiter;
                    accepted_rel_update = trial_rel_update;
                    update_accepted = true;
                }

                if (!update_accepted) { fail_reason = "state_nan_inf"; break; }
                total_limiters += accepted_limiter_added;
                last_rel_update = accepted_rel_update;

                if (params.diag_level != DiagLevel::Off && accepted_limiter_added >= params.diag_clamp_trigger) {
                    std::cout << "    [Limiter] added=" << accepted_limiter_added << " alpha=" << accepted_alpha << "\n";
                }
            }

            if (!converged && fail_reason.empty()) fail_reason = "nonlinear_max_iter";

            if (!converged) {
                total_rollbacks++;
                dt = std::max(dt * 0.5, params.dt_min);
                state = old_state;
                step--;
                std::cout << "    [Rollback] step=" << (step + 1) << " new_dt=" << dt << " reason=" << fail_reason << "\n";
                if (dt <= params.dt_min && total_rollbacks > 20) throw std::runtime_error("[FAIL] dt reached lower bound with repeated rollback.");
                if (total_rollbacks > 80) throw std::runtime_error("[FAIL] Max rollbacks exceeded.");
            }
            else {
                t += dt;
                completed_steps = step;
                final_residual = step_final_residual;
                std::cout << "  [Step Success] step=" << std::setw(3) << step << " dt=" << std::scientific << std::setprecision(3) << dt
                    << " iter_used=" << iter_used << " conv_mode=" << converge_mode
                    << " rollback_count=" << total_rollbacks << " limiter_count=" << total_limiters << "\n";

                if (converge_mode == "abs_res") {
                    if (iter_used <= 3) dt = std::min(dt * 1.2, params.dt_max);
                    else if (iter_used <= 5) dt = std::min(dt * 1.05, params.dt_max);
                    else dt = std::max(dt * 0.8, params.dt_min);
                }
                else if (converge_mode == "rel_res_update") {
                    // rel_res_update often appears with moderate Newton iterations in stiff cases.
                    // Keep dt neutral around 7~10 iterations to avoid artificial collapse to dt_min.
                    if (iter_used <= 6) dt = std::min(dt * 1.08, params.dt_max);
                    else if (iter_used <= 10) dt = dt;
                    else if (iter_used <= 14) dt = std::max(dt * 0.92, params.dt_min);
                    else dt = std::max(dt * 0.8, params.dt_min);
                }
                else if (converge_mode == "stagnation" || converge_mode == "dt_floor_best" || converge_mode == "best_iter_guard" || converge_mode == "dt_floor_hold") {
                    dt = std::max(dt * 0.85, params.dt_min);
                }
                else {
                    dt = std::max(dt * 0.5, params.dt_min);
                }

                if (step % 10 == 0 || step == params.max_steps) {
                    SyncStateToFieldManager(state, fm, mgr);
                    std::string fname = "Test/Transient/Day6/" + caseName + ((step == params.max_steps) ? "/final.vtk" : ("/step_" + std::to_string(step) + ".vtk"));
                    if constexpr (std::is_same_v<MeshMgrType, MeshManager>) PostProcess_2D(mgr, fm).ExportVTK(fname, t);
                    else PostProcess_3D(mgr, fm).ExportVTK(fname, t);
                    VerifyVtkExport(fname, N == 3);
                    vtk_export_count++;
                    std::cout << "    [VTK Export PASS] " << fname << "\n";
                }
            }
        }
        std::cout << "[PASS] Day6 case completed: " << caseName
            << " | steps=" << completed_steps
            << " | rollbacks=" << total_rollbacks
            << " | limiters=" << total_limiters
            << " | final_residual=" << std::scientific << final_residual
            << " | vtk_exports=" << vtk_export_count
            << " | t_end=" << std::scientific << t << " s\n";
    }

} // namespace FIM_Engine
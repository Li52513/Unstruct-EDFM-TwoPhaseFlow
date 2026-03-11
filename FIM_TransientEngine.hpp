/**
 * @file FIM_TransientEngine.hpp
 * @brief 生产级全隐式(FIM)瞬态推进核心引擎 (高度泛化、配置驱动)
 */

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
#include "FIM_TransientSupport.hpp" // [修改] 引用通用支持头文件

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>           // [新增] 引入对病态矩阵极度鲁棒的直接求解器

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#define MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define MKDIR(path) mkdir(path, 0777)
#endif

namespace FIM_Engine { // [修改] 使用通用的 FIM 引擎命名空间

    // =====================================================================
    // 求解器配置与参数模块 (脱离硬编码的核心)
    // =====================================================================

    enum class SolverRoute { FIM, IMPES };
    enum class LinearSolverType { SparseLU, BiCGSTAB };

    /** @brief 物理场初始条件配置 */
    struct InitialConditions {
        double P_init = 2.0e5;   ///< 初始网格压力 (Pa)
        double T_init = 300.0;   ///< 初始网格温度 (K)
        double Sw_init = 0.2;    ///< 初始含水饱和度 (-)
    };

    /** @brief 瞬态牛顿求解器核心控制参数 */
    struct TransientSolverParams {
        int max_steps = 50;                     ///< 模拟最大推进步数

        // 时间步长控制
        double dt_init = 1.0;                   ///< [修复启动冲击] 初始时间步长降至1.0s，应对注水瞬间水锤效应
        double dt_min = 1e-4;                   ///< 最小允许时间步长 (s)
        double dt_max = 86400.0;                ///< 最大允许时间步长 (s)

        // 非线性控制
        int max_newton_iter = 8;                ///< Day6 计划要求：每步最多 8 次非线性迭代
        double abs_res_tol = 1e-6;              ///< Day6 计划要求：绝对残差收敛阈值

        // 防停滞(Stagnation)控制机制
        double stagnation_growth_tol = 0.995;   ///< 停滞容忍生长率
        double stagnation_abs_res_tol = 1.0e6;  ///< 停滞时允许的最大收敛残差 (单相建议 1e4, 两相建议 1e6)
        double stagnation_min_drop = 2.0e-3;    ///< 必须满足的最小相对残差下降率

        // 状态截断与护栏(Limiter & Damping)
        double max_dP = 1.0e4;                  ///< [收紧截断] 牛顿步内最大允许的压力变化量(Pa)
        double max_dT = 2.0;                    ///< [收紧截断] 牛顿步内最大允许的温度变化量(K)
        double max_dSw = 0.05;                  ///< 牛顿步内最大允许的饱和度变化量(-)

        // 线性求解器设置
        LinearSolverType lin_solver = LinearSolverType::SparseLU; ///< [默认SparseLU] 对多物理场极度病态系统免疫
        double bicgstab_droptol = 1e-2;                           ///< 若回退至 BiCGSTAB 时使用的 ILUT 容差
        double well_source_sign = -1.0;                          ///< 井源项符号修正系数: R += sign * Q_well (默认-1用于Qout->Qin对齐)
    };

    /**
     * @brief Day6 可选外部模块注入: 物性初始化与第一/二/三类边界条件
     * @details
     * - property_initializer: 在传输系数计算之前执行，用于显式注入基岩/裂缝物性
     * - pressure_bc / saturation_bc / temperature_bc: 对应方程变量的边界管理器
     */
    template <typename MeshMgrType, typename FieldMgrType>
    struct TransientOptionalModules {
        std::function<void(MeshMgrType&, FieldMgrType&)> property_initializer;
        const BoundarySetting::BoundaryConditionManager* pressure_bc = nullptr;
        const BoundarySetting::BoundaryConditionManager* saturation_bc = nullptr;
        const BoundarySetting::BoundaryConditionManager* temperature_bc = nullptr;
    };

    // =====================================================================
    // Traits 与内部辅助函数
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

        for (int step = 1; step <= params.max_steps; ++step) {
            bool converged = false;
            FIM_StateMap<N> old_state = state;
            std::string fail_reason;
            std::string converge_mode = "none";
            int iter_used = 0;
            double res_iter1 = -1.0;
            double step_final_residual = std::numeric_limits<double>::quiet_NaN();

            for (int iter = 0; iter < params.max_newton_iter; ++iter) {
                iter_used++;
                global_mat.SetZero();

                // 组装积累项
                for (int bi = 0; bi < totalBlocks; ++bi) {
                    const double phi = 0.2, c_pr = 1000.0, rho_r = 2600.0;
                    ADVar<N> P(state.P[bi]); P.grad(0) = 1.0;
                    ADVar<N> T(state.T[bi]); T.grad((N == 2) ? 1 : 2) = 1.0;
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
                        ADVar<N> Sw(state.Sw[bi]); Sw.grad(1) = 1.0;
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
                    FIM_GlobalAssembler<N, ADVar<N>>::AssembleAccumulation(bi, acc_eqs, global_mat);
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
                    FIM_GlobalAssembler<N, ADVar<N>>::AssembleFlux(i, j, evalFlux(true), evalFlux(false), global_mat);
                }

                // 组装井源项（全 Jacobian: d/dP, d/dSw, d/dT）
                SyncStateToFieldManager(state, fm, mgr);
                int totalEq = mgr.getTotalEquationDOFs();
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
                // FIM residual uses Acc + Flux - Qin = 0. Since Qin = -Qout,
                // apply a consistent -1 sign to well residual and Jacobian.
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
                        }

                        std::cout << "    [BC] field=" << fieldLabel
                            << " matrix_faces=" << bc_stats.matrixBCCount
                            << " applied_eq=" << appliedEq << "\n";
                    };

                assembleBoundaryField(modules.pressure_bc, pressureDof, pEqCfg.pressure_field, "P");
                if constexpr (N == 3) {
                    assembleBoundaryField(modules.saturation_bc, saturationDof, sEqCfg.saturation, "Sw");
                }
                assembleBoundaryField(modules.temperature_bc, temperatureDof, tEqCfg.temperatue_field, "T");
                auto A = global_mat.ExportEigenSparseMatrix();
                auto b = global_mat.ExportEigenResidual();

                // 定位最大残差对应的网格/方程，便于调试发散源
                int max_idx = 0;
                double max_res = b.cwiseAbs().maxCoeff(&max_idx);
                if (iter == 0) res_iter1 = max_res;
                step_final_residual = max_res;

                std::cout << "    [NL] step=" << step << " iter=" << iter_used
                    << " dt=" << std::scientific << dt << " res_inf=" << max_res
                    << " (at DOF=" << max_idx << ")\n";

                if (!std::isfinite(max_res)) { fail_reason = "residual_nan_inf"; break; }
                if (max_res < params.abs_res_tol) {
                    converged = true; converge_mode = "abs_res";
                    std::cout << "    [NL-CONVERGED] step=" << step << " mode=" << converge_mode << " res_inf=" << max_res << "\n";
                    break;
                }
                if (iter == params.max_newton_iter - 1) {
                    double growth = (res_iter1 > 0.0) ? (max_res / std::max(res_iter1, 1e-30)) : 1.0;
                    const double rel_drop = 1.0 - growth;
                    const bool stagnation_ok = (growth <= params.stagnation_growth_tol) &&
                        (max_res <= params.stagnation_abs_res_tol) &&
                        (max_res < res_iter1) &&
                        (rel_drop >= params.stagnation_min_drop);
                    if (stagnation_ok) {
                        converged = true; converge_mode = "stagnation";
                        std::cout << "    [NL-CONVERGED] step=" << step << " mode=" << converge_mode << " growth=" << growth << " rel_drop=" << rel_drop << " res_inf=" << max_res << "\n";
                    }
                    else {
                        fail_reason = "nonlinear_diverged";
                        std::cout << "    [NL-DIVERGED] step=" << step << " growth=" << growth << " rel_drop=" << rel_drop << " res_inf=" << max_res << "\n";
                    }
                    break;
                }

                // [修改] 参数驱动的鲁棒线性求解器
                Eigen::VectorXd dx;
                if (params.lin_solver == LinearSolverType::SparseLU) {
                    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
                    solver.compute(A);
                    if (solver.info() != Eigen::Success) { fail_reason = "linear_solve_fail"; break; }
                    dx = solver.solve(-b);
                }
                else {
                    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
                    solver.preconditioner().setDroptol(params.bicgstab_droptol);
                    solver.compute(A);
                    if (solver.info() != Eigen::Success) { fail_reason = "linear_solve_fail"; break; }
                    dx = solver.solve(-b);
                }

                if (!dx.allFinite()) { fail_reason = "dx_nan_inf"; break; }

                // [修改] 配置参数驱动的状态截断与护栏(Damping)
                bool state_valid = true;
                double alpha = 1.0;
                for (int bi = 0; bi < totalBlocks; ++bi) {
                    int eqP = mgr.getEquationIndex(bi, 0), eqT = mgr.getEquationIndex(bi, (N == 3) ? 2 : 1);
                    if (eqP < 0 || eqP >= dx.size() || eqT < 0 || eqT >= dx.size()) { state_valid = false; fail_reason = "invalid_eq_index"; break; }
                    alpha = std::min(alpha, params.max_dP / (std::abs(dx[eqP]) + 1e-14));
                    alpha = std::min(alpha, params.max_dT / (std::abs(dx[eqT]) + 1e-14));
                    if constexpr (N == 3) {
                        int eqSw = mgr.getEquationIndex(bi, 1);
                        if (eqSw < 0 || eqSw >= dx.size()) { state_valid = false; fail_reason = "invalid_eq_index"; break; }
                        alpha = std::min(alpha, params.max_dSw / (std::abs(dx[eqSw]) + 1e-14));
                    }
                }
                if (!state_valid) break;
                alpha = std::max(1e-3, std::min(1.0, alpha));

                for (int bi = 0; bi < totalBlocks; ++bi) {
                    int eqP = mgr.getEquationIndex(bi, 0), eqT = mgr.getEquationIndex(bi, (N == 3) ? 2 : 1);
                    state.P[bi] += alpha * dx[eqP];
                    if (!std::isfinite(state.P[bi])) { state_valid = false; break; }
                    if (state.P[bi] < 1e4) { state.P[bi] = 1e4; total_limiters++; }

                    if constexpr (N == 3) {
                        int eqSw = mgr.getEquationIndex(bi, 1);
                        state.Sw[bi] += alpha * dx[eqSw];
                        if (!std::isfinite(state.Sw[bi])) { state_valid = false; break; }
                        if (state.Sw[bi] < 0.0) { state.Sw[bi] = 0.0; total_limiters++; }
                        if (state.Sw[bi] > 1.0) { state.Sw[bi] = 1.0; total_limiters++; }
                        state.T[bi] += alpha * dx[eqT];
                    }
                    else { state.T[bi] += alpha * dx[eqT]; }

                    if (!std::isfinite(state.T[bi])) { state_valid = false; break; }
                    if (state.T[bi] < 273.15) { state.T[bi] = 273.15; total_limiters++; }
                }
                if (!state_valid) { fail_reason = "state_nan_inf"; break; }
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
                else if (converge_mode == "stagnation") {
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


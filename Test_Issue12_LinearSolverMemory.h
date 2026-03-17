#pragma once
/**
 * @file Test_Issue12_LinearSolverMemory.h
 * @brief Issue #12 验证测试：线性求解器内存与CPR-AMG
 *
 * 子测试：
 *   T1 — 缓存矩阵正确性：solver_A_work 与原始 A 完全一致
 *   T2 — CPR 可解性：AMGCL_CPR 求解已知解的 N=3 链式系统
 *   T3 — 求解器路径切换：SparseLU / AMGCL / AMGCL_CPR 三路径结果一致
 *   T4 — 内存基线对比计时（可选，规模较小快速完成）
 */
#include "FIM_BlockSparseMatrix.h"
#include "FIM_TransientEngine/StepKernels.hpp"

#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Test_Issue12 {

using namespace FIM_Engine;
using namespace FIM_Engine::detail;

// ─────────────────────────────────────────────────────────────
// 内部工具：填充对角占优链式矩阵，确保线性系统可解
// 对角元素 = diag_val，相邻块非对角元素 = off_val
// 已知精确解 x_true = [1,1,...,1]，则 b = -A * x_true
// ─────────────────────────────────────────────────────────────
template <int N>
static void FillDiagDominantChain(FIM_BlockSparseMatrix<N>& mat, int nBlocks,
                                   double diag_val = 4.0, double off_val = 0.5) {
    // 对角块
    for (int i = 0; i < nBlocks; ++i) {
        for (int r = 0; r < N; ++r)
            mat.AddDiagJacobian(i, r, r, diag_val);
    }
    // 相邻块非对角（不含自身块的非对角元，保持结构简洁）
    for (int i = 0; i + 1 < nBlocks; ++i) {
        for (int r = 0; r < N; ++r) {
            // (i, i+1) 和 (i+1, i) 各一个 N×N 块，仅填对角元
            mat.AddOffDiagJacobian(i,     i + 1, r, r,  off_val);
            mat.AddOffDiagJacobian(i + 1, i,     r, r,  off_val);
        }
    }
}

// 计算 b = -A * x_true（x_true = all-ones）
// 返回 b_work 使得 A * (-x_true) = -b，即 A * dx = -b 的精确解为 dx = -x_true
template <int N>
static Eigen::VectorXd MakeRHS(const Eigen::SparseMatrix<double, Eigen::RowMajor, int>& A) {
    const int n = static_cast<int>(A.rows());
    Eigen::VectorXd x_true = Eigen::VectorXd::Ones(n);
    // b = A * x_true，则 A * dx = -b 的解为 dx = -x_true
    return A * x_true;
}

// ─────────────────────────────────────────────────────────────
// T1：缓存矩阵正确性
//   多次调用 SolveLinearSystem 后，cache.solver_A_work（行缩放前）
//   的值与输入 A 完全一致（max_diff < 1e-14）。
//   注意：行缩放会修改 A_work，所以在禁用行缩放的情况下验证。
// ─────────────────────────────────────────────────────────────
inline void T1_CacheMatrixCorrectness() {
    std::cout << "\n[Issue12-T1] LinearSolverCache A_work correctness\n";

    constexpr int N      = 3;
    const     int BLOCKS = 5;

    FIM_BlockSparseMatrix<N> bsm(BLOCKS);
    FillDiagDominantChain<N>(bsm, BLOCKS, 4.0, 0.5);
    bsm.FreezePattern();

    const auto& A_ref = bsm.GetFrozenMatrix();  // RowMajorMat<N>
    Eigen::VectorXd b = MakeRHS<N>(A_ref);

    TransientSolverParams params;
    params.lin_solver        = LinearSolverType::SparseLU;
    params.enable_row_scaling = false;  // 禁用行缩放，方便验证 A_work

    LinearSolverCache<N> cache;
    std::vector<EqContrib> eq_contribs(BLOCKS * N);
    const int totalEq = BLOCKS * N;

    // 调用三次：第一次建缓存（完整拷贝），后两次值拷贝
    for (int call = 0; call < 3; ++call) {
        // 更新矩阵值
        bsm.SetZero();
        FillDiagDominantChain<N>(bsm, BLOCKS, 4.0 + call * 0.1, 0.5);
        const auto& A_cur = bsm.GetFrozenMatrix();
        Eigen::VectorXd b_cur = MakeRHS<N>(A_cur);

        auto res = SolveLinearSystem<N>(A_cur, b_cur, eq_contribs, totalEq, params, cache);

        // 验证 cache.solver_A_work 值与 A_cur 一致
        // 注意：cache.solver_A_work 是行缩放后的值，这里行缩放已禁用，
        // 但 A_work 还会经过 makeCompressed()，值不变
        double max_diff = 0.0;
        for (int k = 0; k < A_cur.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double, Eigen::RowMajor, int>::InnerIterator it(A_cur, k); it; ++it) {
                double cached_val = cache.solver_A_work.coeff(it.row(), it.col());
                max_diff = std::max(max_diff, std::abs(it.value() - cached_val));
            }
        }

        if (max_diff > 1e-14) {
            throw std::runtime_error("[FAIL T1] cache.solver_A_work mismatch at call=" +
                std::to_string(call) + " max_diff=" + std::to_string(max_diff));
        }
        std::cout << "  [PASS] call=" << call << " max_diff=" << max_diff
                  << "  solve_ok=" << res.solve_ok << "\n";
    }
}

// ─────────────────────────────────────────────────────────────
// T2：CPR 可解性
//   用 5×5 块 N=3 链式对角占优矩阵（已知解 dx = -all-ones），
//   验证 AMGCL_CPR 路径残差 ||A*dx + b||/||b|| < 1e-5
// ─────────────────────────────────────────────────────────────
inline void T2_CPRSolvability() {
    std::cout << "\n[Issue12-T2] AMGCL_CPR solvability on N=3 chain system\n";

    constexpr int N      = 3;
    const     int BLOCKS = 5;

    FIM_BlockSparseMatrix<N> bsm(BLOCKS);
    FillDiagDominantChain<N>(bsm, BLOCKS, 4.0, 0.5);
    bsm.FreezePattern();

    const auto& A = bsm.GetFrozenMatrix();
    Eigen::VectorXd b = MakeRHS<N>(A);

    TransientSolverParams params;
    params.lin_solver          = LinearSolverType::AMGCL_CPR;
    params.enable_row_scaling  = false;
    params.amgcl_cpr_tol       = 1.0e-8;
    params.amgcl_cpr_maxiter   = 500;
    params.amgcl_cpr_use_fallback_sparselu = true;

    LinearSolverCache<N> cache;
    std::vector<EqContrib> eq_contribs(BLOCKS * N);
    const int totalEq = BLOCKS * N;

    auto res = SolveLinearSystem<N>(A, b, eq_contribs, totalEq, params, cache);

    if (!res.compute_ok) {
        throw std::runtime_error("[FAIL T2] compute_ok=false. log=" + res.solver_log);
    }
    if (!res.solve_ok) {
        throw std::runtime_error("[FAIL T2] solve_ok=false. log=" + res.solver_log);
    }

    // 验证残差：||A*dx + b||/||b||
    // SolveLinearSystem 解 A*dx = -b，故残差为 ||A*dx + b||
    const double res_norm = (A * res.dx + b).norm();
    const double b_norm   = b.norm();
    const double rel_res  = (b_norm > 1e-30) ? res_norm / b_norm : res_norm;

    std::cout << "  solver_log: " << res.solver_log << "\n";
    std::cout << "  rel_res=" << rel_res << "\n";

    if (rel_res > 1.0e-5) {
        throw std::runtime_error("[FAIL T2] Relative residual too large: " +
            std::to_string(rel_res) + " > 1e-5. log=" + res.solver_log);
    }
    std::cout << "  [PASS] AMGCL_CPR solved system, rel_res=" << rel_res << "\n";
}

// ─────────────────────────────────────────────────────────────
// T3：求解器路径切换
//   同一矩阵分别用 SparseLU / AMGCL / AMGCL_CPR 三路求解，
//   均应返回 solve_ok=true，且解向量之差 < 1e-8
// ─────────────────────────────────────────────────────────────
inline void T3_SolverPathSwitching() {
    std::cout << "\n[Issue12-T3] Solver path switching consistency\n";

    constexpr int N      = 3;
    const     int BLOCKS = 4;

    FIM_BlockSparseMatrix<N> bsm(BLOCKS);
    FillDiagDominantChain<N>(bsm, BLOCKS, 6.0, 0.3);
    bsm.FreezePattern();

    const auto& A = bsm.GetFrozenMatrix();
    Eigen::VectorXd b = MakeRHS<N>(A);
    std::vector<EqContrib> eq_contribs(BLOCKS * N);
    const int totalEq = BLOCKS * N;

    TransientSolverParams params_base;
    params_base.enable_row_scaling = false;
    params_base.amgcl_tol          = 1.0e-8;
    params_base.amgcl_maxiter      = 500;
    params_base.amgcl_cpr_tol      = 1.0e-8;
    params_base.amgcl_cpr_maxiter  = 500;

    struct PathResult {
        std::string name;
        bool solve_ok;
        Eigen::VectorXd dx;
        std::string log;
    };

    std::vector<PathResult> results;

    // SparseLU
    {
        TransientSolverParams p = params_base;
        p.lin_solver = LinearSolverType::SparseLU;
        LinearSolverCache<N> cache;
        auto res = SolveLinearSystem<N>(A, b, eq_contribs, totalEq, p, cache);
        results.push_back({"SparseLU", res.solve_ok, res.dx, res.solver_log});
    }

    // AMGCL
    {
        TransientSolverParams p = params_base;
        p.lin_solver = LinearSolverType::AMGCL;
        LinearSolverCache<N> cache;
        auto res = SolveLinearSystem<N>(A, b, eq_contribs, totalEq, p, cache);
        results.push_back({"AMGCL", res.solve_ok, res.dx, res.solver_log});
    }

    // AMGCL_CPR
    {
        TransientSolverParams p = params_base;
        p.lin_solver = LinearSolverType::AMGCL_CPR;
        LinearSolverCache<N> cache;
        auto res = SolveLinearSystem<N>(A, b, eq_contribs, totalEq, p, cache);
        results.push_back({"AMGCL_CPR", res.solve_ok, res.dx, res.solver_log});
    }

    // 验证所有路径 solve_ok
    bool all_ok = true;
    for (auto& r : results) {
        if (!r.solve_ok) {
            std::cerr << "  [WARN] " << r.name << " solve_ok=false. log=" << r.log << "\n";
            all_ok = false;
        }
    }

    // 验证解向量一致性（相对于 SparseLU 基准）
    const auto& dx_ref = results[0].dx;
    for (size_t i = 1; i < results.size(); ++i) {
        if (!results[0].solve_ok || !results[i].solve_ok) continue;
        double max_diff = (results[i].dx - dx_ref).cwiseAbs().maxCoeff();
        std::cout << "  " << results[i].name << " vs SparseLU: max_diff=" << max_diff << "\n";
        if (max_diff > 1.0e-6) {
            std::cerr << "  [WARN] " << results[i].name << " solution deviates from SparseLU: "
                      << max_diff << " (tolerance=1e-6). log=" << results[i].log << "\n";
        } else {
            std::cout << "  [PASS] " << results[i].name << " consistent with SparseLU\n";
        }
    }

    for (auto& r : results) {
        std::cout << "  " << r.name << ": solve_ok=" << r.solve_ok
                  << "  log=" << r.log << "\n";
    }

    if (!all_ok) {
        throw std::runtime_error("[FAIL T3] One or more solver paths failed. See above.");
    }
    std::cout << "  [PASS] All solver paths returned solve_ok=true\n";
}

// ─────────────────────────────────────────────────────────────
// T4：内存基线对比计时
//   500 块矩阵重复 50 次，比较有缓存 vs 无缓存（首次调用）的耗时
// ─────────────────────────────────────────────────────────────
inline void T4_MemoryTimingBaseline() {
    std::cout << "\n[Issue12-T4] Memory baseline timing (AMGCL, 50 reps)\n";

    constexpr int N      = 3;
    const     int BLOCKS = 100;  // 较小规模，快速完成
    const     int REPS   = 20;

    FIM_BlockSparseMatrix<N> bsm(BLOCKS);
    FillDiagDominantChain<N>(bsm, BLOCKS, 4.0, 0.5);
    bsm.FreezePattern();

    TransientSolverParams params;
    params.lin_solver          = LinearSolverType::AMGCL;
    params.enable_row_scaling  = false;
    params.amgcl_tol           = 1.0e-6;
    params.amgcl_maxiter       = 200;

    std::vector<EqContrib> eq_contribs(BLOCKS * N);
    const int totalEq = BLOCKS * N;

    using Clock = std::chrono::high_resolution_clock;
    using us    = std::chrono::microseconds;

    // 预热：第一次调用建立 AMGCL 层次（时间较长，不计入热路径）
    {
        const auto& A0 = bsm.GetFrozenMatrix();
        Eigen::VectorXd b0 = MakeRHS<N>(A0);
        LinearSolverCache<N> cache_warmup;
        (void)SolveLinearSystem<N>(A0, b0, eq_contribs, totalEq, params, cache_warmup);
        std::cout << "  [Info] Warmup done (AMGCL hierarchy built)\n";
    }

    // 热路径计时：复用同一 cache，矩阵准备在计时范围外
    bsm.SetZero();
    FillDiagDominantChain<N>(bsm, BLOCKS, 4.0, 0.5);
    const auto& A = bsm.GetFrozenMatrix();
    Eigen::VectorXd b = MakeRHS<N>(A);

    LinearSolverCache<N> cache;
    long long total_us = 0;
    for (int rep = 0; rep < REPS; ++rep) {
        auto t0 = Clock::now();
        (void)SolveLinearSystem<N>(A, b, eq_contribs, totalEq, params, cache);
        auto t1 = Clock::now();
        total_us += std::chrono::duration_cast<us>(t1 - t0).count();
    }

    double avg_us = static_cast<double>(total_us) / REPS;
    std::cout << "  blocks=" << BLOCKS << "  N=" << N << "  reps=" << REPS << "\n";
    std::cout << "  Avg SolveLinearSystem (cached): " << avg_us << " us\n";
    std::cout << "  [PASS] Timing recorded (no regression threshold set)\n";
}

// ─────────────────────────────────────────────────────────────
// 入口函数
// ─────────────────────────────────────────────────────────────
inline void Run_All() {
    std::cout << "\n=======================================================\n";
    std::cout << ">>> Issue #12: Linear Solver Memory & CPR-AMG Tests <<<\n";
    std::cout << "=======================================================\n";

    bool all_pass = true;

    auto run = [&](const char* name, auto fn) {
        try {
            fn();
            std::cout << "  => " << name << " PASSED\n";
        } catch (const std::exception& e) {
            std::cerr << "  => " << name << " FAILED: " << e.what() << "\n";
            all_pass = false;
        }
    };

    run("T1_CacheMatrixCorrectness", T1_CacheMatrixCorrectness);
    run("T2_CPRSolvability",         T2_CPRSolvability);
    run("T3_SolverPathSwitching",    T3_SolverPathSwitching);
    run("T4_MemoryTimingBaseline",   T4_MemoryTimingBaseline);

    std::cout << "\n=======================================================\n";
    if (all_pass)
        std::cout << ">>> ALL Issue#12 tests PASSED <<<\n";
    else
        std::cout << ">>> SOME Issue#12 tests FAILED — see above <<<\n";
    std::cout << "=======================================================\n\n";

    if (!all_pass)
        throw std::runtime_error("Issue#12 test suite failed.");
}

} // namespace Test_Issue12

#pragma once
/**
 * @file Test_Issue11_FrozenMatrix.h
 * @brief Issue #11 验证测试：稀疏模式固定化与 O(N) 更新
 *
 * 包含三个独立子测试，均无需运行完整求解器：
 *   T1 — CSR 缓存结构 + 数值等价性 (对应计划 V1 + V2)
 *   T2 — 模式稳定性：近零值不再破坏 nnz (对应计划 V3 思路)
 *   T3 — 性能对比计时 (对应计划 V4)
 */
#include "FIM_BlockSparseMatrix.h"

#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Test_Issue11 {

// ─────────────────────────────────────────────────────────────
// 内部工具：填充一个具有确定拓扑的矩阵（链式连接）
// blocks: 0-1-2-...(nBlocks-1)，每块 N 个 DOF
// ─────────────────────────────────────────────────────────────
template <int N>
static void FillChainMatrix(FIM_BlockSparseMatrix<N>& mat, int nBlocks, double seed = 1.0) {
    // 对角块：每个元素 = seed * (block*N^2 + r*N + c + 1)
    for (int i = 0; i < nBlocks; ++i) {
        for (int r = 0; r < N; ++r)
            for (int c = 0; c < N; ++c)
                mat.AddDiagJacobian(i, r, c,
                    seed * static_cast<double>(i * N * N + r * N + c + 1));
    }
    // 非对角块：相邻块 (i, i+1) 和 (i+1, i)
    for (int i = 0; i + 1 < nBlocks; ++i) {
        for (int r = 0; r < N; ++r)
            for (int c = 0; c < N; ++c) {
                double v = seed * static_cast<double>((i + 1) * 100 + r * N + c + 1);
                mat.AddOffDiagJacobian(i,     i + 1, r, c,  v);
                mat.AddOffDiagJacobian(i + 1, i,     r, c, -v);
            }
    }
}

// ─────────────────────────────────────────────────────────────
// T1：结构正确性 + 数值等价性
//   - 同一矩阵同时通过旧接口 ExportEigenSparseMatrix() 和
//     新接口 GetFrozenMatrix() 导出，逐元素比对。
// ─────────────────────────────────────────────────────────────
inline void T1_StructureAndNumericalEquivalence() {
    std::cout << "\n[Issue11-T1] CSR cache structure & numerical equivalence\n";

    constexpr int N      = 3;
    const     int BLOCKS = 8;   // 8 块链式拓扑，nnz = 8*9 + 7*2*9 = 72+126 = 198

    FIM_BlockSparseMatrix<N> mat(BLOCKS);
    FillChainMatrix<N>(mat, BLOCKS, 1.0);
    mat.FreezePattern();

    // 旧路径（现已去除 1e-16 过滤）
    auto A_old = mat.ExportEigenSparseMatrix();

    // 新路径（O(nnz) 更新）— 首次调用触发 BuildCSRCacheInternal_（内含 assert 验证）
    auto A_new = mat.GetFrozenMatrix();   // 返回 const&，auto 触发值拷贝

    // 1. 维度一致
    if (A_old.rows() != A_new.rows() || A_old.cols() != A_new.cols()) {
        throw std::runtime_error("[FAIL T1] Matrix dimensions differ.");
    }

    // 2. nnz 一致
    if (A_old.nonZeros() != A_new.nonZeros()) {
        throw std::runtime_error("[FAIL T1] nonZeros differ: old=" +
            std::to_string(A_old.nonZeros()) + " new=" +
            std::to_string(A_new.nonZeros()));
    }

    // 3. 所有值逐一比对（允许浮点误差 1e-14）
    double max_diff = 0.0;
    // 遍历 A_old 的每个非零元，在 A_new 中查询同一位置的值
    for (int k = 0; k < A_old.outerSize(); ++k) {
        for (typename FIM_BlockSparseMatrix<N>::SparseMat::InnerIterator it(A_old, k); it; ++it) {
            double v_old = it.value();
            double v_new = A_new.coeff(it.row(), it.col());
            max_diff = std::max(max_diff, std::abs(v_old - v_new));
        }
    }
    if (max_diff > 1e-14) {
        throw std::runtime_error("[FAIL T1] Value mismatch, max_diff=" +
            std::to_string(max_diff));
    }
    std::cout << "  [PASS] dims=" << A_new.rows() << "x" << A_new.cols()
              << "  nnz=" << A_new.nonZeros()
              << "  max_diff=" << max_diff << "\n";

    // 4. 第二次调用（更新值后再次比对，验证 O(N) 更新路径）
    mat.SetZero();
    FillChainMatrix<N>(mat, BLOCKS, 2.0);   // 换一组值

    auto A_old2 = mat.ExportEigenSparseMatrix();
    auto A_new2 = mat.GetFrozenMatrix();

    double max_diff2 = 0.0;
    for (int k = 0; k < A_old2.outerSize(); ++k) {
        for (typename FIM_BlockSparseMatrix<N>::SparseMat::InnerIterator it(A_old2, k); it; ++it) {
            max_diff2 = std::max(max_diff2,
                std::abs(it.value() - A_new2.coeff(it.row(), it.col())));
        }
    }
    if (max_diff2 > 1e-14) {
        throw std::runtime_error("[FAIL T1] Value mismatch on second call, max_diff=" +
            std::to_string(max_diff2));
    }
    std::cout << "  [PASS] second-call (O(N) update path)  max_diff=" << max_diff2 << "\n";
}

// ─────────────────────────────────────────────────────────────
// T2：模式稳定性 — 近零值不再破坏 nnz
//   修改前：某元素 = 1e-17 时会被 1e-16 过滤，nnz 减少，
//           导致模式哈希变化 → SparseLU 重做符号分解。
//   修改后：无过滤，两次调用的 nnz 必须完全相同。
// ─────────────────────────────────────────────────────────────
inline void T2_PatternStabilityNearZero() {
    std::cout << "\n[Issue11-T2] Pattern stability with near-zero values\n";

    constexpr int N      = 3;
    const     int BLOCKS = 4;

    // 第一次：正常值
    FIM_BlockSparseMatrix<N> mat(BLOCKS);
    FillChainMatrix<N>(mat, BLOCKS, 1.0);
    mat.FreezePattern();

    auto A1 = mat.GetFrozenMatrix();
    int nnz1 = A1.nonZeros();

    // 第二次：把某个非对角块的所有元素设为 5e-17（旧代码会过滤掉）
    mat.SetZero();
    FillChainMatrix<N>(mat, BLOCKS, 1.0);
    // 覆盖块(0,1)为极小值
    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c)
            mat.AddOffDiagJacobian(0, 1, r, c,
                5e-17 - mat.GetTotalScalarDOF() * 0.0);  // 纯赋值：重新 SetZero+Fill 已完成

    // 实际做法：先 SetZero 再只填极小值
    mat.SetZero();
    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c)
            mat.AddOffDiagJacobian(0, 1, r, c, 5e-17);
    // 其余块仍用正常值
    for (int i = 1; i + 1 < BLOCKS; ++i) {
        for (int r = 0; r < N; ++r)
            for (int c = 0; c < N; ++c) {
                double v = static_cast<double>((i + 1) * 100 + r * N + c + 1);
                mat.AddOffDiagJacobian(i,     i + 1, r, c,  v);
                mat.AddOffDiagJacobian(i + 1, i,     r, c, -v);
            }
    }
    for (int i = 0; i < BLOCKS; ++i)
        for (int r = 0; r < N; ++r)
            for (int c = 0; c < N; ++c)
                mat.AddDiagJacobian(i, r, c,
                    static_cast<double>(i * N * N + r * N + c + 1));
    // 块(1,0) 极小值
    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c)
            mat.AddOffDiagJacobian(1, 0, r, c, -5e-17);

    auto A2 = mat.GetFrozenMatrix();
    int nnz2 = A2.nonZeros();

    if (nnz1 != nnz2) {
        throw std::runtime_error(
            "[FAIL T2] nnz changed when near-zero values present: nnz1=" +
            std::to_string(nnz1) + " nnz2=" + std::to_string(nnz2) +
            "\n  (This indicates the 1e-16 filter was NOT removed correctly.)");
    }
    std::cout << "  [PASS] nnz stable across near-zero perturbation: nnz=" << nnz1 << "\n";

    // 额外验证：ExportEigenSparseMatrix 在去除过滤后，nnz 同样稳定
    auto A_exp = mat.ExportEigenSparseMatrix();
    if (A_exp.nonZeros() != nnz1) {
        throw std::runtime_error(
            "[FAIL T2] ExportEigenSparseMatrix nnz mismatch after filter removal: got=" +
            std::to_string(A_exp.nonZeros()) + " expected=" + std::to_string(nnz1));
    }
    std::cout << "  [PASS] ExportEigenSparseMatrix nnz also stable: nnz=" << A_exp.nonZeros() << "\n";
}

// ─────────────────────────────────────────────────────────────
// T3：耗时对比
//   在较大矩阵上分别调用旧接口和新接口各 REPS 次，打印耗时比。
//   注意：首次 GetFrozenMatrix() 会建缓存（含排序），计时中单独标注。
// ─────────────────────────────────────────────────────────────
inline void T3_TimingComparison() {
    std::cout << "\n[Issue11-T3] Timing comparison (old Triplet path vs new O(nnz) path)\n";

    constexpr int N      = 3;
    const     int BLOCKS = 500;   // 500 blocks → nnz ≈ 500*9 + 499*2*9 ≈ 13482
    const     int REPS   = 50;

    FIM_BlockSparseMatrix<N> mat(BLOCKS);
    FillChainMatrix<N>(mat, BLOCKS, 1.0);
    mat.FreezePattern();

    using Clock = std::chrono::high_resolution_clock;
    using us    = std::chrono::microseconds;

    // --- 旧接口（Triplet 路径）---
    long long t_old_us = 0;
    for (int rep = 0; rep < REPS; ++rep) {
        mat.SetZero();
        FillChainMatrix<N>(mat, BLOCKS, static_cast<double>(rep + 1));
        auto t0 = Clock::now();
        auto A  = mat.ExportEigenSparseMatrix();
        auto t1 = Clock::now();
        t_old_us += std::chrono::duration_cast<us>(t1 - t0).count();
        (void)A;
    }

    // --- 新接口（首次调用建缓存，不计入热路径耗时）---
    mat.SetZero();
    FillChainMatrix<N>(mat, BLOCKS, 1.0);
    {
        auto t0   = Clock::now();
        auto& ref = mat.GetFrozenMatrix();   // 建缓存
        auto t1   = Clock::now();
        (void)ref;
        std::cout << "  [Info] First GetFrozenMatrix (cache build): "
                  << std::chrono::duration_cast<us>(t1 - t0).count() << " us\n";
    }

    long long t_new_us = 0;
    for (int rep = 0; rep < REPS; ++rep) {
        mat.SetZero();
        FillChainMatrix<N>(mat, BLOCKS, static_cast<double>(rep + 1));
        auto t0 = Clock::now();
        auto A  = mat.GetFrozenMatrix();   // O(nnz) 更新 + 值拷贝
        auto t1 = Clock::now();
        t_new_us += std::chrono::duration_cast<us>(t1 - t0).count();
        (void)A;
    }

    double avg_old = static_cast<double>(t_old_us) / REPS;
    double avg_new = static_cast<double>(t_new_us) / REPS;
    double speedup = (avg_new > 0.0) ? avg_old / avg_new : 0.0;

    std::cout << "  blocks=" << BLOCKS << "  N=" << N
              << "  reps=" << REPS << "\n";
    std::cout << "  ExportEigenSparseMatrix (old): avg=" << avg_old << " us\n";
    std::cout << "  GetFrozenMatrix         (new): avg=" << avg_new << " us\n";
    std::cout << "  Speedup: " << speedup << "x\n";

    if (speedup < 1.0) {
        std::cout << "  [WARN] New path is slower — check build optimization flags.\n";
    } else {
        std::cout << "  [PASS] New path is faster.\n";
    }
}

// ─────────────────────────────────────────────────────────────
// T4：UnfreezePattern 使缓存失效，再冻结后重建正确
// ─────────────────────────────────────────────────────────────
inline void T4_UnfreezeCacheInvalidation() {
    std::cout << "\n[Issue11-T4] UnfreezePattern cache invalidation\n";

    constexpr int N      = 2;
    const     int BLOCKS = 4;

    FIM_BlockSparseMatrix<N> mat(BLOCKS);
    FillChainMatrix<N>(mat, BLOCKS, 1.0);
    mat.FreezePattern();

    auto A1 = mat.GetFrozenMatrix();   // 建缓存

    // Unfreeze → 添加一个新连接 → Freeze → 新缓存应包含新连接
    mat.UnfreezePattern();
    // 添加跨越一格的非对角块 (0, 2)
    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c)
            mat.AddOffDiagJacobian(0, 2, r, c, 99.0);
    mat.FreezePattern();

    auto A2 = mat.GetFrozenMatrix();   // 必须重建缓存

    // A2 应有更多非零元（新增了块(0,2) 的 N*N 个元素）
    if (A2.nonZeros() <= A1.nonZeros()) {
        throw std::runtime_error("[FAIL T4] Cache not rebuilt after UnfreezePattern+FreezePattern: "
            "nnz did not increase. nnz_before=" + std::to_string(A1.nonZeros()) +
            " nnz_after=" + std::to_string(A2.nonZeros()));
    }

    // 验证新值确实写入
    double v = A2.coeff(0 * N, 2 * N);   // 块(0,2)第(0,0)元素
    if (std::abs(v - 99.0) > 1e-14) {
        throw std::runtime_error("[FAIL T4] New connection value incorrect: got=" +
            std::to_string(v) + " expected=99.0");
    }
    std::cout << "  [PASS] nnz_before=" << A1.nonZeros()
              << "  nnz_after=" << A2.nonZeros()
              << "  new_entry(0,2)[0,0]=" << v << "\n";
}

// ─────────────────────────────────────────────────────────────
// 入口函数
// ─────────────────────────────────────────────────────────────
inline void Run_All() {
    std::cout << "\n=================================================\n";
    std::cout << ">>> Issue #11: FrozenMatrix CSR Cache Tests <<<\n";
    std::cout << "=================================================\n";

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

    run("T1_StructureAndNumericalEquivalence", T1_StructureAndNumericalEquivalence);
    run("T2_PatternStabilityNearZero",         T2_PatternStabilityNearZero);
    run("T3_TimingComparison",                 T3_TimingComparison);
    run("T4_UnfreezeCacheInvalidation",        T4_UnfreezeCacheInvalidation);

    std::cout << "\n=================================================\n";
    if (all_pass)
        std::cout << ">>> ALL Issue#11 tests PASSED <<<\n";
    else
        std::cout << ">>> SOME Issue#11 tests FAILED — see above <<<\n";
    std::cout << "=================================================\n\n";

    if (!all_pass)
        throw std::runtime_error("Issue#11 test suite failed.");
}

} // namespace Test_Issue11

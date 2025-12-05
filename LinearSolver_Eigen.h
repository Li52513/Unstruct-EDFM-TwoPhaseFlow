#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include "Solver_AssemblerCOO.h"

// 线性求解配置：默认 BiCGSTAB + ILUT，可切换 SparseLU/LDLT/CG
struct LinearSolverOptions {
    enum class Type { CG, BiCGSTAB, SparseLU, LDLT } type = Type::BiCGSTAB;
    int    maxIters = 2000;
    double tol = 1e-8;
    int    iluFill = 10;     // ILUT fill-factor
    double iluDrop = 1e-4;   // ILUT drop tolerance
    bool   equil = true;   // 是否做对角均衡（强烈建议开启）
    int    reusePreconditioner = 0; // >0 时复用 ILUT 若干步（步计数在内部维护）
};

// 使用 Eigen 求解 SparseSystemCOO：x 为输入初值/输出解
inline bool solveCOO_Eigen
(
    const SparseSystemCOO& sys,
    std::vector<double>& x,
    const LinearSolverOptions& opt,
    int* iters_out = nullptr,
    double* error_out = nullptr
)
{
    using SpMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;
    using Trip = Eigen::Triplet<double>;

    if (iters_out) *iters_out = 0;
    if (error_out) *error_out = 0.0;

    // ---- 空系统防护 ----
    if (sys.n <= 0) {
        x.assign(0, 0.0);
        return true;
    }

    // ---- Triplets → 稀疏矩阵 ----
    std::vector<Trip> triplets;
    triplets.reserve(sys.A.size());
    for (const auto& t : sys.A) triplets.emplace_back(t.r, t.c, t.v);

    SpMat A(sys.n, sys.n);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();
    // 释放峰值 triplet 内存（可选）
    triplets.clear(); triplets.shrink_to_fit();

    // 再次空矩阵防护
    if (A.rows() == 0 || A.cols() == 0) {
        x.assign(sys.n, 0.0);
        return true;
    }

    Eigen::VectorXd b(sys.n);
    for (int i = 0; i < sys.n; ++i)
        b[i] = (i < static_cast<int>(sys.b.size())) ? sys.b[i] : 0.0;

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(sys.n);
    if (x.size() == static_cast<size_t>(sys.n))
        for (int i = 0; i < sys.n; ++i) x0[i] = x[i];

    // ---------- 对角均衡（对称缩放） ----------
    //  Â = D^{-1/2} A D^{-1/2},   b̂ = D^{-1/2} b,   x = D^{-1/2} x̂
    SpMat Aeq = A;
    Eigen::VectorXd beq = b, x0eq = x0;
    Eigen::VectorXd D = Eigen::VectorXd::Ones(sys.n);

    if (opt.equil) {
        D = A.diagonal().cwiseAbs();
        for (int i = 0; i < D.size(); ++i)
            if (D[i] <= 1e-30) D[i] = 1.0;

        const Eigen::VectorXd Dinvsqrt = D.cwiseInverse().cwiseSqrt();
        const Eigen::VectorXd Dsqrt = D.cwiseSqrt();

        // Â = D^{-1/2} A D^{-1/2}  （两次稀疏遍历，保持稀疏）
        for (int k = 0; k < Aeq.outerSize(); ++k)
            for (SpMat::InnerIterator it(Aeq, k); it; ++it)
                it.valueRef() *= Dinvsqrt[it.row()];
        for (int k = 0; k < Aeq.outerSize(); ++k)
            for (SpMat::InnerIterator it(Aeq, k); it; ++it)
                it.valueRef() *= Dinvsqrt[it.col()];

        beq = Dinvsqrt.asDiagonal() * b;     // b̂
        x0eq = Dsqrt.asDiagonal() * x0;    // ★ 注意：x̂0 = D^{+1/2} x0
    }

    auto recoverSolution = [&](const Eigen::VectorXd& sol_eq, Eigen::VectorXd& sol_vec) {
        if (opt.equil) {
            const Eigen::VectorXd Dinvsqrt = D.cwiseInverse().cwiseSqrt();
            sol_vec = Dinvsqrt.asDiagonal() * sol_eq;  // x = D^{-1/2} x̂
        }
        else {
            sol_vec = sol_eq;
        }
        x.assign(sol_vec.data(), sol_vec.data() + sol_vec.size());
        };

    int    iters = 0;
    double errInf = 0.0;
    bool   ok = false;

    switch (opt.type) 
    {
    case LinearSolverOptions::Type::CG:
    {
        // 仅适合 SPD 系统（压力方程若是 TPFA 扩散可近似SPD，温度含对流通常非对称）
        Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper,
            Eigen::DiagonalPreconditioner<double>> solver;
        solver.setMaxIterations(opt.maxIters);
        solver.setTolerance(opt.tol);
        solver.compute(Aeq);
        Eigen::VectorXd sol_eq = solver.solveWithGuess(beq, x0eq);
        iters = solver.iterations();
        ok = (solver.info() == Eigen::Success);

        if (ok) {
            Eigen::VectorXd sol;
            recoverSolution(sol_eq, sol);
            errInf = (A * sol - b).lpNorm<Eigen::Infinity>();
        }
        break;
    }

    case LinearSolverOptions::Type::BiCGSTAB:
    {
        // 预条件器/分析缓存（静态，仅在本 TU 内复用）
        static Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double>> solver;
        static int            lastN = -1;
        static Eigen::Index   lastNNZ = -1;
        static double         lastFill = -1.0;
        static double         lastDrop = -1.0;
        static int            reuseCountdown = 0;
        static bool           haveFactor = false;

        // 复用策略：维度或 nnz 或 ILUT 参数变化 → 重建
        bool recompute = !haveFactor
            || lastN != static_cast<int>(Aeq.rows())
            || lastNNZ != Aeq.nonZeros()
            || lastFill != static_cast<double>(opt.iluFill)
            || lastDrop != opt.iluDrop
            || opt.reusePreconditioner <= 0
            || reuseCountdown <= 0;

        solver.preconditioner().setDroptol(opt.iluDrop);
        solver.preconditioner().setFillfactor(opt.iluFill);
        solver.setMaxIterations(opt.maxIters);
        solver.setTolerance(opt.tol);

        Eigen::VectorXd sol_eq;
        if (recompute) {
            solver.compute(Aeq);
            haveFactor = (solver.info() == Eigen::Success);
            reuseCountdown = std::max(0, opt.reusePreconditioner - 1);
            lastN = static_cast<int>(Aeq.rows());
            lastNNZ = Aeq.nonZeros();
            lastFill = static_cast<double>(opt.iluFill);
            lastDrop = opt.iluDrop;
        }
        else {
            --reuseCountdown;
        }

        sol_eq = solver.solveWithGuess(beq, x0eq);
        iters = solver.iterations();
        ok = (solver.info() == Eigen::Success);

        // 失败一次时：自动重建预条件器再尝试一次
        if (!ok && haveFactor && opt.reusePreconditioner > 0) {
            solver.compute(Aeq);
            haveFactor = (solver.info() == Eigen::Success);
            reuseCountdown = std::max(0, opt.reusePreconditioner - 1);
            lastN = static_cast<int>(Aeq.rows());
            lastNNZ = Aeq.nonZeros();
            lastFill = static_cast<double>(opt.iluFill);
            lastDrop = opt.iluDrop;

            if (haveFactor) {
                sol_eq = solver.solveWithGuess(beq, x0eq);
                iters = solver.iterations();
                ok = (solver.info() == Eigen::Success);
            }
        }

        if (ok) {
            Eigen::VectorXd sol;
            recoverSolution(sol_eq, sol);
            errInf = (A * sol - b).lpNorm<Eigen::Infinity>();
        }
        break;
    }

    case LinearSolverOptions::Type::SparseLU:
    {
        // pattern 复用：结构变化则重建 analyzePattern
        static Eigen::SparseLU<SpMat> solver;
        static Eigen::Index luRows = -1, luCols = -1, luNNZ = -1;

        if (luRows != Aeq.rows() || luCols != Aeq.cols() || luNNZ != Aeq.nonZeros()) {
            solver.analyzePattern(Aeq);
            luRows = Aeq.rows(); luCols = Aeq.cols(); luNNZ = Aeq.nonZeros();
        }
        solver.factorize(Aeq);

        Eigen::VectorXd sol_eq = solver.solve(beq);
        ok = (solver.info() == Eigen::Success);
        iters = 1;

        if (ok) {
            Eigen::VectorXd sol;
            recoverSolution(sol_eq, sol);
            errInf = (A * sol - b).lpNorm<Eigen::Infinity>();
        }
        break;
    }

    case LinearSolverOptions::Type::LDLT:
    {
        // 仅适合 SPD（或至少对称）系统
        Eigen::SimplicialLDLT<SpMat> solver;
        solver.analyzePattern(Aeq);
        solver.factorize(Aeq);

        Eigen::VectorXd sol_eq = solver.solve(beq);
        ok = (solver.info() == Eigen::Success);
        iters = 1;

        if (ok) {
            Eigen::VectorXd sol;
            recoverSolution(sol_eq, sol);
            errInf = (A * sol - b).lpNorm<Eigen::Infinity>();
        }
        break;
    }
    }

    if (iters_out) *iters_out = iters;
    if (error_out) *error_out = errInf;
    return ok;
}


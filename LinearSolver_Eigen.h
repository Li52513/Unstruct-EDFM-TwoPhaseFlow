#pragma once
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include "Solver_AssemblerCOO.h"

// 配置结构：默认 BiCGSTAB+ILUT，可切换 SparseLU/LDLT/CG
struct LinearSolverOptions {
    enum class Type { CG, BiCGSTAB, SparseLU, LDLT } type = Type::BiCGSTAB;
    int    maxIters = 2000;
    double tol = 1e-8;
    int    iluFill = 10;
    double iluDrop = 1e-4;
    bool   equil = true;   // 是否做对角均衡，建议保持 true
    int    reusePreconditioner = 0; // 0->rebuild each call; >0 reuse ILUT for subsequent solves
};

inline bool solveCOO_Eigen(const SparseSystemCOO& sys,
    std::vector<double>& x,
    const LinearSolverOptions& opt,
    int* iters_out = nullptr,
    double* error_out = nullptr)
{
    using SpMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;
    using Trip = Eigen::Triplet<double>;

    std::vector<Trip> triplets;
    triplets.reserve(sys.A.size());
    for (const auto& t : sys.A) triplets.emplace_back(t.r, t.c, t.v);

    SpMat A(sys.n, sys.n);
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd b(sys.n);
    for (int i = 0; i < sys.n; ++i)
        b[i] = (i < static_cast<int>(sys.b.size())) ? sys.b[i] : 0.0;

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(sys.n);
    if (x.size() == static_cast<size_t>(sys.n))
        for (int i = 0; i < sys.n; ++i) x0[i] = x[i];

    // ---------- 可选对角均衡 ----------
    SpMat Aeq = A;
    Eigen::VectorXd beq = b;
    Eigen::VectorXd x0eq = x0;
    Eigen::VectorXd D = Eigen::VectorXd::Ones(sys.n);

    if (opt.equil) {
        D = A.diagonal().cwiseAbs();
        for (int i = 0; i < D.size(); ++i)
            if (D[i] <= 1e-30) D[i] = 1.0;

        Eigen::VectorXd Dinvsqrt = D.cwiseInverse().cwiseSqrt();

        for (int k = 0; k < Aeq.outerSize(); ++k)
            for (SpMat::InnerIterator it(Aeq, k); it; ++it)
                it.valueRef() *= Dinvsqrt[it.row()];
        for (int k = 0; k < Aeq.outerSize(); ++k)
            for (SpMat::InnerIterator it(Aeq, k); it; ++it)
                it.valueRef() *= Dinvsqrt[it.col()];

        beq = Dinvsqrt.asDiagonal() * b;
        x0eq = Dinvsqrt.asDiagonal() * x0;
    }

    auto recoverSolution = [&](const Eigen::VectorXd& sol_eq, Eigen::VectorXd& sol_vec) {
        if (opt.equil) {
            Eigen::VectorXd Dinvsqrt = D.cwiseInverse().cwiseSqrt();
            sol_vec = Dinvsqrt.asDiagonal() * sol_eq;
        }
        else {
            sol_vec = sol_eq;
        }
        x.assign(sol_vec.data(), sol_vec.data() + sol_vec.size());
        };

    int iters = 0;
    double err = 0.0;
    bool ok = false;

    switch (opt.type) {
    case LinearSolverOptions::Type::CG: 
    {
        Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper,
            Eigen::DiagonalPreconditioner<double>> solver;
        solver.setMaxIterations(opt.maxIters);
        solver.setTolerance(opt.tol);
        solver.compute(Aeq);
        Eigen::VectorXd sol_eq = solver.solveWithGuess(beq, x0eq);
        iters = solver.iterations();
        err = solver.error();
        ok = (solver.info() == Eigen::Success);
        if (ok) { Eigen::VectorXd sol; recoverSolution(sol_eq, sol); }
        break;
    }
    case LinearSolverOptions::Type::BiCGSTAB: 
    {
        static Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double>> solver;
        static int lastN = -1;
        static double lastFill = -1.0;
        static double lastDrop = -1.0;
        static int reuseCountdown = 0;
        static bool haveFactor = false;

        bool recompute = !haveFactor
            || lastN != static_cast<int>(Aeq.rows())
            || lastFill != opt.iluFill
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
            lastFill = opt.iluFill;
            lastDrop = opt.iluDrop;
        }
        else {
            --reuseCountdown;
        }

        sol_eq = solver.solveWithGuess(beq, x0eq);
        iters = solver.iterations();
        err = solver.error();
        ok = (solver.info() == Eigen::Success);

        if (!ok && haveFactor && opt.reusePreconditioner > 0) {
            solver.compute(Aeq);
            haveFactor = (solver.info() == Eigen::Success);
            reuseCountdown = std::max(0, opt.reusePreconditioner - 1);
            lastN = static_cast<int>(Aeq.rows());
            lastFill = opt.iluFill;
            lastDrop = opt.iluDrop;
            if (haveFactor) {
                sol_eq = solver.solveWithGuess(beq, x0eq);
                iters = solver.iterations();
                err = solver.error();
                ok = (solver.info() == Eigen::Success);
            }
        }

        if (ok) { Eigen::VectorXd sol; recoverSolution(sol_eq, sol); }
        break;
    }
    case LinearSolverOptions::Type::SparseLU: 
    {
        static bool patternReady = false;
        static Eigen::SparseLU<SpMat> solver;
        if (!patternReady) { solver.analyzePattern(Aeq); patternReady = true; }
        solver.factorize(Aeq);
        Eigen::VectorXd sol_eq = solver.solve(beq);
        err = (Aeq * sol_eq - beq).lpNorm<Eigen::Infinity>();
        iters = 1;
        ok = (solver.info() == Eigen::Success);
        if (ok) { Eigen::VectorXd sol; recoverSolution(sol_eq, sol); }
        break;
    }
    case LinearSolverOptions::Type::LDLT: 
    {
        Eigen::SimplicialLDLT<SpMat> solver;
        solver.analyzePattern(Aeq);
        solver.factorize(Aeq);
        Eigen::VectorXd sol_eq = solver.solve(beq);
        err = (Aeq * sol_eq - beq).lpNorm<Eigen::Infinity>();
        iters = 1;
        ok = (solver.info() == Eigen::Success);
        if (ok) { Eigen::VectorXd sol; recoverSolution(sol_eq, sol); }
        break;
    }
    }

    if (iters_out) *iters_out = iters;
    if (error_out) *error_out = err;
    return ok;
}

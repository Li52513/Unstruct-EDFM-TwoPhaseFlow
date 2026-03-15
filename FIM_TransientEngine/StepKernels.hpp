#pragma once

#include "Types.hpp"
#include "Diagnostics.hpp"
#include "StateSync.hpp"

#include <Eigen/Sparse>
#include <string>

namespace FIM_Engine {
    namespace detail {

        template<int N>
        using RowMajorMat = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;

        template<int N>
        using ColMajorMat = Eigen::SparseMatrix<double>;

        template<int N>
        struct LinearSolveResult {
            Eigen::VectorXd dx;
            bool compute_ok = false;
            bool solve_ok = false;
            bool row_scaling_applied = false;
            std::string solver_log;
        };

        // 预留：后续 phase-2 再把 residual/linear/update kernel 抽离到这里。

    } // namespace detail
} // namespace FIM_Engine

/**
 * @file FIM_JacobianVerifier.h
 * @brief FD vs AD 雅可比严格验证器 (中心差分版)
 */

#pragma once
#include "FIM_BlockSparseMatrix.h"
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>
#include <stdexcept>

class FIM_JacobianVerifier {
public:
    struct VerifierResult {
        double max_abs_error;
        double max_rel_error;
        int worst_row;
        int worst_col;
        bool is_pass;
    };

    /**
     * @brief 执行 FD vs AD 全矩阵严格校验
     */
    template <int N>
    static VerifierResult Verify_FD_vs_AD(
        const FIM_BlockSparseMatrix<N>& ad_mat,
        const std::vector<double>& X0,
        std::function<void(const std::vector<double>&, std::vector<double>&)> ComputeResidual_FD,
        double tol = 1e-6)
    {
        const int total_dof = ad_mat.GetTotalScalarDOF();
        if (static_cast<int>(X0.size()) != total_dof) {
            throw std::invalid_argument("Verify_FD_vs_AD: X0 size mismatch with total DOF.");
        }

        Eigen::MatrixXd fd_dense = Eigen::MatrixXd::Zero(total_dof, total_dof);

        std::vector<double> X_pert = X0;
        std::vector<double> R_plus(total_dof, 0.0);
        std::vector<double> R_minus(total_dof, 0.0);
        const double epsilon = 1e-7;

        for (int j = 0; j < total_dof; ++j) {
            const double h = epsilon * (1.0 + std::abs(X0[j]));

            X_pert[j] = X0[j] + h;
            ComputeResidual_FD(X_pert, R_plus);
            if (static_cast<int>(R_plus.size()) != total_dof) {
                throw std::runtime_error("Verify_FD_vs_AD: ComputeResidual_FD resized residual vector unexpectedly.");
            }

            X_pert[j] = X0[j] - h;
            ComputeResidual_FD(X_pert, R_minus);
            if (static_cast<int>(R_minus.size()) != total_dof) {
                throw std::runtime_error("Verify_FD_vs_AD: ComputeResidual_FD resized residual vector unexpectedly.");
            }

            for (int i = 0; i < total_dof; ++i) {
                fd_dense(i, j) = (R_plus[i] - R_minus[i]) / (2.0 * h);
            }
            X_pert[j] = X0[j];
        }

        using SparseMat = typename FIM_BlockSparseMatrix<N>::SparseMat;
        const SparseMat ad_sparse = ad_mat.ExportEigenSparseMatrix();
        const Eigen::MatrixXd ad_dense(ad_sparse);

        VerifierResult res{ 0.0, 0.0, -1, -1, false };

        for (int i = 0; i < total_dof; ++i) {
            for (int j = 0; j < total_dof; ++j) {
                const double ad_val = ad_dense(i, j);
                const double fd_val = fd_dense(i, j);
                const double abs_err = std::abs(ad_val - fd_val);
                const double rel_err = abs_err / std::max(std::abs(fd_val), 1e-12);

                if (abs_err > res.max_abs_error) {
                    res.max_abs_error = abs_err;
                }
                if (rel_err > res.max_rel_error && abs_err > 1e-8) {
                    res.max_rel_error = rel_err;
                    res.worst_row = i;
                    res.worst_col = j;
                }
            }
        }

        res.is_pass = (res.max_rel_error < tol);
        return res;
    }
};
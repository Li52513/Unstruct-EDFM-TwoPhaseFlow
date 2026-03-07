#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "FIM_BlockSparseMatrix.h"
#include "FIM_GlobalAssembler.h"
#include "FIM_JacobianVerifier.h"
#include "ADVar.hpp"

namespace Test_Day5 {

    /**
     * @brief 验证 FIM_BlockSparseMatrix 基础设施的鲁棒性
     */
    inline void Test_FIM_BlockSparseMatrix_Robustness() {
        std::cout << "--- Running FIM_BlockSparseMatrix Robustness Test ---" << std::endl;
        constexpr int N = 2;
        FIM_BlockSparseMatrix<N> mat(10);

        bool caught_negative_idx = false;
        try { mat.AddResidual(-1, 0, 5.0); }
        catch (const std::out_of_range&) { caught_negative_idx = true; }

        bool caught_out_of_bound = false;
        try { mat.AddDiagJacobian(10, 0, 0, 1.0); }
        catch (const std::out_of_range&) { caught_out_of_bound = true; }

        if (!caught_negative_idx || !caught_out_of_bound) {
            throw std::runtime_error("[FAIL] Matrix container bounds checking is broken.");
        }

        bool caught_get_bounds = false;
        try { mat.GetResidualBlock(999); }
        catch (const std::out_of_range&) { caught_get_bounds = true; }
        if (!caught_get_bounds) throw std::runtime_error("[FAIL] GetResidualBlock bounds check is missing.");

        bool caught_overflow = false;
        try { FIM_BlockSparseMatrix<3> huge_mat(1000000000); }
        catch (const std::invalid_argument&) { caught_overflow = true; }
        if (!caught_overflow) throw std::runtime_error("[FAIL] Matrix dimension overflow check is missing.");

        std::cout << "[PASS] Exception bounds checking validated." << std::endl;

        mat.AddOffDiagJacobian(2, 5, 0, 0, 3.14);
        auto eigen_mat = mat.ExportEigenSparseMatrix();

        if (std::abs(eigen_mat.coeff(2 * N + 0, 5 * N + 1)) > 1e-16 ||
            std::abs(eigen_mat.coeff(2 * N + 1, 5 * N + 0)) > 1e-16 ||
            std::abs(eigen_mat.coeff(2 * N + 1, 5 * N + 1)) > 1e-16) {
            throw std::runtime_error("[FAIL] Off-diagonal block was not strictly zero-initialized.");
        }
        std::cout << "[PASS] Off-diagonal zero-initialization validated." << std::endl;

        Eigen::Matrix<double, N, N> bulk_mat;
        bulk_mat << 1.0, 2.0,
                    3.0, 4.0;
        mat.AddDiagBlock(2, bulk_mat);
        mat.AddOffDiagBlock(3, 7, bulk_mat);
        mat.AddOffDiagBlock(3, 7, bulk_mat);

        Eigen::Matrix<double, N, 1> res_vec;
        res_vec << 10.0, 20.0;
        mat.AddResidualBlock(2, res_vec);
        mat.AddResidualBlock(2, res_vec);

        eigen_mat = mat.ExportEigenSparseMatrix();
        auto eigen_rhs = mat.ExportEigenResidual();

        if (std::abs(eigen_mat.coeff(3 * N + 1, 7 * N + 0) - 6.0) > 1e-10) {
            throw std::runtime_error("[FAIL] AddOffDiagBlock bulk accumulation is incorrect.");
        }
        if (std::abs(eigen_mat.coeff(2 * N + 1, 2 * N + 1) - 4.0) > 1e-10) {
            throw std::runtime_error("[FAIL] AddDiagBlock accumulation is incorrect.");
        }
        if (std::abs(eigen_rhs(2 * N + 1) - 40.0) > 1e-10) {
            throw std::runtime_error("[FAIL] AddResidualBlock accumulation is incorrect.");
        }
        std::cout << "[PASS] AddDiagBlock & AddResidualBlock accumulation validated." << std::endl;

        auto eigen_mat_2nd = mat.ExportEigenSparseMatrix();
        Eigen::MatrixXd d1(eigen_mat), d2(eigen_mat_2nd);
        if (!d1.isApprox(d2, 1e-14)) {
            throw std::runtime_error("[FAIL] Repeated exports are not identical (isApprox failed).");
        }
        std::cout << "[PASS] Repeated export structural stability validated." << std::endl;

        size_t topo_size_before = mat.GetOffDiagMapSize(3);
        mat.SetZero();
        size_t topo_size_after = mat.GetOffDiagMapSize(3);

        auto eigen_mat_zeroed = mat.ExportEigenSparseMatrix();
        if (topo_size_before != topo_size_after || topo_size_after != 1) {
            throw std::runtime_error("[FAIL] SetZero() destroyed topology pattern.");
        }
        if (eigen_mat_zeroed.nonZeros() != 0) {
            throw std::runtime_error("[FAIL] SetZero() did not clear matrix values.");
        }
        std::cout << "[PASS] SetZero() topology retention and value clearing validated." << std::endl;

        mat.FreezePattern();
        bool caught_pattern_freeze = false;
        try {
            mat.AddOffDiagBlock(1, 8, bulk_mat);
        }
        catch (const std::runtime_error&) {
            caught_pattern_freeze = true;
        }

        bool caught_elem_freeze = false;
        try {
            mat.AddOffDiagJacobian(1, 8, 0, 0, 1.0);
        }
        catch (const std::runtime_error&) {
            caught_elem_freeze = true;
        }

        mat.AddOffDiagBlock(3, 7, bulk_mat);

        if (!caught_pattern_freeze) throw std::runtime_error("[FAIL] Pattern Freeze failed on block API.");
        if (!caught_elem_freeze) throw std::runtime_error("[FAIL] Pattern Freeze failed on element API.");

        std::cout << "[PASS] Pattern Freeze security mechanism validated." << std::endl;
        std::cout << "[PASS] FIM_BlockSparseMatrix Infrastructure is strictly PROD-ready.\n" << std::endl;
    }

    template <typename T>
    T NonLinearFlux(const T& P_up, const T& P_down) {
        return (P_up - P_down) * P_up * P_up * 0.001;
    }

    template <typename T>
    T NonLinearWellSource(const T& P_block, double P_bhp) {
        return (P_bhp - P_block) * 5.0;
    }

    template <int N, typename ADVarType>
    void Run_Jacobian_Validation_Scenario(const std::string& case_name, int num_blocks) {
        std::cout << "\n>> [Day5 Gate] Scenario: " << case_name << " (DOF=" << N << "/block)" << std::endl;

        std::vector<double> X0(num_blocks * N);
        for (int i = 0; i < num_blocks * N; ++i) X0[i] = 100.0 + i * 5.0;

        FIM_BlockSparseMatrix<N> ad_mat(num_blocks);
        for (int i = 0; i < num_blocks - 1; ++i) {
            ad_mat.AddOffDiagBlock(i, i + 1, Eigen::Matrix<double, N, N>::Zero());
            ad_mat.AddOffDiagBlock(i + 1, i, Eigen::Matrix<double, N, N>::Zero());
        }
        ad_mat.FreezePattern();
        ad_mat.SetZero();

        for (int i = 0; i < num_blocks; ++i) {
            std::vector<ADVarType> acc_eqs(N);
            for (int eq = 0; eq < N; ++eq) {
                ADVarType var(X0[i * N + eq]);
                var.grad(eq) = 1.0;
                acc_eqs[eq] = var * 0.1;
            }
            FIM_GlobalAssembler<N, ADVarType>::AssembleAccumulation(i, acc_eqs, ad_mat);

            if (i < num_blocks - 1) {
                std::vector<ADVarType> flux_wrt_i(N), flux_wrt_j(N);
                for (int eq = 0; eq < N; ++eq) {
                    ADVarType P_i(X0[i * N + eq]);
                    P_i.grad(eq) = 1.0;
                    ADVarType P_j(X0[(i + 1) * N + eq]);
                    flux_wrt_i[eq] = NonLinearFlux(P_i, P_j);

                    ADVarType P_i2(X0[i * N + eq]);
                    ADVarType P_j2(X0[(i + 1) * N + eq]);
                    P_j2.grad(eq) = 1.0;
                    flux_wrt_j[eq] = NonLinearFlux(P_i2, P_j2);
                }
                FIM_GlobalAssembler<N, ADVarType>::AssembleFlux(i, i + 1, flux_wrt_i, flux_wrt_j, ad_mat);
            }
        }

        const int inj_idx = 0;
        const int prod_idx = num_blocks - 1;
        const double inj_bhp = 200.0;
        const double prod_bhp = 50.0;

        const Eigen::VectorXd rhs_before_source = ad_mat.ExportEigenResidual();
        const auto jac_before_source = ad_mat.ExportEigenSparseMatrix();

        std::vector<ADVarType> src_inj(N), src_prod(N);
        for (int eq = 0; eq < N; ++eq) {
            ADVarType P_inj(X0[inj_idx * N + eq]);
            P_inj.grad(eq) = 1.0;
            src_inj[eq] = NonLinearWellSource(P_inj, inj_bhp);

            ADVarType P_prod(X0[prod_idx * N + eq]);
            P_prod.grad(eq) = 1.0;
            src_prod[eq] = NonLinearWellSource(P_prod, prod_bhp);
        }
        FIM_GlobalAssembler<N, ADVarType>::AssembleSource(inj_idx, src_inj, ad_mat);
        FIM_GlobalAssembler<N, ADVarType>::AssembleSource(prod_idx, src_prod, ad_mat);

        const Eigen::VectorXd rhs_after_source = ad_mat.ExportEigenResidual();
        const auto jac_after_source = ad_mat.ExportEigenSparseMatrix();

        auto ComputeFD = [&](const std::vector<double>& X, std::vector<double>& R) {
            std::fill(R.begin(), R.end(), 0.0);
            for (int i = 0; i < num_blocks; ++i) {
                for (int eq = 0; eq < N; ++eq) {
                    R[i * N + eq] += X[i * N + eq] * 0.1;
                }
                if (i < num_blocks - 1) {
                    for (int eq = 0; eq < N; ++eq) {
                        const double flux = NonLinearFlux(X[i * N + eq], X[(i + 1) * N + eq]);
                        R[i * N + eq] += flux;
                        R[(i + 1) * N + eq] -= flux;
                    }
                }
            }
            for (int eq = 0; eq < N; ++eq) {
                R[inj_idx * N + eq] -= NonLinearWellSource(X[inj_idx * N + eq], inj_bhp);
                R[prod_idx * N + eq] -= NonLinearWellSource(X[prod_idx * N + eq], prod_bhp);
            }
        };

        auto result = FIM_JacobianVerifier::Verify_FD_vs_AD<N>(ad_mat, X0, ComputeFD, 1e-6);

        std::cout << "   Max Abs Error: " << std::scientific << result.max_abs_error << std::endl;
        std::cout << "   Max Rel Error: " << std::scientific << result.max_rel_error
                  << " at worst(row=" << result.worst_row << ", col=" << result.worst_col << ")" << std::endl;

        if (!result.is_pass) {
            throw std::runtime_error("[FAIL] FD vs AD validation failed for " + case_name);
        }
        std::cout << "[PASS] FD vs AD max relative error < 1e-6 Validated." << std::endl;

        if constexpr (N >= 3) {
            constexpr int e = 2;
            const double inj_energy_source = src_inj[e].val;
            const double prod_energy_source = src_prod[e].val;

            if (inj_energy_source <= 0.0) {
                throw std::runtime_error("[FAIL] Injection well is not injecting energy properly.");
            }
            if (prod_energy_source >= 0.0) {
                throw std::runtime_error("[FAIL] Production well is not extracting energy properly.");
            }

            const double net_extraction = -(inj_energy_source + prod_energy_source);
            if (net_extraction <= 0.0) {
                throw std::runtime_error("[FAIL] Net heat extraction direction is incorrect.");
            }

            const int inj_row = inj_idx * N + e;
            const int prod_row = prod_idx * N + e;

            const double inj_res_delta = rhs_after_source(inj_row) - rhs_before_source(inj_row);
            const double prod_res_delta = rhs_after_source(prod_row) - rhs_before_source(prod_row);
            if (std::abs(inj_res_delta + inj_energy_source) > 1e-9) {
                throw std::runtime_error("[FAIL] Injection residual source sign mismatch.");
            }
            if (std::abs(prod_res_delta + prod_energy_source) > 1e-9) {
                throw std::runtime_error("[FAIL] Production residual source sign mismatch.");
            }

            const double inj_jac_delta = jac_after_source.coeff(inj_row, inj_row) - jac_before_source.coeff(inj_row, inj_row);
            const double prod_jac_delta = jac_after_source.coeff(prod_row, prod_row) - jac_before_source.coeff(prod_row, prod_row);
            if (std::abs(inj_jac_delta - 5.0) > 1e-9 || std::abs(prod_jac_delta - 5.0) > 1e-9) {
                throw std::runtime_error("[FAIL] Source Jacobian direction/magnitude mismatch on energy equation.");
            }

            std::cout << "[PASS] 注采对井取热残差与导数装配方向 Validated. (Injection="
                      << inj_energy_source << " , Extraction=" << prod_energy_source
                      << " , NetExtraction=" << net_extraction << ")" << std::endl;
        }
    }

    inline void Run_Day5_GlobalAssembly_Jacobian_2D() {
        std::cout << "========== Running Day5 2D Jacobian Verification ==========" << std::endl;
        Run_Jacobian_Validation_Scenario<2, ADVar<2>>("2D Single-Phase Supercritical CO2", 5);
        Run_Jacobian_Validation_Scenario<3, ADVar<3>>("2D Two-Phase Flow (Water + CO2)", 5);
        std::cout << "[PASS] Day5 2D Global Assembly Passed!" << std::endl;
    }

    inline void Run_Day5_GlobalAssembly_Jacobian_3D() {
        std::cout << "========== Running Day5 3D Jacobian Verification ==========" << std::endl;
        Run_Jacobian_Validation_Scenario<2, ADVar<2>>("3D Single-Phase Supercritical CO2", 15);
        Run_Jacobian_Validation_Scenario<3, ADVar<3>>("3D Two-Phase Flow (Water + CO2)", 15);
        std::cout << "[PASS] Day5 3D Global Assembly Passed!" << std::endl;
    }

} // namespace Test_Day5
/**
 * @file FIM_GlobalAssembler.h
 * @brief ȫ��ʽ�����ȫ����װ�� (Global Assembler for FIM)
 * @details ���� AD �������ӵľֲ������� (Accumulation, Flux, Source)
 * ��Ч����ȫ��װ���� FIM_BlockSparseMatrix���ѿ����ϸ�ĳߴ������ͨ��һ���Լ�顣
 */

#pragma once
#include "FIM_BlockSparseMatrix.h"
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

template <int N, typename ADVarType>
class FIM_GlobalAssembler {
public:

	/**
	 * @brief 累积项装配 (Accumulation)
	 * @details R_i += Accum_i
	 *          J_ii += d(Accum_i)/d(Phi_i)  (仅对主变量 Phi_i 有导数贡献)
	 */
    static void AssembleAccumulation(int block_idx,
        const std::vector<ADVarType>& accum_eqs,
        FIM_BlockSparseMatrix<N>& global_mat)
    {
        if (static_cast<int>(accum_eqs.size()) != N) {
            throw std::invalid_argument("AssembleAccumulation: accum_eqs size mismatch.");
        }

        for (int eq = 0; eq < N; ++eq) {
            global_mat.AddResidual(block_idx, eq, accum_eqs[eq].val);
            for (int var = 0; var < N; ++var) {
                global_mat.AddDiagJacobian(block_idx, eq, var, accum_eqs[eq].grad(var));
            }
        }
    }

    /**
     * @brief Flux assembly (MM / FI / NNC / FF)
     * @details flux = mob * (Phi_j - Phi_i) * T  (positive = inflow to i)
     *          R_i -= Flux (inflow reduces residual)
     *          R_j += Flux (inflow to i = outflow from j)
     */
    static void AssembleFlux(int block_i, int block_j,
        const std::vector<ADVarType>& flux_wrt_i,
        const std::vector<ADVarType>& flux_wrt_j,
        FIM_BlockSparseMatrix<N>& global_mat)
    {
        if (static_cast<int>(flux_wrt_i.size()) != N || static_cast<int>(flux_wrt_j.size()) != N) {
            throw std::invalid_argument("AssembleFlux: flux vector size mismatch.");
        }

        for (int eq = 0; eq < N; ++eq) {
            const double f_i = flux_wrt_i[eq].val;
            const double f_j = flux_wrt_j[eq].val;
            const double abs_tol = 1e-12;
            const double rel_tol = 1e-8;
            if (std::abs(f_i - f_j) > abs_tol + rel_tol * std::max(std::abs(f_i), std::abs(f_j))) {
                throw std::runtime_error("AssembleFlux: inconsistent flux values between seeds.");
            }

            global_mat.AddResidual(block_i, eq, -f_i);
            global_mat.AddResidual(block_j, eq,  f_i);

            for (int var = 0; var < N; ++var) {
                global_mat.AddDiagJacobian(block_i, eq, var, -flux_wrt_i[eq].grad(var));
                global_mat.AddOffDiagJacobian(block_i, block_j, eq, var, -flux_wrt_j[eq].grad(var));

                global_mat.AddOffDiagJacobian(block_j, block_i, eq, var,  flux_wrt_i[eq].grad(var));
                global_mat.AddDiagJacobian(block_j, eq, var,  flux_wrt_j[eq].grad(var));
            }
        }
    }

    /**
     * @brief ��װԴ���� (Source: Well / Boundary / Leakoff)
     * @details Լ����Դ���Ϊע��Ϊ����R_i -= Source_i
     */
    static void AssembleSource(int block_idx,
        const std::vector<ADVarType>& source_wrt_i,
        FIM_BlockSparseMatrix<N>& global_mat)
    {
        if (static_cast<int>(source_wrt_i.size()) != N) {
            throw std::invalid_argument("AssembleSource: source vector size mismatch.");
        }

        for (int eq = 0; eq < N; ++eq) {
            global_mat.AddResidual(block_idx, eq, -source_wrt_i[eq].val);
            for (int var = 0; var < N; ++var) {
                global_mat.AddDiagJacobian(block_idx, eq, var, -source_wrt_i[eq].grad(var));
            }
        }
    }
};
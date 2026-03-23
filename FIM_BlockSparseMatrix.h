#pragma once
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector> 
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <string>
#include <limits> 

template <int N>
class FIM_BlockSparseMatrix {
    static_assert(N > 0, "FIM_BlockSparseMatrix requires N > 0.");

public:
    
    using BlockMat = Eigen::Matrix<double, N, N>;
    using BlockVec = Eigen::Matrix<double, N, 1>;
    using OffDiagMap = std::unordered_map<
        int, BlockMat, std::hash<int>, std::equal_to<int>,
        Eigen::aligned_allocator<std::pair<const int, BlockMat>>
    >;
    using SparseMat = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;

private:
    int total_blocks_;                                                            
    bool is_pattern_frozen_;                                                      

    std::vector<BlockMat, Eigen::aligned_allocator<BlockMat>> diag_blocks_;       
    std::vector<OffDiagMap> off_diag_blocks_;                                     
    std::vector<BlockVec, Eigen::aligned_allocator<BlockVec>> residual_;          

    
    mutable bool     csr_cache_built_ = false;
    mutable SparseMat frozen_mat_;
    mutable std::vector<std::unordered_map<int, int>> col_block_rank_;

    
    inline void CheckBlockIndex(int idx, const char* context) const {
        if (idx < 0 || idx >= total_blocks_) {
            throw std::out_of_range(std::string(context) + ": Block index out of range: " + std::to_string(idx));
        }
    }

    inline void CheckDOFIndex(int dof, const char* context) const {
        if (dof < 0 || dof >= N) {
            throw std::out_of_range(std::string(context) + ": DOF index out of range: " + std::to_string(dof));
        }
    }

    void BuildCSRCacheInternal_() const {
        // Step 1: build sorted column-block rank map for each block-row
        col_block_rank_.resize(total_blocks_);
        for (int i = 0; i < total_blocks_; ++i) {
            std::vector<int> cols;
            cols.push_back(i);
            for (const auto& kv : off_diag_blocks_[i]) cols.push_back(kv.first);
            std::sort(cols.begin(), cols.end());
            col_block_rank_[i].clear();
            for (int rank = 0; rank < (int)cols.size(); ++rank)
                col_block_rank_[i][cols[rank]] = rank;
        }

        // Step 2: build frozen_mat_ with placeholder zeros (one-time Triplet path)
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(static_cast<size_t>(total_blocks_) * N * N * 7);
        for (int i = 0; i < total_blocks_; ++i) {
            int ro = i * N;
            for (int r = 0; r < N; ++r)
                for (int c = 0; c < N; ++c)
                    triplets.emplace_back(ro + r, ro + c, 0.0);
            for (const auto& kv : off_diag_blocks_[i]) {
                int co = kv.first * N;
                for (int r = 0; r < N; ++r)
                    for (int c = 0; c < N; ++c)
                        triplets.emplace_back(ro + r, co + c, 0.0);
            }
        }
        const int M = total_blocks_ * N;
        frozen_mat_.resize(M, M);
        frozen_mat_.setFromTriplets(triplets.begin(), triplets.end());
        frozen_mat_.makeCompressed();

        // Step 3: validate rank offsets match actual innerIndexPtr
        {
            const int*   outer = frozen_mat_.outerIndexPtr();
            const int*   inner = frozen_mat_.innerIndexPtr();
            for (int i = 0; i < total_blocks_; ++i) {
                for (const auto& kv : col_block_rank_[i]) {
                    int j    = kv.first;
                    int rank = kv.second;
                    for (int r = 0; r < N; ++r) {
                        int expected_col = j * N;
                        int actual_col   = inner[outer[i * N + r] + rank * N];
                        assert(actual_col == expected_col &&
                               "CSR cache: rank*N offset does not match innerIndexPtr");
                    }
                }
            }
        }

        csr_cache_built_ = true;
    }

    void UpdateCSRValuesInternal_() const {
        double*    val_ptr = frozen_mat_.valuePtr();
        const int* outer   = frozen_mat_.outerIndexPtr();

        for (int i = 0; i < total_blocks_; ++i) {
            int ro = i * N;

            // diagonal block (i, i)
            const int rank_diag = col_block_rank_[i].at(i);
            for (int r = 0; r < N; ++r) {
                double* dst = val_ptr + outer[ro + r] + rank_diag * N;
                for (int c = 0; c < N; ++c)
                    dst[c] = diag_blocks_[i](r, c);
            }

            // off-diagonal blocks (i, j)
            for (const auto& kv : off_diag_blocks_[i]) {
                const int j    = kv.first;
                const int rank = col_block_rank_[i].at(j);
                for (int r = 0; r < N; ++r) {
                    double* dst = val_ptr + outer[ro + r] + rank * N;
                    for (int c = 0; c < N; ++c)
                        dst[c] = kv.second(r, c);
                }
            }
        }
    }

public:
    explicit FIM_BlockSparseMatrix(int total_blocks) : total_blocks_(total_blocks), is_pattern_frozen_(false) {
        if (total_blocks_ < 0) {
            throw std::invalid_argument("Total blocks cannot be negative.");
        }

        // [�޸� High] ά�������������ֹ������װ�� CSR ��ʽ�±����
        if (total_blocks_ > std::numeric_limits<int>::max() / N) {
            throw std::invalid_argument("total_blocks * N overflows Eigen index range.");
        }

        diag_blocks_.assign(total_blocks_, BlockMat::Zero());
        off_diag_blocks_.resize(total_blocks_);
        residual_.assign(total_blocks_, BlockVec::Zero());
    }

    /**
     * @brief ����ϵͳ���� (Pattern Freeze)�����ú󣬽�ֹ��ʽ�����µķǶԽǿ����ӡ�
     */
    void FreezePattern() { is_pattern_frozen_ = true; }

    /**
     * @brief ���ϵͳ���˶����������ʧЧ CSR ������
     */
    void UnfreezePattern() {
        is_pattern_frozen_ = false;
        csr_cache_built_   = false;
        frozen_mat_.resize(0, 0);
        col_block_rank_.clear();
    }


    void SetZero() {
        for (int i = 0; i < total_blocks_; ++i) {
            diag_blocks_[i].setZero();
            for (auto& pair : off_diag_blocks_[i]) {
                pair.second.setZero();
            }
            residual_[i].setZero();
        }
    }


    inline void AddResidual(int block_idx, int dof_idx, double val) {
        CheckBlockIndex(block_idx, "AddResidual");
        CheckDOFIndex(dof_idx, "AddResidual");
        residual_[block_idx](dof_idx, 0) += val;
    }

    inline void AddDiagJacobian(int block_idx, int row_dof, int col_dof, double val) {
        CheckBlockIndex(block_idx, "AddDiagJacobian");
        CheckDOFIndex(row_dof, "AddDiagJacobian row");
        CheckDOFIndex(col_dof, "AddDiagJacobian col");
        diag_blocks_[block_idx](row_dof, col_dof) += val;
    }

    inline void AddOffDiagJacobian(int row_block, int col_block, int row_dof, int col_dof, double val) {
        CheckBlockIndex(row_block, "AddOffDiagJacobian row_block");
        CheckBlockIndex(col_block, "AddOffDiagJacobian col_block");
        CheckDOFIndex(row_dof, "AddOffDiagJacobian row_dof");
        CheckDOFIndex(col_dof, "AddOffDiagJacobian col_dof");

        if (row_block == col_block) {
            diag_blocks_[row_block](row_dof, col_dof) += val;
            return;
        }

        auto& rowMap = off_diag_blocks_[row_block];
        auto it = rowMap.find(col_block);
        if (it == rowMap.end()) {
            if (is_pattern_frozen_) {
                throw std::runtime_error("Pattern is frozen. Illegal new connection: " +
                    std::to_string(row_block) + " -> " + std::to_string(col_block));
            }
            it = rowMap.emplace(col_block, BlockMat::Zero()).first;
        }
        it->second(row_dof, col_dof) += val;
    }


    inline void AddResidualBlock(int block_idx, const BlockVec& vec) {
        CheckBlockIndex(block_idx, "AddResidualBlock");
        residual_[block_idx] += vec;
    }

    inline void AddDiagBlock(int block_idx, const BlockMat& mat) {
        CheckBlockIndex(block_idx, "AddDiagBlock");
        diag_blocks_[block_idx] += mat;
    }

    inline void AddOffDiagBlock(int row_block, int col_block, const BlockMat& mat) {
        CheckBlockIndex(row_block, "AddOffDiagBlock row_block");
        CheckBlockIndex(col_block, "AddOffDiagBlock col_block");

        if (row_block == col_block) {
            diag_blocks_[row_block] += mat;
            return;
        }

        auto& rowMap = off_diag_blocks_[row_block];
        auto it = rowMap.find(col_block);
        if (it == rowMap.end()) {
            if (is_pattern_frozen_) {
                throw std::runtime_error("Pattern is frozen. Illegal new block connection: " +
                    std::to_string(row_block) + " -> " + std::to_string(col_block));
            }
            it = rowMap.emplace(col_block, BlockMat::Zero()).first;
        }
        it->second += mat;
    }

    inline const BlockVec& GetResidualBlock(int block_idx) const {
        CheckBlockIndex(block_idx, "GetResidualBlock");
        return residual_[block_idx];
    }

    inline int GetTotalScalarDOF() const {
        return total_blocks_ * N;
    }

    inline size_t GetOffDiagMapSize(int block_idx) const {
        CheckBlockIndex(block_idx, "GetOffDiagMapSize");
        return off_diag_blocks_[block_idx].size();
    }

    SparseMat ExportEigenSparseMatrix() const {
        std::vector<Eigen::Triplet<double>> triplets;
        const size_t reserve_nnz =
            static_cast<size_t>(total_blocks_) * static_cast<size_t>(N) *
            static_cast<size_t>(N) * static_cast<size_t>(7);
        triplets.reserve(reserve_nnz);

        for (int i = 0; i < total_blocks_; ++i) {
            int row_offset = i * N;

            for (int r = 0; r < N; ++r) {
                for (int c = 0; c < N; ++c) {
                    triplets.emplace_back(row_offset + r, row_offset + c, diag_blocks_[i](r, c));
                }
            }

            for (const auto& neighbor : off_diag_blocks_[i]) {
                int col_offset = neighbor.first * N;
                for (int r = 0; r < N; ++r) {
                    for (int c = 0; c < N; ++c) {
                        triplets.emplace_back(row_offset + r, col_offset + c, neighbor.second(r, c));
                    }
                }
            }
        }

        std::sort(triplets.begin(), triplets.end(),
            [](const Eigen::Triplet<double>& a, const Eigen::Triplet<double>& b) {
                if (a.row() != b.row()) return a.row() < b.row();
                return a.col() < b.col();
            });

        SparseMat mat(total_blocks_ * N, total_blocks_ * N);
        mat.setFromTriplets(triplets.begin(), triplets.end());
        mat.makeCompressed();
        return mat;
    }

    const SparseMat& GetFrozenMatrix() const {
        if (!is_pattern_frozen_)
            throw std::logic_error("GetFrozenMatrix() requires FreezePattern() to have been called.");
        if (!csr_cache_built_) BuildCSRCacheInternal_();
        UpdateCSRValuesInternal_();
        return frozen_mat_;
    }

    Eigen::VectorXd ExportEigenResidual() const {
        Eigen::VectorXd rhs(total_blocks_ * N);
        for (int i = 0; i < total_blocks_; ++i) {
            rhs.segment<N>(i * N) = residual_[i];
        }
        return rhs;
    }
};
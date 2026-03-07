/**
 * @file FIM_BlockSparseMatrix.h
 * @brief 全隐式求解器块稀疏矩阵与右端项容器 (Block Sparse Matrix for FIM)
 * @details
 * 专为 N 变量系统设计，已硬化：
 * 1. 强制非法索引(如 -1)与越界检查。
 * 2. 保证非对角块按零初始化。
 * 3. 保证导出的 Eigen 矩阵结构确定且排序稳定 (严格 CSR 格式)。
 * 4. 支持 Pattern Freeze (拓扑冻结)，防止非法隐蔽的新连接生成。
 * 5. 严格处理 STL 容器存 Eigen 固定维矩阵的内存对齐问题。
 * 6. 严防极大规模网格组装导致总自由度溢出 32 位整型范围。
 */

#pragma once
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector> // [修复 High] 必须包含此头文件以支持 Eigen STL 对齐
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <string>
#include <limits> // [修复 High] 用于溢出检查

template <int N>
class FIM_BlockSparseMatrix {
    static_assert(N > 0, "FIM_BlockSparseMatrix requires N > 0.");

public:
    // [修复 High] 声明对齐安全的类型别名
    using BlockMat = Eigen::Matrix<double, N, N>;
    using BlockVec = Eigen::Matrix<double, N, 1>;
    using OffDiagMap = std::unordered_map<
        int, BlockMat, std::hash<int>, std::equal_to<int>,
        Eigen::aligned_allocator<std::pair<const int, BlockMat>>
    >;
    using SparseMat = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;

private:
    int total_blocks_;                                                            ///< 系统总块数
    bool is_pattern_frozen_;                                                      ///< 是否冻结拓扑

    // [修复 High] 针对 Eigen 固定大小结构强制使用对齐分配器
    std::vector<BlockMat, Eigen::aligned_allocator<BlockMat>> diag_blocks_;       ///< 主对角块数组
    std::vector<OffDiagMap> off_diag_blocks_;                                     ///< 非对角块数组
    std::vector<BlockVec, Eigen::aligned_allocator<BlockVec>> residual_;          ///< 右端项 (残差) 向量

    // --- 内部安全断言 ---
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

public:
    /**
     * @brief 构造函数，初始化矩阵尺寸
     * @param total_blocks 系统总自由度块数
     */
    explicit FIM_BlockSparseMatrix(int total_blocks) : total_blocks_(total_blocks), is_pattern_frozen_(false) {
        if (total_blocks_ < 0) {
            throw std::invalid_argument("Total blocks cannot be negative.");
        }

        // [修复 High] 维度溢出保护，防止后续组装与 CSR 格式下标溢出
        if (total_blocks_ > std::numeric_limits<int>::max() / N) {
            throw std::invalid_argument("total_blocks * N overflows Eigen index range.");
        }

        diag_blocks_.assign(total_blocks_, BlockMat::Zero());
        off_diag_blocks_.resize(total_blocks_);
        residual_.assign(total_blocks_, BlockVec::Zero());
    }

    /**
     * @brief 冻结系统拓扑 (Pattern Freeze)。调用后，禁止隐式插入新的非对角块连接。
     */
    void FreezePattern() { is_pattern_frozen_ = true; }

    /**
     * @brief 解除系统拓扑冻结。
     */
    void UnfreezePattern() { is_pattern_frozen_ = false; }

    /**
     * @brief 清空矩阵和残差向量（保留拓扑内存以备下一牛顿步使用）
     */
    void SetZero() {
        for (int i = 0; i < total_blocks_; ++i) {
            diag_blocks_[i].setZero();
            for (auto& pair : off_diag_blocks_[i]) {
                pair.second.setZero();
            }
            residual_[i].setZero();
        }
    }

    // --- 逐元素累加接口 ---

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

    // --- 整块累加接口 ---

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

    // --- Getter 接口 ---

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

    // --- 导出接口 ---

    /**
     * @brief 导出为 Eigen::SparseMatrix (CSR 格式)，供线性求解器使用
     */
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
                    double val = diag_blocks_[i](r, c);
                    if (std::abs(val) > 1e-16) {
                        triplets.emplace_back(row_offset + r, row_offset + c, val);
                    }
                }
            }

            for (const auto& neighbor : off_diag_blocks_[i]) {
                int col_offset = neighbor.first * N;
                for (int r = 0; r < N; ++r) {
                    for (int c = 0; c < N; ++c) {
                        double val = neighbor.second(r, c);
                        if (std::abs(val) > 1e-16) {
                            triplets.emplace_back(row_offset + r, col_offset + c, val);
                        }
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

    Eigen::VectorXd ExportEigenResidual() const {
        Eigen::VectorXd rhs(total_blocks_ * N);
        for (int i = 0; i < total_blocks_; ++i) {
            rhs.segment<N>(i * N) = residual_[i];
        }
        return rhs;
    }
};
#pragma once

#include "Types.hpp"
#include "Diagnostics.hpp"
#include "StateSync.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/eigen.hpp>
#include <amgcl/preconditioner/cpr.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>

#include <cstdint>
#include <memory>
#include <tuple>
#include <string>
#include <vector>


namespace FIM_Engine {
    namespace detail {

        template<int N>
        using RowMajorMat = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;

        template<int N>
        using ColMajorMat = Eigen::SparseMatrix<double>;

        template<int N>
        using AMGCLBackend = amgcl::backend::builtin<double>;

        template<int N>
        using AMGCLSolver = amgcl::make_solver<
            amgcl::amg<
            AMGCLBackend<N>,
            amgcl::coarsening::smoothed_aggregation,
            amgcl::relaxation::spai0
            >,
            amgcl::solver::bicgstab<AMGCLBackend<N>>
        >;

        template<int N>
        using AMGCLCPRSolver = amgcl::make_solver<
            amgcl::preconditioner::cpr<
                amgcl::amg<
                    AMGCLBackend<N>,
                    amgcl::coarsening::smoothed_aggregation,
                    amgcl::relaxation::spai0
                >,
                amgcl::relaxation::as_preconditioner<
                    AMGCLBackend<N>,
                    amgcl::relaxation::ilu0
                >
            >,
            amgcl::solver::bicgstab<AMGCLBackend<N>>
        >;

        template<int N>
        struct LinearSolverCache {
            Eigen::SparseLU<ColMajorMat<N>> sparse_lu_solver;
            bool sparse_lu_pattern_ready = false;
            bool sparse_lu_initialized = false;
            int sparse_lu_rows = -1;
            int sparse_lu_cols = -1;
            Eigen::Index sparse_lu_nnz = -1;
            std::uint64_t sparse_lu_pattern_hash = 0;

            Eigen::BiCGSTAB<RowMajorMat<N>, Eigen::IncompleteLUT<double>> bicgstab_solver;

            typename AMGCLSolver<N>::params amgcl_prm{};
            std::unique_ptr<AMGCLSolver<N>> amgcl_solver;
            bool configured = false;

            // ── 雷区27修复：工作矩阵/向量缓存（避免每次 Newton 迭代堆分配）──
            RowMajorMat<N>      solver_A_work;       // rows()==0 → not yet built
            Eigen::VectorXd     solver_b_work;
            std::vector<double> amgcl_rhs_cache;
            std::vector<double> amgcl_dx_cache;

            // ── 雷区28修复：CPR-AMG 求解器 ──
            typename AMGCLCPRSolver<N>::params amgcl_cpr_prm{};
            std::unique_ptr<AMGCLCPRSolver<N>> amgcl_cpr_solver; // nullptr → not yet built

            void Configure(const TransientSolverParams& params) {
                bicgstab_solver.preconditioner().setDroptol(params.bicgstab_droptol);
                amgcl_prm.solver.tol = params.amgcl_tol;
                amgcl_prm.solver.maxiter = params.amgcl_maxiter;
                amgcl_prm.precond.allow_rebuild = true;
                // CPR 参数初始化
                amgcl_cpr_prm.precond.block_size = N;
                amgcl_cpr_prm.solver.tol         = params.amgcl_cpr_tol;
                amgcl_cpr_prm.solver.maxiter      = params.amgcl_cpr_maxiter;
                configured = true;
            }
        };

        template<int N>
        struct LinearSolveResult {
            Eigen::VectorXd dx;
            bool compute_ok = false;
            bool solve_ok = false;
            bool row_scaling_applied = false;
            std::string solver_log;
        };

        template<typename SparseMatType>
        inline std::uint64_t HashSparsePattern(const SparseMatType& mat);

        template<int N, typename SparseMatType>
        LinearSolveResult<N> SolveLinearSystem(
            const SparseMatType& A,
            const Eigen::VectorXd& b,
            const std::vector<EqContrib>& eq_contribs,
            int totalEq,
            const TransientSolverParams& params,
            LinearSolverCache<N>& cache);

        template<int N, typename MeshMgrType>
        bool ApplyTrialUpdate(
            const FIM_StateMap<N>& base_state,
            double alpha_try,
            const Eigen::VectorXd& dx,
            const MeshMgrType& mgr,
            int totalBlocks,
            const TransientSolverParams& params,
            double p_floor,
            double p_ceil,
            double t_floor,
            double t_ceil,
            FIM_StateMap<N>& out_state,
            int& limiter_added_local,
            double& rel_update_inf,
            int num_well_blocks = 0);

    } // namespace detail
} // namespace FIM_Engine

#include "StepKernels_impl.hpp"

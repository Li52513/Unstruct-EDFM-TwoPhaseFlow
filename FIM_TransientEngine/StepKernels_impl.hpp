#pragma once

namespace FIM_Engine {
    namespace detail {

        template<typename SparseMatType>
        inline std::uint64_t HashSparsePattern(const SparseMatType& mat) {
            constexpr std::uint64_t kOffset = 1469598103934665603ull;
            constexpr std::uint64_t kPrime = 1099511628211ull;
            std::uint64_t h = kOffset;
            const auto* outer = mat.outerIndexPtr();
            const auto* inner = mat.innerIndexPtr();
            for (int i = 0; i <= mat.outerSize(); ++i) {
                h ^= static_cast<std::uint64_t>(static_cast<unsigned int>(outer[i]));
                h *= kPrime;
            }
            const Eigen::Index nnz = mat.nonZeros();
            for (Eigen::Index k = 0; k < nnz; ++k) {
                h ^= static_cast<std::uint64_t>(static_cast<unsigned int>(inner[k]));
                h *= kPrime;
            }
            return h;
        }

        template<int N, typename SparseMatType>
        inline LinearSolveResult<N> SolveLinearSystem(
            const SparseMatType& A,
            const Eigen::VectorXd& b,
            const std::vector<EqContrib>& eq_contribs,
            int totalEq,
            const TransientSolverParams& params,
            LinearSolverCache<N>& cache) {

            LinearSolveResult<N> out;
            if (!cache.configured) {
                cache.Configure(params);
            }

            // ── 雷区27修复：首次调用完整拷贝（含结构），后续仅 memcpy values ──
            if (cache.solver_A_work.rows() == 0 || cache.solver_A_work.rows() != A.rows()) {
                cache.solver_A_work = A;
                cache.solver_b_work.resize(b.size());
            } else {
                // 稀疏模式已冻结（Issue #11）：只需拷贝值数组，O(nnz) 无堆分配
                std::copy(A.valuePtr(),
                          A.valuePtr() + A.nonZeros(),
                          cache.solver_A_work.valuePtr());
            }
            cache.solver_b_work = b;
            auto& A_work = cache.solver_A_work;
            auto& b_work = cache.solver_b_work;

            if (params.enable_row_scaling) {
                Eigen::VectorXd row_scale(totalEq);
                for (int r = 0; r < totalEq; ++r) {
                    const double diag_acc_abs = (r >= 0 && r < static_cast<int>(eq_contribs.size()))
                        ? std::abs(eq_contribs[r].D_acc) : 0.0;
                    const double diag_abs = std::abs(A_work.coeff(r, r));
                    const double denom = std::max({ diag_acc_abs, diag_abs, params.row_scale_floor });
                    double s = 1.0 / std::max(denom, 1e-30);
                    s = std::max(params.row_scale_min, std::min(params.row_scale_max, s));
                    row_scale[r] = s;
                }

                for (int r = 0; r < A_work.rows(); ++r) {
                    const double s = row_scale[r];
                    for (typename RowMajorMat<N>::InnerIterator it(A_work, r); it; ++it) {
                        it.valueRef() *= s;
                    }
                }
                b_work.array() *= row_scale.array();
                out.row_scaling_applied = true;
            }

            A_work.makeCompressed();

            // ── 共用辅助 lambda：AMGCL 向量准备（消除临时 Eigen 堆分配）──
            auto prepare_amgcl_rhs = [&](int rows) {
                if ((int)cache.amgcl_rhs_cache.size() != rows) {
                    cache.amgcl_rhs_cache.resize(rows);
                    cache.amgcl_dx_cache.resize(rows);
                }
                for (int i = 0; i < rows; ++i) cache.amgcl_rhs_cache[i] = -b_work[i];
                std::fill(cache.amgcl_dx_cache.begin(), cache.amgcl_dx_cache.end(), 0.0);
            };

            // ── 共用辅助 lambda：迭代求解失败时 fallback 到 SparseLU ──
            auto sparselu_fallback = [&](bool should_fallback) -> std::string {
                if (!out.solve_ok && should_fallback) {
                    ColMajorMat<N> A_lu = A_work;
                    A_lu.makeCompressed();
                    cache.sparse_lu_solver.compute(A_lu);
                    out.compute_ok = (cache.sparse_lu_solver.info() == Eigen::Success);
                    if (out.compute_ok) {
                        out.dx = cache.sparse_lu_solver.solve(-b_work);
                        out.solve_ok = (cache.sparse_lu_solver.info() == Eigen::Success);
                    }
                    return "SparseLU";
                }
                return "none";
            };

            if (params.lin_solver == LinearSolverType::AMGCL_CPR) {
                const int rows = static_cast<int>(A_work.rows());

                // 首次建立 CPR 层次，后续用 partial_update 仅重建值
                if (!cache.amgcl_cpr_solver) {
                    cache.amgcl_cpr_solver = std::make_unique<AMGCLCPRSolver<N>>(
                        A_work, cache.amgcl_cpr_prm);
                } else {
                    cache.amgcl_cpr_solver->precond().partial_update(A_work);
                }

                prepare_amgcl_rhs(rows);
                std::size_t iters = 0;
                double error = 0.0;
                std::tie(iters, error) = (*cache.amgcl_cpr_solver)(
                    cache.amgcl_rhs_cache, cache.amgcl_dx_cache);

                out.dx = Eigen::VectorXd::Map(cache.amgcl_dx_cache.data(), rows);
                out.compute_ok = true;
                out.solve_ok = std::isfinite(error) && (error <= params.amgcl_cpr_tol);
                const std::string fallback_str = sparselu_fallback(params.amgcl_cpr_use_fallback_sparselu);

                out.solver_log = "solver=AMGCL_CPR"
                    " compute_ok=" + std::string(out.compute_ok ? "true" : "false") +
                    " solve_ok="   + std::string(out.solve_ok   ? "true" : "false") +
                    " iters="      + std::to_string(iters) +
                    " error="      + std::to_string(error) +
                    " fallback="   + fallback_str +
                    " scaled="     + std::string(out.row_scaling_applied ? "true" : "false");
                return out;
            }

            if (params.lin_solver == LinearSolverType::AMGCL) {
                const int rows = static_cast<int>(A_work.rows());

                if (!cache.amgcl_solver) {
                    cache.amgcl_solver = std::make_unique<AMGCLSolver<N>>(A_work, cache.amgcl_prm);
                } else {
                    cache.amgcl_solver->precond().rebuild(A_work);
                }

                prepare_amgcl_rhs(rows);
                std::size_t iters = 0;
                double error = 0.0;
                std::tie(iters, error) = (*cache.amgcl_solver)(cache.amgcl_rhs_cache, cache.amgcl_dx_cache);

                out.dx = Eigen::VectorXd::Map(cache.amgcl_dx_cache.data(), rows);
                out.compute_ok = true;
                out.solve_ok = std::isfinite(error) && (error <= params.amgcl_tol);
                const std::string fallback_str = sparselu_fallback(params.amgcl_use_fallback_sparselu);

                out.solver_log = "solver=AMGCL compute_ok=" + std::string(out.compute_ok ? "true" : "false") +
                    " solve_ok=" + std::string(out.solve_ok ? "true" : "false") +
                    " iters=" + std::to_string(iters) +
                    " error=" + std::to_string(error) +
                    " fallback=" + fallback_str +
                    " scaled=" + std::string(out.row_scaling_applied ? "true" : "false");
                return out;
            }

            if (params.lin_solver == LinearSolverType::SparseLU) {
                ColMajorMat<N> A_lu = A_work;
                A_lu.makeCompressed();

                const std::uint64_t pattern_hash = HashSparsePattern(A_lu);
                const bool shape_changed =
                    (!cache.sparse_lu_pattern_ready) ||
                    (cache.sparse_lu_rows != A_lu.rows()) ||
                    (cache.sparse_lu_cols != A_lu.cols()) ||
                    (cache.sparse_lu_nnz != A_lu.nonZeros()) ||
                    (cache.sparse_lu_pattern_hash != pattern_hash);
                const bool pattern_reused = cache.sparse_lu_initialized && !shape_changed;

                if (!cache.sparse_lu_initialized || shape_changed) {
                    cache.sparse_lu_solver.compute(A_lu);
                    out.compute_ok = (cache.sparse_lu_solver.info() == Eigen::Success);
                    cache.sparse_lu_pattern_ready = out.compute_ok;
                    cache.sparse_lu_initialized = out.compute_ok;
                }
                else {
                    cache.sparse_lu_solver.factorize(A_lu);
                    out.compute_ok = (cache.sparse_lu_solver.info() == Eigen::Success);
                    if (!out.compute_ok) {
                        cache.sparse_lu_solver.compute(A_lu);
                        out.compute_ok = (cache.sparse_lu_solver.info() == Eigen::Success);
                        cache.sparse_lu_pattern_ready = out.compute_ok;
                        cache.sparse_lu_initialized = out.compute_ok;
                    }
                }

                if (out.compute_ok) {
                    cache.sparse_lu_rows = static_cast<int>(A_lu.rows());
                    cache.sparse_lu_cols = static_cast<int>(A_lu.cols());
                    cache.sparse_lu_nnz = A_lu.nonZeros();
                    cache.sparse_lu_pattern_hash = pattern_hash;
                    out.dx = cache.sparse_lu_solver.solve(-b_work);
                    out.solve_ok = (cache.sparse_lu_solver.info() == Eigen::Success);
                }

                out.solver_log = "solver=SparseLU compute_ok=" + std::string(out.compute_ok ? "true" : "false") +
                    " solve_ok=" + std::string(out.solve_ok ? "true" : "false") +
                    " nnzA=" + std::to_string(A_lu.nonZeros()) +
                    " scaled=" + std::string(out.row_scaling_applied ? "true" : "false") +
                    " pattern_reused=" + std::string(pattern_reused ? "true" : "false") +
                    " info=" + std::to_string(cache.sparse_lu_solver.info());
                return out;
            }

            cache.bicgstab_solver.compute(A_work);
            out.compute_ok = (cache.bicgstab_solver.info() == Eigen::Success);
            if (out.compute_ok) {
                out.dx = cache.bicgstab_solver.solve(-b_work);
                out.solve_ok = (cache.bicgstab_solver.info() == Eigen::Success);
            }
            const std::string iters = out.compute_ok ? std::to_string(cache.bicgstab_solver.iterations()) : "NA";
            const std::string error = out.compute_ok ? std::to_string(cache.bicgstab_solver.error()) : "NA";
            out.solver_log = "solver=BiCGSTAB compute_ok=" + std::string(out.compute_ok ? "true" : "false") +
                " solve_ok=" + std::string(out.solve_ok ? "true" : "false") +
                " iters=" + iters +
                " error=" + error +
                " scaled=" + std::string(out.row_scaling_applied ? "true" : "false") +
                " info=" + std::to_string(cache.bicgstab_solver.info());

            return out;
        }

        template<int N, typename MeshMgrType>
        inline bool ApplyTrialUpdate(
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
            double& rel_update_inf) {

            out_state = base_state;
            limiter_added_local = 0;
            rel_update_inf = 0.0;

            for (int bi = 0; bi < totalBlocks; ++bi) {
                const int eqP = mgr.getEquationIndex(bi, 0);
                const int eqT = mgr.getEquationIndex(bi, (N == 3) ? 2 : 1);
                if (eqP < 0 || eqP >= dx.size() || eqT < 0 || eqT >= dx.size()) return false;

                const double dP = alpha_try * dx[eqP];
                const double p_ref = std::max(std::abs(base_state.P[bi]), 1.0e5);
                rel_update_inf = std::max(rel_update_inf, std::abs(dP) / p_ref);
                out_state.P[bi] += dP;
                if (!std::isfinite(out_state.P[bi])) return false;
                if (params.clamp_state_to_eos_bounds) {
                    if (out_state.P[bi] < p_floor) { out_state.P[bi] = p_floor; limiter_added_local++; }
                    if (out_state.P[bi] > p_ceil) { out_state.P[bi] = p_ceil; limiter_added_local++; }
                }

                const double dT = alpha_try * dx[eqT];
                const double t_ref = std::max(std::abs(base_state.T[bi]), 300.0);
                rel_update_inf = std::max(rel_update_inf, std::abs(dT) / t_ref);

                if constexpr (N == 3) {
                    const int eqSw = mgr.getEquationIndex(bi, 1);
                    if (eqSw < 0 || eqSw >= dx.size()) return false;
                    const double dSw = alpha_try * dx[eqSw];
                    rel_update_inf = std::max(rel_update_inf, std::abs(dSw));
                    out_state.Sw[bi] += dSw;
                    if (!std::isfinite(out_state.Sw[bi])) return false;
                    if (out_state.Sw[bi] < 0.0) { out_state.Sw[bi] = 0.0; limiter_added_local++; }
                    if (out_state.Sw[bi] > 1.0) { out_state.Sw[bi] = 1.0; limiter_added_local++; }
                }

                out_state.T[bi] += dT;
                if (!std::isfinite(out_state.T[bi])) return false;
                if (params.clamp_state_to_eos_bounds) {
                    if (out_state.T[bi] < t_floor) { out_state.T[bi] = t_floor; limiter_added_local++; }
                    if (out_state.T[bi] > t_ceil) { out_state.T[bi] = t_ceil; limiter_added_local++; }
                }
            }

            return true;
        }

    } // namespace detail
} // namespace FIM_Engine

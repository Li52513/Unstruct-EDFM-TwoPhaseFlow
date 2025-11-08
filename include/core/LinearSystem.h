#pragma once

#include <vector>
#include <cstddef>

namespace bndf {

struct SparseCOO {
  int n{0};
  std::vector<int> rows;
  std::vector<int> cols;
  std::vector<double> vals;

  void reset(int n_, std::size_t reserve_nz) {
    n = n_;
    rows.clear(); cols.clear(); vals.clear();
    rows.reserve(reserve_nz); cols.reserve(reserve_nz); vals.reserve(reserve_nz);
  }
  void add(int r, int c, double v) { rows.push_back(r); cols.push_back(c); vals.push_back(v); }
};

struct SparseCSR {
  int n{0};
  std::vector<int> row_ptr; // size n+1
  std::vector<int> col_idx;
  std::vector<double> val;

  void matvec(const std::vector<double>& x, std::vector<double>& y) const;
};

SparseCSR coo_to_csr_sum_duplicates(const SparseCOO& coo);

struct DiagJacobi {
  std::vector<double> invD; // size n; zero means no preconditioning
  void build(const SparseCSR& A);
  void apply(const std::vector<double>& r, std::vector<double>& z) const;
};

struct BiCGSTABResult {
  int iters{0};
  double residual{0};
  bool converged{false};
};

BiCGSTABResult bicgstab(const SparseCSR& A,
                        const std::vector<double>& b,
                        std::vector<double>& x,
                        int maxIters,
                        double tol,
                        const DiagJacobi* M = nullptr);

} // namespace bndf


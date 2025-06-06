#ifndef DGD_BENCHMARK_BENCHMARK_RESULT_H_
#define DGD_BENCHMARK_BENCHMARK_RESULT_H_

#include <algorithm>
#include <numeric>
#include <vector>

// Benchmark result.
struct BenchmarkResult {
  double solve_time;
  double prim_dual_gap;  // Relative primal-dual gaps.
  double prim_feas_err;
  double dual_feas_err;
  int iter;
  bool optimal_flag;
};

// Benchmark result array.
struct BenchmarkResultArray {
  std::vector<double> solve_time;
  std::vector<double> prim_dual_gap;
  std::vector<double> prim_feas_err;
  std::vector<double> dual_feas_err;
  std::vector<int> iter;
  std::vector<bool> optimal_flag;
  int nrun;
  int idx;

  explicit BenchmarkResultArray(int nrun);

  void AddResult(const BenchmarkResult& res);

  double AvgSolveTime() const;

  double MaxPrimalDualGap() const;

  double MaxPrimalFeasibilityError() const;

  double MaxDualFeasibilityError() const;

  double AvgIterations() const;

  int NumSuboptimalResults() const;
};

inline BenchmarkResultArray::BenchmarkResultArray(int nrun)
    : solve_time(nrun),
      prim_dual_gap(nrun),
      prim_feas_err(nrun),
      dual_feas_err(nrun),
      iter(nrun),
      optimal_flag(nrun),
      nrun(nrun),
      idx(0) {}

inline void BenchmarkResultArray::AddResult(const BenchmarkResult& res) {
  if (idx < nrun) {
    solve_time[idx] = res.solve_time;
    prim_dual_gap[idx] = res.prim_dual_gap;
    prim_feas_err[idx] = res.prim_feas_err;
    dual_feas_err[idx] = res.dual_feas_err;
    iter[idx] = res.iter;
    optimal_flag[idx] = res.optimal_flag;
    ++idx;
  }
}

inline double BenchmarkResultArray::AvgSolveTime() const {
  if (idx == 0) return 0.0;
  return std::accumulate(solve_time.begin(), solve_time.end(), 0.0) / idx;
}

inline double BenchmarkResultArray::MaxPrimalDualGap() const {
  if (idx == 0) return 0.0;
  return *std::max_element(prim_dual_gap.begin(), prim_dual_gap.end());
}

inline double BenchmarkResultArray::MaxPrimalFeasibilityError() const {
  if (idx == 0) return 0.0;
  return *std::max_element(prim_feas_err.begin(), prim_feas_err.end());
}

inline double BenchmarkResultArray::MaxDualFeasibilityError() const {
  if (idx == 0) return 0.0;
  return *std::max_element(dual_feas_err.begin(), dual_feas_err.end());
}

inline double BenchmarkResultArray::AvgIterations() const {
  if (idx == 0) return 0.0;
  return std::accumulate(iter.begin(), iter.end(), 0.0) / idx;
}

inline int BenchmarkResultArray::NumSuboptimalResults() const {
  if (idx == 0) return 0;
  return idx - std::accumulate(optimal_flag.begin(), optimal_flag.end(), 0);
}

#endif  // DGD_BENCHMARK_BENCHMARK_RESULT_H_

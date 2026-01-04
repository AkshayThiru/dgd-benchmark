#ifndef DGD_BENCHMARK_HELPERS_BENCHMARK_RESULT_H_
#define DGD_BENCHMARK_HELPERS_BENCHMARK_RESULT_H_

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

namespace bench {

// Benchmark result.
struct BenchmarkResult {
  double solve_time;
  double prim_dual_gap;  // Relative primal-dual gap.
  double prim_infeas_err;
  double dual_infeas_err;
  int iter;
  bool optimal_flag;
};

// Benchmark result array.
struct BenchmarkResultArray {
  std::vector<double> solve_times;
  std::vector<double> prim_dual_gaps;
  std::vector<double> prim_infeas_errs;
  std::vector<double> dual_infeas_errs;
  std::vector<int> iters;
  std::vector<bool> optimal_flags;
  int size;
  int idx;

  explicit BenchmarkResultArray(int size);

  void AddResult(const BenchmarkResult& res);

  double AvgSolveTime() const;

  double MaxPrimalDualGap() const;

  double MaxPrimalInfeasibilityError() const;

  double MaxDualInfeasibilityError() const;

  double AvgIterations() const;

  int NumSuboptimalResults() const;

  void PrintStatistics() const;

  // Writes benchmark output to a Feather file.
  bool SaveToFile(const std::string& filename);
};

inline BenchmarkResultArray::BenchmarkResultArray(int size)
    : solve_times(size),
      prim_dual_gaps(size),
      prim_infeas_errs(size),
      dual_infeas_errs(size),
      iters(size),
      optimal_flags(size),
      size(size),
      idx(0) {}

inline void BenchmarkResultArray::AddResult(const BenchmarkResult& res) {
  if (idx < size) {
    solve_times[idx] = res.solve_time;
    prim_dual_gaps[idx] = res.prim_dual_gap;
    prim_infeas_errs[idx] = res.prim_infeas_err;
    dual_infeas_errs[idx] = res.dual_infeas_err;
    iters[idx] = res.iter;
    optimal_flags[idx] = res.optimal_flag;
    ++idx;
  }
}

inline double BenchmarkResultArray::AvgSolveTime() const {
  if (idx == 0) return 0.0;
  return std::accumulate(solve_times.begin(), solve_times.end(), 0.0) / idx;
}

inline double BenchmarkResultArray::MaxPrimalDualGap() const {
  if (idx == 0) return 0.0;
  return *std::max_element(prim_dual_gaps.begin(), prim_dual_gaps.end());
}

inline double BenchmarkResultArray::MaxPrimalInfeasibilityError() const {
  if (idx == 0) return 0.0;
  return *std::max_element(prim_infeas_errs.begin(), prim_infeas_errs.end());
}

inline double BenchmarkResultArray::MaxDualInfeasibilityError() const {
  if (idx == 0) return 0.0;
  return *std::max_element(dual_infeas_errs.begin(), dual_infeas_errs.end());
}

inline double BenchmarkResultArray::AvgIterations() const {
  if (idx == 0) return 0.0;
  return std::accumulate(iters.begin(), iters.end(), 0.0) / idx;
}

inline int BenchmarkResultArray::NumSuboptimalResults() const {
  if (idx == 0) return 0;
  return idx - std::accumulate(optimal_flags.begin(), optimal_flags.end(), 0);
}

inline void BenchmarkResultArray::PrintStatistics() const {
  std::cout << "Avg. solve time (us)    : " << AvgSolveTime() << std::endl;
  std::cout << std::scientific;
  std::cout << "Max. prim dual gap      : " << MaxPrimalDualGap() << std::endl;
  std::cout << "Max. prim infeas err (m): " << MaxPrimalInfeasibilityError()
            << std::endl;
  std::cout << "Max. dual infeas err    : " << MaxDualInfeasibilityError()
            << std::endl;
  std::cout.unsetf(std::ios::fixed | std::ios::scientific);
  std::cout << "Avg. iterations         : " << AvgIterations() << std::endl;
  std::cout << "Num. suboptimal runs    : " << NumSuboptimalResults()
            << std::endl;
  std::cout << std::endl;
}

}  // namespace bench

#endif  // DGD_BENCHMARK_HELPERS_BENCHMARK_RESULT_H_

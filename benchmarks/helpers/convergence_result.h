#ifndef DGD_BENCHMARK_HELPERS_CONVERGENCE_RESULT_H_
#define DGD_BENCHMARK_HELPERS_CONVERGENCE_RESULT_H_

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace bench {

// Convergence statistics for each iteration.
struct IterationStats {
  explicit IterationStats(int max_iter) : max_iter_(max_iter) {
    histograms_.resize(max_iter, std::vector<long long>(nbins_, 0));
  }

  int AddResult(int iter, double ub, double lb);

  void AddResults(int last_iter, const std::vector<double>& ubs,
                  const std::vector<double>& lbs);

  double GetQuantile(int iter, double q) const;

  std::vector<double> GetQuantiles(int iter,
                                   const std::vector<double>& qs) const;

 private:
  std::vector<std::vector<long long>> histograms_;

  const int max_iter_;
  const int nbins_ = 1000;
  const double min_gap_ = 1e-12;
  const double max_gap_ = 1e3;
  const double log_min_ = std::log2(min_gap_);
  const double log_max_ = std::log2(max_gap_);
  const double inv_log_range_ = 1.0 / (log_max_ - log_min_);
};

inline int IterationStats::AddResult(int iter, double ub, double lb) {
  if (iter >= static_cast<int>(histograms_.size())) return -1;

  double log_gap;
  if (std::isinf(ub) || (lb <= 0.0)) {
    log_gap = log_max_;
  } else {
    const double gap = std::max(min_gap_, std::min(max_gap_, (ub / lb) - 1.0));
    log_gap = std::log2(gap);
  }

  int bin_idx =
      static_cast<int>((log_gap - log_min_) * inv_log_range_ * nbins_);
  bin_idx = std::max(0, std::min(nbins_ - 1, bin_idx));

  histograms_[iter][bin_idx]++;
  return bin_idx;
}

inline void IterationStats::AddResults(int last_iter,
                                       const std::vector<double>& ubs,
                                       const std::vector<double>& lbs) {
  int iter = 0;
  for (; iter < last_iter; ++iter) {
    AddResult(iter, ubs[iter], lbs[iter]);
  }
  const int bin_idx = AddResult(iter, ubs[iter], lbs[iter]);
  for (iter = last_iter + 1; iter < max_iter_; ++iter)
    histograms_[iter][bin_idx]++;
}

inline double IterationStats::GetQuantile(int iter, double q) const {
  const auto& hist = histograms_[iter];

  long long total = 0LL;
  for (const auto& count : hist) total += count;

  if (total == 0LL) return max_gap_;

  q = std::max(0.0, std::min(1.0, q));
  long long target = static_cast<long long>(q * static_cast<double>(total));
  long long cumulative = 0LL;
  for (int i = 0; i < nbins_; ++i) {
    cumulative += hist[i];
    if (cumulative >= target) {
      return std::pow(2.0,
                      log_min_ + (i + 0.5) * (log_max_ - log_min_) / nbins_);
    }
  }
  return max_gap_;
}

inline std::vector<double> IterationStats::GetQuantiles(
    int iter, const std::vector<double>& qs) const {
  std::vector<double> res;
  res.reserve(qs.size());
  for (const double q : qs) res.push_back(GetQuantile(iter, q));
  return res;
}

// Convergence result array.
struct ConvergenceResultArray {
  const int max_iter;
  std::vector<double> min, q01, q25, q50, q75, q99, max;

  explicit ConvergenceResultArray(int max_iter);

  void ComputeStatistics(const IterationStats& stats);

  // Writes benchmark output to a Feather file.
  bool SaveToFile(const std::string& filename) const;
};

inline ConvergenceResultArray::ConvergenceResultArray(int max_iter)
    : max_iter(max_iter),
      min(max_iter),
      q01(max_iter),
      q25(max_iter),
      q50(max_iter),
      q75(max_iter),
      q99(max_iter),
      max(max_iter) {}

inline void ConvergenceResultArray::ComputeStatistics(
    const IterationStats& stats) {
  for (int iter = 0; iter < max_iter; ++iter) {
    min[iter] = stats.GetQuantile(iter, 0.0);
    q01[iter] = stats.GetQuantile(iter, 0.01);
    q25[iter] = stats.GetQuantile(iter, 0.25);
    q50[iter] = stats.GetQuantile(iter, 0.5);
    q75[iter] = stats.GetQuantile(iter, 0.75);
    q99[iter] = stats.GetQuantile(iter, 0.99);
    max[iter] = stats.GetQuantile(iter, 1.0);
  }
}

}  // namespace bench

#endif  // DGD_BENCHMARK_HELPERS_CONVERGENCE_RESULT_H_

#ifndef DGD_BENCHMARK_BENCHMARK_INTERFACE_H_
#define DGD_BENCHMARK_BENCHMARK_INTERFACE_H_

#include <memory>
#include <vector>

#include "dgd/geometry/convex_set.h"
#include "dgd/output.h"
#include "dgd/settings.h"
#include "dgd_benchmark/benchmark_result.h"
#include "dgd_benchmark/data_types.h"
#include "dgd_benchmark/dsf/dsf.h"
#include "dgd_benchmark/dsf/dsf_collision.h"
#include "dgd_benchmark/ie/internal_expanding.h"
#include "dgd_benchmark/inc/incremental.h"
#include "dgd_benchmark/inc/polyhedron.h"
#include "dgd_benchmark/timer.h"

// Interface class for calling distance functions and retrieving results.
class BenchmarkInterface {
 public:
  explicit BenchmarkInterface(int ncold, int nwarm);

  ~BenchmarkInterface() = default;

  void SetDsfSettings(const dsf::Settings& settings);

  void SetIeSettings(const ie::Settings& settings);

  void SetIncSettings(const inc::Settings& settings);

  void SetDgdSettings(const dgd::Settings& settings);

  // DSF benchmarks.
  void DsfColdStart(dsf::DSF* set1, const Transform3& tf1, dsf::DSF* set2,
                    const Transform3& tf2, BenchmarkResult& res);

  // IE benchmarks.
  void IeColdStart(const dgd::ConvexSet<3>* set1, const Transform3& tf1,
                   const dgd::ConvexSet<3>* set2, const Transform3& tf2,
                   BenchmarkResult& res);

  // Incremental benchmarks.
  void IncColdStart(const inc::Polyhedron* set1,
                    const dgd::ConvexSet<3>* set1_v, const Transform3& tf1,
                    const inc::Polyhedron* set2,
                    const dgd::ConvexSet<3>* set2_v, const Transform3& tf2,
                    BenchmarkResult& res);

  void IncWarmStart(const inc::Polyhedron* set1,
                    const dgd::ConvexSet<3>* set1_v, const Transform3& tf1,
                    const inc::Polyhedron* set2,
                    const dgd::ConvexSet<3>* set2_v, const Transform3& tf2,
                    const Vec3& dx, const Rotation3& drot,
                    BenchmarkResultArray& resa);

  // DGD benchmarks.
  void DgdColdStart(const dgd::ConvexSet<3>* set1, const Transform3& tf1,
                    const dgd::ConvexSet<3>* set2, const Transform3& tf2,
                    BenchmarkResult& res);

  void DgdWarmStart(const dgd::ConvexSet<3>* set1, const Transform3& tf1,
                    const dgd::ConvexSet<3>* set2, const Transform3& tf2,
                    const Vec3& dx, const Rotation3& drot,
                    BenchmarkResultArray& resa);

  int ncold() const { return ncold_; }

  int nwarm() const { return nwarm_; }

 private:
  // Differential support function variables.
  struct {
    dsf::Settings settings;
  } dsf_;

  // Internal expanding variables.
  struct {
    ie::Settings settings;
  } ie_;

  // Incremental variables.
  struct {
    std::shared_ptr<inc::GrowthDistanceSolver> solver;
    std::vector<inc::SolutionStatus> status_ws;  // Temporary.
  } inc_;

  // Differentiable growth distance variables.
  struct {
    dgd::Settings settings;
    std::vector<dgd::SolutionStatus> status_ws;  // Temporary.
  } dgd_;

  // Temporary storage variables.
  std::vector<Vec3> z1_ws_, z2_ws_, normal_ws_;
  std::vector<double> gd_ws_;
  std::vector<int> iter_ws_;

  Timer timer_;
  // Number of calls per distance function call.
  const int ncold_;
  // Number of timesteps after the initial solve for which warm-start is used.
  const int nwarm_;
};

inline void BenchmarkInterface::SetDsfSettings(const dsf::Settings& settings) {
  dsf_.settings = settings;
}

inline void BenchmarkInterface::SetIeSettings(const ie::Settings& settings) {
  ie_.settings = settings;
}

inline void BenchmarkInterface::SetIncSettings(const inc::Settings& settings) {
  inc_.solver->SolverSettings() = settings;
}

inline void BenchmarkInterface::SetDgdSettings(const dgd::Settings& settings) {
  dgd_.settings = settings;
}

#endif  // DGD_BENCHMARK_BENCHMARK_INTERFACE_H_

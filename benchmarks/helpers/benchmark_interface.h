#ifndef DGD_BENCHMARK_HELPERS_BENCHMARK_INTERFACE_H_
#define DGD_BENCHMARK_HELPERS_BENCHMARK_INTERFACE_H_

#include <memory>
#include <vector>

#include "dcf/dcf_collision.h"
#include "dcf/dsf.h"
#include "dcf/dsf_interface.h"
#include "dcf/precompiled.h"
#include "dgd/data_types.h"
#include "dgd/geometry/convex_set.h"
#include "dgd/output.h"
#include "dgd/settings.h"
#include "dgd/utils/random.h"
#include "dgd/utils/timer.h"
#include "helpers/benchmark_result.h"
#include "ie/internal_expanding.h"
#include "inc/data_types.h"
#include "inc/incremental.h"
#include "inc/polyhedron.h"
#include "internal_helpers/set_generator.h"

namespace bench {

enum class DgdSolverType { CuttingPlane, TrustRegionNewton };

enum class DgdBcSolverType {
  Cramer,
  LU,
};

struct OptimalSolution {
  dgd::Vec3r z1, z2, normal;
  dgd::Real gd;
  int iter;
  int status;
};

// Interface class for calling distance functions and retrieving results.
class BenchmarkInterface {
 public:
  static constexpr unsigned int kVdsfExp = 8;

  using DsfPtr = std::shared_ptr<dcf::VDSFInterface<kVdsfExp>>;

  explicit BenchmarkInterface(int ncold, int nwarm);

  ~BenchmarkInterface() = default;

  void SetDcfSettings(const dcf::Settings& settings) {
    dcf_.settings = settings;
  }

  void SetIeSettings(const ie::Settings& settings) { ie_.settings = settings; }

  void SetIncSettings(const inc::Settings& settings) {
    inc_.solver->SolverSettings() = settings;
  }

  void SetDgdSettings(const dgd::Settings& settings) {
    dgd_.settings = settings;
  }

  // Load meshes.
  void LoadMeshesFromObjFiles(const std::vector<std::string>& filenames);

  // DCF benchmarks.
  void DcfColdStart(int set1_idx, const dcf::Transform3& tf1, int set2_idx,
                    const dcf::Transform3& tf2, BenchmarkResultArray& res_arr);

  // IE benchmarks.
  void IeColdStart(const dgd::ConvexSet<3>* set1, const dgd::Transform3r& tf1,
                   const dgd::ConvexSet<3>* set2, const dgd::Transform3r& tf2,
                   BenchmarkResultArray& res_arr);

  // Incremental benchmarks.
  void IncColdStart(int set1_idx, const inc::Transform3& tf1, int set2_idx,
                    const inc::Transform3& tf2, BenchmarkResultArray& res_arr);

  void IncWarmStart(int set1_idx, const inc::Transform3& tf1, int set2_idx,
                    const inc::Transform3& tf2, const inc::Vec3& dx,
                    const inc::Rotation3& drot, BenchmarkResultArray& res_arr);

  // DGD benchmarks.
  template <DgdSolverType S, DgdBcSolverType BST>
  void DgdColdStart(const dgd::ConvexSet<3>* set1, const dgd::Transform3r& tf1,
                    const dgd::ConvexSet<3>* set2, const dgd::Transform3r& tf2,
                    BenchmarkResultArray& res_arr);

  template <DgdSolverType S, DgdBcSolverType BST>
  void DgdWarmStart(const dgd::ConvexSet<3>* set1, const dgd::Transform3r& tf1,
                    const dgd::ConvexSet<3>* set2, const dgd::Transform3r& tf2,
                    const dgd::Vec3r& dx, const dgd::Rotation3r& drot,
                    BenchmarkResultArray& res_arr, dgd::WarmStartType ws_type);

  void SetRngSeed(unsigned int seed = 5489u) {
    rng_.SetSeed(seed);
    generator_->SetRngSeed(seed);
  }

  void SetRandomRngSeed() {
    rng_.SetRandomSeed();
    generator_->SetRandomRngSeed();
  }

  int RandomMeshIndex() { return rng_.RandomInt(0, nmeshes_ - 1); }

  dgd::Rng& rng() { return rng_; }

  const std::vector<DsfPtr>& vdsfs() const { return vdsfs_; }

  const std::vector<std::shared_ptr<inc::Polyhedron>>& polyhedra() const {
    return polyhedra_;
  }

  const std::vector<std::shared_ptr<dgd::Mesh>>& meshes() const {
    return generator_->meshes();
  }

  const std::shared_ptr<dgd::ConvexSet<3>> RandomCurvedPrimitiveSet() const {
    return generator_->GetRandomCurvedPrimitive3DSet();
  }

  int ncold() const { return ncold_; }

  int nwarm() const { return nwarm_; }

  int nmeshes() const { return nmeshes_; }

 private:
  // Differential contact features variables.
  struct {
    dcf::Settings settings;
  } dcf_;

  // Internal expanding variables.
  struct {
    ie::Settings settings;
  } ie_;

  // Incremental variables.
  struct {
    std::shared_ptr<inc::GrowthDistanceSolver> solver;
  } inc_;

  // Differentiable growth distance variables.
  struct {
    dgd::Settings settings;
  } dgd_;

  std::vector<DsfPtr> vdsfs_;
  std::vector<std::shared_ptr<inc::Polyhedron>> polyhedra_;
  std::shared_ptr<dgd::bench::ConvexSetGenerator> generator_;

  // Temporary storage variables.
  std::vector<OptimalSolution> opt_sols_;

  dgd::Rng rng_;
  dgd::Timer timer_;
  // Number of calls per distance function call.
  const int ncold_;
  // Number of timesteps after the initial solve for which warm-start is used.
  const int nwarm_;
  // Number of meshes.
  int nmeshes_;
};

}  // namespace bench

#endif  // DGD_BENCHMARK_HELPERS_BENCHMARK_INTERFACE_H_

#include "helpers/benchmark_interface.h"

#include <stdexcept>

#include "dgd/dgd.h"
#include "dgd/error_metrics.h"
#include "dgd/geometry/3d/mesh.h"
#include "dgd/mesh_loader.h"
#include "inc/solution_error.h"
#include "internal_helpers/math_utils.h"
#include "internal_helpers/mesh_loader.h"

namespace bench {

BenchmarkInterface::BenchmarkInterface(int ncold, int nwarm)
    : ncold_(ncold), nwarm_(nwarm), nmeshes_(0) {
  inc_.solver = std::make_shared<inc::GrowthDistanceSolver>();

  vdsfs_.clear();
  polyhedra_.clear();
  const dgd::bench::ConvexSetFeatureRange fr{};
  generator_ = std::make_shared<dgd::bench::ConvexSetGenerator>(fr);

  opt_sols_.clear();
  opt_sols_.resize(nwarm);

  timer_.Start();
  timer_.Stop();
}

void BenchmarkInterface::LoadMeshesFromObjFiles(
    const std::vector<std::string>& filenames) {
  generator_->LoadMeshesFromObjFiles(filenames);
  const auto meshes = generator_->meshes();
  nmeshes_ = generator_->nmeshes();

  vdsfs_.resize(nmeshes_);
  polyhedra_.resize(nmeshes_);
  dgd::bench::MeshProperties mp;
  for (int i = 0; i < nmeshes_; ++i) {
    vdsfs_[i] = std::make_shared<dcf::VDSFInterface<kVdsfExp>>(
        meshes[i]->vertices(), meshes[i]->inradius(), 0.0);
    mp.SetFacetMeshFromVertices(meshes[i]->vertices());
    polyhedra_[i] =
        std::make_shared<inc::Polyhedron>(mp.normal, mp.offset, mp.fgraph);
  }
}

void BenchmarkInterface::DcfColdStart(int set1_idx, const dcf::Transform3& tf1,
                                      int set2_idx, const dcf::Transform3& tf2,
                                      BenchmarkResultArray& res_arr) {
  if ((set1_idx >= nmeshes_) || (set2_idx >= nmeshes_)) {
    throw std::range_error("Set indices are out of range");
  }
  dcf::Output out;
  dcf::DSF* set1 = vdsfs_[set1_idx]->VDSFPtr();
  dcf::DSF* set2 = vdsfs_[set2_idx]->VDSFPtr();
  // Initial run to account for cache misses.
  dcf::GrowthDistance(set1, tf1, set2, tf2, dcf_.settings, out);

  timer_.Stop();
  timer_.Start();
  for (int i = 0; i < ncold_; ++i) {
    dcf::GrowthDistance(set1, tf1, set2, tf2, dcf_.settings, out);
  }
  timer_.Stop();
  const auto err = dcf::ComputeSolutionError(set1, tf1, set2, tf2, out);

  BenchmarkResult res;
  res.solve_time = timer_.Elapsed() / double(ncold_);
  res.prim_dual_gap = err.prim_dual_gap;
  res.prim_infeas_err = err.prim_infeas_err;
  res.dual_infeas_err = err.dual_infeas_err;
  res.iter = out.iter;
  res.optimal_flag = (out.status == dcf::SolutionStatus::CoincidentCenters) ||
                     (out.status == dcf::SolutionStatus::Optimal);
  res_arr.AddResult(res);
}

void BenchmarkInterface::IeColdStart(const dgd::ConvexSet<3>* set1,
                                     const dgd::Transform3r& tf1,
                                     const dgd::ConvexSet<3>* set2,
                                     const dgd::Transform3r& tf2,
                                     BenchmarkResultArray& res_arr) {
  ie::Output out;
  ie::GrowthDistance(set1, tf1, set2, tf2, ie_.settings, out);

  timer_.Stop();
  timer_.Start();
  for (int i = 0; i < ncold_; ++i) {
    ie::GrowthDistance(set1, tf1, set2, tf2, ie_.settings, out);
  }
  timer_.Stop();
  const auto err = ie::ComputeSolutionError(set1, tf1, set2, tf2, out);

  BenchmarkResult res;
  res.solve_time = timer_.Elapsed() / double(ncold_);
  res.prim_dual_gap = err.prim_dual_gap;
  res.prim_infeas_err = err.prim_infeas_err;
  res.dual_infeas_err = err.dual_infeas_err;
  res.iter = out.iter;
  res.optimal_flag = (out.status == ie::SolutionStatus::CoincidentCenters) ||
                     (out.status == ie::SolutionStatus::Optimal);
  res_arr.AddResult(res);
}

void BenchmarkInterface::IncColdStart(int set1_idx, const inc::Transform3& tf1,
                                      int set2_idx, const inc::Transform3& tf2,
                                      BenchmarkResultArray& res_arr) {
  if ((set1_idx >= nmeshes_) || (set2_idx >= nmeshes_)) {
    throw std::range_error("Set indices are out of range");
  }
  inc::Output out;
  const inc::Polyhedron* set1 = polyhedra_[set1_idx].get();
  const inc::Polyhedron* set2 = polyhedra_[set2_idx].get();
  inc_.solver->GrowthDistance(set1, tf1, set2, tf2, out);

  timer_.Stop();
  timer_.Start();
  for (int i = 0; i < ncold_; ++i) {
    inc_.solver->GrowthDistance(set1, tf1, set2, tf2, out);
  }
  timer_.Stop();
  const dgd::ConvexSet<3>* set1_v = generator_->meshes()[set1_idx].get();
  const dgd::ConvexSet<3>* set2_v = generator_->meshes()[set2_idx].get();
  const auto err = inc::ComputeSolutionError(set1_v, tf1, set2_v, tf2, out);

  BenchmarkResult res;
  res.solve_time = timer_.Elapsed() / double(ncold_);
  res.prim_dual_gap = err.prim_dual_gap;
  res.prim_infeas_err = err.prim_infeas_err;
  res.dual_infeas_err = err.dual_infeas_err;
  res.iter = out.iter;
  res.optimal_flag = (out.status == inc::SolutionStatus::CoincidentCenters) ||
                     (out.status == inc::SolutionStatus::Optimal);
  res_arr.AddResult(res);
}

namespace {

inline void SetOptimalSolution(const inc::Output& out,
                               OptimalSolution& opt_sol) {
  opt_sol.z1 = out.z1;
  opt_sol.z2 = out.z2;
  opt_sol.normal = out.normal;
  opt_sol.gd = out.growth_dist;
  opt_sol.iter = out.iter;
  opt_sol.status = static_cast<int>(out.status);
}

inline void SetOutput(const OptimalSolution& opt_sol, inc::Output& out) {
  out.z1 = opt_sol.z1;
  out.z2 = opt_sol.z2;
  out.normal = opt_sol.normal;
  out.growth_dist = opt_sol.gd;
  out.status = static_cast<inc::SolutionStatus>(opt_sol.status);
}

}  // namespace

void BenchmarkInterface::IncWarmStart(int set1_idx, const inc::Transform3& tf1,
                                      int set2_idx, const inc::Transform3& tf2,
                                      const inc::Vec3& dx,
                                      const inc::Rotation3& drot,
                                      BenchmarkResultArray& res_arr) {
  if ((set1_idx >= nmeshes_) || (set2_idx >= nmeshes_)) {
    throw std::range_error("Set indices are out of range");
  }
  inc::Output out;
  // Cold start.
  const inc::Polyhedron* set1 = polyhedra_[set1_idx].get();
  const inc::Polyhedron* set2 = polyhedra_[set2_idx].get();
  inc_.solver->GrowthDistance(set1, tf1, set2, tf2, out, false);

  inc::Transform3 tf1_t{tf1};
  opt_sols_.clear();
  opt_sols_.resize(nwarm_);
  timer_.Stop();
  timer_.Start();
  for (int i = 0; i < nwarm_; ++i) {
    dgd::bench::UpdateTransform(tf1_t, dx, drot);
    inc_.solver->GrowthDistance(set1, tf1_t, set2, tf2, out, true);
    SetOptimalSolution(out, opt_sols_[i]);
  }
  timer_.Stop();
  const double solve_time_avg = timer_.Elapsed() / double(nwarm_);

  const dgd::ConvexSet<3>* set1_v = generator_->meshes()[set1_idx].get();
  const dgd::ConvexSet<3>* set2_v = generator_->meshes()[set2_idx].get();
  tf1_t = tf1;
  for (int i = 0; i < nwarm_; ++i) {
    dgd::bench::UpdateTransform(tf1_t, dx, drot);
    SetOutput(opt_sols_[i], out);
    const auto err = inc::ComputeSolutionError(set1_v, tf1_t, set2_v, tf2, out);
    BenchmarkResult res;
    res.solve_time = solve_time_avg;
    res.prim_dual_gap = err.prim_dual_gap;
    res.prim_infeas_err = err.prim_infeas_err;
    res.dual_infeas_err = err.dual_infeas_err;
    res.iter = opt_sols_[i].iter;
    res.optimal_flag = (out.status == inc::SolutionStatus::CoincidentCenters) ||
                       (out.status == inc::SolutionStatus::Optimal);
    res_arr.AddResult(res);
  }
}

template <DgdSolverType S, DgdBcSolverType BST>
void BenchmarkInterface::DgdColdStart(const dgd::ConvexSet<3>* set1,
                                      const dgd::Transform3r& tf1,
                                      const dgd::ConvexSet<3>* set2,
                                      const dgd::Transform3r& tf2,
                                      BenchmarkResultArray& res_arr) {
  using dgd::detail::BcSolverType;
  [[maybe_unused]] constexpr BcSolverType bst = (BST == DgdBcSolverType::Cramer)
                                                    ? BcSolverType::kCramer
                                                    : BcSolverType::kLU;

  using C = dgd::ConvexSet<3>;
  dgd::Output<3> out;
  if constexpr (S == DgdSolverType::CuttingPlane) {
    dgd::GrowthDistanceCp<3, C, C, bst>(set1, tf1, set2, tf2, dgd_.settings,
                                        out);
  } else {
    dgd::GrowthDistanceTrn(set1, tf1, set2, tf2, dgd_.settings, out);
  }

  timer_.Stop();
  timer_.Start();
  for (int i = 0; i < ncold_; ++i) {
    if constexpr (S == DgdSolverType::CuttingPlane) {
      dgd::GrowthDistanceCp<3, C, C, bst>(set1, tf1, set2, tf2, dgd_.settings,
                                          out);
    } else {
      dgd::GrowthDistanceTrn(set1, tf1, set2, tf2, dgd_.settings, out);
    }
  }
  timer_.Stop();
  const auto err = dgd::ComputeSolutionError(set1, tf1, set2, tf2, out);

  // TODO:
  if (out.status != dgd::SolutionStatus::Optimal) {
    exit(EXIT_FAILURE);
  }

  BenchmarkResult res;
  res.solve_time = timer_.Elapsed() / double(ncold_);
  res.prim_dual_gap = err.prim_dual_gap;
  res.prim_infeas_err = err.prim_infeas_err;
  res.dual_infeas_err = err.dual_infeas_err;
  res.iter = out.iter;
  res.optimal_flag = (out.status == dgd::SolutionStatus::CoincidentCenters) ||
                     (out.status == dgd::SolutionStatus::Optimal);
  res_arr.AddResult(res);
}

template void BenchmarkInterface::DgdColdStart<DgdSolverType::CuttingPlane,
                                               DgdBcSolverType::Cramer>(
    const dgd::ConvexSet<3>* set1, const dgd::Transform3r& tf1,
    const dgd::ConvexSet<3>* set2, const dgd::Transform3r& tf2,
    BenchmarkResultArray& res_arr);

template void BenchmarkInterface::DgdColdStart<DgdSolverType::CuttingPlane,
                                               DgdBcSolverType::LU>(
    const dgd::ConvexSet<3>* set1, const dgd::Transform3r& tf1,
    const dgd::ConvexSet<3>* set2, const dgd::Transform3r& tf2,
    BenchmarkResultArray& res_arr);

template void BenchmarkInterface::DgdColdStart<DgdSolverType::TrustRegionNewton,
                                               DgdBcSolverType::LU>(
    const dgd::ConvexSet<3>* set1, const dgd::Transform3r& tf1,
    const dgd::ConvexSet<3>* set2, const dgd::Transform3r& tf2,
    BenchmarkResultArray& res_arr);

namespace {

inline void SetOptimalSolution(const dgd::Output<3>& out,
                               OptimalSolution& opt_sol) {
  opt_sol.z1 = out.z1;
  opt_sol.z2 = out.z2;
  opt_sol.normal = out.normal;
  opt_sol.gd = out.growth_dist_ub;
  opt_sol.iter = out.iter;
  opt_sol.status = static_cast<int>(out.status);
}

inline void SetOutput(const OptimalSolution& opt_sol, dgd::Output<3>& out) {
  out.z1 = opt_sol.z1;
  out.z2 = opt_sol.z2;
  out.normal = opt_sol.normal;
  out.growth_dist_ub = opt_sol.gd;
  out.status = static_cast<dgd::SolutionStatus>(opt_sol.status);
}

}  // namespace

template <DgdSolverType S, DgdBcSolverType BST>
void BenchmarkInterface::DgdWarmStart(
    const dgd::ConvexSet<3>* set1, const dgd::Transform3r& tf1,
    const dgd::ConvexSet<3>* set2, const dgd::Transform3r& tf2,
    const dgd::Vec3r& dx, const dgd::Rotation3r& drot,
    BenchmarkResultArray& res_arr, dgd::WarmStartType ws_type) {
  static_assert(S != DgdSolverType::TrustRegionNewton,
                "The TRN solver does not support warm starting.");
  using dgd::detail::BcSolverType;
  constexpr BcSolverType bst = (BST == DgdBcSolverType::Cramer)
                                   ? BcSolverType::kCramer
                                   : BcSolverType::kLU;

  using C = dgd::ConvexSet<3>;
  dgd_.settings.ws_type = ws_type;
  dgd::Output<3> out;
  if constexpr (S == DgdSolverType::CuttingPlane) {
    dgd::GrowthDistanceCp<3, C, C, bst>(set1, tf1, set2, tf2, dgd_.settings,
                                        out);
  }

  dgd::Transform3r tf1_t{tf1};
  opt_sols_.clear();
  opt_sols_.resize(nwarm_);
  timer_.Stop();
  timer_.Start();
  for (int i = 0; i < nwarm_; ++i) {
    dgd::bench::UpdateTransform(tf1_t, dx, drot);
    if constexpr (S == DgdSolverType::CuttingPlane) {
      dgd::GrowthDistanceCp<3, C, C, bst>(set1, tf1_t, set2, tf2, dgd_.settings,
                                          out, true);
    }
    SetOptimalSolution(out, opt_sols_[i]);
  }
  timer_.Stop();
  const double solve_time_avg = timer_.Elapsed() / double(nwarm_);

  tf1_t = tf1;
  for (int i = 0; i < nwarm_; ++i) {
    dgd::bench::UpdateTransform(tf1_t, dx, drot);
    SetOutput(opt_sols_[i], out);
    const auto err = dgd::ComputeSolutionError(set1, tf1_t, set2, tf2, out);
    BenchmarkResult res;
    res.solve_time = solve_time_avg;
    res.prim_dual_gap = err.prim_dual_gap;
    res.prim_infeas_err = err.prim_infeas_err;
    res.dual_infeas_err = err.dual_infeas_err;
    res.iter = opt_sols_[i].iter;
    res.optimal_flag = (out.status == dgd::SolutionStatus::CoincidentCenters) ||
                       (out.status == dgd::SolutionStatus::Optimal);
    res_arr.AddResult(res);
  }
}

template void BenchmarkInterface::DgdWarmStart<DgdSolverType::CuttingPlane,
                                               DgdBcSolverType::Cramer>(
    const dgd::ConvexSet<3>* set1, const dgd::Transform3r& tf1,
    const dgd::ConvexSet<3>* set2, const dgd::Transform3r& tf2,
    const dgd::Vec3r& dx, const dgd::Rotation3r& drot,
    BenchmarkResultArray& res_arr, dgd::WarmStartType ws_type);

template void BenchmarkInterface::DgdWarmStart<DgdSolverType::CuttingPlane,
                                               DgdBcSolverType::LU>(
    const dgd::ConvexSet<3>* set1, const dgd::Transform3r& tf1,
    const dgd::ConvexSet<3>* set2, const dgd::Transform3r& tf2,
    const dgd::Vec3r& dx, const dgd::Rotation3r& drot,
    BenchmarkResultArray& res_arr, dgd::WarmStartType ws_type);

}  // namespace bench

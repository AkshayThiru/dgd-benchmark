#include "dgd_benchmark/benchmark_interface.h"

#include "dgd/error_metrics.h"
#include "dgd/growth_distance.h"
#include "dgd/output.h"
#include "dgd_benchmark/benchmark_result.h"
#include "dgd_benchmark/data_types.h"
#include "dgd_benchmark/inc/solution_error.h"

BenchmarkInterface::BenchmarkInterface(int ncold, int nwarm)
    : ncold_(ncold), nwarm_(nwarm) {
  inc_.solver = std::make_shared<inc::GrowthDistanceSolver>();

  z1_ws_.resize(nwarm);
  z2_ws_.resize(nwarm);
  normal_ws_.resize(nwarm);
  gd_ws_.resize(nwarm);
  iter_ws_.resize(nwarm);
  inc_.status_ws.resize(nwarm);
  dgd_.status_ws.resize(nwarm);

  timer_.Start();
  timer_.Stop();
}

void BenchmarkInterface::DsfColdStart(dsf::DSF* set1, const Transform3& tf1,
                                      dsf::DSF* set2, const Transform3& tf2,
                                      BenchmarkResult& res) {
  dsf::Output out;
  // Initial run to account for CPU cache misses.
  dsf::GrowthDistance(set1, tf1, set2, tf2, dsf_.settings, out);

  timer_.Stop();
  timer_.Start();
  for (int i = 0; i < ncold_; ++i) {
    dsf::GrowthDistance(set1, tf1, set2, tf2, dsf_.settings, out);
  }
  timer_.Stop();
  const auto err = dsf::ComputeSolutionError(set1, tf1, set2, tf2, out);

  res.solve_time = timer_.Elapsed() / double(ncold_);
  res.prim_dual_gap = err.prim_dual_gap;
  res.prim_feas_err = err.prim_feas_err;
  res.dual_feas_err = err.dual_feas_err;
  res.iter = out.iter;
  res.optimal_flag = (out.status == dsf::SolutionStatus::CoincidentCenters) ||
                     (out.status == dsf::SolutionStatus::Optimal);
}

void BenchmarkInterface::IeColdStart(const dgd::ConvexSet<3>* set1,
                                     const Transform3& tf1,
                                     const dgd::ConvexSet<3>* set2,
                                     const Transform3& tf2,
                                     BenchmarkResult& res) {
  ie::Output out;
  ie::GrowthDistance(set1, tf1, set2, tf2, ie_.settings, out);

  timer_.Stop();
  timer_.Start();
  for (int i = 0; i < ncold_; ++i) {
    ie::GrowthDistance(set1, tf1, set2, tf2, ie_.settings, out);
  }
  timer_.Stop();
  const auto err = ie::ComputeSolutionError(set1, tf1, set2, tf2, out);

  res.solve_time = timer_.Elapsed() / double(ncold_);
  res.prim_dual_gap = err.prim_dual_gap;
  res.prim_feas_err = err.prim_feas_err;
  res.dual_feas_err = err.dual_feas_err;
  res.iter = out.iter;
  res.optimal_flag = (out.status == ie::SolutionStatus::CoincidentCenters) ||
                     (out.status == ie::SolutionStatus::Optimal);
}

void BenchmarkInterface::IncColdStart(const inc::Polyhedron* set1,
                                      const dgd::ConvexSet<3>* set1_v,
                                      const Transform3& tf1,
                                      const inc::Polyhedron* set2,
                                      const dgd::ConvexSet<3>* set2_v,
                                      const Transform3& tf2,
                                      BenchmarkResult& res) {
  inc::Output out;
  inc_.solver->GrowthDistance(set1, tf1, set2, tf2, out);

  timer_.Stop();
  timer_.Start();
  for (int i = 0; i < ncold_; ++i) {
    inc_.solver->GrowthDistance(set1, tf1, set2, tf2, out);
  }
  timer_.Stop();
  const auto err = inc::ComputeSolutionError(set1_v, tf1, set2_v, tf2, out);

  res.solve_time = timer_.Elapsed() / double(ncold_);
  res.prim_dual_gap = err.prim_dual_gap;
  res.prim_feas_err = err.prim_feas_err;
  res.dual_feas_err = err.dual_feas_err;
  res.iter = out.iter;
  res.optimal_flag = (out.status == inc::SolutionStatus::CoincidentCenters) ||
                     (out.status == inc::SolutionStatus::Optimal);
}

void BenchmarkInterface::IncWarmStart(
    const inc::Polyhedron* set1, const dgd::ConvexSet<3>* set1_v,
    const Transform3& tf1, const inc::Polyhedron* set2,
    const dgd::ConvexSet<3>* set2_v, const Transform3& tf2, const Vec3& dx,
    const Rotation3& drot, BenchmarkResultArray& resa) {
  inc::Output out;
  // Cold start.
  inc_.solver->GrowthDistance(set1, tf1, set2, tf2, out, false);

  Transform3 tf1_t{tf1}, tf2_t{tf2};
  timer_.Stop();
  timer_.Start();
  for (int i = 0; i < nwarm_; ++i) {
    tf1_t.block<3, 1>(0, 3) += dx;
    tf1_t.block<3, 3>(0, 0) *= drot;
    inc_.solver->GrowthDistance(set1, tf1_t, set2, tf2_t, out, true);
    z1_ws_[i] = out.z1;
    z2_ws_[i] = out.z2;
    normal_ws_[i] = out.normal;
    gd_ws_[i] = out.growth_dist;
    iter_ws_[i] = out.iter;
    inc_.status_ws[i] = out.status;
  }
  timer_.Stop();
  const double solve_time_avg = timer_.Elapsed() / double(nwarm_);

  inc::SolutionError err;
  BenchmarkResult res;
  tf1_t = tf1;
  tf2_t = tf2;
  for (int i = 0; i < nwarm_; ++i) {
    tf1_t.block<3, 1>(0, 3) += dx;
    tf1_t.block<3, 3>(0, 0) *= drot;
    out.z1 = z1_ws_[i];
    out.z2 = z2_ws_[i];
    out.normal = normal_ws_[i];
    out.growth_dist = gd_ws_[i];
    out.status = inc_.status_ws[i];
    err = inc::ComputeSolutionError(set1_v, tf1_t, set2_v, tf2_t, out);
    res.solve_time = solve_time_avg;
    res.prim_dual_gap = err.prim_dual_gap;
    res.prim_feas_err = err.prim_feas_err;
    res.dual_feas_err = err.dual_feas_err;
    res.iter = iter_ws_[i];
    res.optimal_flag = (out.status == inc::SolutionStatus::CoincidentCenters) ||
                       (out.status == inc::SolutionStatus::Optimal);
    resa.AddResult(res);
  }
}

void BenchmarkInterface::DgdColdStart(const dgd::ConvexSet<3>* set1,
                                      const Transform3& tf1,
                                      const dgd::ConvexSet<3>* set2,
                                      const Transform3& tf2,
                                      BenchmarkResult& res) {
  dgd::Output<3> out;
  dgd::GrowthDistance(set1, tf1, set2, tf2, dgd_.settings, out);

  timer_.Stop();
  timer_.Start();
  for (int i = 0; i < ncold_; ++i) {
    dgd::GrowthDistance(set1, tf1, set2, tf2, dgd_.settings, out);
  }
  timer_.Stop();
  const auto err = dgd::ComputeSolutionError(set1, tf1, set2, tf2, out);

  res.solve_time = timer_.Elapsed() / double(ncold_);
  res.prim_dual_gap = err.prim_dual_gap;
  res.prim_feas_err = err.prim_feas_err;
  res.dual_feas_err = err.dual_feas_err;
  res.iter = out.iter;
  res.optimal_flag = (out.status == dgd::SolutionStatus::CoincidentCenters) ||
                     (out.status == dgd::SolutionStatus::Optimal);
}

void BenchmarkInterface::DgdWarmStart(const dgd::ConvexSet<3>* set1,
                                      const Transform3& tf1,
                                      const dgd::ConvexSet<3>* set2,
                                      const Transform3& tf2, const Vec3& dx,
                                      const Rotation3& drot,
                                      BenchmarkResultArray& resa) {
  dgd::Output<3> out;
  dgd::GrowthDistance(set1, tf1, set2, tf2, dgd_.settings, out);

  Transform3 tf1_t{tf1}, tf2_t{tf2};
  timer_.Stop();
  timer_.Start();
  for (int i = 0; i < nwarm_; ++i) {
    tf1_t.block<3, 1>(0, 3) += dx;
    tf1_t.block<3, 3>(0, 0) *= drot;
    dgd::GrowthDistance(set1, tf1_t, set2, tf2_t, dgd_.settings, out, true);
    z1_ws_[i] = out.z1;
    z2_ws_[i] = out.z2;
    normal_ws_[i] = out.normal;
    gd_ws_[i] = out.growth_dist_ub;
    iter_ws_[i] = out.iter;
    dgd_.status_ws[i] = out.status;
  }
  timer_.Stop();
  const double solve_time_avg = timer_.Elapsed() / double(nwarm_);

  dgd::SolutionError err;
  BenchmarkResult res;
  tf1_t = tf1;
  tf2_t = tf2;
  for (int i = 0; i < nwarm_; ++i) {
    tf1_t.block<3, 1>(0, 3) += dx;
    tf1_t.block<3, 3>(0, 0) *= drot;
    out.z1 = z1_ws_[i];
    out.z2 = z2_ws_[i];
    out.normal = normal_ws_[i];
    out.growth_dist_ub = gd_ws_[i];
    out.status = dgd_.status_ws[i];
    err = dgd::ComputeSolutionError(set1, tf1_t, set2, tf2_t, out);
    res.solve_time = solve_time_avg;
    res.prim_dual_gap = err.prim_dual_gap;
    res.prim_feas_err = err.prim_feas_err;
    res.dual_feas_err = err.dual_feas_err;
    res.iter = iter_ws_[i];
    res.optimal_flag = (out.status == dgd::SolutionStatus::CoincidentCenters) ||
                       (out.status == dgd::SolutionStatus::Optimal);
    resa.AddResult(res);
  }
}

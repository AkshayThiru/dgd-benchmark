#ifndef DGD_BENCHMARK_INC_SOLUTION_ERROR_H_
#define DGD_BENCHMARK_INC_SOLUTION_ERROR_H_

#include <cmath>

#include "dgd/geometry/convex_set.h"
#include "dgd_benchmark/inc/data_types.h"
#include "dgd_benchmark/inc/incremental.h"

namespace inc {

// Solution error.
struct SolutionError {
  double prim_dual_gap;
  double prim_feas_err = 0.0;
  double dual_feas_err;
};

// Computes solution error.
inline SolutionError ComputeSolutionError(const dgd::ConvexSet<3>* set1,
                                          const Transform3& tf1,
                                          const dgd::ConvexSet<3>* set2,
                                          const Transform3& tf2,
                                          const Output& out) {
  SolutionError err{};
  err.dual_feas_err = 0.0;
  if (out.status == SolutionStatus::CoincidentCenters) {
    err.prim_dual_gap = 0.0;
    return err;
  } else if (out.status != SolutionStatus::Optimal) {
    err.prim_dual_gap = kInf;
    return err;
  }

  const Vec3 p1 = tf1.block<3, 1>(0, 3), p2 = tf2.block<3, 1>(0, 3);
  const Rotation3 rot1 = tf1.block<3, 3>(0, 0), rot2 = tf2.block<3, 3>(0, 0);

  Vec3 sp;
  const double sv1 = set1->SupportFunction(rot1.transpose() * out.normal, sp);
  const double sv2 = set2->SupportFunction(-rot2.transpose() * out.normal, sp);
  const double lb = (p2 - p1).dot(out.normal) / (sv1 + sv2);

  err.prim_dual_gap = std::abs(out.growth_dist / lb - 1.0);
  return err;
}

}  // namespace inc

#endif  // DGD_BENCHMARK_INC_SOLUTION_ERROR_H_

#ifndef DGD_BENCHMARK_DCF_DCF_COLLISION_H_
#define DGD_BENCHMARK_DCF_DCF_COLLISION_H_

#include <cstdint>

#include "dcf/dsf.h"
#include "dcf/precompiled.h"
#include "dcf/utils.h"

namespace dcf {

// Solver settings.
struct Settings {
  double min_center_dist = kSqrtEps;
  double tol = kSqrtEps;  // Relative tolerance.
  int max_iter = 50;
  int ie_iter = 4;
};

// Differentiable contact feature.
struct DCF {
  Mat<3, 12> dnormal;
  Mat<3, 12> dp1, dp2;
  Mat<1, 12> dgap;
  // normal points along pos2 - pos1.
  Vec3 normal;
  Vec3 p1, p2;
  double gap;
};

// Solution status.
enum class SolutionStatus {
  Optimal,
  MaxIterReached,
  CoincidentCenters,
};

// Solver output.
struct Output {
  DCF dcf;
  int iter = 0;
  SolutionStatus status = SolutionStatus::MaxIterReached;
};

// Computes the growth distance.
double GrowthDistance(DSF* dsf1, const Transform3& tf1, DSF* dsf2,
                      const Transform3& tf2, const Settings& settings,
                      Output& out);

inline double GrowthDistance(DSF* dsf1, const Vec<7>& tfq1, DSF* dsf2,
                             const Vec<7>& tfq2, const Settings& settings,
                             Output& out) {
  Transform3 tf1, tf2;
  Quaternion2RotationSe3(tfq1, tf1);
  Quaternion2RotationSe3(tfq2, tf2);
  return GrowthDistance(dsf1, tf1, dsf2, tf2, settings, out);
}

// Solution error.
struct SolutionError {
  double prim_dual_gap;
  double prim_infeas_err;
  double dual_infeas_err = 0.0;
};

SolutionError ComputeSolutionError(DSF* dsf1, const Transform3& tf1, DSF* dsf2,
                                   const Transform3& tf2, Output& out);

}  // namespace dcf

#endif  // DGD_BENCHMARK_DCF_DCF_COLLISION_H_

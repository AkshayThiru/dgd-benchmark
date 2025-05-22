#ifndef DSF_DSF_COLLISION_H_
#define DSF_DSF_COLLISION_H_

#include <cstdint>

#include "dsf/dsf.h"
#include "dsf/precompiled.h"
#include "dsf/utils.h"

namespace dsf {

// Solver settings.
struct Settings {
  int max_iter{50};
  double min_center_dist{kEpsSqrt};
  double tol{kEpsSqrt};
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
enum class SolutionStatus : uint8_t {
  Optimal,
  MaxIterReached,
  CoincidentCenters,
};

// Solver output.
struct Output {
  DCF dcf;
  int iter{0};
  SolutionStatus status{SolutionStatus::MaxIterReached};
};

// Computes the growth distance.
double GrowthDistance(DSF* dsf1, const Transform3& tf1, DSF* dsf2,
                      const Transform3& tf2, Output& out,
                      const Settings& settings);

inline double GrowthDistance(DSF* dsf1, const Vec<7>& tfq1, DSF* dsf2,
                             const Vec<7>& tfq2, Output& out,
                             const Settings& settings) {
  Transform3 tf1, tf2;
  Quaternion2RotationSe3(tfq1, tf1);
  Quaternion2RotationSe3(tfq2, tf2);
  return GrowthDistance(dsf1, tf1, dsf2, tf2, out, settings);
}

// Solution error.
struct SolutionError {
  double prim_dual_gap;
  double prim_feas_err;
  double dual_feas_err{0.0};
};

SolutionError ComputeSolutionError(DSF* dsf1, const Transform3& tf1, DSF* dsf2,
                                   const Transform3& tf2, Output& out);

}  // namespace dsf

#endif  // DSF_DSF_COLLISION_H_

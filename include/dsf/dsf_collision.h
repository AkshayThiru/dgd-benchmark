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
  Matrix<double, 3, 12> dnormal;
  Matrix<double, 3, 12> dp1, dp2;
  Matrix<double, 1, 12> dgap;
  Vector3d normal;
  Vector3d p1, p2;
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
  SolutionStatus status;
};

// Computes the growth distance.
double GrowthDistance(DSF* dsf1, const Matrix4d& tf1, DSF* dsf2,
                      const Matrix4d& tf2, Output& out,
                      const Settings& settings);

inline double GrowthDistance(DSF* dsf1, const Vector<double, 7>& tfq1,
                             DSF* dsf2, const Vector<double, 7>& tfq2,
                             Output& out, const Settings& settings) {
  Matrix4d tf1, tf2;
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

SolutionError ComputeSolutionError(DSF* dsf1, const Matrix4d& tf1, DSF* dsf2,
                                   const Matrix4d& tf2, Output& out);

}  // namespace dsf

#endif  // DSF_DSF_COLLISION_H_

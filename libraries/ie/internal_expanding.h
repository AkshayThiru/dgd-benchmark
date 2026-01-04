#ifndef DGD_BENCHMARK_IE_INTERNAL_EXPANDING_H_
#define DGD_BENCHMARK_IE_INTERNAL_EXPANDING_H_

#include <Eigen/LU>
#include <cassert>
#include <cmath>
#include <cstdint>

#include "dgd/data_types.h"
#include "dgd/geometry/convex_set.h"

namespace ie {

// Constants.
using dgd::kSqrtEps;

// Types.
using dgd::Matr;
using dgd::Vec3r;
using dgd::Vecr;

using dgd::Affine;
using dgd::Linear;
using dgd::Rotation3r;
using dgd::Transform3r;

// Other types.
using dgd::ConvexSet;
using dgd::SupportFunctionHint;
using dgd::detail::ConvexSetValidator;

// Solver settings.
struct Settings {
  double min_center_dist = kSqrtEps;
  double tol = kSqrtEps;
  double min_simplex_size = 1e-3;
  int max_iter = 100;
};

// Solution status.
enum class SolutionStatus {
  Optimal,
  MaxIterReached,
  CoincidentCenters,
};

// Solver output.
struct Output {
  Matr<3, 4> s1_ = Matr<3, 4>::Zero();
  Matr<3, 4> s2_ = Matr<3, 4>::Zero();
  Vecr<4> bc_ = Vecr<4>::Zero();

  SupportFunctionHint<3> hint1_{};
  SupportFunctionHint<3> hint2_{};

  // normal points along p2 - p1.
  Vec3r normal = Vec3r::Zero();
  Vec3r z1 = Vec3r::Zero();
  Vec3r z2 = Vec3r::Zero();
  double growth_dist = 0.0;

  int iter = 0;
  SolutionStatus status = SolutionStatus::MaxIterReached;
};

// Growth distance function.
template <class C1, class C2>
double GrowthDistance(const C1* set1, const Transform3r& tf1, const C2* set2,
                      const Transform3r& tf2, const Settings& settings,
                      Output& out) {
  static_assert(ConvexSetValidator<3, C1>::valid, "Incompatible set C1");
  static_assert(ConvexSetValidator<3, C2>::valid, "Incompatible set C2");

  out.hint2_.n_prev = out.hint1_.n_prev = Vec3r::Zero();
  out.iter = 0;

  // Check center distance.
  const Vec3r p1 = Affine(tf1), p2 = Affine(tf2);
  const Vec3r p12 = p1 - p2;
  const double cdist = p12.norm();
  if (cdist < settings.min_center_dist) {
    out.normal = Vec3r::Zero();
    out.z1 = p1;
    out.z2 = p2;
    out.growth_dist = 0.0;
    out.status = SolutionStatus::CoincidentCenters;
    return 0.0;
  }

  const Rotation3r rot1 = Linear(tf1), rot2 = Linear(tf2);

  // Simplex variables.
  const Vec3r u0 = Vec3r::Ones();
  Matr<3, 3> W, Winv;
  Vec3r cone_coeff;
  Vec3r sp1, sp2, sp;

  Matr<3, 4> V;
#ifdef IE_SMALL_SIMPLEX
  const double simplex_size = settings.min_simplex_size;
#else
  const double simplex_size = set1->inradius() + set2->inradius();
#endif
  V.col(0) = simplex_size * Vec3r(0.8, 0.0, -0.5);
  V.col(1) = simplex_size * Vec3r(-0.5, 0.5, -0.5);
  V.col(2) = simplex_size * Vec3r(-0.5, -0.5, -0.5);
  V.col(3) = simplex_size * Vec3r(0.0, 0.0, 1.0);

  // Internal expanding algorithm.
  int idx_prev = -1;
  double gd = 0.0;
  while (true) {
    int idx = 0;
    for (; idx < 4; ++idx) {
      if (idx == idx_prev) continue;
      W.col(0) = V.col((idx > 0) ? 0 : 1);
      W.col(1) = V.col((idx > 1) ? 1 : 2);
      W.col(2) = V.col((idx > 2) ? 2 : 3);
      Winv = W.inverse();
      cone_coeff = -Winv * p12;
      if (cone_coeff.minCoeff() >= 0) break;
    }
    assert(idx < 4);
    out.normal = Winv.transpose() * u0;
    out.normal.normalize();

    set1->SupportFunction(rot1.transpose() * out.normal, sp1, &out.hint1_);
    set2->SupportFunction(-rot2.transpose() * out.normal, sp2, &out.hint2_);
    sp.noalias() = rot1 * sp1 - rot2 * sp2;
    ++out.iter;

    V.col(idx) = sp;
    out.s1_.col(idx) = sp1;
    out.s2_.col(idx) = sp2;
    gd = cone_coeff.sum();
    cone_coeff /= gd;
    out.bc_((idx > 0) ? 0 : 1) = cone_coeff(0);
    out.bc_((idx > 1) ? 1 : 2) = cone_coeff(1);
    out.bc_((idx > 2) ? 2 : 3) = cone_coeff(2);
    out.bc_(idx) = 0.0;

#ifdef IE_DISTANCE_TC
    if (out.normal.dot(sp - W.col(0)) < settings.tol) {
      out.status = SolutionStatus::Optimal;
      break;
    }
#else
    const double lb = -out.normal.dot(p12) / out.normal.dot(sp);
    if (gd < (1.0 + settings.tol) * lb) {
      out.status = SolutionStatus::Optimal;
      break;
    }
#endif
    if (out.iter >= settings.max_iter) {
      out.status = SolutionStatus::MaxIterReached;
      break;
    }

    idx_prev = idx;
  }

  out.growth_dist = gd;
  out.z1 = p1 + rot1 * out.s1_ * out.bc_;
  out.z2 = p2 + rot2 * out.s2_ * out.bc_;
  return out.growth_dist;
}

// Solution error.
struct SolutionError {
  double prim_dual_gap;
  double prim_infeas_err;
  double dual_infeas_err = 0.0;
};

// Computes solution error.
inline SolutionError ComputeSolutionError(const ConvexSet<3>* set1,
                                          const Transform3r& tf1,
                                          const ConvexSet<3>* set2,
                                          const Transform3r& tf2,
                                          const Output& out) {
  SolutionError err{};
  if (out.status == SolutionStatus::CoincidentCenters) {
    err.prim_dual_gap = err.prim_infeas_err = 0.0;
    return err;
  }

  const Vec3r p1 = Affine(tf1), p2 = Affine(tf2);
  const Rotation3r rot1 = Linear(tf1), rot2 = Linear(tf2);
  const Vec3r cp1 = p1 + out.growth_dist * (out.z1 - p1);
  const Vec3r cp2 = p2 + out.growth_dist * (out.z2 - p2);

  Vec3r sp;
  const double sv1 = set1->SupportFunction(rot1.transpose() * out.normal, sp);
  const double sv2 = set2->SupportFunction(-rot2.transpose() * out.normal, sp);
  const double lb = (p2 - p1).dot(out.normal) / (sv1 + sv2);

  err.prim_dual_gap = std::abs(out.growth_dist / lb - 1.0);
  err.prim_infeas_err = (cp1 - cp2).norm() / (p1 - p2).norm();
  return err;
}

}  // namespace ie

#endif  // DGD_BENCHMARK_IE_INTERNAL_EXPANDING_H_

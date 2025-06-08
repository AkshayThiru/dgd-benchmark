#ifndef DGD_BENCHMARK_IE_INTERNAL_EXPANDING_H_
#define DGD_BENCHMARK_IE_INTERNAL_EXPANDING_H_

#include <cassert>
#include <cmath>
#include <cstdint>

#include "dgd/data_types.h"
#include "dgd/geometry/convex_set.h"

namespace ie {

// Solver settings.
struct Settings {
  int max_iter = 100;
  double min_center_dist = dgd::kSqrtEps;
  double tol = dgd::kSqrtEps;
  double min_simplex_size = 1e-3;
};

// Solution status.
enum class SolutionStatus : uint8_t {
  Optimal,
  MaxIterReached,
  CoincidentCenters,
};

// Solver output.
struct Output {
  dgd::Matr<3, 4> s1_ = dgd::Matr<3, 4>::Zero();
  dgd::Matr<3, 4> s2_ = dgd::Matr<3, 4>::Zero();
  dgd::Vecr<4> bc_ = dgd::Vecr<4>::Zero();

  dgd::SupportFunctionHint<3> hint1_{};
  dgd::SupportFunctionHint<3> hint2_{};

  // normal points along p2 - p1.
  dgd::Vec3r normal = dgd::Vec3r::Zero();
  dgd::Vec3r z1 = dgd::Vec3r::Zero();
  dgd::Vec3r z2 = dgd::Vec3r::Zero();
  double growth_dist = 0.0;

  int iter = 0;
  SolutionStatus status = SolutionStatus::MaxIterReached;
};

// Growth distance function.
template <class C1, class C2>
double GrowthDistance(const C1* set1, const dgd::Transform3r& tf1,
                      const C2* set2, const dgd::Transform3r& tf2,
                      const Settings& settings, Output& out) {
  static_assert((C1::dimension() == 3) && (C2::dimension() == 3),
                "Convex sets are not three-dimensional");

  out.hint2_.n_prev = out.hint1_.n_prev = dgd::Vec3r::Zero();
  out.iter = 0;

  // Check center distance.
  const dgd::Vec3r p1 = tf1.block<3, 1>(0, 3), p2 = tf2.block<3, 1>(0, 3);
  const dgd::Vec3r p12 = p1 - p2;
  const double cdist = p12.norm();
  if (cdist < settings.min_center_dist) {
    out.normal = dgd::Vec3r::Zero();
    out.z1 = p1;
    out.z2 = p2;
    out.growth_dist = 0.0;
    out.status = SolutionStatus::CoincidentCenters;
    return 0.0;
  }

  const dgd::Rotation3r rot1 = tf1.block<3, 3>(0, 0),
                        rot2 = tf2.block<3, 3>(0, 0);

  // Simplex variables.
  const dgd::Vec3r u0 = dgd::Vec3r::Ones();
  dgd::Matr<3, 3> W, Winv;
  dgd::Vec3r cone_coeff;
  dgd::Vec3r sp1, sp2, sp;

  dgd::Matr<3, 4> V;
#ifdef IE_SMALL_SIMPLEX
  const double simplex_size = settings.min_simplex_size;
#else
  const double simplex_size = set1->inradius() + set2->inradius();
#endif
  V.col(0) = simplex_size * dgd::Vec3r(0.8, 0.0, -0.5);
  V.col(1) = simplex_size * dgd::Vec3r(-0.5, 0.5, -0.5);
  V.col(2) = simplex_size * dgd::Vec3r(-0.5, -0.5, -0.5);
  V.col(3) = simplex_size * dgd::Vec3r(0.0, 0.0, 1.0);

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
  double prim_feas_err;
  double dual_feas_err = 0.0;
};

// Computes solution error.
inline SolutionError ComputeSolutionError(const dgd::ConvexSet<3>* set1,
                                          const dgd::Transform3r& tf1,
                                          const dgd::ConvexSet<3>* set2,
                                          const dgd::Transform3r& tf2,
                                          const Output& out) {
  SolutionError err{};
  if (out.status == SolutionStatus::CoincidentCenters) {
    err.prim_dual_gap = err.prim_feas_err = 0.0;
    return err;
  } else if (out.status != SolutionStatus::Optimal) {
    err.prim_dual_gap = err.prim_feas_err = dgd::kInf;
    return err;
  }

  const dgd::Vec3r p1 = tf1.block<3, 1>(0, 3), p2 = tf2.block<3, 1>(0, 3);
  const dgd::Rotation3r rot1 = tf1.block<3, 3>(0, 0),
                        rot2 = tf2.block<3, 3>(0, 0);
  const dgd::Vec3r cp1 = p1 + out.growth_dist * (out.z1 - p1);
  const dgd::Vec3r cp2 = p2 + out.growth_dist * (out.z2 - p2);

  dgd::Vec3r sp;
  const double sv1 = set1->SupportFunction(rot1.transpose() * out.normal, sp);
  const double sv2 = set2->SupportFunction(-rot2.transpose() * out.normal, sp);
  const double lb = (p2 - p1).dot(out.normal) / (sv1 + sv2);

  err.prim_dual_gap = std::abs(out.growth_dist / lb - 1.0);
  err.prim_feas_err = (cp1 - cp2).norm();
  return err;
}

}  // namespace ie

#endif  // DGD_BENCHMARK_IE_INTERNAL_EXPANDING_H_

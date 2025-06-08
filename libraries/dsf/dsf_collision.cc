#include "dsf/dsf_collision.h"

#include <cmath>

#include "dsf/dsf.h"
#include "dsf/precompiled.h"

namespace dsf {

namespace {

// Computes the differentiable contact feature.
inline void DcfImpDiff(DSF* dsf1, const Vec3& pos1, const Rotation3& R1,
                       DSF* dsf2, const Vec3& pos2, const Rotation3& R2,
                       const Vec4& var, DCF& dcf) {
  Vec3 s1, s2;
  Mat3 ds1dx, ds2dx;
  Mat<3, 6> ds1dq, ds2dq;
  dsf1->SupportFunction(var.head<3>(), pos1, R1, s1, ds1dx, ds1dq);
  dsf2->SupportFunction(-var.head<3>(), pos2, R2, s2, ds2dx, ds2dq);

  const Vec3 q_bar = pos1 - pos2;
  Mat4 dedvar;
  dedvar.block<3, 3>(0, 0) = var(3) * (ds1dx + ds2dx);
  dedvar.block<3, 1>(0, 3) = s1 - s2 - q_bar;
  dedvar.block<1, 3>(3, 0) = 2.0 * var.head<3>().transpose();
  dedvar(3, 3) = 0.0;
  Mat<4, 12> dedq = Mat<4, 12>::Zero();
  dedq.block<3, 6>(0, 0) = var(3) * ds1dq;
  dedq.block<3, 3>(0, 0) += (1.0 - var(3)) * Mat3::Identity();
  dedq.block<3, 6>(0, 6) = -var(3) * ds2dq;
  dedq.block<3, 3>(0, 6) -= (1.0 - var(3)) * Mat3::Identity();
  const Mat<3, 12> dxdq =
      -(dedvar.colPivHouseholderQr().solve(dedq).topRows<3>());

  dcf.p1 = s1;
  dcf.p2 = s2;
  dcf.dp1 = ds1dx * dxdq;
  dcf.dp1.block<3, 6>(0, 0) += ds1dq;
  dcf.dp2 = -ds2dx * dxdq;
  dcf.dp2.block<3, 6>(0, 6) += ds2dq;
  dcf.normal = var.head<3>();
  dcf.dnormal = (Mat3::Identity() - dcf.normal * dcf.normal.transpose()) * dxdq;
#ifdef DSF_SIGNED_DISTANCE
  // Signed distance.
  const double sign = var(3) > 1.0 ? 1.0 : -1.0;
  dcf.gap = sign * (s1 - s2).norm();
  dcf.dgap = sign * (s1 - s2).normalized().transpose() * (dcf.dp1 - dcf.dp2);
#else
  // Growth distance.
  const Vec3 sp = s1 - s2 - q_bar;
  const double spnorm = sp.norm(), q_barnorm = q_bar.norm();
  dcf.gap = q_barnorm / spnorm;
  dcf.dgap.setZero();
  dcf.dgap.leftCols<3>() =
      (q_bar / q_barnorm + dcf.gap / spnorm * sp).transpose() / spnorm;
  dcf.dgap.middleCols<3>(6) = -dcf.dgap.leftCols<3>();
  dcf.dgap -=
      dcf.gap * sp.transpose() * (dcf.dp1 - dcf.dp2) / (spnorm * spnorm);
#endif
}

}  // namespace

double GrowthDistance(DSF* dsf1, const Transform3& tf1, DSF* dsf2,
                      const Transform3& tf2, const Settings& settings,
                      Output& out) {
  const Vec3 pos1 = tf1.block<3, 1>(0, 3), pos2 = tf2.block<3, 1>(0, 3);
  const Rotation3 R1 = tf1.block<3, 3>(0, 0), R2 = tf2.block<3, 3>(0, 0);

  const Vec3 q_bar = pos1 - pos2;

  // Check if the centers are coincident.
  if (q_bar.norm() < settings.min_center_dist) {
    out.dcf.gap = 0.0;
    out.iter = 0;
    out.status = SolutionStatus::CoincidentCenters;
    return 0.0;
  }

  // Simplex-Newton hybrid algorithm.
  const Vec3 u0(1.0, 1.0, 1.0);
  Mat3 W_sim, Winv_sim;
  Vec3 c_sim, x;
  Mat<3, 4> V_sim;
  const double small = 1e-3;
  V_sim.col(0) = small * Vec3(1.0, 0.0, -0.5);
  V_sim.col(1) = small * Vec3(-0.5, 0.5, -0.5);
  V_sim.col(2) = small * Vec3(-0.5, -0.5, -0.5);
  V_sim.col(3) = small * Vec3(0.0, 0.0, 1.0);

  Vec3 s1, s2;
  Mat3 ds1dx, ds2dx;

  Vec4 var, dir, e;
  Mat4 dedvar;
  dedvar(3, 3) = 0.0;

  out.iter = 0;
  while (true) {
    // Simplex intersection.
    for (int i = 0; i < 4; ++i) {
      W_sim.col(0) = V_sim.col((i > 0) ? 0 : 1);
      W_sim.col(1) = V_sim.col((i > 1) ? 1 : 2);
      W_sim.col(2) = V_sim.col((i > 2) ? 2 : 3);
      Winv_sim = W_sim.inverse();
      c_sim = -Winv_sim * q_bar;
      if (c_sim.minCoeff() >= 0.0) break;
    }
    V_sim.block<3, 3>(0, 0) = W_sim;

    if (out.iter <= 4) {
      x = Winv_sim.transpose() * u0;
      if (out.iter == 4) {
        dsf1->SupportFunction(x, pos1, R1, s1, ds1dx);
        dsf2->SupportFunction(-x, pos2, R2, s2, ds2dx);
      } else {
        dsf1->SupportFunction(x, pos1, R1, s1);
        dsf2->SupportFunction(-x, pos2, R2, s2);
      }
      V_sim.col(3) = s1 - s2 - q_bar;
    } else {
      // Newton step.
      var.head<3>() = x.normalized();
      var(3) = q_bar.norm() / (s1 - s2 - q_bar).norm();
      e.head<3>() = q_bar + var(3) * (s1 - s2 - q_bar);
      e(3) = 0.0;
      // COUT_SCALAR(e.norm());
      if (e.norm() < settings.tol) {
        out.status = SolutionStatus::Optimal;
        break;
      }

      dedvar.topLeftCorner<3, 3>() = var(3) * (ds1dx + ds2dx);
      dedvar.block<3, 1>(0, 3) = s1 - pos1 - (s2 - pos2);
      dedvar.block<1, 3>(3, 0) = 2.0 * var.head<3>().transpose();
      dir = -dedvar.colPivHouseholderQr().solve(e);
      x = var.head<3>() + dir.head<3>();
      dsf1->SupportFunction(x, pos1, R1, s1, ds1dx);
      dsf2->SupportFunction(-x, pos2, R2, s2, ds2dx);

      // Safety.
      c_sim = Winv_sim * (s1 - s2 - q_bar);
      if (c_sim.minCoeff() < 0.0) {
        x = Winv_sim.transpose() * u0;
        dsf1->SupportFunction(x, pos1, R1, s1, ds1dx);
        dsf2->SupportFunction(-x, pos2, R2, s2, ds2dx);
      }
      V_sim.col(3) = s1 - s2 - q_bar;
    }

    out.iter++;
    if (out.iter > settings.max_iter) {
      out.status = SolutionStatus::MaxIterReached;
      break;
    }
  }

  var.head<3>() = x.normalized();
  var(3) = q_bar.norm() / (s1 - s2 - q_bar).norm();
  DcfImpDiff(dsf1, pos1, R1, dsf2, pos2, R2, var, out.dcf);
  return out.dcf.gap;
}

SolutionError ComputeSolutionError(DSF* dsf1, const Transform3& tf1, DSF* dsf2,
                                   const Transform3& tf2, Output& out) {
  SolutionError err = SolutionError{};
  if (out.status == SolutionStatus::CoincidentCenters) {
    err.prim_dual_gap = err.prim_feas_err = 0.0;
    return err;
  } else if (out.status != SolutionStatus::Optimal) {
    err.prim_dual_gap = err.prim_feas_err = kInf;
    return err;
  }

  const Vec3 pos1 = tf1.block<3, 1>(0, 3), pos2 = tf2.block<3, 1>(0, 3);
  const Rotation3 R1 = tf1.block<3, 3>(0, 0), R2 = tf2.block<3, 3>(0, 0);
  const DCF& dcf = out.dcf;
  const Vec3 normal = dcf.normal, q_bar = pos1 - pos2;

  Vec3 s1, s2;
  dsf1->SupportFunction(normal, pos1, R1, s1);
  dsf2->SupportFunction(-normal, pos2, R2, s2);
  const double sv1 = (s1 - pos1).dot(normal), sv2 = (s2 - pos2).dot(-normal);
  const double gap_prim = q_bar.norm() / (s1 - s2 - q_bar).norm();
  const double gap_dual = -q_bar.dot(normal) / (sv1 + sv2);

  err.prim_dual_gap = std::abs(gap_prim / gap_dual - 1.0);
  err.prim_feas_err = (q_bar + dcf.gap * (dcf.p1 - dcf.p2 - q_bar)).norm();
  return err;
}

}  // namespace dsf

#include "inc/incremental.h"

#include <array>
#include <cassert>
#include <iostream>

#include "inc/data_types.h"
#include "inc/polyhedron.h"

#define PRINT_SUBALGORITHM(ep) std::cout << ep << std::endl;
#define PRINT_GROWTH_DISTANCE std::cout << "sigma: " << res_(3) << std::endl;

namespace inc {

namespace {

inline int OtherFacetIdx(const std::array<int, 3>& fids, int fid, int& ofid) {
  if (fids[1] == fid) {
    ofid = fids[0];
    return 1;
  } else {
    ofid = fids[1];
    return 0;
  }
}

inline int OtherFacetIdx(const std::array<int, 3>& fids, int fid, int& ofid1,
                         int& ofid2) {
  if (fids[2] == fid) {
    ofid1 = fids[0];
    ofid2 = fids[1];
    return 2;
  } else {
    ofid2 = fids[2];
    return OtherFacetIdx(fids, fid, ofid1);
  }
}

inline int BestIdx(const Polyhedron* set, const Vec3& ray,
                   const std::array<int, 3>& fids, int nact) {
  int fid = fids[0];

  if (nact > 1) {
    double val, max = set->Gamma()[fid].dot(ray);
    if ((val = set->Gamma()[fids[1]].dot(ray)) > max) {
      max = val;
      fid = fids[1];
    }

    if (nact > 2) {
      if ((val = set->Gamma()[fids[2]].dot(ray)) > max) {
        fid = fids[2];
      }
    }
  }
  return fid;
}

inline std::ostream& operator<<(std::ostream& os, EntryPoint ep) {
  switch (ep) {
    case EntryPoint::IZ:
      os << "[IZ ]";
      break;
    case EntryPoint::FFD:
      os << "[FFD]";
      break;
    case EntryPoint::FED:
      os << "[FED]";
      break;
    case EntryPoint::FD:
      os << "[FD ]";
  }
  return os;
}

}  // namespace

// Initial state:
//   idx_sf_[0][0], idx_sf_[1][0] (= 0 for cold start),
//    0 <= idx_sf_[s][0] < sets[s]->nfacet().
inline SolutionStatus GrowthDistanceSolver::Initialize(
    const Polyhedron* sets[2], int& iter) {
  // PRINT_SUBALGORITHM(EntryPoint::IZ)
  const Vec3 p12 = p_[0] - p_[1];

  const double rv1 = sets[0]->RayCast(-R_[0].transpose() * p12, idx_sf_[0][0]);
  const double rv2 = sets[1]->RayCast(R_[1].transpose() * p12, idx_sf_[1][0]);
  nact_[0] = nact_[1] = 1;
  idx_lp_[0] = {idx_sf_[0][0], 0};
  idx_lp_[1] = {idx_sf_[1][0], 1};

  res_(3) = (rv1 * rv2) / (rv1 + rv2);
  res_.head<3>() = p_[1] + res_(3) / rv2 * p12;
  ++iter;
  // PRINT_GROWTH_DISTANCE
  return SolutionStatus::NoError_;
}

// Initial state:
//   idx_sf_[0][0], idx_sf_[1][0], idx_lp_[0], idx_lp_[1], nact_ = {1, 1}, res_.
inline SolutionStatus GrowthDistanceSolver::FaceFaceDescent(
    const Polyhedron* sets[2], int& iter) {
  // PRINT_SUBALGORITHM(EntryPoint::FFD)
  A_.block<3, 1>(0, 0) = R_[0] * sets[0]->Gamma()[idx_sf_[0][0]];
  A_.block<3, 1>(0, 1) = R_[1] * sets[1]->Gamma()[idx_sf_[1][0]];
  if (!ls_.ComputeFaceFaceGradient(A_, dx_)) {
    idx_lp_[2] = idx_lp_[3] = {-1, -1};
    return SolutionStatus::Optimal;
  }

  double step[2];
  int idx_f[2];
  UpdateLocalXDx();
  for (int s = 0; s < 2; ++s) {
    step[s] = sets[s]->LocalFaceStepSize(Rtx_p_[s], Rtdx_[s], res_(3),
                                         idx_sf_[s][0], idx_f[s]);
  }

  const int si = (step[0] > step[1]);  // Set with the incoming facet.
  idx_sf_[si][1] = idx_f[si];
  ++nact_[si];
  idx_lp_[2] = {idx_f[si], si};

  res_.head<3>() += step[si] * dx_;
  res_(3) -= step[si];
  ++iter;
  // PRINT_GROWTH_DISTANCE
  return SolutionStatus::NoError_;
}

// Initial state:
//   idx_sf_[0][0 ... nact_[0]-1], idx_sf_[1][0 ... nact_[1]-1],
//   nact_ (nact_[0] + nact_[1] = 3), idx_lp_[0 ... 2], A_.leftCols(2), res_.
inline SolutionStatus GrowthDistanceSolver::FaceEdgeDescent(
    const Polyhedron* sets[2], int& iter) {
  // PRINT_SUBALGORITHM(EntryPoint::FED);
  const int sn = idx_lp_[2][1];  // The set corresponding to the new facet.
  const int snc = 1 - sn;        // The other set.
  A_.block<3, 1>(0, 2) = R_[sn] * sets[sn]->Gamma()[idx_lp_[2][0]];
  if (!ls_.ComputeFaceEdgeGradient(A_, dx_)) {
    idx_lp_[3] = {-1, -1};
    return SolutionStatus::Optimal;
  }

  double step[2];
  int idx_f[2];
  UpdateLocalXDx();
  step[sn] =
      sets[sn]->LocalEdgeStepSize(Rtx_p_[sn], Rtdx_[sn], res_(3),
                                  idx_sf_[sn][0], idx_sf_[sn][1], idx_f[sn]);
  step[snc] = sets[snc]->LocalFaceStepSize(Rtx_p_[snc], Rtdx_[snc], res_(3),
                                           idx_sf_[snc][0], idx_f[snc]);

  const int si = (step[0] > step[1]);  // Set with the incoming facet.
  A_.block<3, 1>(0, 3) = R_[si] * sets[si]->Gamma()[idx_f[si]];
  idx_sf_[si][nact_[si]++] = idx_f[si];
  idx_lp_[3] = {idx_f[si], si};

  res_.head<3>() += step[si] * dx_;
  res_(3) -= step[si];
  ++iter;
  // PRINT_GROWTH_DISTANCE
  return SolutionStatus::NoError_;
}

// Initial state:
//   idx_sf_[0][0 ... nact_[0]-1], idx_sf_[1][0 ... nact_[1]-1],
//   nact_ (nact_[0] + nact_[1] = 4), idx_lp_, A_, res_.
inline SolutionStatus GrowthDistanceSolver::FinalDescent(
    const Polyhedron* sets[2], int& iter) {
  // PRINT_SUBALGORITHM(EntryPoint::FD);
  double step[2];
  int idx_f[2];
  int lpo;  // Outgoing index from the LP active constraint set (0 <= lpo < 4).
  // Sets corresponding to the outgoing and incoming active facets.
  int so, soc, si;

  do {
    if (!ls_.ComputeFinalDual(A_, lambda_)) {
      return SolutionStatus::SingularityError;
    }
    const double lambda_min = lambda_.minCoeff(&lpo);
    if (lambda_min >= settings_.lambda_min) {
      return SolutionStatus::Optimal;
    }
    ls_.ComputeFinalGradient(lpo, lambda_min, dx_);

    so = idx_lp_[lpo][1];
    if (nact_[so] == 1) {
      step[0] = sets[so]->GlobalStepSize(
          R_[so].transpose() * (res_.head<3>() - p_[so]),
          R_[so].transpose() * dx_, res_(3), idx_lp_[lpo][0], idx_f[0]);

      A_.block<3, 1>(0, lpo) = R_[so] * sets[so]->Gamma()[idx_f[0]];
      idx_sf_[so][0] = idx_f[0];
      idx_lp_[lpo][0] = idx_f[0];

      res_.head<3>() += step[0] * dx_;
      res_(3) -= step[0];
    } else {
      soc = 1 - so;
      UpdateLocalXDx();
      assert(((nact_[so] == 2) && (nact_[soc] == 2)) ||
             ((nact_[so] == 3) && (nact_[soc] == 1)));
      int so_lpo_idx, f_other1, f_other2;
      if (nact_[so] == 2) {
        so_lpo_idx = OtherFacetIdx(idx_sf_[so], idx_lp_[lpo][0], f_other1);
        step[so] = sets[so]->LocalFaceStepSize(Rtx_p_[so], Rtdx_[so], res_(3),
                                               f_other1, idx_f[so]);
        step[soc] = sets[soc]->LocalEdgeStepSize(Rtx_p_[soc], Rtdx_[soc],
                                                 res_(3), idx_sf_[soc][0],
                                                 idx_sf_[soc][1], idx_f[soc]);
      } else {
        so_lpo_idx =
            OtherFacetIdx(idx_sf_[so], idx_lp_[lpo][0], f_other1, f_other2);
        step[so] = sets[so]->LocalEdgeStepSize(Rtx_p_[so], Rtdx_[so], res_(3),
                                               f_other1, f_other2, idx_f[so]);
        step[soc] = sets[soc]->LocalFaceStepSize(
            Rtx_p_[soc], Rtdx_[soc], res_(3), idx_sf_[soc][0], idx_f[soc]);
      }
      si = (step[0] > step[1]);
      A_.block<3, 1>(0, lpo) = R_[si] * sets[si]->Gamma()[idx_f[si]];
      if (si == so) {
        idx_sf_[so][so_lpo_idx] = idx_f[so];
      } else {
        idx_sf_[so][so_lpo_idx] = idx_sf_[so][--nact_[so]];
        idx_sf_[soc][nact_[soc]++] = idx_f[soc];
      }
      idx_lp_[lpo] = {idx_f[si], si};

      res_.head<3>() += step[si] * dx_;
      res_(3) -= step[si];
    }
    assert((nact_[0] > 0) && (nact_[1] > 0));
    ++iter;
    // PRINT_GROWTH_DISTANCE

    if (iter >= settings_.max_iter) {
      return SolutionStatus::MaxIterReached;
    }
  } while (true);
}

inline bool GrowthDistanceSolver::FeasibilityCheck(
    const Polyhedron* sets[2], std::array<bool, 2>& feas) const {
  if (nact_[0] == 1) {
    feas[0] = sets[0]->LocalFaceFeasible(Rtx_p_[0], res_(3), idx_sf_[0][0]);
    feas[1] = true;
  } else if (nact_[0] == 2) {
    feas[0] = sets[0]->LocalEdgeFeasible(Rtx_p_[0], res_(3), idx_sf_[0][0],
                                         idx_sf_[0][1]);
    feas[1] = sets[1]->LocalEdgeFeasible(Rtx_p_[1], res_(3), idx_sf_[1][0],
                                         idx_sf_[1][1]);
  } else {
    feas[0] = true;
    feas[1] = sets[1]->LocalFaceFeasible(Rtx_p_[1], res_(3), idx_sf_[1][0]);
  }
  return feas[0] & feas[1];
}

inline void GrowthDistanceSolver::FeasibilityRecovery(
    const Polyhedron* sets[2], std::array<bool, 2>& feas) {
  double sigma[2] = {res_(3), res_(3)};
  int bfid[2] = {-1, -1};
  for (int s = 0; s < 2; ++s) {
    if (!feas[s]) {
      bfid[s] = BestIdx(sets[s], Rtx_p_[s], idx_sf_[s], nact_[s]);
      sigma[s] = sets[s]->RayCast(Rtx_p_[s], bfid[s]);
    }
  }

  const int sgc = (sigma[0] < sigma[1]);
  const int sg = 1 - sgc;  // Set for which GlobalStepSize will be called.
  assert(bfid[sgc] != -1);
  assert(sigma[sgc] > res_(3));
  assert(nact_[sgc] < 3);
  res_(3) = sigma[sgc];
  dx_ = p_[sgc] - res_.head<3>();
  const double mdsigma =
      -sets[sgc]->Gamma()[bfid[sgc]].dot(R_[sgc].transpose() * dx_);
  dx_ /= mdsigma;

  const double step = sets[sg]->GlobalStepSize(
      Rtx_p_[sg], R_[sg].transpose() * dx_, res_(3), idx_sf_[sg][0], bfid[sg]);

  idx_sf_[0][0] = bfid[0];
  idx_sf_[1][0] = bfid[1];
  idx_lp_[0] = {bfid[0], 0};
  idx_lp_[1] = {bfid[1], 1};
  nact_[0] = nact_[1] = 1;

  res_.head<3>() += step * dx_;
  res_(3) -= step;
}

inline EntryPoint GrowthDistanceSolver::Incremental(const Polyhedron* sets[2],
                                                    int& iter) {
  assert((nact_[0] > 0) && (nact_[1] > 0));
  const Vec3 p12 = p_[0] - p_[1];
  const int nact = nact_[0] + nact_[1];

  auto best_idx = [this, &sets, &p12]() -> void {
    idx_sf_[0][0] =
        BestIdx(sets[0], -R_[0].transpose() * p12, idx_sf_[0], nact_[0]);
    idx_sf_[1][0] =
        BestIdx(sets[1], R_[1].transpose() * p12, idx_sf_[1], nact_[1]);
  };

  if (nact < 4) {
    best_idx();
    return EntryPoint::IZ;
  } else {
    Vec4 b;
    int s;
    for (int i = 0; i < 4; ++i) {
      s = idx_lp_[i][1];
      A_.block<3, 1>(0, i) = R_[s] * sets[s]->Gamma()[idx_lp_[i][0]];
      b(i) = A_.block<3, 1>(0, i).dot(p_[s]);
    }
    if (!ls_.ComputeWarmStartSolution(A_, b, res_) || res_(3) < 0) {
      best_idx();
      return EntryPoint::IZ;
    }

    UpdateLocalX();
    std::array<bool, 2> feas;
    if (FeasibilityCheck(sets, feas)) return EntryPoint::FD;

    ++iter;
    FeasibilityRecovery(sets, feas);
    return EntryPoint::FFD;
  }
}

inline void GrowthDistanceSolver::ComputeDualOptimalSolution(
    const Polyhedron* sets[2]) {
  if (nact_[0] == 1) {
    normal_ = R_[0] * sets[0]->Gamma()[idx_sf_[0][0]].normalized();
  } else if (nact_[1] == 1) {
    normal_ = -R_[1] * sets[1]->Gamma()[idx_sf_[1][0]].normalized();
  } else {
    Vec4 lambda1 = Vec4::Zero();
    for (int i = 0; i < 4; ++i) {
      if (idx_lp_[i][1] == 0) lambda1(i) = lambda_(i);
    }
    normal_ = (A_.topRows<3>() * lambda1).normalized();
  }
}

double GrowthDistanceSolver::GrowthDistance(const Polyhedron* set1,
                                            const Transform3& tf1,
                                            const Polyhedron* set2,
                                            const Transform3& tf2, Output& out,
                                            bool warm_start) {
  const Polyhedron* sets[2] = {set1, set2};
  p_[0] = tf1.block<3, 1>(0, 3);
  p_[1] = tf2.block<3, 1>(0, 3);
  R_[0] = tf1.block<3, 3>(0, 0);
  R_[1] = tf2.block<3, 3>(0, 0);
  out.iter = 0;

  // Check center distance.
  if ((p_[0] - p_[1]).norm() < settings_.min_center_dist) {
    out.normal = normal_ = Vec3::Zero();
    out.z1 = res_.head<3>() = p_[0];
    out.z2 = p_[1];
    out.status = status_ = SolutionStatus::CoincidentCenters;
    return (out.growth_dist = res_(3) = 0.0);
  }

  EntryPoint ep{EntryPoint::IZ};
  if (warm_start && (status_ == SolutionStatus::Optimal)) {
    ep = Incremental(sets, out.iter);
  } else {
    idx_sf_[0][0] = idx_sf_[1][0] = 0;
  }

  switch (ep) {
    case EntryPoint::IZ:
      out.status = status_ = Initialize(sets, out.iter);
      [[fallthrough]];

    case EntryPoint::FFD:
      out.status = status_ = FaceFaceDescent(sets, out.iter);
      if (status_ == SolutionStatus::Optimal) break;
      [[fallthrough]];

    case EntryPoint::FED:
      out.status = status_ = FaceEdgeDescent(sets, out.iter);
      if (status_ == SolutionStatus::Optimal) break;
      [[fallthrough]];

    case EntryPoint::FD:
      out.status = status_ = FinalDescent(sets, out.iter);
  }

  switch (status_) {
    case SolutionStatus::Optimal:
      ComputeDualOptimalSolution(sets);
      out.normal = normal_;
      [[fallthrough]];

    case SolutionStatus::MaxIterReached:
      // assert(res_(3) >= 0.0);
      out.z1 = p_[0] + (res_.head<3>() - p_[0]) / res_(3);
      out.z2 = p_[1] + (res_.head<3>() - p_[1]) / res_(3);
      return (out.growth_dist = res_(3));

    default:
      return -1.0;
  }
}

}  // namespace inc

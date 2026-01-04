#ifndef DGD_BENCHMARK_INC_INCREMENTAL_H_
#define DGD_BENCHMARK_INC_INCREMENTAL_H_

#include <Eigen/Dense>
#include <array>

#include "inc/data_types.h"
#include "inc/linear_solver.h"

namespace inc {

class Polyhedron;

template <int outer, int inner>
using MdArray = std::array<std::array<int, inner>, outer>;

const double kLsTol = kEps;  // Linear solver tolerance.

// Solver settings.
struct Settings {
  double min_center_dist = kSqrtEps;
  double lambda_min = -kEps;
  int max_iter = 1000;
  bool lp_act_cons_inv = false;
};

// Algorithm entry point.
enum class EntryPoint {
  IZ,   // Initialization.
  FFD,  // Face-face descent.
  FED,  // Face-edge descent.
  FD,   // Final descent.
};

// Solution status.
enum class SolutionStatus {
  NoError_,  // Internal.
  Optimal,
  CoincidentCenters,
  SingularityError,
  MaxIterReached,
};

// Solver output.
struct Output {
  // normal points along p2 - p1.
  Vec3 normal = Vec3::Zero();
  Vec3 z1 = Vec3::Zero();
  Vec3 z2 = Vec3::Zero();
  double growth_dist = 0.0;

  int iter = 0;
  SolutionStatus status = SolutionStatus::MaxIterReached;
};

// Growth distance solver.
class GrowthDistanceSolver {
 public:
  explicit GrowthDistanceSolver();

  ~GrowthDistanceSolver() = default;

  double GrowthDistance(const Polyhedron* set1, const Transform3& tf1,
                        const Polyhedron* set2, const Transform3& tf2,
                        Output& out, bool warm_start = false);

  Settings& SolverSettings();

 private:
  // Linear system variables.
  LinearSolver ls_;
  Mat4 A_;
  Vec4 lambda_;
  Vec3 dx_;

  // Settings.
  Settings settings_;

  // Rigid body transformation variables.
  Mat3 R_[2];
  Vec3 p_[2];
  Vec3 Rtx_p_[2];
  Vec3 Rtdx_[2];

  // Results.
  Vec4 res_;
  // Normal vector (dual optimal solution).
  Vec3 normal_;
  // For each set, the list of active facet indices (at most 3).
  MdArray<2, 3> idx_sf_;
  // For each set, the number of active facets.
  int nact_[2];
  // For each active LP constraint (at most 4), [active facet index, set].
  MdArray<4, 2> idx_lp_;
  // Solution status.
  SolutionStatus status_;

  void UpdateLocalX();

  void UpdateLocalDx();

  void UpdateLocalXDx();

  SolutionStatus Initialize(const Polyhedron* sets[2], int& iter);

  SolutionStatus FaceFaceDescent(const Polyhedron* sets[2], int& iter);

  SolutionStatus FaceEdgeDescent(const Polyhedron* sets[2], int& iter);

  SolutionStatus FinalDescent(const Polyhedron* sets[2], int& iter);

  bool FeasibilityCheck(const Polyhedron* sets[2],
                        std::array<bool, 2>& feas) const;

  void FeasibilityRecovery(const Polyhedron* sets[2],
                           std::array<bool, 2>& feas);

  EntryPoint Incremental(const Polyhedron* sets[2], int& iter);

  void ComputeDualOptimalSolution(const Polyhedron* sets[2]);
};

inline GrowthDistanceSolver::GrowthDistanceSolver() : ls_(kLsTol), settings_{} {
  A_.setConstant(-1.0);
}

inline Settings& GrowthDistanceSolver::SolverSettings() { return settings_; }

inline void GrowthDistanceSolver::UpdateLocalX() {
  Rtx_p_[0] = R_[0].transpose() * (res_.head<3>() - p_[0]);
  Rtx_p_[1] = R_[1].transpose() * (res_.head<3>() - p_[1]);
}

inline void GrowthDistanceSolver::UpdateLocalDx() {
  Rtdx_[0] = R_[0].transpose() * dx_;
  Rtdx_[1] = R_[1].transpose() * dx_;
}

inline void GrowthDistanceSolver::UpdateLocalXDx() {
  UpdateLocalX();
  UpdateLocalDx();
}

}  // namespace inc

#endif  // DGD_BENCHMARK_INC_INCREMENTAL_H_

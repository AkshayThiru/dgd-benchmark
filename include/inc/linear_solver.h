#ifndef INC_LINEAR_SOLVER_H_
#define INC_LINEAR_SOLVER_H_

#include <Eigen/Dense>
#include <cmath>

#include "inc/data_types.h"

namespace inc {

class LinearSolver {
 public:
  explicit LinearSolver(double tolerance);

  ~LinearSolver() = default;

  // Solves Gamma * dx = (-1, -1), where
  // Gamma = A_{0:3, 0:2}^T * A_{0:3, 0:2}.
  bool ComputeFaceFaceGradient(const Mat4& A, Vec3& dx);

  // Solves A_{0:3, 0:3}^T * dx = (-1, -1, -1).
  bool ComputeFaceEdgeGradient(const Mat4& A, Vec3& dx);

  bool ComputeFinalDual(const Mat4& A, Vec4& lambda);

  void ComputeFinalGradient(int i, double lambdai, Vec3& dx);

  void ComputeActiveConstraintSetInverse(Mat4& inv);

  bool ComputeWarmStartSolution(const Mat4& A, const Vec4& b, Vec4& res);

 private:
  struct {
    Eigen::PartialPivLU<Mat4> lu;
  } f_;  // Variables for final descent and warm start.

  struct {
    Eigen::PartialPivLU<Mat3> lu;
  } fe_;  // Variables for face-edge descent.

  struct {
    Mat2 mat;
    Vec2 alpha;
  } ff_;  // Variables for face-face descent.

  const double tol_;
};

inline LinearSolver::LinearSolver(double tolerance) : tol_(tolerance) {}

inline bool LinearSolver::ComputeFaceFaceGradient(const Mat4& A, Vec3& dx) {
  ff_.mat = A.block<3, 2>(0, 0).transpose() * A.block<3, 2>(0, 0);
  const double det{ff_.mat.determinant()};

  if (det < tol_ * ff_.mat(0, 0) * ff_.mat(1, 1)) return false;

  ff_.alpha =
      Vec2(ff_.mat(0, 1) - ff_.mat(1, 1), ff_.mat(0, 1) - ff_.mat(0, 0)) / det;
  dx = A.block<3, 2>(0, 0) * ff_.alpha;
  return true;
}

inline bool LinearSolver::ComputeFaceEdgeGradient(const Mat4& A, Vec3& dx) {
  fe_.lu.compute(A.block<3, 3>(0, 0).transpose());

  if (std::abs(fe_.lu.determinant()) < tol_) return false;

  dx = fe_.lu.solve(-Vec3::Ones());
  return true;
}

inline bool LinearSolver::ComputeFinalDual(const Mat4& A, Vec4& lambda) {
  f_.lu.compute(A);

  if (std::abs(f_.lu.determinant()) < tol_) return false;

  lambda = f_.lu.solve(-Vec4::UnitW());
  return true;
}

inline void LinearSolver::ComputeFinalGradient(int i, double lambdai,
                                               Vec3& dx) {
  Vec4 ei{Vec4::Zero()};
  ei(i) = 1.0;
  dx = f_.lu.transpose().solve(ei).eval().head<3>() / lambdai;
}

inline void LinearSolver::ComputeActiveConstraintSetInverse(Mat4& inv) {
  inv = f_.lu.inverse();
}

inline bool LinearSolver::ComputeWarmStartSolution(const Mat4& A, const Vec4& b,
                                                   Vec4& res) {
  f_.lu.compute(A.transpose());

  if (std::abs(f_.lu.determinant()) < tol_) return false;

  res = f_.lu.solve(b);
  return true;
}

}  // namespace inc

#endif  // INC_LINEAR_SOLVER_H_

#ifndef DGD_BENCHMARK_INC_DATA_TYPES_H_
#define DGD_BENCHMARK_INC_DATA_TYPES_H_

#include <Eigen/Core>
#include <cmath>
#include <limits>

namespace inc {

// Constants.
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kEps = std::numeric_limits<double>::epsilon();
inline const double kSqrtEps = std::sqrt(kEps);

// Types.
template <int dim>
using Vec = Eigen::Matrix<double, dim, 1>;
using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Vec4 = Eigen::Vector4d;

template <int row, int col>
using Mat = Eigen::Matrix<double, row, col>;
using Mat2 = Eigen::Matrix2d;
using Mat3 = Eigen::Matrix3d;
using Mat4 = Eigen::Matrix4d;

using Rotation3 = Mat3;
using Transform3 = Mat4;

inline auto Affine(Mat4& tf) -> decltype(tf.topRightCorner<3, 1>()) {
  return tf.topRightCorner<3, 1>();
}

inline auto Affine(const Mat4& tf) -> decltype(tf.topRightCorner<3, 1>()) {
  return tf.topRightCorner<3, 1>();
}

inline auto Linear(Mat4& tf) -> decltype(tf.topLeftCorner<3, 3>()) {
  return tf.topLeftCorner<3, 3>();
}

inline auto Linear(const Mat4& tf) -> decltype(tf.topLeftCorner<3, 3>()) {
  return tf.topLeftCorner<3, 3>();
}

}  // namespace inc

#endif  // DGD_BENCHMARK_INC_DATA_TYPES_H_

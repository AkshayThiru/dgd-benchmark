#ifndef DGD_BENCHMARK_DATA_TYPES_H_
#define DGD_BENCHMARK_DATA_TYPES_H_

#include <Eigen/Core>

constexpr int kDsfExp = 16;  // Exponent for VDSF.

using Vec3 = Eigen::Vector3d;

using Rotation3 = Eigen::Matrix3d;
using Transform3 = Eigen::Matrix4d;

#endif  // DGD_BENCHMARK_DATA_TYPES_H_

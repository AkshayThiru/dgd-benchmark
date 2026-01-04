#ifndef DGD_BENCHMARK_DCF_PRECOMPILED_H_
#define DGD_BENCHMARK_DCF_PRECOMPILED_H_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

namespace dcf {

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kEps = std::numeric_limits<double>::epsilon();
inline const double kSqrtEps = std::sqrt(kEps);

#define COUT_CHECK std::cout << "CHECK" << std::endl;
#define COUT_SCALAR(val)                                                 \
  std::cout << std::fixed << std::setprecision(9) << #val << ": " << val \
            << std::setprecision(5) << std::endl;
#define COUT_VECTOR(vec)                                          \
  std::cout << std::fixed << std::setprecision(5) << #vec << ": " \
            << (vec).transpose() << std::endl;
#define COUT_MATRIX(mat)                                                       \
  std::cout << std::fixed << std::setprecision(5) << #mat << ": " << std::endl \
            << mat << std::endl;

#define ERROR_EXIT(msg)                         \
  do {                                          \
    std::cerr << "Error: " << msg << std::endl; \
    std::exit(EXIT_FAILURE);                    \
  } while (0)

template <int dim>
using Vec = Eigen::Matrix<double, dim, 1>;
using Vec3 = Eigen::Vector3d;
using Vec4 = Eigen::Vector4d;

template <int row, int col>
using Mat = Eigen::Matrix<double, row, col>;
using Mat3 = Eigen::Matrix3d;
using Mat4 = Eigen::Matrix4d;

using Rotation3 = Mat3;
using Transform3 = Mat4;

using Quaternion = Eigen::Quaterniond;

}  // namespace dcf

#endif  // DGD_BENCHMARK_DCF_PRECOMPILED_H_

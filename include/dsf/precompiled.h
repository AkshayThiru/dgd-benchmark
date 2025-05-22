#ifndef DSF_PRECOMPILED_H_
#define DSF_PRECOMPILED_H_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

namespace dsf {

const double kEps{std::numeric_limits<double>::epsilon()};
const double kEpsSqrt{std::sqrt(kEps)};

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
using Vec = Eigen::Vector<double, dim>;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector4d Vec4;

template <int row, int col>
using Mat = Eigen::Matrix<double, row, col>;
typedef Eigen::Matrix3d Mat3;
typedef Eigen::Matrix4d Mat4;

typedef Mat3 Rotation3;
typedef Mat4 Transform3;

typedef Eigen::Quaterniond Quaternion;

}  // namespace dsf

#endif  // DSF_PRECOMPILED_H_

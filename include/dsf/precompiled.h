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

using Eigen::Vector;
using Eigen::Vector3d;
using Eigen::Vector4d;

using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::Matrix4d;

using Eigen::Quaterniond;

}  // namespace dsf

#endif  // DSF_PRECOMPILED_H_

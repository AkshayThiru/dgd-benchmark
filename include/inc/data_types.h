#ifndef INC_DATA_TYPES_H_
#define INC_DATA_TYPES_H_

#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <limits>

namespace inc {

// Constants.
const double kInf{std::numeric_limits<double>::infinity()};
const double kEps{std::numeric_limits<double>::epsilon()};
const double kEpsSqrt{std::sqrt(kEps)};

// Types.
template <int dim>
using Vec = Eigen::Vector<double, dim>;
typedef Eigen::Vector2d Vec2;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector4d Vec4;

template <int row, int col>
using Mat = Eigen::Matrix<double, row, col>;
typedef Eigen::Matrix2d Mat2;
typedef Eigen::Matrix3d Mat3;
typedef Eigen::Matrix4d Mat4;

typedef Mat3 Rotation3;
typedef Mat4 Transform3;

}  // namespace inc

#endif  // INC_DATA_TYPES_H_

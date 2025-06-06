#ifndef DGD_BENCHMARK_DSF_UTILS_H_
#define DGD_BENCHMARK_DSF_UTILS_H_

#include "dgd_benchmark/dsf/precompiled.h"

namespace dsf {

// Quaternion SE(3) to Matrix SE(3).
inline void Quaternion2RotationSe3(const Vec<7>& tfq, Transform3& tf) {
  tf.block<3, 3>(0, 0) = Quaternion(tfq.tail<4>()).toRotationMatrix();
  tf.block<3, 1>(0, 3) = tfq.head<3>();
  tf.row(3) = Vec4::UnitW().transpose();
}

// Hat map.
inline void Skew(const Vec3& w, Mat3& W) {
  W << 0.0, -w(2), w(1), w(2), 0.0, -w(0), -w(1), w(0), 0.0;
}

// Exponentiation function
template <unsigned int exp>
inline double Power(double base) {
  const double half_power = Power<exp / 2>(base);
  if constexpr (exp % 2 == 0) {
    return half_power * half_power;
  } else {
    return base * half_power * half_power;
  }
}

template <>
inline double Power<0>(double /*base*/) {
  return 1.0;
}

template <>
inline double Power<1>(double base) {
  return base;
}

}  // namespace dsf

#endif  // DGD_BENCHMARK_DSF_UTILS_H_

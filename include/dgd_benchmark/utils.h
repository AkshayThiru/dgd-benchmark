#ifndef DGD_BENCHMARK_UTILS_H_
#define DGD_BENCHMARK_UTILS_H_

#include "dgd/utils.h"
#include "dgd_benchmark/data_types.h"

inline void SetRandomTransforms(Transform3& tf1, Transform3& tf2,
                                double position_range) {
  dgd::RandomRigidBodyTransform<3>(-position_range, position_range, tf1);
  dgd::RandomRigidBodyTransform<3>(-position_range, position_range, tf2);
}

#endif  // DGD_BENCHMARK_UTILS_H_

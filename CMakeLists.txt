cmake_minimum_required(VERSION 3.15)

project(
  GrowthDistanceBenchmark
  VERSION 1.0.0
  DESCRIPTION "Benchmarks for growth distance algorithms"
)

message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID}")


# Options
option(DCF_SIGNED_DISTANCE "[DCF]: Return signed distance output" OFF)
option(IE_SMALL_SIMPLEX    "[IE] : Use small initial simplex" OFF)
option(IE_DISTANCE_TC      "[IE] : Use the original distance termination criterion" OFF)

option(BENCHMARK_BUILD_TESTS "Build tests" ON)


# Dependencies
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Eigen3 5.0 REQUIRED)
find_package(Arrow REQUIRED) # For data storage

set(DGD_LIB dgd)
add_subdirectory(extern/differentiable-growth-distance)
add_subdirectory(extern/differentiable-growth-distance/benchmarks/internal_helpers)

# Benchmark libraries
add_subdirectory(libraries/dcf)
add_subdirectory(libraries/ie)
add_subdirectory(libraries/inc)

# Tests
if (BENCHMARK_BUILD_TESTS)
  add_subdirectory(tests)
endif()

# Benchmarks
add_subdirectory(benchmarks)

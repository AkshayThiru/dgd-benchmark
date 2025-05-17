cmake_minimum_required(VERSION 3.15)

project(
  GrowthDistanceBenchmark
  VERSION 0.1.0
  DESCRIPTION "Benchmarks for growth distance algorithms"
)

message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID}")


# Options
option(DSF_SIGNED_DISTANCE "[DSF]: Return signed distance output" OFF)
option(IE_SMALL_SIMPLEX    "[IE] : Use small initial simplex" OFF)
option(IE_DISTANCE_TC      "[IE] : Use the original distance termination criterion" ON)


# Dependencies
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(Eigen3 3.4.90 REQUIRED NO_MODULE)
find_package(Qhull REQUIRED)
find_package(dgd 1.0.0 REQUIRED)

# Benchmark library
add_subdirectory(src)


# Tests
add_subdirectory(tests)

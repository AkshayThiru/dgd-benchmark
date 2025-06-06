# Source files
set(DSF_SOURCES
  "${PROJECT_SOURCE_DIR}/src/dsf/dsf_collision.cc"
  "${PROJECT_SOURCE_DIR}/src/dsf/dsf_interface.cc"
)
set(INC_SOURCES
  "${PROJECT_SOURCE_DIR}/src/inc/incremental.cc"
  "${PROJECT_SOURCE_DIR}/src/inc/polyhedron.cc"
)

set(BENCHMARK_SOURCES
  "${PROJECT_SOURCE_DIR}/src/benchmark_interface.cc"
  "${PROJECT_SOURCE_DIR}/src/convex_set_generator.cc"
  "${PROJECT_SOURCE_DIR}/src/file_io.cc"
)


# Benchmark library
add_library(benchmark STATIC
  ${DSF_SOURCES}
  ${INC_SOURCES}
  ${BENCHMARK_SOURCES}
)

target_include_directories(benchmark
  PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
)
target_link_libraries(benchmark
  PUBLIC
    Eigen3::Eigen
    dgd::dgd
    Arrow::arrow_shared
)

target_compile_features(benchmark PUBLIC cxx_std_17)
set_target_properties(benchmark PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(benchmark
  PUBLIC
    -Wall -Wextra -Werror -Wconversion -Wpedantic
    $<$<CONFIG:Debug>:-O0 -g>
    $<$<CONFIG:Release>:-O3 -DNDEBUG>
)

if (DSF_SIGNED_DISTANCE)
  target_compile_definitions(benchmark PRIVATE DSF_SIGNED_DISTANCE)
endif()
if (IE_SMALL_SIMPLEX)
  target_compile_definitions(benchmark PRIVATE IE_SMALL_SIMPLEX)
endif()
if (IE_DISTANCE_TC)
  target_compile_definitions(benchmark PRIVATE IE_DISTANCE_TC)
endif()

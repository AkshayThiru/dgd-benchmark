add_library(benchmark_helper STATIC
  "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_result.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_interface.cc"
)

target_include_directories(benchmark_helper
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/../"
)
target_link_libraries(benchmark_helper
  PUBLIC
    Eigen3::Eigen
    dgd::dgd
    dgd_bench
    dcf
    ie
    inc
  PRIVATE
    Arrow::arrow_shared
)

target_compile_features(benchmark_helper PUBLIC cxx_std_17)
target_compile_options(benchmark_helper
  PUBLIC
    -Wall -Wextra -Werror -Wconversion -Wpedantic
    $<$<CONFIG:Debug>:-O0 -g>
    $<$<CONFIG:Release>:-O3 -DNDEBUG>
    # $<$<CONFIG:Release>:-funroll-loops -march=native -flto>
)

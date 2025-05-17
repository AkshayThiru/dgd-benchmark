# Source files
set(DSF_SOURCES
  "${PROJECT_SOURCE_DIR}/src/dsf/dsf_collision.cc"
  "${PROJECT_SOURCE_DIR}/src/dsf/dsf_interface.cc"
)

# Benchmark library
add_library(bench STATIC
  ${DSF_SOURCES}
)

target_include_directories(bench
  PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
)
target_link_libraries(bench
  PUBLIC
    Eigen3::Eigen
    dgd::dgd
)

target_compile_features(bench PUBLIC cxx_std_17)
set_target_properties(bench PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(bench
  PUBLIC
    -Wall -Wextra -Werror -Wconversion -Wpedantic
    $<$<CONFIG:Debug>:-O0 -g>
    $<$<CONFIG:Release>:-O3 -DNDEBUG>
)

if (DSF_SIGNED_DISTANCE)
  target_compile_definitions(bench PRIVATE DSF_SIGNED_DISTANCE)
endif()
if (IE_SMALL_SIMPLEX)
  target_compile_definitions(bench PRIVATE IE_SMALL_SIMPLEX)
endif()
if (IE_DISTANCE_TC)
  target_compile_definitions(bench PRIVATE IE_DISTANCE_TC)
endif()

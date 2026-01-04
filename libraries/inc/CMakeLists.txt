add_library(inc STATIC
  "${CMAKE_CURRENT_SOURCE_DIR}/incremental.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyhedron.cc"
)

target_include_directories(inc
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/../"
)
target_link_libraries(inc
  PUBLIC
    Eigen3::Eigen
)

target_compile_features(inc PUBLIC cxx_std_17)
target_compile_options(inc
  PRIVATE
    -Wall -Wextra -Werror -Wconversion -Wpedantic
    $<$<CONFIG:Debug>:-O0 -g>
    $<$<CONFIG:Release>:-O3 -DNDEBUG>
    # $<$<CONFIG:Release>:-funroll-loops -march=native -flto>
)

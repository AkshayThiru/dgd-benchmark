add_library(dsf STATIC
  "${CMAKE_CURRENT_SOURCE_DIR}/dsf_collision.cc"
)

target_include_directories(dsf
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/../"
)
target_link_libraries(dsf
  PUBLIC
    Eigen3::Eigen
    dgd::dgd
)

target_compile_features(dsf PUBLIC cxx_std_17)
target_compile_options(dsf
  PRIVATE
    -Wall -Wextra -Werror -Wconversion -Wpedantic
    $<$<CONFIG:Debug>:-O0 -g>
    $<$<CONFIG:Release>:-O3 -DNDEBUG>
)

if (DSF_SIGNED_DISTANCE)
  target_compile_definitions(dsf PRIVATE DSF_SIGNED_DISTANCE)
endif()

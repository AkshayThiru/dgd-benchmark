add_library(dcf STATIC
  "${CMAKE_CURRENT_SOURCE_DIR}/dcf_collision.cc"
)

target_include_directories(dcf
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/../"
)
target_link_libraries(dcf
  PUBLIC
    Eigen3::Eigen
    dgd::dgd
)

target_compile_features(dcf PUBLIC cxx_std_17)
target_compile_options(dcf
  PRIVATE
    -Wall -Wextra -Werror -Wconversion -Wpedantic
    $<$<CONFIG:Debug>:-O0 -g>
    $<$<CONFIG:Release>:-O3 -DNDEBUG>
)

if (DCF_SIGNED_DISTANCE)
  target_compile_definitions(dcf PRIVATE DCF_SIGNED_DISTANCE)
endif()

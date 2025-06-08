add_library(ie INTERFACE)

target_include_directories(ie
  INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/../"
)
target_link_libraries(ie
  INTERFACE
    Eigen3::Eigen
    dgd::dgd
)

target_compile_features(ie INTERFACE cxx_std_17)

if (IE_SMALL_SIMPLEX)
  target_compile_definitions(ie INTERFACE IE_SMALL_SIMPLEX)
endif()
if (IE_DISTANCE_TC)
  target_compile_definitions(ie INTERFACE IE_DISTANCE_TC)
endif()

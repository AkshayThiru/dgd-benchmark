add_subdirectory(helpers)

add_executable(dsf_bm
  "${CMAKE_CURRENT_SOURCE_DIR}/dsf_bm.cc"
)
target_link_libraries(dsf_bm PRIVATE benchmark_helper)

add_executable(mesh_bm
  "${CMAKE_CURRENT_SOURCE_DIR}/mesh_bm.cc"
)
target_link_libraries(mesh_bm PRIVATE benchmark_helper)

add_executable(primitive_bm
  "${CMAKE_CURRENT_SOURCE_DIR}/primitive_bm.cc"
)
target_link_libraries(primitive_bm PRIVATE benchmark_helper)

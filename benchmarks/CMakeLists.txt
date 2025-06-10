add_subdirectory(helpers)


add_executable(dsf_bm
  "${CMAKE_CURRENT_SOURCE_DIR}/dsf_bm.cc"
)
target_link_libraries(dsf_bm PRIVATE benchmark_helper)

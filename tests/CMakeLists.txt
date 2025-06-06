add_executable(dsf_test "${PROJECT_SOURCE_DIR}/tests/dsf_test.cc")
target_link_libraries(dsf_test PRIVATE benchmark)

add_executable(ie_test "${PROJECT_SOURCE_DIR}/tests/ie_test.cc")
target_link_libraries(ie_test PRIVATE benchmark)

add_executable(inc_test "${PROJECT_SOURCE_DIR}/tests/inc_test.cc")
target_link_libraries(inc_test PRIVATE benchmark)

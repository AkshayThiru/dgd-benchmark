add_executable(dsf_test "${PROJECT_SOURCE_DIR}/tests/dsf_test.cc")
target_link_libraries(dsf_test PRIVATE bench)

add_executable(ie_test "${PROJECT_SOURCE_DIR}/tests/ie_test.cc")
target_link_libraries(ie_test PRIVATE bench)

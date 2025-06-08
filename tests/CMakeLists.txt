add_executable(dsf_test "${CMAKE_CURRENT_SOURCE_DIR}/dsf_test.cc")
target_link_libraries(dsf_test PRIVATE dsf)

add_executable(ie_test "${CMAKE_CURRENT_SOURCE_DIR}/ie_test.cc")
target_link_libraries(ie_test PRIVATE ie)

add_executable(inc_test "${CMAKE_CURRENT_SOURCE_DIR}/inc_test.cc")
target_link_libraries(inc_test PRIVATE inc dgd::dgd)

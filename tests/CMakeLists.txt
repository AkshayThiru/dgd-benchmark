add_executable(dcf_test "${CMAKE_CURRENT_SOURCE_DIR}/dcf_test.cc")
target_link_libraries(dcf_test PRIVATE dcf)

add_executable(ie_test "${CMAKE_CURRENT_SOURCE_DIR}/ie_test.cc")
target_link_libraries(ie_test PRIVATE ie)

add_executable(inc_test "${CMAKE_CURRENT_SOURCE_DIR}/inc_test.cc")
target_link_libraries(inc_test PRIVATE inc dgd::dgd)

add_executable(test_fmempool test_fmempool.F90)
target_link_libraries(test_fmempool fmempool)
add_test(NAME fmempool COMMAND test_fmempool)
set_tests_properties(fmempool PROPERTIES TIMEOUT 20)

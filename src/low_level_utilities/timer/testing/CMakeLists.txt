add_executable(test_timer test_timer.F90)
target_link_libraries(test_timer sll_timer)
add_test(NAME timer COMMAND test_timer)
set_tests_properties(timer PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

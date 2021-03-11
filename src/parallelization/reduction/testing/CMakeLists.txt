add_executable(test_reduction test_reduction.F90)
target_link_libraries(test_reduction sll_reduction)
add_test(NAME reduction COMMAND test_reduction)
set_tests_properties(reduction PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

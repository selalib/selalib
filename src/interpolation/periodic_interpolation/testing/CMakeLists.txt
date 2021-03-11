add_executable(test_periodic_interpolation test_periodic_interpolation.F90)
target_link_libraries(test_periodic_interpolation sll_periodic_interpolation)
add_test(NAME periodic_interpolation COMMAND test_periodic_interpolation)
set_tests_properties(periodic_interpolation PROPERTIES PASS_REGULAR_EXPRESSION
                                                       "PASSED")

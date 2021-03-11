add_executable(test_dg_fields test_dg_fields.F90)
target_link_libraries(test_dg_fields sll_fields)
add_executable(test_scalar_field_1d test_scalar_field_1d.F90)
target_link_libraries(test_scalar_field_1d sll_fields)
add_test(NAME scalar_field_1d COMMAND test_scalar_field_1d)
set_tests_properties(scalar_field_1d PROPERTIES PASS_REGULAR_EXPRESSION
                                                "PASSED")

add_executable(test_scalar_field_2d test_scalar_field_2d.F90
                                    helper_functions.F90)
target_link_libraries(test_scalar_field_2d sll_fields)
add_test(NAME scalar_field_2d COMMAND test_scalar_field_2d)
set_tests_properties(scalar_field_2d PROPERTIES PASS_REGULAR_EXPRESSION
                                                "PASSED")

add_executable(test_plot_array_4d test_plot_array_4d.F90)
target_link_libraries(test_plot_array_4d sll_fields)
add_test(NAME plot_array_4d COMMAND test_plot_array_4d)
set_tests_properties(plot_array_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

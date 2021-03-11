add_executable(test_control_variate test_control_variate.F90)
target_link_libraries(test_control_variate sll_pic_weights sll_particle_groups)

add_test(NAME control_variate COMMAND test_control_variate)
set_tests_properties(control_variate PROPERTIES PASS_REGULAR_EXPRESSION
                                                "PASSED")

if(HDF5_ENABLED)

  add_executable(test_linear_pendulum test_linear_pendulum.F90)
  target_link_libraries(test_linear_pendulum sll_splitting_methods)
  add_test(NAME splitting_methods_linear_pendulum COMMAND test_linear_pendulum)
  set_tests_properties(splitting_methods_linear_pendulum
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

  add_executable(test_const_coef_adv_2d test_const_coef_adv_2d.F90)
  target_link_libraries(test_const_coef_adv_2d sll_splitting_methods)
  # ADD_TEST( NAME splitting_methods_const_coef_advection_2d COMMAND
  # test_const_coef_adv_2d ) SET_TESTS_PROPERTIES(
  # splitting_methods_const_coef_advection_2d PROPERTIES PASS_REGULAR_EXPRESSION
  # "PASSED" )

  add_executable(test_adv_2d test_adv_2d.F90)
  target_link_libraries(test_adv_2d sll_splitting_methods)
  # ADD_TEST( NAME splitting_methods_advection_2d COMMAND test_adv_2d )
  # SET_TESTS_PROPERTIES( splitting_methods_advection_2d PROPERTIES
  # PASS_REGULAR_EXPRESSION "PASSED" )

endif(HDF5_ENABLED)

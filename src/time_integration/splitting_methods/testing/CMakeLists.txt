IF(HDF5_ENABLED)

  ADD_EXECUTABLE( test_linear_pendulum test_linear_pendulum.F90 )
  TARGET_LINK_LIBRARIES( test_linear_pendulum sll_splitting_methods )
  ADD_TEST( NAME splitting_methods_linear_pendulum COMMAND test_linear_pendulum )
  SET_TESTS_PROPERTIES( splitting_methods_linear_pendulum PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  
  ADD_EXECUTABLE( test_const_coef_adv_2d test_const_coef_adv_2d.F90 )
  TARGET_LINK_LIBRARIES( test_const_coef_adv_2d sll_splitting_methods )
  #ADD_TEST( NAME splitting_methods_const_coef_advection_2d COMMAND test_const_coef_adv_2d )
  #SET_TESTS_PROPERTIES( splitting_methods_const_coef_advection_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  ADD_EXECUTABLE( test_adv_2d test_adv_2d.F90 )
  TARGET_LINK_LIBRARIES( test_adv_2d sll_splitting_methods )
  #ADD_TEST( NAME splitting_methods_advection_2d COMMAND test_adv_2d )
  #SET_TESTS_PROPERTIES( splitting_methods_advection_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ENDIF(HDF5_ENABLED)

ADD_LIBRARY( sll_operator_splitting STATIC
  sll_time_splitting_base.F90 
  sll_operator_splitting_base.F90
  sll_linear_pendulum_operators.F90
  sll_const_coef_adv_2d.F90
  sll_adv_2d.F90
  sll_vlasov_poisson_adv_cart2d.F90 )

TARGET_LINK_LIBRARIES( sll_operator_splitting
  sll_distribution_function
  sll_interpolators
  sll_utilities	
  sll_file_io
  sll_advection
  sll_poisson_solvers
  ode_solvers )

ADD_LIBRARY( sll_time_splitting STATIC
  sll_time_splitting_coeff.F90 )

ADD_LIBRARY( sll_time_composition STATIC
  sll_time_composition_coeff.F90 )

  
TARGET_LINK_LIBRARIES( sll_time_composition 
                       sll_working_precision
                       sll_memory)

  
# Tests
IF(BUILD_TESTING)

  #ADD_EXECUTABLE( vp1d_bsl_operator_split vp1d_bsl_operator_split.F90 )
  #TARGET_LINK_LIBRARIES( vp1d_bsl_operator_split sll_operator_splitting )
  ADD_EXECUTABLE( test_linear_pendulum unit_test_linear_pendulum.F90 )
  TARGET_LINK_LIBRARIES( test_linear_pendulum sll_operator_splitting )
  ADD_TEST( NAME splitting_linear_pendulum COMMAND test_linear_pendulum )
  # Ctest
  IF(HDF5_ENABLED)
    ADD_EXECUTABLE( test_const_coef_adv_2d unit_test_const_coef_adv_2d.F90 )
    TARGET_LINK_LIBRARIES( test_const_coef_adv_2d sll_operator_splitting )
    ADD_TEST( NAME splitting_const_coef_adv_2d COMMAND test_const_coef_adv_2d )

    ADD_EXECUTABLE( test_adv_2d unit_test_adv_2d.F90 )
    TARGET_LINK_LIBRARIES( test_adv_2d sll_operator_splitting )
    ADD_TEST( NAME splitting_adv_2d COMMAND test_adv_2d )
  ENDIF(HDF5_ENABLED)

SET_TESTS_PROPERTIES( splitting_linear_pendulum
                      #splitting_const_coef_adv_2d
                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(BUILD_TESTING)

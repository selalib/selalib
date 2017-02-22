  IF(HDF5_ENABLED)
    IF(HDF5_IS_PARALLEL)
      SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
    ENDIF(HDF5_IS_PARALLEL)
    #ADD_EXECUTABLE( vp1d_bsl_operator_split vp1d_bsl_operator_split.F90 )
    #TARGET_LINK_LIBRARIES( vp1d_bsl_operator_split sll_splitting_methods )
    ADD_EXECUTABLE( test_linear_pendulum test_linear_pendulum.F90 )
    TARGET_LINK_LIBRARIES( test_linear_pendulum sll_splitting_methods )
    ADD_TEST( NAME linear_pendulum COMMAND test_linear_pendulum )
    # Ctest
    ADD_EXECUTABLE( test_const_coef_adv_2d test_const_coef_adv_2d.F90 )
    TARGET_LINK_LIBRARIES( test_const_coef_adv_2d sll_splitting_methods )
    ADD_TEST( NAME const_coef_adv_2d COMMAND test_const_coef_adv_2d )

    ADD_EXECUTABLE( test_adv_2d test_adv_2d.F90 )
    TARGET_LINK_LIBRARIES( test_adv_2d sll_splitting_methods )
    ADD_TEST( NAME adv_2d COMMAND test_adv_2d )

    SET_TESTS_PROPERTIES( linear_pendulum
                      #splitting_const_coef_adv_2d
                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  ENDIF(HDF5_ENABLED)

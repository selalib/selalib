SET( GYROAVERAGE_SRC_FILES 
  sll_gyroaverage_2d_base.F90
  sll_gyroaverage_2d_polar.F90
  #sll_gyroaverage_2d_polar_computation.F90
  sll_gyroaverage_2d_polar_pade_solver.F90
  sll_gyroaverage_2d_polar_splines_solver.F90
  sll_gyroaverage_2d_polar_hermite_solver.F90 )

ADD_LIBRARY( sll_gyroaverage STATIC ${GYROAVERAGE_SRC_FILES} )

TARGET_LINK_LIBRARIES( sll_gyroaverage
  sll_characteristics
  sll_interpolators
  sll_periodic_interpolation
  sll_splines )

# Ctests
IF(BUILD_TESTING)
  
  ADD_EXECUTABLE( test_gyroaverage_2d_polar_hermite
    unit_test_gyroaverage_2d_polar_hermite.F90 )

  ADD_EXECUTABLE( test_gyroaverage_2d_polar_splines
    unit_test_gyroaverage_2d_polar_splines.F90 )

  ADD_EXECUTABLE( test_gyroaverage_2d_polar_pade
    unit_test_gyroaverage_2d_polar_pade.F90 )

  TARGET_LINK_LIBRARIES( test_gyroaverage_2d_polar_hermite sll_gyroaverage )
  TARGET_LINK_LIBRARIES( test_gyroaverage_2d_polar_splines sll_gyroaverage )
  TARGET_LINK_LIBRARIES( test_gyroaverage_2d_polar_pade    sll_gyroaverage )

  ADD_TEST( NAME gyroaverage_polar_hermite COMMAND test_gyroaverage_2d_polar_hermite )
  ADD_TEST( NAME gyroaverage_polar_splines COMMAND test_gyroaverage_2d_polar_splines )
  ADD_TEST( NAME gyroaverage_polar_pade    COMMAND test_gyroaverage_2d_polar_pade    )

  SET_TESTS_PROPERTIES( gyroaverage_polar_hermite PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  SET_TESTS_PROPERTIES( gyroaverage_polar_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  SET_TESTS_PROPERTIES( gyroaverage_polar_pade    PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ENDIF(BUILD_TESTING)

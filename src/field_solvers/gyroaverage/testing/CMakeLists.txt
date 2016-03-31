  SET( INPUT_FILES zeros_bessel.txt )
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

  ADD_EXECUTABLE( test_gyroaverage_2d_polar_hermite
    test_gyroaverage_2d_polar_hermite.F90 )

  ADD_EXECUTABLE( test_gyroaverage_2d_polar_splines
    test_gyroaverage_2d_polar_splines.F90 )

  ADD_EXECUTABLE( test_gyroaverage_2d_polar_pade
    test_gyroaverage_2d_polar_pade.F90 )

  ADD_EXECUTABLE(  test_qn_2d_polar_splines
    test_qn_2d_polar_splines.F90)  
 
  ADD_EXECUTABLE(  test_qn_2d_polar_precompute
    test_qn_2d_polar_precompute.F90)

  TARGET_LINK_LIBRARIES( test_gyroaverage_2d_polar_hermite sll_gyroaverage )
  TARGET_LINK_LIBRARIES( test_gyroaverage_2d_polar_splines sll_gyroaverage )
  TARGET_LINK_LIBRARIES( test_gyroaverage_2d_polar_pade    sll_gyroaverage )
  TARGET_LINK_LIBRARIES( test_qn_2d_polar_splines          sll_gyroaverage )  
  TARGET_LINK_LIBRARIES( test_qn_2d_polar_precompute       sll_gyroaverage )

  ADD_TEST( NAME gyroaverage_2d_polar_hermite COMMAND test_gyroaverage_2d_polar_hermite )
  ADD_TEST( NAME gyroaverage_2d_polar_splines COMMAND test_gyroaverage_2d_polar_splines )
  ADD_TEST( NAME gyroaverage_2d_polar_pade    COMMAND test_gyroaverage_2d_polar_pade    )
  ADD_TEST( NAME qn_2d_polar_splines          COMMAND test_qn_2d_polar_splines          )
  ADD_TEST( NAME qn_2d_polar_precompute       COMMAND test_qn_2d_polar_precompute       )

  SET_TESTS_PROPERTIES( gyroaverage_2d_polar_hermite PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  SET_TESTS_PROPERTIES( gyroaverage_2d_polar_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  SET_TESTS_PROPERTIES( gyroaverage_2d_polar_pade    PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  SET_TESTS_PROPERTIES( qn_2d_polar_splines          PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  SET_TESTS_PROPERTIES( qn_2d_polar_precompute       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

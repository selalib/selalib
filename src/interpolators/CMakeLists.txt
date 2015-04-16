SET( INTERPOLATORS_SRC_FILES 
  sll_interpolators_1d_base.F90
  sll_interpolators_2d_base.F90 
  sll_cubic_spline_interpolator_1d.F90
  sll_quintic_spline_interpolator_1d.F90
  sll_cubic_spline_interpolator_1d_nonuniform.F90
  sll_lagrange_interpolator_1d.F90 
  sll_periodic_interpolator_1d.F90
  sll_cubic_spline_interpolator_2d.F90
  sll_hermite_interpolator_2d.F90
  sll_arbitrary_degree_spline_interpolator_1d.F90
  sll_arbitrary_degree_spline_interpolator_2d.F90 ) 

ADD_LIBRARY( sll_interpolators STATIC ${INTERPOLATORS_SRC_FILES} )

TARGET_LINK_LIBRARIES( sll_interpolators sll_splines
                                         sll_periodic_interpolation
                                         sll_lagrange_interpolation
                                         sll_hermite_interpolation
                                         sll_fft
                                         sll_utilities 
                                         sll_timer
                                         deboor_splines
                                         ${LAPACK_LIBRARIES} )

IF(BUILD_TESTING)

  SLL_TEST( cubic_spline_interpolator_1d sll_interpolators )
  SLL_TEST( cubic_spline_interpolator_2d sll_interpolators )
  SLL_TEST( arbitrary_degree_spline_interpolator_1d sll_interpolators sll_file_io )
  SLL_TEST( arbitrary_degree_spline_interpolator_2d sll_interpolators sll_file_io )
  SLL_TEST( periodic_interpolator_1d sll_interpolators )
  SLL_TEST( arbitrary_degree_splines_1d_periodic  sll_interpolators )
  SLL_TEST( arbitrary_degree_splines_1d_dirichlet sll_interpolators )
  SLL_TEST( arbitrary_degree_splines_1d_neumann sll_interpolators )
  SLL_TEST( quintic_spline_interpolator_1d sll_interpolators )

  SET_TESTS_PROPERTIES( cubic_spline_interpolator_1d 
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
  SET_TESTS_PROPERTIES( cubic_spline_interpolator_2d 
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
  SET_TESTS_PROPERTIES( arbitrary_degree_spline_interpolator_1d 
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
  SET_TESTS_PROPERTIES( arbitrary_degree_spline_interpolator_2d 
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20 )
  SET_TESTS_PROPERTIES( periodic_interpolator_1d 
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
  SET_TESTS_PROPERTIES( arbitrary_degree_splines_1d_periodic
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
  SET_TESTS_PROPERTIES( arbitrary_degree_splines_1d_dirichlet
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
  SET_TESTS_PROPERTIES( arbitrary_degree_splines_1d_neumann 
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
  SET_TESTS_PROPERTIES( quintic_spline_interpolator_1d 
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )

ENDIF(BUILD_TESTING)

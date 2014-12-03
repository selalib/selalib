SET( INTERPOLATORS_SRC_FILES 
  sll_interpolators_1d_base.F90
  sll_interpolators_2d_base.F90 
  sll_cubic_spline_interpolator_1d.F90
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
                                         deboor_splines 
                                         sll_timer
                                         ${LAPACK_LIBRARIES} )

# Ctests
IF(BUILD_TESTING)
  SLL_TEST( cubic_spline_interpolator_1d sll_interpolators )
  SLL_TEST( cubic_spline_interpolator_2d sll_interpolators )
  SLL_TEST( arbitrary_degree_spline_interpolator_1d sll_interpolators sll_file_io )
  SLL_TEST( arbitrary_degree_spline_interpolator_2d sll_interpolators sll_file_io )
  SLL_TEST( periodic_interpolator_1d sll_interpolators )
ENDIF(BUILD_TESTING)

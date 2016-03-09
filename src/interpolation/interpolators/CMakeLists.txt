SET( INTERPOLATORS_SRC_FILES 
  sll_m_interpolators_1d_base.F90
  sll_m_interpolators_2d_base.F90 
  sll_m_cubic_spline_interpolator_1d.F90
  sll_m_quintic_spline_interpolator_1d.F90
  sll_m_cubic_spline_interpolator_1d_nonuniform.F90
  sll_m_lagrange_interpolator_1d.F90 
  sll_m_hermite_interpolator_1d.F90 
  sll_m_periodic_interpolator_1d.F90
  sll_m_cubic_spline_interpolator_2d.F90
  sll_m_hermite_interpolator_2d.F90
  sll_m_arbitrary_degree_spline_interpolator_1d.F90
  sll_m_arbitrary_degree_spline_interpolator_2d.F90)
  #sll_m_bspline_interpolator_1d.F90
  #sll_m_bspline_interpolator_2d.F90)

ADD_LIBRARY( sll_interpolators STATIC ${INTERPOLATORS_SRC_FILES} )

TARGET_LINK_LIBRARIES( sll_interpolators sll_splines
                                         sll_periodic_interpolation
                                         sll_lagrange_interpolation
                                         sll_hermite_interpolation
                                         sll_fft
                                         sll_utilities 
                                         sll_timer
                                         ${LAPACK_LIBRARIES} )

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)

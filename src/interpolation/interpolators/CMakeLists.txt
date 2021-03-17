set(INTERPOLATORS_SRC_FILES
    sll_m_interpolators_1d_base.F90
    sll_m_interpolators_2d_base.F90
    sll_m_lagrange_interpolator_1d.F90
    sll_m_periodic_interpolator_1d.F90
    sll_m_cubic_spline_interpolator_1d.F90
    sll_m_cubic_spline_interpolator_1d_nonuniform.F90
    sll_m_cubic_spline_interpolator_2d.F90
    sll_m_hermite_interpolator_1d.F90
    sll_m_hermite_interpolator_2d.F90
    sll_m_nufft_interpolator_2d.F90
    sll_m_cubic_spline_halo_1d.F90)

add_library(sll_interpolators STATIC ${INTERPOLATORS_SRC_FILES})

target_link_libraries(
  sll_interpolators
  sll_splines
  sll_periodic_interpolation
  sll_lagrange_interpolation
  sll_hermite_interpolation
  sll_nufft_interpolation
  sll_utilities
  sll_timer
  ${LAPACK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)

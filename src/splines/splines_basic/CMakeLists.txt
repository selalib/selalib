set(SPLINES_FILES
    sll_m_bsplines.F90
    sll_m_bsplines_base.F90
    sll_m_bsplines_uniform.F90
    sll_m_bsplines_non_uniform.F90
    sll_m_spline_1d.F90
    sll_m_spline_2d.F90
    sll_m_spline_interpolator_1d.F90
    sll_m_spline_interpolator_2d.F90
    sll_m_polar_spline_interpolator_2d.F90
    sll_m_spline_matrix.F90
    sll_m_spline_matrix_base.F90
    sll_m_spline_matrix_dense.F90
    sll_m_spline_matrix_banded.F90
    sll_m_spline_matrix_periodic_banded.F90
    sll_m_low_level_bsplines.F90
    sll_m_splines_pp.F90
    sll_m_cubic_non_uniform_splines.F90
    sll_m_hex_pre_filters.F90
    sll_m_box_splines.F90)

list(APPEND SPLINES_FILES sll_m_cubic_splines.F90)

add_library(sll_splines STATIC ${SPLINES_FILES})

target_link_libraries(
  sll_splines
  sll_working_precision
  sll_assert
  sll_errors
  sll_constants
  sll_utilities
  sll_linear_solvers
  sll_meshes
  sll_timer
  sll_quadrature
  pppack
  ${LAPACK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)

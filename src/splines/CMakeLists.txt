SET( SPLINES_FILES
  #
  sll_m_bsplines.F90
  sll_m_bsplines_base.F90
  sll_m_bsplines_uniform.F90
  sll_m_bsplines_non_uniform.F90
  sll_m_spline_1d.F90
  sll_m_spline_2d.F90
  sll_m_spline_interpolator_1d.F90
  sll_m_spline_interpolator_2d.F90
  #
  sll_m_spline_matrix.F90
  sll_m_spline_matrix_base.F90
  sll_m_spline_matrix_dense.F90
  sll_m_spline_matrix_banded.F90
  sll_m_spline_matrix_periodic_banded.F90
  #
  sll_m_low_level_bsplines.F90
  #
  sll_m_deboor_splines_1d.F90 
  sll_m_splines_pp.F90
  sll_m_cubic_non_uniform_splines.F90
  sll_m_quintic_splines.F90
  sll_m_hex_pre_filters.F90
  sll_m_box_splines.F90 )

IF(PGI)
  LIST(APPEND SPLINES_FILES sll_m_cubic_splines_pgi.F90)
ELSE(PGI)
  LIST(APPEND SPLINES_FILES sll_m_cubic_splines.F90)
ENDIF(PGI)

ADD_LIBRARY( sll_splines STATIC ${SPLINES_FILES} )

TARGET_LINK_LIBRARIES( sll_splines
  sll_working_precision
  sll_assert
  sll_errors
  sll_constants
  sll_utilities
  sll_linear_solvers
  sll_meshes
  sll_timer
  pppack
  ${LAPACK_LIBRARIES} )

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)

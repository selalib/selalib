ADD_LIBRARY( sll_splines_smooth STATIC
  sll_m_polar_bsplines_2d.F90 )

TARGET_LINK_LIBRARIES( sll_splines_smooth
  sll_splines
  sll_coordinate_transformations
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
  ${LAPACK_LIBRARIES} )

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)

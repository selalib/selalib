SET( SPLINES_FILES
  sll_m_polar_bsplines_2d.F90
  sll_m_singular_mapping_base.F90
  sll_m_singular_mapping_discrete.F90
  sll_m_singular_mapping_analytic.F90
  sll_m_singular_mapping_analytic_target.F90
  sll_m_singular_mapping_analytic_czarny.F90
  sll_m_singular_mapping_advector_base.F90
  sll_m_singular_mapping_advector_constant.F90
  sll_m_singular_mapping_advector_rotating.F90
  sll_m_jacobian_2d_pseudo_cartesian.F90 )

ADD_LIBRARY( sll_splines_smooth STATIC ${SPLINES_FILES} )

TARGET_LINK_LIBRARIES( sll_splines_smooth
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
  ${LAPACK_LIBRARIES} )

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)

ADD_LIBRARY( sll_splines STATIC
  sll_m_cubic_splines.F90
  sll_m_quintic_splines.F90
  sll_m_bsplines.F90
  sll_m_arbitrary_degree_splines.F90
  sll_m_hex_pre_filters.F90
  sll_m_box_splines.F90
  sll_m_cubic_non_uniform_splines.F90 )

TARGET_LINK_LIBRARIES( sll_splines
  sll_linear_solvers
  sll_boundary_condition_descriptors
  sll_meshes
  sll_integration
  sll_linear_solvers
  ${LAPACK_LIBRARIES} )

# Tests

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)

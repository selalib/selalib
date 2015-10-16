ADD_LIBRARY( sll_splines STATIC
  sll_cubic_splines.F90
  sll_quintic_splines.F90
  sll_bsplines.F90
  sll_arbitrary_degree_splines.F90
  sll_hex_prefilters.F90
  sll_box_splines.F90
  cubic_nonuniform_splines.F90 )

TARGET_LINK_LIBRARIES( sll_splines
  sll_linear_solvers
  sll_boundary_condition_descriptors
  sll_meshes
  sll_integration
  ${LAPACK_LIBRARIES} )

# Tests

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)

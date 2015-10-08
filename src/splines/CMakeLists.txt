ADD_LIBRARY( sll_splines STATIC
  sll_cubic_splines.F90
  sll_quintic_splines.F90
  sll_bsplines.F90
  sll_arbitrary_degree_splines.F90
  cubic_nonuniform_splines.F90 )

TARGET_LINK_LIBRARIES( sll_splines
  sll_utilities
  sll_boundary_condition_descriptors
  ${LAPACK_LIBRARIES} )

ADD_LIBRARY( hex_pre_filters STATIC hex_pre_filters.F90 )
TARGET_LINK_LIBRARIES( hex_pre_filters sll_meshes )

ADD_LIBRARY( sll_box_splines STATIC sll_box_splines.F90 )
TARGET_LINK_LIBRARIES( sll_box_splines
  sll_splines
  sll_utilities
  hex_pre_filters
  sll_meshes
  sll_integration
  sll_boundary_condition_descriptors
  ${LAPACK_LIBRARIES})


# Tests
IF(BUILD_TESTING)

  ADD_EXECUTABLE( test_cubic_splines unit_test.F90 test_func_module.F90
                  util_constants.F90 test_processes_module.F90 )
  TARGET_LINK_LIBRARIES( test_cubic_splines sll_splines )

  ADD_EXECUTABLE( test_quintic_splines unit_test_quintic_splines.F90 )
  TARGET_LINK_LIBRARIES( test_quintic_splines sll_splines )

  ADD_EXECUTABLE( test_bsplines_1d unit_test_bsplines_1d.F90 )
  TARGET_LINK_LIBRARIES( test_bsplines_1d sll_splines )

  ADD_EXECUTABLE( test_bsplines_2d unit_test_bsplines_2d.F90 )
  TARGET_LINK_LIBRARIES( test_bsplines_2d sll_splines )

  ADD_EXECUTABLE( test_non_unif_splines unit_test_nonuniform.F90 )
  TARGET_LINK_LIBRARIES( test_non_unif_splines sll_splines )

  ADD_EXECUTABLE( test_arbitrary_degree_splines unit_test_arbitrary_degree.F90 )
  TARGET_LINK_LIBRARIES( test_arbitrary_degree_splines sll_splines )

  ADD_EXECUTABLE( test_box_splines_derivatives test_box_splines_derivatives.F90 )
  TARGET_LINK_LIBRARIES( test_box_splines_derivatives sll_box_splines )

  # Ctests
  ADD_TEST( NAME cubic_splines             COMMAND test_cubic_splines            )
  ADD_TEST( NAME quintic_splines           COMMAND test_quintic_splines          )
  ADD_TEST( NAME bsplines_1d               COMMAND test_bsplines_1d              )
  ADD_TEST( NAME bsplines_2d               COMMAND test_bsplines_2d              )
  ADD_TEST( NAME cubic_non_uniform_splines COMMAND test_non_unif_splines         )
  ADD_TEST( NAME splines_arbitrary_degree  COMMAND test_arbitrary_degree_splines )

  SET_TESTS_PROPERTIES( cubic_splines
                        bsplines_1d
                        bsplines_2d
                        cubic_non_uniform_splines
                        splines_arbitrary_degree
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )


ENDIF(BUILD_TESTING)


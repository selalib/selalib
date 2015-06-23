ADD_LIBRARY( sll_general_coordinate_elliptic_solvers_multipatch STATIC
  connectivity_qn.f90
  sll_general_coordinate_elliptic_solver_multipatch.F90
)

TARGET_LINK_LIBRARIES( sll_general_coordinate_elliptic_solvers_multipatch 
  sll_scalar_field_multipatch
  sll_coordinate_transformation_multipatch
  sll_sparse_matrix
  sll_integration
  sll_general_coordinate_elliptic_solvers )

IF(BUILD_TESTING AND MPI_ENABLED)
  ADD_EXECUTABLE( test_general_coordinate_elliptic_solver_multipatch 
                  test_functions.F90 unit_test.F90)
  TARGET_LINK_LIBRARIES( test_general_coordinate_elliptic_solver_multipatch
    sll_general_coordinate_elliptic_solvers_multipatch)

  IF(SUITESPARSE_ENABLED)
  ELSE()
  IF(PYTHON3_FOUND)
    ADD_TEST( NAME general_coordinate_elliptic_solver_multipatch 
    COMMAND test_general_coordinate_elliptic_solver_multipatch )
    SET_TESTS_PROPERTIES( general_coordinate_elliptic_solver_multipatch
    PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  ENDIF(PYTHON3_FOUND)
  ENDIF(SUITESPARSE_ENABLED)

  SET( INPUT_FILES_DIR ${CMAKE_SOURCE_DIR}/CAID_coordinate_transformations/circle_mp5_pts12/ )
  FILE( GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml" )
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

ENDIF(BUILD_TESTING AND MPI_ENABLED)

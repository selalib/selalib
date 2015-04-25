ADD_LIBRARY( sll_general_coordinate_elliptic_solvers_multipatch STATIC
  connectivity_qn.f90
  sll_general_coordinate_elliptic_solver_multipatch.F90
)
# profilmat.f90
# GC.f90
# sparse_matrix.F90
# sparse_matrix_con.F90
# sparse_matrix_operations.F90 
# connectivity.F90
# sparse_matrix_composed.F90
# tracelog.F90
# connectivities_def.F90
# sparse_matrix_basic.F90
# sparse_matrix_LM.F90  
# quicksort.f90
# sparse_matrix_def.F90
# precision.f90
# coocsr.f90 
# linpack.f90)

TARGET_LINK_LIBRARIES( sll_general_coordinate_elliptic_solvers_multipatch 
  sll_scalar_field_multipatch
  sll_coordinate_transformation_multipatch
  sll_sparse_matrix
  sll_integration
  original_deboor_splines)

  #sll_general_coordinate_elliptic_solvers
  #ll_sparse_matrix 
  #sll_fields
  #sll_integration)

IF(BUILD_TESTING AND MPI_ENABLED)
  ADD_EXECUTABLE( test_general_coordinate_elliptic_solver_multipatch unit_test.F90)
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

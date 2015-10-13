IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_ENABLED)


   ADD_LIBRARY(sll_distribution_function_multipatch
     STATIC sll_distribution_function_4d_multipatch.F90)

   TARGET_LINK_LIBRARIES(sll_distribution_function_multipatch
     sll_scalar_field_multipatch
     sll_coordinate_transformation_multipatch
     sll_remap
     sll_parallel_array_utilities
     sll_timer
     sll_meshes_multipatch
     sll_memory)

# Ctest
IF(BUILD_TESTING)

  ADD_EXECUTABLE( test_distribution_function_4d_multipatch unit_test.F90 )
  TARGET_LINK_LIBRARIES( test_distribution_function_4d_multipatch
    sll_distribution_function_multipatch
    sll_scalar_field_multipatch
    sll_common_array_initializers
    sll_timer
    sll_file_io
    sll_memory )

  IF(PYTHON3_FOUND)
    SET( ARGS " " )
     ADD_MPI_TEST( distribution_function_4d_multipatch 
       test_distribution_function_4d_multipatch ${PROCS} ${ARGS} )
     SET_TESTS_PROPERTIES( distribution_function_4d_multipatch 
       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  ENDIF(PYTHON3_FOUND)
  SET( INPUT_FILES_DIR ${CMAKE_SOURCE_DIR}/CAID_coordinate_transformations/square_4p_n10/ )
  FILE( GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml" )
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

ENDIF(BUILD_TESTING)

ADD_LIBRARY(sll_meshes_multipatch STATIC 
  sll_meshes_multipatch.F90)

TARGET_LINK_LIBRARIES( sll_meshes_multipatch
 sll_meshes
 sll_memory
 sll_assert
 sll_working_precision )


# Ctest
IF(BUILD_TESTING)
  ADD_EXECUTABLE( test_meshes_multipatch unit_test.F90)
  TARGET_LINK_LIBRARIES( test_meshes_multipatch sll_meshes_multipatch sll_assert )
  ADD_TEST( NAME meshes_multipatch COMMAND test_meshes_multipatch )
  SET_TESTS_PROPERTIES( meshes_multipatch PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(BUILD_TESTING)

ADD_LIBRARY( sll_coordinate_transformation_multipatch STATIC
  sll_coordinate_transformation_multipatch.F90 )

TARGET_LINK_LIBRARIES( sll_coordinate_transformation_multipatch
  sll_coordinate_transformations
  sll_meshes_multipatch
  sll_memory )

# Ctest
IF(BUILD_TESTING)

  ADD_EXECUTABLE( test_coordinate_transformation_multipatch_2d unit_test_2d.F90 )
  TARGET_LINK_LIBRARIES(test_coordinate_transformation_multipatch_2d
     sll_coordinate_transformation_multipatch
     sll_file_io
     sll_memory)

  IF(PYTHON3_FOUND)
    ADD_TEST( NAME coordinate_transformation_multipatch_2d
    COMMAND test_coordinate_transformation_multipatch_2d )
    SET_TESTS_PROPERTIES( coordinate_transformation_multipatch_2d 
    PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  ENDIF(PYTHON3_FOUND)

  SET( INPUT_FILES_DIR ${CMAKE_SOURCE_DIR}/CAID_coordinate_transformations/circle_mp5_pts12/ )
  FILE( GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml" )
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

ENDIF(BUILD_TESTING)

ADD_LIBRARY( sll_scalar_field_multipatch STATIC
             sll_scalar_field_2d_multipatch.F90 )

TARGET_LINK_LIBRARIES( sll_scalar_field_multipatch
                       sll_coordinate_transformation_multipatch
                       sll_fields
                       sll_meshes_multipatch )

# Ctest
IF(BUILD_TESTING)

  ADD_EXECUTABLE( test_scalar_field_multipatch_2d unit_test_2d.F90 )
  TARGET_LINK_LIBRARIES( test_scalar_field_multipatch_2d sll_scalar_field_multipatch )

  IF(PYTHON3_FOUND)
    ADD_TEST( NAME scalar_field_multipatch_2d COMMAND test_scalar_field_multipatch_2d )
    SET_TESTS_PROPERTIES( scalar_field_multipatch_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  ENDIF(PYTHON3_FOUND)
  
  SET( INPUT_FILES_DIR ${CMAKE_SOURCE_DIR}/CAID_coordinate_transformations/square_4p_n10/ )
  FILE( GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml" )
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

ENDIF(BUILD_TESTING)


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


ENDIF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_ENABLED)

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

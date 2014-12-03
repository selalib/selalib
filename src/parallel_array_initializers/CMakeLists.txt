ADD_LIBRARY(sll_common_array_initializers STATIC 
            sll_common_array_initializers_module.F90 )

TARGET_LINK_LIBRARIES(sll_common_array_initializers
                      sll_assert
                      sll_working_precision
                      sll_constants)

IF(MPI_MODULE_ENABLED)

   ADD_LIBRARY(sll_parallel_array_initializers STATIC 
               sll_parallel_array_initializer_module.F90
               sll_common_array_initializers_module.F90 )

   TARGET_LINK_LIBRARIES(sll_parallel_array_initializers 
                         sll_remap
                         sll_coordinate_transformations
                         sll_meshes)

  # Ctest
  IF(BUILD_TESTING AND (PROCESSOR_COUNT GREATER 1))

    ADD_EXECUTABLE( test_parallel_array_initializer unit_test_4d.F90 )
    TARGET_LINK_LIBRARIES( test_parallel_array_initializer 
      sll_parallel_array_initializers )

    SET( ARGS " " )
    ADD_MPI_TEST( parallel_array_initializers 
      test_parallel_array_initializer ${PROCS} ${ARGS} )
    SET_TESTS_PROPERTIES( parallel_array_initializers
      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  ENDIF(BUILD_TESTING AND (PROCESSOR_COUNT GREATER 1))

ENDIF(MPI_MODULE_ENABLED)

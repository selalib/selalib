ADD_LIBRARY(sll_common_array_initializers STATIC 
            sll_common_array_initializers_module.F90 )

TARGET_LINK_LIBRARIES(sll_common_array_initializers
                      sll_assert
                      sll_working_precision
                      sll_constants)

IF(MPI_ENABLED)

  ADD_LIBRARY( sll_parallel_utilities STATIC buffer_loader_utilities.F90 )
  
  TARGET_LINK_LIBRARIES( sll_parallel_utilities sll_remap sll_collective ) 
   
  IF(BUILD_TESTING)
    ADD_EXECUTABLE( test_buffer_loader_utilities
      unit_test_buffer_loader_utilities.F90 )
    TARGET_LINK_LIBRARIES( test_buffer_loader_utilities
      sll_parallel_utilities )
  ENDIF(BUILD_TESTING)

  ADD_LIBRARY(sll_parallel_array_utilities STATIC 
              sll_parallel_array_initializer_module.F90
              sll_common_array_initializers_module.F90 )

  TARGET_LINK_LIBRARIES(sll_parallel_array_utilities 
                        sll_remap
                        sll_coordinate_transformations
                        sll_file_io_parallel
                        sll_file_io
                        sll_utilities
                        sll_meshes)

  IF(BUILD_TESTING AND HDF5_ENABLED AND HDF5_PARALLEL_ENABLED)
    ADD_EXECUTABLE( test_parallel_array_initializer unit_test_4d.F90
              sll_parallel_array_output_module.F90)

    TARGET_LINK_LIBRARIES( test_parallel_array_initializer 
                           sll_parallel_array_utilities )

    SET( ARGS " " )
    ADD_MPI_TEST( parallel_array_initializers 
      test_parallel_array_initializer ${PROCS} ${ARGS} )
    SET_TESTS_PROPERTIES( parallel_array_initializers
      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  ENDIF()

  #ADD_EXECUTABLE(local_to_global local_to_global.F90)
  #ADD_EXECUTABLE(scatter scatter.F90)

ENDIF()

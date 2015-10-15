ADD_EXECUTABLE( test_buffer_loader_utilities test_buffer_loader_utilities.F90 )
TARGET_LINK_LIBRARIES( test_buffer_loader_utilities sll_parallel_utilities )

IF(HDF5_ENABLED AND HDF5_PARALLEL_ENABLED)
  ADD_EXECUTABLE( test_parallel_array_initializer 
    test_parallel_array_initializer.F90 
    sll_parallel_array_output_module.F90 )

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

ADD_LIBRARY( sll_remap STATIC sll_remap.F90 )
TARGET_LINK_LIBRARIES( sll_remap sll_collective ) 


# Ctests
IF(BUILD_TESTING)
  SET( ARGS " " )
  
  # 2D test
  ADD_EXECUTABLE( test_remap_2d unit_test_2d.F90 )
  TARGET_LINK_LIBRARIES( test_remap_2d sll_remap )
  ADD_MPI_TEST( remap_2d test_remap_2d ${PROCS} ${ARGS} )
  SET_TESTS_PROPERTIES( remap_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  # 3D test
  ADD_EXECUTABLE( test_remap_3d unit_test_3d.F90 )
  TARGET_LINK_LIBRARIES( test_remap_3d sll_remap )
  ADD_MPI_TEST( remap_3d test_remap_3d ${PROCS} ${ARGS} )
  SET_TESTS_PROPERTIES( remap_3d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  IF(PROCESSOR_COUNT GREATER 1)
    # 4D test
    ADD_EXECUTABLE(test_remap_4d unit_test_4d.F90 )
    TARGET_LINK_LIBRARIES(test_remap_4d sll_remap )
    ADD_MPI_TEST( remap_4d test_remap_4d ${PROCS} ${ARGS} )
    SET_TESTS_PROPERTIES( remap_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

    # 5D test
    ADD_EXECUTABLE( test_remap_5d unit_test_5d.F90 )
    TARGET_LINK_LIBRARIES( test_remap_5d sll_remap )
    ADD_MPI_TEST( remap_5d test_remap_5d ${PROCS} ${ARGS} )
    SET_TESTS_PROPERTIES( remap_5d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

    # 6D test
    ADD_EXECUTABLE( test_remap_6d unit_test_6d.F90 )
    TARGET_LINK_LIBRARIES( test_remap_6d sll_remap )
    ADD_MPI_TEST( remap_6d test_remap_6d ${PROCS} ${ARGS} )
    SET_TESTS_PROPERTIES( remap_6d PROPERTIES 
      PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 200 )

  ENDIF(PROCESSOR_COUNT GREATER 1)

  IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)
    # Layout output (not in Ctest)
    ADD_EXECUTABLE( test_layout_output unit_test_output.F90 )
    TARGET_LINK_LIBRARIES( test_layout_output sll_remap sll_file_io_parallel )
    
    # Parallel I/O
    ADD_EXECUTABLE( test_io_parallel unit_test_parallel.F90 )
    TARGET_LINK_LIBRARIES( test_io_parallel sll_file_io_parallel sll_remap )
    SET(PROCS 4)
    ADD_MPI_TEST( io_parallel test_io_parallel ${PROCS} ${ARGS} )
    SET_TESTS_PROPERTIES( io_parallel PROPERTIES
      PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20 )

  ENDIF()

ENDIF(BUILD_TESTING)

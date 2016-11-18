add_executable( test_io_serial_blocks test_serial_blocks.F90 )
target_link_libraries( test_io_serial_blocks sll_file_io_parallel )

set( ARGS " " )
set( PROCS  4 )
add_mpi_test( io_serial_blocks test_io_serial_blocks ${PROCS} ${ARGS} )
set_tests_properties( io_serial_blocks PROPERTIES
  FAIL_REGULAR_EXPRESSION "NOT PASS" TIMEOUT 20 )

if( HDF5_PARALLEL_ENABLED )
  add_executable( test_hdf5_io_parallel test_hdf5_io_parallel.F90)
  target_link_libraries( test_hdf5_io_parallel sll_file_io_parallel )
  set( ARGS " " )
  set( PROCS  4 )
  add_mpi_test( hdf5_io_parallel test_hdf5_io_parallel ${PROCS} ${ARGS} )
  set_tests_properties( hdf5_io_parallel PROPERTIES
    FAIL_REGULAR_EXPRESSION "ERROR" TIMEOUT 10 )
endif()

#ADD_EXECUTABLE(test_io_xdmf_parallel test_xdmf_parallel.F90)
#TARGET_LINK_LIBRARIES(test_io_xdmf_parallel sll_file_io_parallel )
#ADD_MPI_TEST( io_xdmf_parallel test_xdmf_parallel ${PROCS} ${ARGS} )
#SET_TESTS_PROPERTIES( io_xdmf_parallel PROPERTIES FAIL_REGULAR_EXPRESSION "NOT PASS" )

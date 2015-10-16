ADD_EXECUTABLE(test_io_serial_blocks test_serial_blocks.F90)
TARGET_LINK_LIBRARIES(test_io_serial_blocks sll_file_io_parallel )

SET(ARGS " ")
SET(PROCS  4)
ADD_MPI_TEST( io_serial_blocks test_io_serial_blocks ${PROCS} ${ARGS} )
SET_TESTS_PROPERTIES( io_serial_blocks PROPERTIES FAIL_REGULAR_EXPRESSION "NOT PASS" )

#ADD_EXECUTABLE(test_io_xdmf_parallel test_xdmf_parallel.F90)
#TARGET_LINK_LIBRARIES(test_io_xdmf_parallel sll_file_io_parallel )
#ADD_MPI_TEST( io_xdmf_parallel test_xdmf_parallel ${PROCS} ${ARGS} )
#SET_TESTS_PROPERTIES( io_xdmf_parallel PROPERTIES FAIL_REGULAR_EXPRESSION "NOT PASS" )

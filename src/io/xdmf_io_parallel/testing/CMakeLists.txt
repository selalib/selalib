# XDMF file creation using xdmf parallel type
ADD_EXECUTABLE( test_xdmf_parallel unit_test_xdmf_parallel.F90 )
TARGET_LINK_LIBRARIES( test_xdmf_parallel
  sll_xdmf_io_parallel
  sll_collective )

SET( XML_REFERENCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/ref1.xml )
ADD_MPI_TEST( xdmf_parallel test_xdmf_parallel ${PROCS} ${XML_REFERENCE_FILE} )
SET_TESTS_PROPERTIES( xdmf_parallel PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
  TIMEOUT 10 )

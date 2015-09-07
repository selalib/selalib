# Test 1.) XML file creation (using xml type)
ADD_EXECUTABLE( test_xml unit_test_xml.F90 )
TARGET_LINK_LIBRARIES( test_xml sll_xdmf )

SET( XML_REFERENCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/ref1.xml )
ADD_TEST( NAME xml COMMAND test_xml ${XML_REFERENCE_FILE} )
SET_TESTS_PROPERTIES( xml PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
  TIMEOUT 1 )


# Test 2.) XDMF file creation (identical to test 1, but uses xdmf type)
SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ADD_EXECUTABLE( test_xdmf_parallel unit_test_xdmf_parallel.F90 )
TARGET_LINK_LIBRARIES( test_xdmf_parallel
  sll_working_precision
  sll_collective
  sll_xdmf )

ADD_TEST( NAME xdmf_parallel COMMAND test_xdmf_parallel ${XML_REFERENCE_FILE} )
SET_TESTS_PROPERTIES( xdmf_parallel PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
  TIMEOUT 1 )

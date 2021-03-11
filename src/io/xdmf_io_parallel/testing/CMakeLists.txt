# XDMF file creation using xdmf parallel type
add_executable(test_xdmf_parallel unit_test_xdmf_parallel.F90)
target_link_libraries(test_xdmf_parallel sll_xdmf_io_parallel sll_collective)

set(XML_REFERENCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/ref1.xml)
add_mpi_test(xdmf_parallel test_xdmf_parallel ${PROCS} ${XML_REFERENCE_FILE})
set_tests_properties(xdmf_parallel PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
                                              TIMEOUT 10)

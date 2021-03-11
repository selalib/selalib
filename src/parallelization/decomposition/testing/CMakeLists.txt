set(ARGS " ")

add_executable(test_decomposition test_decomposition.F90)
target_link_libraries(test_decomposition sll_decomposition sll_collective)
add_mpi_test(decomposition test_decomposition ${PROCS} ${ARGS})
set_tests_properties(decomposition PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

# IF(PROCESSOR_COUNT GREATER 1) ENDIF(PROCESSOR_COUNT GREATER 1)
#
# IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED) ENDIF()

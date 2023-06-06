SET(ARGS " ")

ADD_EXECUTABLE( test_decomposition test_decomposition.F90 )
TARGET_LINK_LIBRARIES( test_decomposition sll_decomposition sll_collective sll_utilities )
ADD_MPI_TEST( decomposition test_decomposition ${PROCS} ${ARGS} )
SET_TESTS_PROPERTIES( decomposition PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )


ADD_EXECUTABLE( test_decomposition_6d3d test_decomposition_6d3d.F90 )
TARGET_LINK_LIBRARIES( test_decomposition_6d3d sll_decomposition sll_collective sll_utilities )
#ADD_MPI_TEST( decomposition_6d3d test_decomposition_6d3d ${PROCS} ${ARGS} )
#SET_TESTS_PROPERTIES( decomposition_6d3d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )


# --- temporarily disabled due to linking problems with poisson
#ADD_EXECUTABLE( test_decomposition_6d3d_remap test_decomposition_6d3d_remap.F90 )
#TARGET_LINK_LIBRARIES( test_decomposition_6d3d_remap sll_decomposition sll_collective sll_remap sll_poisson_solvers sll_poisson_parallel )
#ADD_MPI_TEST( decomposition_6d3d_remap test_decomposition_6d3d_remap ${PROCS} ${ARGS} )
#SET_TESTS_PROPERTIES( decomposition_6d3d_remap PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )


ADD_EXECUTABLE( test_decomposition_slim test_decomposition_slim.F90 )
TARGET_LINK_LIBRARIES( test_decomposition_slim sll_decomposition sll_collective sll_utilities )


ADD_EXECUTABLE( test_decomposition_advanced test_decomposition_advanced.F90 )
TARGET_LINK_LIBRARIES( test_decomposition_advanced sll_decomposition sll_collective sll_utilities )


#  IF(PROCESSOR_COUNT GREATER 1)
#  ENDIF(PROCESSOR_COUNT GREATER 1)
#
#  IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)
#  ENDIF()

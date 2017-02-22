IF(HDF5_IS_PARALLEL)
SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL)
ADD_EXECUTABLE(test_distribution_function test_distribution_function.F90)
TARGET_LINK_LIBRARIES(test_distribution_function sll_distribution_function)
ADD_TEST( NAME distribution_function COMMAND test_distribution_function )


ADD_EXECUTABLE( test_initial_distribution test_initial_distribution.F90)
TARGET_LINK_LIBRARIES( test_initial_distribution sll_distribution_function sll_xdmf_io)

SET( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/initial_distribution.nml )
ADD_TEST( NAME initial_distribution COMMAND test_initial_distribution ${NML_FILE})
SET_TESTS_PROPERTIES(initial_distribution PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

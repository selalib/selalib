ADD_EXECUTABLE( test_lagrange_fast test_lagrange_fast.F90 )
TARGET_LINK_LIBRARIES( test_lagrange_fast sll_lagrange_fast )
ADD_TEST( NAME lagrange_fast COMMAND test_lagrange_fast )
SET_TESTS_PROPERTIES( lagrange_fast PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )


IF( MPI_ENABLED )
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   ADD_EXECUTABLE( test_lagrange_fast_parallel test_lagrange_fast_parallel.F90 )
   TARGET_LINK_LIBRARIES( test_lagrange_fast_parallel sll_lagrange_fast sll_decomposition sll_collective )
   ADD_TEST( NAME lagrange_fast_parallel COMMAND test_lagrange_fast_parallel )
   SET_TESTS_PROPERTIES( lagrange_fast_parallel PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF( MPI_ENABLED )

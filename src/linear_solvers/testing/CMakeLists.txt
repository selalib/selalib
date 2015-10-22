ADD_EXECUTABLE( test_tridiagonal test_tridiagonal.F90 )
TARGET_LINK_LIBRARIES( test_tridiagonal sll_linear_solvers )
ADD_TEST( NAME tridiagonal COMMAND test_tridiagonal )

ADD_EXECUTABLE( test_toeplitz_penta_diagonal test_toeplitz_penta_diagonal.F90 )
TARGET_LINK_LIBRARIES( test_toeplitz_penta_diagonal sll_linear_solvers )
ADD_TEST( NAME toeplitz_penta_diagonal COMMAND test_toeplitz_penta_diagonal )
SET_TESTS_PROPERTIES( toeplitz_penta_diagonal PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( example_tridiagonal example_tridiagonal.F90 )
TARGET_LINK_LIBRARIES( example_tridiagonal sll_linear_solvers )

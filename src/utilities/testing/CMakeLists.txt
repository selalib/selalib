# General test
ADD_EXECUTABLE( test_utilities test_utilities.F90 )
TARGET_LINK_LIBRARIES( test_utilities sll_utilities )
ADD_TEST( NAME utilities COMMAND test_utilities )

# Tri-diagonal test
ADD_EXECUTABLE( test_tridiagonal test_tridiagonal.F90 )
TARGET_LINK_LIBRARIES( test_tridiagonal sll_utilities )
ADD_TEST( NAME tridiagonal COMMAND test_tridiagonal )

# Penta-diagonal test
ADD_EXECUTABLE( test_toeplitz_penta_diagonal test_toeplitz_penta_diagonal.F90 )
TARGET_LINK_LIBRARIES( test_toeplitz_penta_diagonal sll_utilities )
ADD_TEST( NAME toeplitz_penta_diagonal COMMAND test_toeplitz_penta_diagonal )
SET_TESTS_PROPERTIES( toeplitz_penta_diagonal PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_fornberg test_fornberg.F90 )
TARGET_LINK_LIBRARIES( test_fornberg sll_utilities )
ADD_TEST( NAME fornberg COMMAND test_fornberg )

# Example (tridiagonal solver)
ADD_EXECUTABLE( example_tridiagonal example_tridiagonal.F90 )
TARGET_LINK_LIBRARIES( example_tridiagonal sll_utilities )

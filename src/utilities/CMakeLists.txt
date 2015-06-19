SET( SRC_UTILITIES
  sll_utilities.F90
  tridiagonal.F90
  sll_toeplitz_penta_diagonal.F90
  sll_primitives.F90
  choleski.F90
  fornberg.F90
  kiss.F90
  sol.f
  muls.f )

IF(Fortran_COMPILER MATCHES IBM)
   SET_SOURCE_FILES_PROPERTIES( sol.f muls.f 
     PROPERTIES COMPILE_FLAGS '-qfixed' )
ENDIF()

ADD_LIBRARY( sll_utilities STATIC ${SRC_UTILITIES} )
TARGET_LINK_LIBRARIES( sll_utilities
  sll_errors sll_memory sll_assert sll_working_precision )

# Example (tridiagonal solver)
ADD_EXECUTABLE( example_tridiagonal example.F90 )
TARGET_LINK_LIBRARIES( example_tridiagonal sll_utilities )

# Example (tridiagonal solver)
ADD_EXECUTABLE( test_fornberg test_fornberg.F90 )
TARGET_LINK_LIBRARIES( test_fornberg sll_utilities sll_constants )

# Ctests
IF(BUILD_TESTING)

  # General test
  ADD_EXECUTABLE( test_utilities unit_test.F90 )
  TARGET_LINK_LIBRARIES( test_utilities sll_utilities )
  ADD_TEST( NAME utilities COMMAND test_utilities )

  # Tri-diagonal test
  ADD_EXECUTABLE( test_tridiagonal unit_test_tridiag.F90 )
  TARGET_LINK_LIBRARIES( test_tridiagonal sll_utilities )
  ADD_TEST( NAME tridiagonal COMMAND test_tridiagonal )
  
  # Penta-diagonal test
  ADD_EXECUTABLE( test_toeplitz_penta_diagonal test_sll_toeplitz_penta_diagonal.F90 )
  TARGET_LINK_LIBRARIES( test_toeplitz_penta_diagonal sll_utilities )
  ADD_TEST( NAME toeplitz_penta_diagonal COMMAND test_toeplitz_penta_diagonal )
  SET_TESTS_PROPERTIES( toeplitz_penta_diagonal
    PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ENDIF(BUILD_TESTING)

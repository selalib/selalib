SET(SRC_UTILITIES sll_utilities.F90 
                  tridiagonal.F90
                  sll_toeplitz_penta_diagonal.F90
                  lagrange_interpolation.F90)

IF(Fortran_COMPILER MATCHES IBM)
ELSE()
   SET(SRC_UTILITIES ${SRC_UTILITIES} sol.f muls.f )
ENDIF()

ADD_LIBRARY(sll_utilities STATIC ${SRC_UTILITIES})

TARGET_LINK_LIBRARIES(sll_utilities sll_memory sll_assertion sll_working_precision)

IF(NOT STDF95)
   ADD_EXECUTABLE(test_utilities unit_test.F90)
   TARGET_LINK_LIBRARIES(test_utilities sll_utilities )
ENDIF()

#  penta diagonal test
ADD_EXECUTABLE(test_toeplitz_penta_diagonal test_sll_toeplitz_penta_diagonal.F90)
TARGET_LINK_LIBRARIES(test_toeplitz_penta_diagonal sll_utilities) 

#  tridiagonal test
ADD_EXECUTABLE(test_tridiagonal unit_test_tridiag.F90)
TARGET_LINK_LIBRARIES(test_tridiagonal sll_utilities)

#  tridiagonal example
ADD_EXECUTABLE(example_tridiagonal example.F90)
TARGET_LINK_LIBRARIES(example_tridiagonal sll_utilities)

# lagrange interpolator
ADD_EXECUTABLE( test_lagrange unit_test_lagrange.F90 )
TARGET_LINK_LIBRARIES( test_lagrange sll_utilities )

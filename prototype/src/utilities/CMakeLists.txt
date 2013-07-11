ADD_LIBRARY(sll_utilities STATIC sll_utilities.F90 
                                 tridiagonal.F90
                                 sll_toeplitz_penta_diagonal.F90
                                 lagrange_interpolation.F90
				         sol.f
				         muls.f )

TARGET_LINK_LIBRARIES(sll_utilities sll_memory sll_assertion sll_working_precision)

IF(NOT STDF95)
  add_executable(test_utilities unit_test.F90)
  target_link_libraries(test_utilities sll_utilities )
ENDIF()

#  penta diagonal test
add_executable(test_toeplitz_penta_diagonal test_sll_toeplitz_penta_diagonal.F90)
target_link_libraries(test_toeplitz_penta_diagonal sll_utilities) 

#  tridiagonal test
add_executable(test_tridiagonal unit_test_tridiag.F90)
target_link_libraries(test_tridiagonal sll_utilities)

#  tridiagonal example
add_executable(example_tridiagonal example.F90)
target_link_libraries(example_tridiagonal sll_utilities)

# lagrange interpolator
add_executable( test_lagrange unit_test_lagrange.F90 )
target_link_libraries( test_lagrange sll_utilities )

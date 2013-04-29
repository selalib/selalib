ADD_LIBRARY(sll_utilities STATIC sll_utilities.F90 
                                 tridiagonal.F90
                                 sll_toeplitz_penta_diagonal.F90
                                 lagrange_interpolation.F90
				 sol.f
				 muls.f )

ADD_DEPENDENCIES(sll_utilities sll_memory
                               sll_assertion
                               sll_working_precision)

IF(NOT STDF95)
  add_executable(test_utilities unit_test.F90)
  target_link_libraries(test_utilities sll_utilities sll_working_precision)
ENDIF()

#  penta diagonal test
add_library(sll_toeplitz_penta_diagonal STATIC sll_toeplitz_penta_diagonal.F90)
add_executable(test_toeplitz_penta_diagonal test_sll_toeplitz_penta_diagonal.F90)
target_link_libraries(test_toeplitz_penta_diagonal sll_toeplitz_penta_diagonal 
                                         sll_memory sll_working_precision)

#  tridiagonal test
add_library(sll_tridiagonal STATIC tridiagonal.F90)
add_dependencies(sll_tridiagonal sll_working_precision)
add_executable(test_tridiagonal unit_test_tridiag.F90)
target_link_libraries(test_tridiagonal sll_tridiagonal
	                               sll_memory
	                               sll_working_precision)

#  tridiagonal example
add_executable(example_tridiagonal example.F90)
target_link_libraries(example_tridiagonal sll_tridiagonal
	                                  sll_memory
	                                  sll_working_precision)

# lagrange interpolator
add_library( lagrange_interpolation STATIC lagrange_interpolation.F90 )
add_dependencies( lagrange_interpolation sll_working_precision sll_assertion)
add_executable( test_lagrange unit_test_lagrange.F90 )
target_link_libraries( test_lagrange lagrange_interpolation sll_working_precision sll_assertion )

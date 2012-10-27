SET(NUMERICAL_UTILITIES_SRC tridiagonal.F90
                            sll_toeplitz_penta_diagonal.F90
                            lagrange_interpolation.F90 )

ADD_LIBRARY(sll_numerical_utilities STATIC ${NUMERICAL_UTILITIES_SRC})

ADD_DEPENDENCIES(sll_numerical_utilities sll_memory
                                         sll_assertion
                                         sll_working_precision
                                         sll_utilities)



ADD_LIBRARY(sll_tridiagonal STATIC tridiagonal.F90)
add_dependencies(sll_tridiagonal sll_working_precision)

add_executable(test_tridiagonal unit_test.F90)

target_link_libraries(test_tridiagonal sll_tridiagonal 
                     sll_memory sll_working_precision)

set(SPD_SRC sll_toeplitz_penta_diagonal.F90)
add_library(sll_toeplitz_penta_diagonal STATIC ${SPD_SRC})

add_executable(test_toeplitz_penta_diagonal test_sll_toeplitz_penta_diagonal.F90)
target_link_libraries(test_toeplitz_penta_diagonal sll_toeplitz_penta_diagonal 
                                         sll_memory sll_working_precision)

target_link_libraries(test_tridiagonal
	sll_tridiagonal
	sll_memory
	sll_working_precision
)

add_executable(example_tridiagonal example.F90)
target_link_libraries(example_tridiagonal
	sll_tridiagonal
	sll_memory
	sll_working_precision
)

# lagrange interpolator
add_library( lagrange_interpolation STATIC lagrange_interpolation.F90 )
add_dependencies( lagrange_interpolation
  sll_working_precision
  sll_assertion
  )

add_executable( test_lagrange unit_test_lagrange.F90 )
target_link_libraries( test_lagrange 
  lagrange_interpolation 
  sll_working_precision 
  sll_assertion 
  )


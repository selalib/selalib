set(TRIDIAGONAL_SRC tridiagonal.F90)
add_library(sll_tridiagonal STATIC tridiagonal.F90)
add_dependencies(sll_tridiagonal sll_working_precision)

add_executable(test_tridiagonal unit_test.F90)

target_link_libraries(test_tridiagonal sll_tridiagonal sll_memory sll_working_precision)

set(SPD_SRC toep_penta_diagonal.c)
add_library(toep_penta_diagonal STATIC ${SPD_SRC})

add_executable(toep_penta_diagonal_tester toep_penta_diagonal_tester.c)
target_link_libraries(toep_penta_diagonal_tester toep_penta_diagonal)

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


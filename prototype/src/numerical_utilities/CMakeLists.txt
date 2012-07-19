set(TRIDIAGONAL_SRC tridiagonal.F90)
add_library(sll_tridiagonal STATIC ${TRIDIAGONAL_SRC})

add_executable(test_tridiagonal unit_test.F90)
target_link_libraries(test_tridiagonal sll_tridiagonal sll_memory)

add_executable(example_tridiagonal example.F90)
target_link_libraries(example_tridiagonal sll_tridiagonal sll_memory)

set(SPD_SRC sym_penta_diagonal.c)
add_library(sym_penta_diagonal STATIC ${SPD_SRC})

add_executable(sym_penta_diagonal_tester sym_penta_diagonal_tester.c)
target_link_libraries(sym_penta_diagonal_tester sym_penta_diagonal)

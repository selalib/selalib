# use all sll_*F90 files
add_library(sll_splines STATIC sll_splines.F90)

set(SPL_LIBS sll_splines sll_tridiagonal sll_memory sll_assertion sll_constants)

add_library(util_constants STATIC util_constants.F90)

# add first test
set(SRC_TEST1 unit_test.F90 test_func_module.F90 test_processes_module.F90 util_constants.F90)
add_executable(test_splines ${SRC_TEST1})
target_link_libraries(test_splines ${SPL_LIBS})

# add second test
set(SRC_TEST2 unit_test_nonuniform.F90 cubic_nonuniform_splines.F90)
add_executable(test_non_unif_splines ${SRC_TEST2})
target_link_libraries(test_non_unif_splines ${SPL_LIBS})

# third test: arbitrary order splines
set(SRC_TEST3 unit_test_arbitrary_degree.F90 sll_arbitrary_degree_splines.F90)
add_executable(test_arbitrary_degree_splines ${SRC_TEST3})
target_link_libraries(test_arbitrary_degree_splines ${SPL_LIBS})

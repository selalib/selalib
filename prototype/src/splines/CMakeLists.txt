# use all sll_*F90 files
add_library(sll_splines STATIC sll_splines.F90)
add_library(sll_quintic_splines STATIC quintic_spline.c)
add_library(sll_odd_degree_splines STATIC sll_odd_degree_splines.F90)
add_library(sll_arbitrary_degree_splines STATIC sll_arbitrary_degree_splines.F90)

set(SPL_LIBS sll_splines sll_quintic_splines sll_odd_degree_splines 
sll_arbitrary_degree_splines sll_toeplitz_penta_diagonal sll_tridiagonal 
sll_memory sll_assertion sll_constants util_constants ${LAPACK_LIBRARIES})

add_dependencies(sll_splines sll_arbitrary_degree_splines sll_memory sll_assertion sll_working_precision)

add_library(cubic_non_unif_splines STATIC cubic_nonuniform_splines.F90)
add_dependencies(cubic_non_unif_splines sll_memory sll_assertion sll_working_precision
sll_tridiagonal sll_splines)


add_library(util_constants STATIC util_constants.F90)
add_dependencies(util_constants sll_constants sll_working_precision)


# add uniform splines test
set(SRC_TEST1 unit_test.F90 test_func_module.F90 test_processes_module.F90 )
add_executable(test_splines ${SRC_TEST1})
target_link_libraries(test_splines ${SPL_LIBS} )

# add non uniform splines test
set(SRC_TEST2 unit_test_nonuniform.F90 cubic_nonuniform_splines.F90)
add_executable(test_non_unif_splines ${SRC_TEST2})
target_link_libraries(test_non_unif_splines ${SPL_LIBS})

# third test: arbitrary degree splines
set(SRC_TEST3 unit_test_arbitrary_degree.F90 sll_arbitrary_degree_splines.F90)
add_executable(test_arbitrary_degree_splines ${SRC_TEST3})
target_link_libraries(test_arbitrary_degree_splines ${SPL_LIBS})

# add quintic splines test
set(SRC_TEST4 test_sll_quintic_splines.F90 sll_quintic_splines.F90)
add_executable(test_quintic_splines ${SRC_TEST4})
target_link_libraries(test_quintic_splines ${SPL_LIBS})

# add any odd degree splines test
set(SRC_TEST5 test_sll_odd_degree_splines.F90 sll_odd_degree_splines.F90)
add_executable(test_odd_degree_splines ${SRC_TEST5})
target_link_libraries(test_odd_degree_splines ${SPL_LIBS} -llapack)


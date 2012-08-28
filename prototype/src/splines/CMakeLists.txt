# use all sll_*F90 files
add_library(sll_splines STATIC sll_splines.F90)
add_library(quintic_spline STATIC quintic_spline.c)
add_library(sll_odd_order_splines STATIC sll_odd_order_splines.F90)
add_library(util_constants STATIC util_constants.F90)


set(SPL_LIBS sll_splines quintic_spline sll_odd_order_splines
toep_penta_diagonal sll_tridiagonal sll_memory sll_assertion sll_constants util_constants )


# add uniform splines test
set(SRC_TEST1 unit_test.F90 test_func_module.F90 test_processes_module.F90 )
add_executable(test_splines ${SRC_TEST1})
target_link_libraries(test_splines ${SPL_LIBS} )

# add non uniform splines test
set(SRC_TEST2 unit_test_nonuniform.F90 cubic_nonuniform_splines.F90)
add_executable(test_non_unif_splines ${SRC_TEST2})
target_link_libraries(test_non_unif_splines ${SPL_LIBS})

# add quintic splines test
set(SRC_TEST3 quintic_spline_tester.c quintic_spline.c)
add_executable(quintic_spline_tester ${SRC_TEST3})
target_link_libraries(quintic_spline_tester ${SPL_LIBS})

# add any odd order splines test
set(SRC_TEST4 test_sll_odd_order_splines.F90 sll_odd_order_splines.F90)
add_executable(test_sll_odd_order_splines ${SRC_TEST4})
target_link_libraries(test_sll_odd_order_splines ${SPL_LIBS})

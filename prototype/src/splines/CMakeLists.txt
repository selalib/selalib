# use all sll_*F90 files
set(SPLINES_SRC "sll_splines.F90")
add_library(sll_splines STATIC ${SPLINES_SRC})
set(SPL_LIBS sll_memory sll_tridiagonal sll_working_precision sll_assertion sll_constants)
# add first test
set(SRC_TEST1 unit_test.F90 test_func_module.F90 test_processes_module.F90 util_constants.F90)
add_executable(test_splines ${SRC_TEST1})
target_link_libraries(test_splines sll_splines ${SPL_LIBS})
# add first test
set(SRC_TEST2 unit_test_nonuniform.F90 cubic_nonuniform_splines.F90)
add_executable(test_non_unif_splines ${SRC_TEST2})
target_link_libraries(test_non_unif_splines sll_splines ${SPL_LIBS})


# lagrange
add_executable(test_lagrange test_lagrange.F90)
target_link_libraries(test_lagrange sll_lagrange_interpolation)
add_test(NAME lagrange COMMAND test_lagrange)
# SET_TESTS_PROPERTIES( lagrange PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

# lagrange_interpolation
add_executable(test_lagrange_interpolation test_lagrange_interpolation.F90)
target_link_libraries(test_lagrange_interpolation sll_lagrange_interpolation)
add_test(NAME lagrange_interpolation COMMAND test_lagrange_interpolation)
set_tests_properties(lagrange_interpolation PROPERTIES PASS_REGULAR_EXPRESSION
                                                       "PASSED")

# lagrange_fast (serial)
add_executable(test_lagrange_interpolation_1d_fast
               test_lagrange_interpolation_1d_fast.F90)
target_link_libraries(test_lagrange_interpolation_1d_fast
                      sll_lagrange_interpolation)
add_test(NAME lagrange_interpolation_1d_fast
         COMMAND test_lagrange_interpolation_1d_fast)
set_tests_properties(lagrange_interpolation_1d_fast
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

# lagrange fast (parallel)
if(MPI_ENABLED)
  add_executable(test_lagrange_interpolation_1d_fast_parallel
                 test_lagrange_interpolation_1d_fast_parallel.F90)
  target_link_libraries(
    test_lagrange_interpolation_1d_fast_parallel sll_lagrange_interpolation
    sll_decomposition sll_collective)
  add_mpi_test(lagrange_interpolation_1d_fast_parallel
               test_lagrange_interpolation_1d_fast_parallel "1" "")
  set_tests_properties(lagrange_interpolation_1d_fast_parallel
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
endif(MPI_ENABLED)

add_executable(bsl_1d_cubic_uniform_periodic test_1d_cubic_uniform_periodic.F90)
target_link_libraries(bsl_1d_cubic_uniform_periodic sll_interpolators)

add_executable(bsl_1d_cubic_uniform_compact test_1d_cubic_uniform_compact.F90)
target_link_libraries(bsl_1d_cubic_uniform_compact sll_interpolators)

add_executable(bsl_1d_cubic_nonuniform_periodic
               test_1d_cubic_nonuniform_periodic.F90)
target_link_libraries(bsl_1d_cubic_nonuniform_periodic sll_interpolators)

add_executable(bsl_1d_cubic_nonuniform_compact
               test_1d_cubic_nonuniform_compact.F90)
target_link_libraries(bsl_1d_cubic_nonuniform_compact sll_interpolators)

add_executable(bsl_2d_omp unit_test_2d_omp.F90)
target_link_libraries(bsl_2d_omp sll_interpolators)

add_executable(bsl_landau_1d landau_1d.F90)
target_link_libraries(bsl_landau_1d sll_interpolators sll_poisson_solvers)

add_executable(bsl_landau_2d landau_2d.F90)
target_link_libraries(bsl_landau_2d sll_interpolators sll_poisson_solvers
                      sll_maxwell_solvers)

add_executable(bsl_advection_2d_nufft bsl_advection_2d_nufft.F90)
target_link_libraries(bsl_advection_2d_nufft sll_file_io sll_fft
                      sll_nufft_interpolation)

if(BUILD_MUDPACK)
  add_executable(bsl_2d_multigrid landau_2d_multigrid.F90)
  target_link_libraries(bsl_2d_multigrid sll_interpolators sll_poisson_solvers
                        mudpack)
endif(BUILD_MUDPACK)

add_executable(bsl_2d_maxwell landau_2d_maxwell.F90)
target_link_libraries(bsl_2d_maxwell sll_interpolators sll_poisson_solvers
                      sll_maxwell_solvers)

add_test(NAME bsl_1d_uniform_periodic COMMAND bsl_1d_cubic_uniform_periodic)
add_test(NAME bsl_1d_uniform_compact COMMAND bsl_1d_cubic_uniform_compact)
# ADD_TEST( NAME bsl_1d_nonuniform_periodic COMMAND
# bsl_1d_cubic_nonuniform_periodic ) # ISSUE
add_test(NAME bsl_1d_nonuniform_compact COMMAND bsl_1d_cubic_nonuniform_compact)

set_tests_properties(bsl_1d_uniform_periodic PROPERTIES PASS_REGULAR_EXPRESSION
                                                        "PASSED")
set_tests_properties(bsl_1d_uniform_compact PROPERTIES PASS_REGULAR_EXPRESSION
                                                       "PASSED")
# SET_TESTS_PROPERTIES( bsl_1d_nonuniform_periodic PROPERTIES
# PASS_REGULAR_EXPRESSION "PASSED" ) # ISSUE
set_tests_properties(bsl_1d_nonuniform_compact
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

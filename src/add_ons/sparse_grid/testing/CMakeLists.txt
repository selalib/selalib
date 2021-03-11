add_executable(test_sparse_grid_2d test_sparse_grid_2d.F90)
target_link_libraries(test_sparse_grid_2d sll_sparse_grid sll_interpolators
                      sll_periodic_interpolation)
add_test(NAME sparse_grid_2d COMMAND test_sparse_grid_2d)
set_tests_properties(sparse_grid_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_executable(test_sparse_grid_3d test_sparse_grid_3d.F90)
target_link_libraries(test_sparse_grid_3d sll_sparse_grid sll_interpolators
                      sll_periodic_interpolation)
add_test(NAME sparse_grid_3d COMMAND test_sparse_grid_3d)
set_tests_properties(sparse_grid_3d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_executable(test_sparse_grid_4d test_sparse_grid_4d.F90)
target_link_libraries(test_sparse_grid_4d sll_sparse_grid sll_interpolators
                      sll_periodic_interpolation)
add_test(NAME sparse_grid_4d COMMAND test_sparse_grid_4d)
set_tests_properties(sparse_grid_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_executable(test_poisson_2d_sparse_grid_fft
               test_poisson_2d_sparse_grid_fft.F90)
target_link_libraries(test_poisson_2d_sparse_grid_fft sll_poisson_solvers
                      sll_sparse_grid ${FFTW_LIBRARIES})
add_test(NAME poisson_2d_sparse_grid_fft
         COMMAND test_poisson_2d_sparse_grid_fft)
set_tests_properties(poisson_2d_sparse_grid_fft
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_executable(test_poisson_3d_sparse_grid_fft
               test_poisson_3d_sparse_grid_fft.F90)
target_link_libraries(test_poisson_3d_sparse_grid_fft sll_poisson_solvers
                      sll_sparse_grid ${FFTW_LIBRARIES})
add_test(NAME poisson_3d_sparse_grid_fft
         COMMAND test_poisson_3d_sparse_grid_fft)
set_tests_properties(poisson_3d_sparse_grid_fft
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

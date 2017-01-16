# Test 1
add_executable( test_poisson_2d_polar_par
  test_poisson_2d_polar_par.F90
  m_test_case_poisson_par_2d_base.F90
  m_test_case_poisson_par_2d_dirichlet.F90
  m_test_case_poisson_par_2d_neumann_mode0.F90 )
target_link_libraries( test_poisson_2d_polar_par sll_poisson_solvers_parallel )

ADD_MPI_TEST( poisson_2d_polar_par test_poisson_2d_polar_par ${PROCS} " " )
set_tests_properties( poisson_2d_polar_par
  PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20 )


if( HDF5_PARALLEL_ENABLED AND HDF5_IS_PARALLEL )

  # Test 2
  add_executable( test_poisson_2d_periodic_par test_poisson_2d_periodic_par.F90 )
  target_link_libraries( test_poisson_2d_periodic_par
    sll_poisson_solvers_parallel
    sll_file_io_parallel )

  ADD_MPI_TEST( poisson_2d_periodic_par test_poisson_2d_periodic_par ${PROCS} " " )
  set_tests_properties( poisson_2d_periodic_par
    PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20 )

  # Test 3
  add_executable( test_poisson_3d_periodic_par test_poisson_3d_periodic_par.F90 )
  target_link_libraries( test_poisson_3d_periodic_par
    sll_poisson_solvers_parallel
    sll_file_io_parallel )

  ADD_MPI_TEST( poisson_3d_periodic_par test_poisson_3d_periodic_par ${PROCS} " " )
  set_tests_properties( poisson_3d_periodic_par
    PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20 )

endif()

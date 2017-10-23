# Save .mod files locally (this does not pollute the main "modules" directory)
set( CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )

# Test 1
add_executable( test_poisson_2d_polar_par
  test_poisson_2d_polar_par.F90
  m_test_poisson_2d_polar_base.F90
  m_test_poisson_2d_polar_dirichlet.F90
  m_test_poisson_2d_polar_neumann_mode0.F90
  m_test_poisson_2d_polar_circle_dirichlet.F90 )
target_link_libraries( test_poisson_2d_polar_par sll_poisson_solvers_parallel )

ADD_MPI_TEST( poisson_2d_polar_par_np1 test_poisson_2d_polar_par 1 " " )
ADD_MPI_TEST( poisson_2d_polar_par_np4 test_poisson_2d_polar_par 4 " " )
set_tests_properties(
  poisson_2d_polar_par_np1
  poisson_2d_polar_par_np4
  PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20 )


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

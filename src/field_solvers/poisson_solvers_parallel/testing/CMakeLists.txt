add_executable( test_poisson_2d_polar_par test_poisson_2d_polar_par.F90 )
target_link_libraries( test_poisson_2d_polar_par sll_poisson_parallel )

add_executable( test_poisson_2d_periodic_par test_poisson_2d_periodic_par.F90 )
target_link_libraries( test_poisson_2d_periodic_par sll_poisson_parallel )

add_executable( test_poisson_3d_periodic_par test_poisson_3d_periodic_par.F90 )
target_link_libraries( test_poisson_3d_periodic_par sll_poisson_parallel )

# CTest
if(HDF5_PARALLEL_ENABLED AND HDF5_IS_PARALLEL)
  set( ARGS " " )
  ADD_MPI_TEST( poisson_2d_polar_par    test_poisson_2d_polar_par    ${PROCS} ${ARGS} )
  ADD_MPI_TEST( poisson_2d_periodic_par test_poisson_2d_periodic_par ${PROCS} ${ARGS} )
  ADD_MPI_TEST( poisson_3d_periodic_par test_poisson_3d_periodic_par ${PROCS} ${ARGS} )
  set_tests_properties(
    poisson_2d_polar_par
    poisson_2d_periodic_par
    poisson_3d_periodic_par 
    PROPERTIES
    PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20 )
endif()

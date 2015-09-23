# ****************************************************************************
#
#                           PARALLEL POISSON SOLVERS
#
# ****************************************************************************

# This test has parallel diagnostics. Maybe these should only be turned on
# when debugging this specific module.
IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)

  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

  SET( POISSON_PARALLEL_SRCS
    sll_poisson_polar_parallel.F90 
    sll_poisson_3d_periodic_par.F90 
    sll_poisson_2d_periodic_cartesian_par.F90
    sll_quasi_neutral_3d_solver.F90 )

  SET( POISSON_PARALLEL_SRCS ${POISSON_PARALLEL_SRCS} 
    sll_poisson_2d_polar_parallel_solver.F90 )

  ADD_LIBRARY( sll_poisson_parallel STATIC ${POISSON_PARALLEL_SRCS} )
  TARGET_LINK_LIBRARIES( sll_poisson_parallel
    sll_remap 
    sll_fft
    sll_file_io_parallel
    sll_collective
    sll_utilities
    sll_poisson_solvers )

  # Tests
  IF(BUILD_TESTING)

    ADD_EXECUTABLE( test_poisson_3d_periodic_par test_poisson_3d_periodic_par.F90 )
    TARGET_LINK_LIBRARIES( test_poisson_3d_periodic_par sll_poisson_parallel )

    ADD_EXECUTABLE( test_poisson_2d_per_cart_par 
      test_poisson_2d_periodic_cartesian_par.F90 )
    TARGET_LINK_LIBRARIES( test_poisson_2d_per_cart_par sll_poisson_parallel )

    ADD_EXECUTABLE( test_poisson_polar_parallel test_poisson_polar_parallel.F90 )
    TARGET_LINK_LIBRARIES( test_poisson_polar_parallel sll_poisson_parallel )

    ADD_EXECUTABLE( test_poisson_polar_parallel_solver 
      unit_test_poisson_polar_parallel_solver.F90 )
    TARGET_LINK_LIBRARIES( test_poisson_polar_parallel_solver sll_poisson_parallel )

  ENDIF(BUILD_TESTING)

ENDIF()

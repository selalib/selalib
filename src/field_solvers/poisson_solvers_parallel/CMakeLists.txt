# This test has parallel diagnostics. Maybe these should only be turned on
# when debugging this specific module.
IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)

  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

  SET( POISSON_PARALLEL_SRCS
    sll_m_poisson_2d_polar_par.F90 
    sll_m_poisson_3d_periodic_par.F90 
    sll_m_poisson_2d_periodic_cartesian_par.F90 )

  SET( POISSON_PARALLEL_SRCS ${POISSON_PARALLEL_SRCS} 
    sll_m_poisson_2d_polar_parallel_solver.F90 )

  ADD_LIBRARY( sll_poisson_parallel STATIC ${POISSON_PARALLEL_SRCS} )
  TARGET_LINK_LIBRARIES( sll_poisson_parallel
    sll_remap 
    sll_fft
    sll_file_io_parallel
    sll_collective
    sll_utilities
    sll_poisson_solvers )

  IF(BUILD_TESTING)
    ADD_SUBDIRECTORY(testing)
  ENDIF(BUILD_TESTING)

ENDIF()

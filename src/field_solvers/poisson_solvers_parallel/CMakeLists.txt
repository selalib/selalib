# This test has parallel diagnostics. Maybe these should only be turned on
# when debugging this specific module.
if( HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED )

  set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )

  set( POISSON_PARALLEL_SRCS
    sll_m_poisson_2d_periodic_par.F90
    sll_m_poisson_2d_polar_par.F90
    sll_m_poisson_2d_polar_parallel_solver.F90
    sll_m_poisson_3d_periodic_par.F90 )

  add_library( sll_poisson_parallel STATIC ${POISSON_PARALLEL_SRCS} )
  target_link_libraries( sll_poisson_parallel
    sll_remap 
    sll_fft
    sll_file_io_parallel
    sll_collective
    sll_utilities
    sll_poisson_solvers )

  if( BUILD_TESTING )
    add_subdirectory( testing )
  endif( BUILD_TESTING )

endif()

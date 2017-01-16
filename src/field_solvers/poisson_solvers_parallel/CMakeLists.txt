set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )

set( POISSON_PARALLEL_SRCS
  sll_m_poisson_2d_periodic_par.F90
  sll_m_poisson_2d_polar_par.F90
  sll_m_poisson_3d_periodic_par.F90 )

add_library( sll_poisson_solvers_parallel STATIC ${POISSON_PARALLEL_SRCS} )
target_link_libraries( sll_poisson_solvers_parallel
  sll_assert
  sll_errors
  sll_working_precision
  sll_boundary_condition_descriptors
  sll_fft
  sll_linear_solvers
  sll_utilities
  sll_collective
  sll_remap )

# Ctests
if( BUILD_TESTING )
  add_subdirectory( testing )
endif( BUILD_TESTING )

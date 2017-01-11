# Parallel library
set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )

add_library( sll_quasi_neutral_solvers_parallel STATIC
  sll_m_qn_solver_2d_polar_par.F90
  sll_m_qn_solver_3d_polar_par.F90 )
  
target_link_libraries( sll_quasi_neutral_solvers_parallel
  sll_constants
  sll_boundary_condition_descriptors
  sll_fft
  sll_linear_solvers
  sll_collective
  sll_remap )

# Ctests
if( BUILD_TESTING )
  add_subdirectory( testing )
endif( BUILD_TESTING )

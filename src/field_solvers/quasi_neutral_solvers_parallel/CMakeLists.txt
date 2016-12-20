# Parallel library
SET( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )
ADD_LIBRARY( sll_quasi_neutral_solvers_parallel STATIC
  sll_m_qn_solver_2d_polar_par.F90
  sll_m_qn_solver_2d_parallel.F90
  sll_m_qn_solver_3d_polar_parallel_x1.F90
  sll_m_qn_solver_3d_polar_parallel_x1_wrapper.F90 )
  
TARGET_LINK_LIBRARIES( sll_quasi_neutral_solvers_parallel
  sll_remap sll_fft sll_poisson_solvers sll_quasi_neutral_solvers)


# Ctests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)

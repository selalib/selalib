add_library( sll_quasi_neutral_solvers STATIC sll_m_qn_solver_2d_polar.F90 )

target_link_libraries( sll_quasi_neutral_solvers
  sll_assert
  sll_errors
  sll_working_precision
  sll_boundary_condition_descriptors
  sll_constants
  sll_fft
  sll_linear_solvers )

# Ctests
if( BUILD_TESTING )
  add_subdirectory( testing )
endif( BUILD_TESTING )

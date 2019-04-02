add_library( sll_quasi_neutral_solvers STATIC
  sll_m_qn_solver_2d_polar.F90
  sll_m_qn_solver_2d_fem_sps_weak_form.F90
  sll_m_qn_solver_2d_fem_sps_stencil_new_assembler.F90
  sll_m_qn_solver_2d_fem_sps_stencil_new.F90 )

target_link_libraries( sll_quasi_neutral_solvers
  sll_assert
  sll_errors
  sll_working_precision
  sll_boundary_condition_descriptors
  sll_constants
  sll_fft
  sll_linear_solvers
  sll_general_coordinate_elliptic_solvers )

# Ctests
if( BUILD_TESTING )
  add_subdirectory( testing )
endif( BUILD_TESTING )

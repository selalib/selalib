add_library(
  sll_splitting_methods STATIC
  sll_m_hamiltonian_splitting_base.F90
  sll_m_time_splitting_coeff.F90
  sll_m_time_splitting.F90
  sll_m_operator_splitting.F90
  sll_m_linear_pendulum_operators.F90
  sll_m_const_coef_advection_2d.F90
  sll_m_time_composition.F90
  sll_m_split_advection_2d.F90
  sll_m_vp_cartesian_2d.F90)

target_link_libraries(
  sll_splitting_methods
  sll_distribution_function
  sll_interpolators
  sll_advection
  sll_poisson_solvers
  sll_memory
  sll_ode_solvers)

# Tests
if(BUILD_TESTING)

  add_subdirectory(testing)

endif(BUILD_TESTING)

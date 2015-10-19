ADD_LIBRARY( sll_operator_splitting STATIC
  sll_m_hamiltonian_splitting_base.F90
  sll_m_time_splitting.F90 
  sll_m_operator_splitting.F90
  sll_m_linear_pendulum_operators.F90
  sll_m_const_coef_advection_2d.F90
  sll_m_split_advection_2d.F90
  sll_m_vp_cartesian_2d.F90 )

TARGET_LINK_LIBRARIES( sll_operator_splitting
  sll_distribution_function
  sll_interpolators
  sll_advection
  sll_poisson_solvers
  ode_solvers )

ADD_LIBRARY( sll_time_splitting STATIC
  sll_m_time_splitting_coeff.F90 )

ADD_LIBRARY( sll_time_composition STATIC
  sll_m_time_composition.F90 )

  
TARGET_LINK_LIBRARIES( sll_time_composition 
                       sll_working_precision
                       sll_memory)

  
# Tests
IF(BUILD_TESTING)

  ADD_SUBDIRECTORY(testing)


ENDIF(BUILD_TESTING)

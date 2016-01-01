add_library( sll_ode_integrators STATIC
             sll_m_ode_integrator_base.F90
             sll_m_rk_explicit.F90 )

target_link_libraries( sll_ode_integrators
                       sll_vector_space
                       sll_errors )
# Tests
if( BUILD_TESTING )
  add_subdirectory( testing )
endif( BUILD_TESTING )

add_library( sll_ode_integrators STATIC
             sll_ode_integrator_base.F90
             rk_explicit.F90 )

target_link_libraries( sll_ode_integrators
                       sll_vector_space
                       sll_working_precision
                       sll_errors )
# Tests
if( BUILD_TESTING )
  add_subdirectory( testing )
endif( BUILD_TESTING )

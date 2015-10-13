  ADD_EXECUTABLE( test_ode_solvers 
    test_ode_solvers.F90 )
  TARGET_LINK_LIBRARIES( test_ode_solvers ode_solvers sll_memory )
  ADD_TEST( NAME ode_solvers COMMAND test_ode_solvers )

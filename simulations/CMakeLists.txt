set( SLL_LIB selalib )
set( INSTALL_DIR "${CMAKE_INSTALL_PREFIX_SIMS}/bin/selalib_simulations" )
set( SIM_BASE_CLASS
  "${CMAKE_SOURCE_DIR}/src/simulation_tools/simulation_base/sll_m_sim_base_class.F90" )

# Parallel simulations
if( MPI_ENABLED )
  set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )
  add_subdirectory( parallel )
endif()

# Serial simulations
add_subdirectory( serial )

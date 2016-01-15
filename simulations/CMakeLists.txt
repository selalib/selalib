SET( INSTALL_DIR "${CMAKE_INSTALL_PREFIX_SIMS}/simulations" )
SET( SIM_BASE_CLASS
  "${CMAKE_SOURCE_DIR}/src/simulation_tools/simulation_base/sll_m_sim_base_class.F90" )

ADD_DEFINITIONS( -DMUDPACK )

# Parallel simulations
IF( MPI_ENABLED )
  SET( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )
  ADD_SUBDIRECTORY( parallel )
ENDIF()

# Serial simulations
ADD_SUBDIRECTORY( serial )

set( SLL_LIB selalib )
set( INSTALL_DIR "${CMAKE_INSTALL_PREFIX_SIMS}/bin/selalib_simulations" )
set( SIM_BASE_CLASS
  "${CMAKE_SOURCE_DIR}/src/simulation_tools/simulation_base/sll_m_sim_base_class.F90" )

include_directories(BEFORE ${HDF5_INCLUDE_DIR_FORTRAN})

# Define MUDPACK macro if MUDPACK is built
IF(BUILD_MUDPACK)
  ADD_DEFINITIONS(-DMUDPACK)
ENDIF()

# Serial simulations
add_subdirectory( serial )

# Parallel simulations
if( MPI_ENABLED )
  add_subdirectory( parallel )
endif()


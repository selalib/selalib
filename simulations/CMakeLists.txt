set( SLL_LIB selalib )
set( INSTALL_DIR "${CMAKE_INSTALL_PREFIX_SIMS}/bin/selalib_simulations" )
set( SIM_BASE_CLASS
  "${CMAKE_SOURCE_DIR}/src/simulation_tools/simulation_base/sll_m_sim_base_class.F90" )

if(HDF5_IS_PARALLEL)
  set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )
endif(HDF5_IS_PARALLEL)
include_directories(BEFORE ${HDF5_INCLUDE_DIR_FORTRAN})

# Serial simulations
add_subdirectory( serial )

# Parallel simulations
if( MPI_ENABLED )
  add_subdirectory( parallel )
endif()


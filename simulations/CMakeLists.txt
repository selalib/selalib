set(SLL_LIB selalib)
set(INSTALL_DIR "${CMAKE_INSTALL_PREFIX_SIMS}/bin/selalib_simulations")
set(SIM_BASE_CLASS
    "${CMAKE_SOURCE_DIR}/src/simulation_tools/simulation_base/sll_m_sim_base_class.F90"
)

include_directories(BEFORE ${HDF5_Fortran_INCLUDE_DIRS})

# Define MUDPACK macro if MUDPACK is built
if(BUILD_MUDPACK)
  add_definitions(-DMUDPACK)
endif()

# Serial simulations
add_subdirectory(serial)

# Parallel simulations
if(MPI_ENABLED)
  add_subdirectory(parallel)
endif()

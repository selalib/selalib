SET( INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/simulations" )
SET( SIM_BASE_CLASS
  "${CMAKE_SOURCE_DIR}/src/simulation_tools/simulation_base/sll_m_sim_base_class.F90" )

ADD_DEFINITIONS( -DMUDPACK )

ADD_SUBDIRECTORY( serial )
ADD_SUBDIRECTORY( parallel )

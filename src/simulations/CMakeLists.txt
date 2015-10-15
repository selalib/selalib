SET( INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/simulations" )

ADD_DEFINITIONS( -DMUDPACK )

ADD_SUBDIRECTORY(simulations_sequential)
ADD_SUBDIRECTORY(simulations_parallel  )
ADD_SUBDIRECTORY(eulerian-finite-volume)

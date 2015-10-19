ADD_SUBDIRECTORY( quadrature )
ADD_SUBDIRECTORY( time_integration )

ADD_SUBDIRECTORY(low_level_utilities)
ADD_SUBDIRECTORY(interfaces)
ADD_SUBDIRECTORY(utilities)
ADD_SUBDIRECTORY(linear_solvers)
ADD_SUBDIRECTORY(splines)
ADD_SUBDIRECTORY(mesh)
#ADD_SUBDIRECTORY(pic_particle_types)
#ADD_SUBDIRECTORY(pic_particle_initializers)
#ADD_SUBDIRECTORY(pic_accumulators)
#ADD_SUBDIRECTORY(pic_utilities)
#ADD_SUBDIRECTORY(field_accumulators)
ADD_SUBDIRECTORY(interpolation)
ADD_SUBDIRECTORY(field_solvers)
ADD_SUBDIRECTORY(io)

# Container for all particle methods (mostly PIC)
ADD_SUBDIRECTORY( particle_methods )

IF(MPI_ENABLED)
   
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   ADD_SUBDIRECTORY(parallelization)

ENDIF()

ADD_SUBDIRECTORY(data_structures)

ADD_SUBDIRECTORY(add_ons)

ADD_SUBDIRECTORY(namelist_initialization)    



IF(BUILD_SIMULATIONS)
  GET_FILENAME_COMPONENT( INSTALL_DIR_SIMS "${CMAKE_SOURCE_DIR}/usr" ABSOLUTE )
  SET(CMAKE_INSTALL_PREFIX_SIMS ${INSTALL_DIR_SIMS} CACHE PATH 
    "simulations install dir")
  ADD_SUBDIRECTORY(simulations)
ENDIF(BUILD_SIMULATIONS)

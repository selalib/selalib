ADD_SUBDIRECTORY(low_level_utilities)
ADD_SUBDIRECTORY(interfaces)
ADD_SUBDIRECTORY(utilities)
ADD_SUBDIRECTORY(linear_solvers)
ADD_SUBDIRECTORY(splines)
ADD_SUBDIRECTORY(mesh)
ADD_SUBDIRECTORY(meshes)
ADD_SUBDIRECTORY(integration)
ADD_SUBDIRECTORY(reduction)
#ADD_SUBDIRECTORY(pic_particle_types)
ADD_SUBDIRECTORY(random_deviate_generators)
#ADD_SUBDIRECTORY(pic_particle_initializers)
#ADD_SUBDIRECTORY(pic_accumulators)
#ADD_SUBDIRECTORY(pic_utilities)
#ADD_SUBDIRECTORY(field_accumulators)
ADD_SUBDIRECTORY(interpolation)
ADD_SUBDIRECTORY(characteristics)
ADD_SUBDIRECTORY(advection)
ADD_SUBDIRECTORY(fcisl)
ADD_SUBDIRECTORY(gyroaverage)
ADD_SUBDIRECTORY(field_solvers)
ADD_SUBDIRECTORY(descriptors)
ADD_SUBDIRECTORY(io)

#Fortran version is checked in CMakeModules/FortranConfig.cmake
IF(F2008)
  ADD_SUBDIRECTORY( vector_space )
  ADD_SUBDIRECTORY( ode_integrators )
ENDIF(F2008)

# Container for all particle methods (mostly PIC)
ADD_SUBDIRECTORY( particle_methods )

IF(MPI_ENABLED)
   
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   ADD_SUBDIRECTORY(parallelization)
   ADD_SUBDIRECTORY(pif)

ENDIF()

ADD_SUBDIRECTORY(data_structures)
ADD_SUBDIRECTORY(time_solvers)
ADD_SUBDIRECTORY(ode_solvers)

ADD_SUBDIRECTORY(add_ons)

ADD_SUBDIRECTORY(namelist_initialization)    



IF(BUILD_SIMULATIONS)
  GET_FILENAME_COMPONENT( INSTALL_DIR_SIMS "${CMAKE_SOURCE_DIR}/usr" ABSOLUTE )
  SET(CMAKE_INSTALL_PREFIX_SIMS ${INSTALL_DIR_SIMS} CACHE PATH 
    "simulations install dir")
  ADD_SUBDIRECTORY(simulations)
ENDIF(BUILD_SIMULATIONS)

IF( MPI_ENABLED )
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF( MPI_ENABLED )

ADD_SUBDIRECTORY( add_ons )
ADD_SUBDIRECTORY( data_structures )
ADD_SUBDIRECTORY( field_solvers )
ADD_SUBDIRECTORY( initialization_tools )
ADD_SUBDIRECTORY( interfaces )
ADD_SUBDIRECTORY( interpolation )
ADD_SUBDIRECTORY( io )
ADD_SUBDIRECTORY( linear_solvers )
ADD_SUBDIRECTORY( low_level_utilities )
ADD_SUBDIRECTORY( mesh )
ADD_SUBDIRECTORY( particle_methods ) # Container for all particle methods (mostly PIC)
ADD_SUBDIRECTORY( quadrature )
ADD_SUBDIRECTORY( semi_lagrangian )
ADD_SUBDIRECTORY( time_integration )
ADD_SUBDIRECTORY( splines )

# Parallelization tools
IF( MPI_ENABLED )
   ADD_SUBDIRECTORY( parallelization )
ENDIF( MPI_ENABLED )

# Simulations
IF( BUILD_SIMULATIONS )
  GET_FILENAME_COMPONENT( INSTALL_DIR_SIMS "${CMAKE_SOURCE_DIR}/usr" ABSOLUTE )
  SET(CMAKE_INSTALL_PREFIX_SIMS ${INSTALL_DIR_SIMS} CACHE PATH 
    "simulations install dir")
  ADD_SUBDIRECTORY( simulations )
ENDIF( BUILD_SIMULATIONS )

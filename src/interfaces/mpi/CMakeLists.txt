# Detect MPI library. A similar if-elseif-else block is placed in "sll_mpi.F90"
# TODO: use CMake variables to detect Intel-MPI, MPICH, OpenMPI
if( DEFINED __INTEL_COMPILER )
  set( MPI_LIB "intel-mpi" )
else()
  set( MPI_LIB "openmpi" )
endif()

# Specific operations for OpenMPI
if( MPI_LIB EQUAL "openmpi" )
  execute_process( COMMAND ompi_info OUTPUT_VARIABLE OMPI_INFO )
  string( REGEX MATCH "Open MPI: ([0-9])\\.([0-9])\\.([0-9])" OMPI_VERSION ${OMPI_INFO} )
  set( OMPI_V_MAJOR ${CMAKE_MATCH_1} )
  set( OMPI_V_MINOR ${CMAKE_MATCH_2} )
  set( OMPI_V_REVSN ${CMAKE_MATCH_3} )
  add_definitions( -DOMPI_V_MINOR=${OMPI_V_MINOR} )
endif()

#message( STATUS ${OMPI_INFO} )
#message( STATUS "OpenMPI version is: ${OMPI_V_MAJOR}.${OMPI_V_MINOR}.${OMPI_V_REVSN}" )
#message( STATUS "=======================================" )
#message( STATUS ${CMAKE_MATCH_0} )
#message( STATUS ${CMAKE_MATCH_1} )
#message( STATUS ${CMAKE_MATCH_2} )
#message( STATUS ${CMAKE_MATCH_3} )
#message( STATUS "=======================================" )

add_library( sll_mpi OBJECT sll_mpi.F90 )

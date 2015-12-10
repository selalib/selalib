message( STATUS "=========================================================")
message( STATUS "MPI_Fortran_FOUND         = ${MPI_Fortran_FOUND}"         )
message( STATUS "MPI_Fortran_COMPILER      = ${MPI_Fortran_COMPILER}"      )
message( STATUS "MPI_Fortran_COMPILE_FLAGS = ${MPI_Fortran_COMPILE_FLAGS}" )
message( STATUS "MPI_Fortran_INCLUDE_PATH  = ${MPI_Fortran_INCLUDE_PATH}"  )
message( STATUS "MPI_Fortran_LINK_FLAGS    = ${MPI_Fortran_LINK_FLAGS}"    )
message( STATUS "MPI_Fortran_LIBRARIES     = ${MPI_Fortran_LIBRARIES}"     )
message( STATUS "=========================================================")

#-------------------------------------------------------------------------------
# Detect MPI library. A similar if-elseif-else block is placed in "sll_mpi.F90"
#-------------------------------------------------------------------------------
# TODO: use CMake variables to detect Intel-MPI, MPICH, OpenMPI
if( Fortran_COMPILER_NAME MATCHES ifort )
  set( MPI_LIB "intel" )
else()
#  set( TEST_STRING ${MPI_Fortran_INCLUDE_PATH} )
#  set( TEST_STRING "/usr/include/mpich2/")
#  set( TEST_STRING "/usr/local/Cellar/open-mpi/1.10.1/include;/usr/local/Cellar/open-mpi/1.10.1/lib" )
#  string( REGEX MATCH "open\\-?mpi|mpich" MPI_LIB ${TEST_STRING} )
  string( REGEX MATCH "open\\-?mpi|mpich" MPI_LIB ${MPI_Fortran_INCLUDE_PATH} )
  string( REPLACE "-" "" MPI_LIB ${MPI_LIB} )
endif()

message( STATUS "MPI_LIB = ${MPI_LIB}" )

#-------------------------------------------------------------------------------
# Specific operations for Intel MPI
#-------------------------------------------------------------------------------
if( ${MPI_LIB} STREQUAL "intel" )
  add_definitions( -DINTEL_MPI )
endif()

#-------------------------------------------------------------------------------
# Specific operations for OpenMPI
#-------------------------------------------------------------------------------
if( ${MPI_LIB} STREQUAL "openmpi" )
  find_program(OMPI_INFO_PROG NAMES ompi_info HINTS /opt/local/libexec/openmpi-mp/ /usr/local/bin )
  execute_process( COMMAND ${OMPI_INFO_PROG} OUTPUT_VARIABLE OMPI_INFO )
  string( REGEX MATCH "Open MPI: ([0-9]+)\\.([0-9]+)\\.([0-9]+)" OMPI_VERSION ${OMPI_INFO} )
  message( STATUS ${CMAKE_MATCH_0} )
  set( OMPI_V_MAJOR ${CMAKE_MATCH_1} )
  set( OMPI_V_MINOR ${CMAKE_MATCH_2} )
  set( OMPI_V_REVSN ${CMAKE_MATCH_3} )
  add_definitions( -DOMPI )
  add_definitions( -DOMPI_V_MAJOR=${OMPI_V_MAJOR} )
  add_definitions( -DOMPI_V_MINOR=${OMPI_V_MINOR} )
  add_definitions( -DOMPI_V_REVSN=${OMPI_V_REVSN} )
endif()

#-------------------------------------------------------------------------------
# Specific operations for MPICH
#-------------------------------------------------------------------------------
if( ${MPI_LIB} STREQUAL "mpich" )
  execute_process(
    COMMAND mpichversion
    COMMAND grep "MPICH Version:\\s*[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+"
    OUTPUT_VARIABLE MPICH_INFO
    )
  message( STATUS ${MPICH_INFO} )
  string( REGEX MATCH "([0-9]+)\\.([0-9]+)\\.([0-9]+)" MPICH_VERSION ${MPICH_INFO} )
#  execute_process( COMMAND mpichversion OUTPUT_VARIABLE MPICH_INFO )
#  string( REGEX MATCH "MPICH Version:    	([0-9]+)\\.([0-9]+)\\.([0-9]+)" MPICH_VERSION ${MPICH_INFO} )
  set( MPICH_V_MAJOR ${CMAKE_MATCH_1} )
  set( MPICH_V_MINOR ${CMAKE_MATCH_2} )
  set( MPICH_V_REVSN ${CMAKE_MATCH_3} )
  add_definitions( -DMPICH )
  add_definitions( -DMPICH_V_MAJOR=${MPICH_V_MAJOR} )
  add_definitions( -DMPICH_V_MINOR=${MPICH_V_MINOR} )
  add_definitions( -DMPICH_V_REVSN=${MPICH_V_REVSN} )
endif()

#-------------------------------------------------------------------------------
#message( STATUS ${OMPI_INFO} )
#message( STATUS "OpenMPI version is: ${OMPI_V_MAJOR}.${OMPI_V_MINOR}.${OMPI_V_REVSN}" )
#message( STATUS "=======================================" )
#message( STATUS ${CMAKE_MATCH_0} )
#message( STATUS ${CMAKE_MATCH_1} )
#message( STATUS ${CMAKE_MATCH_2} )
#message( STATUS ${CMAKE_MATCH_3} )
#message( STATUS "=======================================" )

#-------------------------------------------------------------------------------
# ADD LIBRARY
#-------------------------------------------------------------------------------
add_library( sll_mpi OBJECT sll_mpi.F90 )

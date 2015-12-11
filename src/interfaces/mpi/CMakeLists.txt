message( STATUS "=========================================================")
message( STATUS "MPI_Fortran_FOUND         = ${MPI_Fortran_FOUND}"         )
message( STATUS "MPI_Fortran_COMPILER      = ${MPI_Fortran_COMPILER}"      )
message( STATUS "MPI_Fortran_COMPILE_FLAGS = ${MPI_Fortran_COMPILE_FLAGS}" )
message( STATUS "MPI_Fortran_INCLUDE_PATH  = ${MPI_Fortran_INCLUDE_PATH}"  )
message( STATUS "MPI_Fortran_LINK_FLAGS    = ${MPI_Fortran_LINK_FLAGS}"    )
message( STATUS "MPI_Fortran_LIBRARIES     = ${MPI_Fortran_LIBRARIES}"     )
message( STATUS "=========================================================")

#-------------------------------------------------------------------------------
# Detect MPI library: OpenMPI | MPICH | Intel MPI | Bull-X MPI
#-------------------------------------------------------------------------------
set( TEST_STRING ${MPI_Fortran_INCLUDE_PATH} )
string( REGEX MATCH "open\\-?mpi|mpich|impi|bullxmpi" MPI_LIB ${TEST_STRING} )
if( MPI_LIB )
  string( REPLACE "-" "" MPI_LIB ${MPI_LIB} )
else()
  set( MPI_LIB "unknown" )
endif()
message( STATUS "MPI_LIB = ${MPI_LIB}" )

#if( Fortran_COMPILER_NAME MATCHES ifort )
#  set( MPI_LIB "intel" )
#else()
#  set( TEST_STRING ${MPI_Fortran_INCLUDE_PATH} )
##  set( TEST_STRING "/usr/include/mpich2/")
##  set( TEST_STRING "/usr/local/Cellar/open-mpi/1.10.1/include;/usr/local/Cellar/open-mpi/1.10.1/lib" )
##  set( TEST_STRING  "/opt/mpi/bullxmpi/1.2.8.4/include;/opt/mpi/bullxmpi/1.2.8.4/lib" )
#  string( REGEX MATCH "open\\-?mpi|mpich" MPI_LIB ${TEST_STRING} )
#  if( MPI_LIB )
#    string( REPLACE "-" "" MPI_LIB ${MPI_LIB} )
#  else()
#    set( MPI_LIB "unknown" )
#  endif()
#endif()
#
#message( STATUS "MPI_LIB = ${MPI_LIB}" )

#-------------------------------------------------------------------------------
# Specific operations for Intel MPI
#-------------------------------------------------------------------------------
#if( ${MPI_LIB} STREQUAL "intel" )
if( ${MPI_LIB} STREQUAL "impi" )
  set( INTEL_MPI_INFO $ENV{I_MPI_ROOT} )
  string( REGEX MATCH "impi/([0-9]+)\\.([0-9]+)\\.([0-9]+)" INTEL_MPI_VERSION ${INTEL_MPI_INFO} )
  message( STATUS ${CMAKE_MATCH_0} )
  math( EXPR INTEL_MPI_VERSION
    "${CMAKE_MATCH_1}*10000 + ${CMAKE_MATCH_2}*100 + ${CMAKE_MATCH_3}" )
  add_definitions( -DINTEL_MPI )
  add_definitions( -DINTEL_MPI_VERSION=${INTEL_MPI_VERSION} )
  message( STATUS "INTEL_MPI_VERSION = ${INTEL_MPI_VERSION}" )
endif()

#-------------------------------------------------------------------------------
# Specific operations for OpenMPI
#-------------------------------------------------------------------------------
if( ${MPI_LIB} STREQUAL "openmpi" )
  find_program(OMPI_INFO_PROG NAMES ompi_info HINTS /opt/local/libexec/openmpi-mp/ /usr/local/bin )
  execute_process( COMMAND ${OMPI_INFO_PROG} OUTPUT_VARIABLE OMPI_INFO )
  string( REGEX MATCH "Open MPI: ([0-9]+)\\.([0-9]+)\\.([0-9]+)" OMPI_VERSION ${OMPI_INFO} )
  message( STATUS ${CMAKE_MATCH_0} )
  math( EXPR OMPI_VERSION
    "${CMAKE_MATCH_1}*10000 + ${CMAKE_MATCH_2}*100 + ${CMAKE_MATCH_3}" )
  add_definitions( -DOMPI )
  add_definitions( -DOMPI_VERSION=${OMPI_VERSION} )
  message( STATUS "OMPI_VERSION = ${OMPI_VERSION}" )
endif()

#-------------------------------------------------------------------------------
# Specific operations for MPICH
#-------------------------------------------------------------------------------
if( ${MPI_LIB} STREQUAL "mpich" )
#  execute_process( COMMAND mpichversion OUTPUT_VARIABLE MPICH_INFO )
#  string( REGEX MATCH "MPICH Version:    	([0-9]+)\\.([0-9]+)\\.([0-9]+)" MPICH_VERSION ${MPICH_INFO} )
  execute_process(
    COMMAND mpichversion
    COMMAND grep "MPICH Version:\\s*[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+"
    OUTPUT_VARIABLE MPICH_INFO
    )
  message( STATUS ${MPICH_INFO} )
  string( REGEX MATCH "([0-9]+)\\.([0-9]+)\\.([0-9]+)" MPICH_VERSION ${MPICH_INFO} )
  math( EXPR MPICH_VERSION
    "${CMAKE_MATCH_1}*10000 + ${CMAKE_MATCH_2}*100 + ${CMAKE_MATCH_3}" )
  add_definitions( -DMPICH )
  add_definitions( -DMPICH_VERSION=${MPICH_VERSION} )
  message( STATUS "MPICH_VERSION = ${MPICH_VERSION}" )
endif()

#-------------------------------------------------------------------------------
# ADD LIBRARY
#-------------------------------------------------------------------------------
add_library( sll_mpi OBJECT sll_mpi.F90 )

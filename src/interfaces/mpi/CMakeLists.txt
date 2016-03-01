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
#  set( TEST_STRING "/usr/include/mpich2/")
#  set( TEST_STRING "/usr/local/Cellar/open-mpi/1.10.1/include;/usr/local/Cellar/open-mpi/1.10.1/lib" )
#  set( TEST_STRING  "/opt/mpi/bullxmpi/1.2.8.4/include;/opt/mpi/bullxmpi/1.2.8.4/lib" )
set( TEST_STRING ${MPI_Fortran_INCLUDE_PATH} )
string( REGEX MATCH "open\\-?mpi|mpich|impi|bullxmpi" MPI_LIB ${TEST_STRING} )
if( MPI_LIB )
  string( REPLACE "-" "" MPI_LIB ${MPI_LIB} )
else()
  set( MPI_LIB "unknown" )
endif()
message( STATUS "MPI_LIB = ${MPI_LIB}" )

#-------------------------------------------------------------------------------
# Specific operations for Intel MPI
#-------------------------------------------------------------------------------
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
  find_program(MPICH_INFO_PROG NAMES mpichversion HINTS /opt/local/bin /usr/local/bin )
  message( STATUS "MPICH_INFO_PROG:${MPICH_INFO_PROG}" )
  if( MPICH_INFO_PROG )
    execute_process(
      COMMAND ${MPICH_INFO_PROG}
      COMMAND grep "MPICH Version:\\s*[0-9]\\+\\.[0-9]"#\\+\\.[0-9]\\+"
      OUTPUT_VARIABLE MPICH_INFO
      )
    #message( STATUS "MPICH_INFO:${MPICH_INFO}" )
    string( REGEX MATCH "([0-9]+)\\.([0-9])" MPICH_VERSION ${MPICH_INFO} )
    math( EXPR MPICH_VERSION
      "${CMAKE_MATCH_1}*10000 + ${CMAKE_MATCH_2}*100 ") #+ ${CMAKE_MATCH_3}" )
    add_definitions( -DMPICH )
    add_definitions( -DMPICH_VERSION=${MPICH_VERSION} )
    message( STATUS "MPICH_VERSION = ${MPICH_VERSION}" )
  else()
    message( "WARNING: cannot find 'mpichversion' executable." )
  endif()
endif()

#-------------------------------------------------------------------------------
# Specific operations for Bull-X MPI
#-------------------------------------------------------------------------------
if( ${MPI_LIB} STREQUAL "bullxmpi" )
  # TODO: extract version number
  add_definitions( -DBULLX_MPI )
endif()

#-------------------------------------------------------------------------------
# AUTOMATIC APPROACH: used if the above does not work!
#-------------------------------------------------------------------------------

# # TEST
# add_definitions( -DTEST -DTEST_VAL=2 )

# Get the new MPI compiler definitions
get_property( cur_defs DIRECTORY .  PROPERTY COMPILE_DEFINITIONS )
get_property( par_defs DIRECTORY .. PROPERTY COMPILE_DEFINITIONS )
list( REMOVE_ITEM cur_defs ${par_defs} )
foreach( def ${cur_defs} )
  list( APPEND MPI_DEFS -D${def} )
endforeach()
string( REPLACE ";" " " MPI_DEFS_STR "${MPI_DEFS}" )
message( STATUS ${MPI_DEFS_STR} )

# Variables
set( SRC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/sll_mpi.F90 )
set( ERR_FILE ${CMAKE_CURRENT_SOURCE_DIR}/mpi_errors.txt )
set( EXT_FILE ${CMAKE_CURRENT_SOURCE_DIR}/external_auto.F90 )
set( PYSCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/create_external_list.py )

# Test if 'sll_mpi' is correctly compiled
set( CMD ${MPI_Fortran_COMPILER} ${SRC_FILE} -I${MPI_Fortran_MOD_DIR} -c ${MPI_DEFS} )
string( REPLACE ";" " " CMD_MSG "${CMD}" )
message( STATUS ${CMD_MSG} )
execute_process( COMMAND ${CMD} ERROR_VARIABLE ERR_STAT )

# If there is a compilation problem, generate a list of external symbols
if( ERR_STAT )
  message( "" )
  message( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
  message( "WARNING: file 'sll_mpi.F90' cannot be compiled!"                   )
  message( "         fall back to automatic generation of external symbols..." )
  message( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
  message( "" )
  set( CMD ${MPI_Fortran_COMPILER} ${SRC_FILE} -c )
  string( REPLACE ";" " " CMD_MSG "${CMD}" )
  message( STATUS ${CMD_MSG} )
  execute_process( COMMAND ${CMD} ERROR_FILE ${ERR_FILE} )
  execute_process( COMMAND ${PYTHON_EXECUTABLE} ${PYSCRIPT} ${ERR_FILE} ${EXT_FILE} )
  remove_definitions( ${MPI_DEFS} )
  add_definitions( -DMPI_AUTO )
  get_property( cur_defs DIRECTORY .  PROPERTY COMPILE_DEFINITIONS )
  message( STATUS "${cur_defs}" )
else()
  message( STATUS "'sll_mpi.F90' is correctly compiled.. proceed with configuration" )
endif()

#-------------------------------------------------------------------------------
# ADD LIBRARY
#-------------------------------------------------------------------------------
add_library( sll_mpi STATIC sll_mpi.F90 )

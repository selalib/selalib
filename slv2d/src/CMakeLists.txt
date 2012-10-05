project(slv2d)
enable_language(Fortran)
cmake_minimum_required(VERSION 2.8)

# Determine how-to install the modules. CMAKE_BINARY_DIR is the directory
# in which the make command is invoked.
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# Add the modules directory to the list of include directories
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# Specify the output for binary
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Add new directory to find the FIND<PACKAGE>.cmake files. CMAKE_SOURCE_DIR is
# the directory where the project sources are located. CMAKE_MODULE_PATH is 
# the list of directories in which to search for CMake modules.
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules ${CMAKE_MODULE_PATH} )

# make sure that the default is a RELEASE. Here is where we need to set
# the parameters for the DEBUG or RELEASE versions. For instance, DEBUG should
# turn on assertions.
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# DGFORTRAN is used to set the way in which numbers get converted to strings
# by the preprocessor. gfotran uses "x" while other preprocessors use the
# cpp #x. Defining the GFORTRAN flag chooses the first.
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
IF(Fortran_COMPILER_NAME STREQUAL "gfortran")
  message(STATUS "gfortran compiler")
  add_definitions(-DGFORTRAN)
ENDIF()
add_definitions(-DDEBUG -DMPIF90)

# Define folder to be deleted by make clean
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/include/)
add_custom_target( clear
                   COMMAND find . -name *.mod -delete
		   COMMAND find . -name *.a -delete
		   COMMAND find . -name *.o -delete
		   COMMAND find . -name *.cmake -delete
		   COMMAND rm -rf bin
		   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
		   COMMENT "Delete all .a, .mod, .o and folders /modules, /bin" VERBATIM)

##########################################################
# Try to determine the compiler
try_run( RUN_RESULT_VAR
	 COMPILE_RESULT_VAR
         ${CMAKE_BINARY_DIR}
         ${CMAKE_SOURCE_DIR}/check_compiler.F90
)


set(CMAKE_Fortran_FLAGS "-llapack -lblas -g -Wall -cpp -pedantic -ffree-line-length-none -std=f2003 -fall-intrinsics")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
set(CMAKE_C_FLAGS "-lm")
set(CMAKE_C_FLAGS "-llapack -lblas")

set(HDF5_ENABLED ON CACHE BOOL "Use HDF5 format for data output ")
set(HDF5_PARALLEL_ENABLED OFF CACHE BOOL "Use Parallel HDF5")

IF(HDF5_ENABLED)
	IF( DEFINED ENV{HDF5_ROOT} )
		SEt(HDF5_ROOT $ENV{HDF5_ROOT})
	ENDIF()
	FIND_PACKAGE(HDF5 QUIET)
	IF(HDF5_FOUND)
		MESSAGE(STATUS "HDF5 FOUND")
		INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
	ENDIF()
ENDIF(HDF5_ENABLED)

# add the cache entry FFTPACK_ENABLED for enable/disable fftpack
SET(FFTPACK_ENABLED OFF CACHE BOOL " ")

FIND_PACKAGE(LAPACK REQUIRED)

GET_FILENAME_COMPONENT(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
IF(Fortran_COMPILER STREQUAL "GFORTRAN")
	set(CMAKE_Fortran_FLAGS "-g -Wall -cpp -pedantic -ffree-line-length-none -std=f2003 -fall-intrinsics")
ELSEIF(Fortran_COMPILER STREQUAL "INTEL")
	set(CMAKE_Fortran_FLAGS " ")
ELSEIF(Fortran_COMPILER_NAME STREQUAL "xlf")
	set(CMAKE_Fortran_FLAGS "-qextname=flush -qthreaded -qhalt=e")
ELSE()
	message(STATUS "NO KNOWN FORTRAN COMPILER FOUND")
ENDIF()

#add flag to disable MPI stuff for debug
SET(MPI_MODULE_ENABLED ON CACHE BOOL " ")

FIND_PACKAGE(MPI)
IF(MPI_FOUND)
   MESSAGE(STATUS "MPI FOUND")
   INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ELSE(MPI_FOUND)
   MESSAGE(STATUS "MPI NOT FOUND")
   SET(MPI_MODULE_ENABLED OFF CACHE BOOL " " FORCE)
ENDIF(MPI_FOUND)

ENABLE_TESTING()

ADD_TEST(NAME vm2d COMMAND test_vm2d)

IF(HDF5_PARALLEL_ENABLED AND MPI_MODULE_ENABLED)
   SET(PROCS 4)
   SET(ARGS "")
   ADD_TEST(NAME vp2d COMMAND 
            ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_vp2d
            ${MPIEXEC_POSTFLAGS} ${ARGS})
   SET_TESTS_PROPERTIES(vp2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF()

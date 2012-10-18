project(sll)
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
#set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/modules/)
#set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/bin/)
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


set(CMAKE_Fortran_FLAGS "-g -Wall -cpp -pedantic -ffree-line-length-none -std=f2003 -fall-intrinsics")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# COMPILE_RESULT_VAR is set to true if try_run succeed
# RUN_RESULT_VAR is a string that represent the exit status
# message(STATUS "TRY_RUN_STATUS : ${COMPILE_RESULT_VAR}, EXIT_STATUS : ${RUN_RESULT_VAR}")

set(NO_FORTRAN_2003 NO)
IF(COMPILE_RESULT_VAR)
	IF(${RUN_RESULT_VAR} STREQUAL 20 )
		message(STATUS "COMPILER IS GNU, FORTRAN 2003 NOT SUPPORTED")
		set(Fortran_COMPILER "GFORTRAN")
		set(NO_FORTRAN_2003 YES)
		add_definitions(-DNOF03SUPPORT)
	ELSEIF(${RUN_RESULT_VAR} STREQUAL 21 )
		message(STATUS "COMPILER IS GNU, FORTRAN 2003 SUPPORTED")
		set(Fortran_COMPILER "GFORTRAN")
	ELSEIF(${RUN_RESULT_VAR} STREQUAL 30 )
		message(STATUS "COMPILER IS INTEL, FORTRAN 2003 NOT SUPPORTED")
		set(NO_FORTRAN_2003 YES)
		set(Fortran_COMPILER "INTEL")
	ELSEIF(${RUN_RESULT_VAR} STREQUAL 31 )
		message(STATUS "COMPILER IS INTEL, FORTRAN 2003 SUPPORTED")
		set(Fortran_COMPILER "INTEL")
	ENDIF()
ELSE()
	message(STATUS "UNABLE TO DETERMINE WHICH COMPILER IS USED")
ENDIF()
##########################################################

# add the cache entry HDF5_ENABLED for enable/disable hdf5
set(HDF5_ENABLED ON CACHE BOOL "Use HDF5 format for data output ")
set(HDF5_PARALLEL_ENABLED OFF CACHE BOOL "Use Parallel HDF5")

IF (APPLE)
   FIND_LIBRARY(LAPACK_LIBRARIES  lapack)
   FIND_LIBRARY(BLAS_LIBRARIES    blas)
ELSE()
   FIND_PACKAGE(LAPACK REQUIRED)
   FIND_PACKAGE(BLAS   REQUIRED)
ENDIF()

SET(LINK_LIBRARIES  ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# search for HDF5
IF(HDF5_ENABLED)
	IF( DEFINED ENV{HDF5_ROOT} )
		set(HDF5_ROOT $ENV{HDF5_ROOT})
	ENDIF()
        #IF( DEFINED ENV{HDF5_HOME} )
	#	set(HDF5_ROOT $ENV{HDF5_HOME})
	#	message(STATUS "HDF5_HOME FOUND")
	#ENDIF()
	find_package(HDF5 QUIET)
	IF(HDF5_FOUND)
		message(STATUS "HDF5 FOUND")
	ELSE()
		add_definitions(-DNOHDF5)
		set(HDF5_ENABLED OFF CACHE BOOL " " FORCE)
		set(HDF5_LIBRARIES "")
	ENDIF()
ELSE(HDF5_ENABLED)
	add_definitions(-DNOHDF5)
ENDIF(HDF5_ENABLED)
##########################################################

##########################################################
# add the cache entry FFT_DEFAULT_LIBRARY to define the default library use by sll_fft
set(FFT_DEFAULT_LIBRARY SLLFFT CACHE STRING "specify the library use in sll_fft")


# add the cache entry FFTPACK_ENABLED for enable/disable fftpack
set(FFTPACK_ENABLED OFF CACHE BOOL " ")

# recherche le package FFTW
set(FFTW_ENABLED ON CACHE BOOL " ")

IF(FFTW_ENABLED)
	IF( DEFINED ENV{FFTW_ROOT} )
		set(FFTW_ROOT $ENV{FFTW_ROOT})
	ENDIF()
	find_package(FFTW QUIET)
	IF(FFTW_FOUND)
		message(STATUS "FFTW FOUND")
  		include_directories(${FFTW_INCLUDE_DIRS})
	ELSE(FFTW_FOUND)
		message(STATUS "FFTW NOT FOUND")
		set(FFTW_ENABLED NO CACHE BOOL " " FORCE)
	ENDIF(FFTW_FOUND)
ENDIF(FFTW_ENABLED)


find_package(LAPACK REQUIRED)

# recherche de MPI
#IF( DEFINED ENV{MPI_ROOT} )
#  set(MPI_ROOT $ENV{MPI_ROOT})
#  set(_MPI_BASE_DIR ${MPI_ROOT})

IF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTPACK" AND NOT FFTPACK_ENABLED)
	message("Please put on FFTPACK_ENABLED to use fftpack")

ENDIF()
IF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTW" AND NOT FFTW_FOUND)
	message("You can't use fftw library because it's not installed")
ELSEIF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTW" AND NOT FFTW_ENABLED)
	message("Please put on FFTW_ENABLED to use fftw")
ENDIF()
##########################################################
# find out which compiler we are using.
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
IF(Fortran_COMPILER STREQUAL "GFORTRAN")
	set(CMAKE_Fortran_FLAGS "-g -Wall -cpp -pedantic -ffree-line-length-none -std=f2003 -fall-intrinsics -fbounds-check")
ELSEIF(Fortran_COMPILER STREQUAL "INTEL")
	set(CMAKE_Fortran_FLAGS " ")
ELSEIF(Fortran_COMPILER_NAME STREQUAL "xlf")
	set(CMAKE_Fortran_FLAGS "-qextname=flush -qthreaded -qhalt=e")
ELSE()
	message(STATUS "NO KNOWN FORTRAN COMPILER FOUND")
ENDIF()

# this checks the existing collection of compiler flags to find if the
# fortran95 standard was requested. 
STRING(FIND ${CMAKE_Fortran_FLAGS} "-std=f95" VAR)
IF("${VAR}" STREQUAL "-1")
  set(STDF95 NO)
ELSE()
  set(STDF95 YES)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fmax-identifier-length=63")
  add_definitions(-DSTDF95)
ENDIF()

find_package(FFTW QUIET)
if (FFTW_FOUND)
	message(STATUS "FFTW FOUND")
	message(STATUS "libraries = " ${FFTW_LIBRARIES})
	message(STATUS "include = " ${FFTW_INCLUDE_DIRS})
	include_directories(${FFTW_INCLUDE_DIRS})
	set(FFT_ADD_MODULE fftpack_module fftw_module)
else (FFTW_FOUND)
	message(STATUS "FFTW NOT FOUND")
	add_definitions(-D_NOFFTW)
	set(FFT_ADD_MODULE fftpack_module)
endif(FFTW_FOUND)

# include's
include_directories(remap poisson_solvers utilities memory assert 
splines numerical_utilities precision fields) 


# pour la doc 


#add flag to disable MPI stuff for debug
set(MPI_MODULE_ENABLED ON CACHE BOOL " ")

# search for MPI
# Samuel DeSantis: added the line "HINTS $ENV{MPI_ROOT}" to the 
# find_program(MPIEXEC) in the FindMPI.cmake file (src/CMakeModules). 
# If specifying MPI_ROOT in your environment doesn't work, then try changing
# your PATH variable.

find_package(MPI)
IF(MPI_FOUND)
   message(STATUS "MPI FOUND")
   include_directories(${MPI_Fortran_INCLUDE_PATH})
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ELSE(MPI_FOUND)
   message(STATUS "MPI NOT FOUND")
   set(MPI_MODULE_ENABLED OFF CACHE BOOL " " FORCE)
ENDIF(MPI_FOUND)

# include's. Here we include the folders which contain '.h' files.
include_directories(fft remap poisson_solvers utilities memory assert precision fields) 

# for the documentation

find_package(Doxygen)
if(DOXYGEN_FOUND)
   get_filename_component(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR} PATH)
   set(DOXYGEN_OUTPUT_DIR "${DOXYGEN_OUTPUT_DIR}/documentation/")
   message(STATUS "The documentation is in ${DOXYGEN_OUTPUT_DIR}")
   configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in 
   ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
   add_custom_target(doc ${DOXYGEN_EXECUTABLE} 
   ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY 
   ${CMAKE_CURRENT_BINARY_DIR} 
   COMMENT "Generating API documentation with Doxygen" VERBATIM)
else(DOXYGEN_FOUND)
   MESSAGE(STATUS "DOXYGEN NOT FOUND")
endif(DOXYGEN_FOUND)

set(FISHPACK_ENABLED OFF CACHE BOOL " ")
set(MUDPACK_ENABLED OFF CACHE BOOL " ")

######################################################################
#  III. SUBDIRECTORY                                                 #
######################################################################

#low-level
add_subdirectory(memory)
add_subdirectory(precision)
add_subdirectory(assert)
add_subdirectory(constants)
add_subdirectory(utilities)
add_subdirectory(low_level_file_io)
add_subdirectory(timer)

#mid-level
add_subdirectory(numerical_utilities)
add_subdirectory(non_linear_solvers)
add_subdirectory(splines)
add_subdirectory(WENO)
add_subdirectory(integration)
IF(FFTPACK_ENABLED)
   add_subdirectory(fftpack)
ENDIF()

add_subdirectory(fft)

IF(MPI_MODULE_ENABLED)
   add_subdirectory(collective)
   add_subdirectory(remap)
ENDIF()

IF(NOT STDF95)
   add_subdirectory(poisson_solvers)
ENDIF()
add_subdirectory(electric_field_accumulators)

#high-level
add_subdirectory(interpolators)
add_subdirectory(mapped_meshes)
add_subdirectory(fields)
add_subdirectory(time_solvers)
add_subdirectory(ode_solvers)
add_subdirectory(distribution_function)
add_subdirectory(advection_field)
add_subdirectory(simulation)

IF(FISHPACK_ENABLED)
   add_subdirectory(fishpack)
ENDIF()

IF(MUDPACK_ENABLED)
   IF(MPI_MODULE_ENABLED)
	add_subdirectory(multigrid)
   ELSE(MPI_MODULE_ENABLED)
	add_subdirectory(mudpack)
   ENDIF(MPI_MODULE_ENABLED)
ENDIF(MUDPACK_ENABLED)
add_subdirectory(LinRood)
#add_subdirectory(maxwell_solvers) #use old mesh_types 
#add_subdirectory(BSL) top
add_subdirectory(CSL)
#add_subdirectory(quasi_neutral_solvers) don't compile problem with fftw3_thread! top
add_subdirectory(selalib_level EXCLUDE_FROM_ALL)
#add_custom_target(install make selalib_low_level)



######################################################################
#  IV. TEST                                                          #
######################################################################

enable_testing()

add_test(NAME memory COMMAND test_memory)
set_tests_properties(memory PROPERTIES TIMEOUT 20)

#add_test(precision precision/test_precision) #too long
set(passRegex "Assertion error triggered in file")

add_test(NAME assert COMMAND test_assert)
set_tests_properties(assert PROPERTIES PASS_REGULAR_EXPRESSION "${passRegex}")

add_test(NAME constants COMMAND test_constants)

add_test(NAME utilities COMMAND test_utilities)

add_test(NAME low_level_file_io COMMAND test_io)
set_tests_properties(low_level_file_io PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

IF(HDF5_PARALLEL_ENABLED AND MPI_MODULE_ENABLED)
   set(PROCS 4)
   set(ARGS "")
   add_test(NAME low_level_file_io_parallel COMMAND 
            ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_io_parallel
            ${MPIEXEC_POSTFLAGS} ${ARGS})
   set_tests_properties(low_level_file_io_parallel PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF()

add_test(NAME timer COMMAND test_timer)
add_test(NAME tridiagonal COMMAND test_tridiagonal)
add_test(NAME lagrange COMMAND test_lagrange)

add_test(NAME toeplitz_penta_diagonal COMMAND test_toeplitz_penta_diagonal)
set_tests_properties(toeplitz_penta_diagonal PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_test(NAME newton_raphson COMMAND test_newton_raphson)
add_test(NAME splines COMMAND test_splines #test_non_unif_splines 
  test_arbitrary_degree_splines)
add_test(NAME quintic_splines COMMAND test_quintic_splines)
add_test(NAME odd_degree_splines COMMAND test_odd_degree_splines)
set_tests_properties(splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
add_test(NAME integration COMMAND test_integration)
add_test(NAME fft COMMAND test_fft)
add_test(NAME electric_field_accumulators COMMAND test_e_field_accumulator_2d)
IF(FISHPACK_ENABLED)
	add_test(fishpack fishpack/thwscyl fishpack/thwsplr fishpack/thwsssp fishpack/tpois3d)
ENDIF()
IF(MUDPACK_ENABLED)
        IF(MPI_MODULE_ENABLED)
           set(PROCS 4)
           set(ARGS "")
	   add_test(NAME multigrid COMMAND 
	   ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	   ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_multigrid
	   ${MPIEXEC_POSTFLAGS} ${ARGS})
           set_tests_properties(multigrid PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
        ELSE()
	   add_test(NAME mudpack COMMAND tmud34sp tmud24sp )
           set_tests_properties(mudpack PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
        ENDIF()
ENDIF()

# we can call ${MPIEXEC_MAX_NUMPROCS} to have the number of procs
set(PROCS 2)
set(ARGS "")
add_test(NAME collective
	COMMAND
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_collective
	${MPIEXEC_POSTFLAGS} ${ARGS})
set_tests_properties(collective PROPERTIES FAIL_REGULAR_EXPRESSION "NOT PASS")

set(PROCS 8)
set(ARGS "")
add_test(NAME remap
	COMMAND
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap
	${MPIEXEC_POSTFLAGS} ${ARGS})
set_tests_properties(remap PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_test(NAME remap_2d
	COMMAND
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap_2d
	${MPIEXEC_POSTFLAGS} ${ARGS})
set_tests_properties(remap_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
# to make the test meaningful at higher dimensions, we need more processes...
set(PROCS 16)
add_test(NAME remap_4d
	COMMAND
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap_4d
	${MPIEXEC_POSTFLAGS} ${ARGS})
set_tests_properties(remap_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

set(PROCS 16)
add_test(NAME remap_6d
	COMMAND
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap_6d
	${MPIEXEC_POSTFLAGS} ${ARGS})
set_tests_properties(remap_6d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")


set(PROCS 16)
add_test(NAME poisson_periodic_cartesian_par_2d
         COMMAND
	 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_poisson_2d_per_cart_par
	 ${MPIEXEC_POSTFLAGS} ${ARGS})
set_tests_properties(remap_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_test(NAME WENO COMMAND test_WENO_interp test_WENO_recon)
add_test(NAME interpolators COMMAND test_interpolators_1d
				    test_interpolators_2d)
add_test(NAME mapped_meshes COMMAND test_mapped_meshes_1d
				    test_mapped_meshes_2d)
add_test(NAME fields COMMAND test_scalar_field)
add_test(NAME time_splitting COMMAND test_time_splitting)
add_test(NAME ode_solvers COMMAND test_implicit_ode_nonuniform)
add_test(NAME distribution_function COMMAND test_distribution_function)
add_test(NAME advection_field COMMAND test_advection_field)
IF(HDF5_ENABLE_PARALLEL AND HDF5_PARALLEL_ENABLED)
  set(PROCS 8)
  add_test(NAME test_vp4d_sim
         COMMAND
	 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_vp4d_sim
	 ${MPIEXEC_POSTFLAGS} ${ARGS})
  set_tests_properties(test_vp4d_sim PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF()
#add_test(maxwell_solvers maxwell_solvers/test_maxwell_2d)
   #add_test(quasi_neutral_solvers quasi_neutral_solvers/test_quasi_neutral_solver)
IF(FFTPACK_ENABLED)
	add_test(NAME poisson_solvers COMMAND test_poisson_1d)
ENDIF(FFTPACK_ENABLED)
#add_test(poisson_solvers poisson_solvers/test_poisson_2d)
#add_test(poisson_solvers poisson_solvers/test_poisson_3d_seq)
#add_test(poisson_solvers mpirun -np 2 poisson_solvers/test_poisson_3d_par)
   #add_test(poisson_solvers poisson_solvers/unit_test)
#add_test(BSL BSL/test_bsl)
#add_test(BSL BSL/test_landau_4d)
#add_test(BSL BSL/test_landau_2d)
   #add_test(CSL CSL/test_csl)
#add_test(simulation_1d simulation/test_simulation_1d)


set(FORTRANCL_ENABLED NO CACHE BOOL "Enable OpenCL test programs")
# recherche de FORTRANCL
IF(FORTRANCL_ENABLED)
	find_package(FORTRANCL QUIET)
	if (FORTRANCL_FOUND)
	   message(STATUS "FORTRANCL FOUND")
	   include_directories(${FORTRANCL_INCLUDE_DIRS})
           add_subdirectory(opencl)
           add_test(NAME opencl COMMAND test_opencl)
           set_tests_properties(opencl PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
	endif(FORTRANCL_FOUND)
ENDIF(FORTRANCL_ENABLED)

##########################################################
# add the cache entry PASTIX_ENABLED for enable/disable pastix
set(PASTIX_ENABLED NO CACHE BOOL " ")

IF(PASTIX_ENABLED)
   IF( DEFINED ENV{SCOTCH_ROOT} )
      set(SCOTCH_ROOT $ENV{SCOTCH_ROOT})
      set(PTSCOTCH_ROOT $ENV{SCOTCH_ROOT})
   ENDIF()
   IF( DEFINED ENV{PASTIX_ROOT} )
      set(PASTIX_ROOT $ENV{PASTIX_ROOT})
      set(MURGE_ROOT $ENV{PASTIX_ROOT})
   ENDIF()
   find_package(RT QUIET)
   find_package(SCOTCH QUIET)
   find_package(PTSCOTCH QUIET)
   find_package(MURGE QUIET)
   find_package(PASTIX QUIET)
   IF(PASTIX_FOUND AND PTSCOTCH_FOUND AND MURGE_FOUND AND SCOTCH_FOUND)
      message(STATUS "PASTIX FOUND")
      include_directories(${PASTIX_INCLUDE_DIRS})
      add_subdirectory(pastix)
      set(PROCS 4)
      set(ARGS 1000 3)
      add_test(NAME pastix 
               COMMAND
               ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_pastix ${ARGS})
      set_tests_properties(pastix PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ELSE()
      message(STATUS "PASTIX NOT FOUND")
   ENDIF()
ENDIF(PASTIX_ENABLED)
##########################################################


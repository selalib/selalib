project(sll)
enable_language(Fortran)
cmake_minimum_required(VERSION 2.8)

# Determine how-to install the modules
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
# Add the modules directory to the list of include directories
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# Specify the output for binary
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Add new directory to find the FIND<PACKAGE>.cmake files
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules ${CMAKE_MODULE_PATH} )

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Define folder to be delete by make clean
#set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/modules/)
#set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/bin/)
add_custom_target(clear
		COMMAND find . -name *.mod -delete
		COMMAND find . -name *.a -delete
		COMMAND find . -name *.o -delete
		COMMAND find . -name *.cmake -delete
		COMMAND rm -rf bin
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
		COMMENT "Delete all .a, .mod, .o and folders /modules, /bin" VERBATIM)

##########################################################
# Try to determine the compiler
try_run(RUN_RESULT_VAR
	COMPILE_RESULT_VAR
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/check_compiler.F90
)

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
	message(STATUS "IMPOSSIBLE TO DETERMINE WHICH COMPILER IS USED")
ENDIF()
##########################################################

# add the cache entry HDF5_ENABLED for enable/disable hdf5
set(HDF5_ENABLED ON CACHE BOOL " ")

# recherche de HDF5
IF(HDF5_ENABLED)
	IF( DEFINED ENV{HDF5_ROOT} )
		set(HDF5_ROOT $ENV{HDF5_ROOT})
	ENDIF()
	find_package(HDF5 QUIET)
	IF(HDF5_FOUND)
		message(STATUS "HDF5 FOUND")
		include_directories(${HDF5_INCLUDE_DIRS})
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

IF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTPACK" AND NOT FFTPACK_ENABLED)
	message("Please put on FFTPACK_ENABLED to use fftpack")
ENDIF()
IF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTW" AND NOT FFTW_FOUND)
	message("You can't use fftw library because it's not installed")
ELSEIF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTW" AND NOT FFTW_ENABLED)
	message("Please put on FFTW_ENABLED to use fftw")
ENDIF()
##########################################################

##########################################################
# add the cache entry PASTIX_ENABLED for enable/disable pastix
set(PASTIX_ENABLED ON CACHE BOOL " ")

IF(PASTIX_ENABLED)
   IF( DEFINED ENV{SCOTCH_ROOT} )
      set(SCOTCH_ROOT $ENV{SCOTCH_ROOT})
      set(PTSCOTCH_ROOT $ENV{SCOTCH_ROOT})
   ENDIF()
   IF( DEFINED ENV{PASTIX_ROOT} )
      set(PASTIX_ROOT $ENV{PASTIX_ROOT})
      set(MURGE_ROOT $ENV{PASTIX_ROOT})
   ENDIF()
   find_package(BLAS QUIET)
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

add_definitions(-DGFORTRAN -DDEBUG -DMPIF90)

get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
IF(Fortran_COMPILER STREQUAL "GFORTRAN")
	set(CMAKE_Fortran_FLAGS "-g -Wall -cpp -pedantic -ffree-line-length-none -std=f2003 -fall-intrinsics")
ELSEIF(Fortran_COMPILER STREQUAL "INTEL")
	set(CMAKE_Fortran_FLAGS " ")
ELSEIF(Fortran_COMPILER_NAME STREQUAL "xlf")
	set(CMAKE_Fortran_FLAGS "-qextname=flush -qthreaded -qhalt=e")
ELSE()
	message(STATUS "NO FORTRAN COMPILER FOUND")
ENDIF()

STRING(FIND ${CMAKE_Fortran_FLAGS} "-std=f95" VAR)
IF("${VAR}" STREQUAL "-1")
  set(STDF95 NO)
ELSE()
  set(STDF95 YES)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fmax-identifier-length=63")
  add_definitions(-DSTDF95)
ENDIF()

#add flag to disabled MPI stuff for debug
set(DISABLED_MPI_MODULE OFF CACHE BOOL " ")

# recherche de MPI
# I add HINTS $ENV{MPI_ROOT} to the find_program(MPIEXEC) in the FindMPI.cmake file (src/CMakeModules)
# If specify MPI_ROOT doesn't work try to change your PATH variable
find_package(MPI)
IF(MPI_FOUND)
	message(STATUS "MPI FOUND")
	include_directories(${MPI_Fortran_INCLUDE_PATH})
ELSE(MPI_FOUND)
	message(STATUS "MPI NOT FOUND")
	#set(DISABLED_MPI_MODULE ON CACHE BOOL " " FORCE)
ENDIF(MPI_FOUND)

# include's
include_directories(fft remap poisson_solvers utilities memory assert precision fields) 

# pour la doc 
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

#############
#  RESUME
#############
#
#message("DISABLED_MPI_MODULE = ${DISABLED_MPI_MODULE}")
#
#############

#low-level
set(LOW_LEVEL_LIST_DIR memory precision assert constants utilities low_level_file_io timer)
#set(LOW_LEVEL_MODULES sll_memory sll_working_precision sll_assertion sll_constants sll_utilities sll_low_level_file_io sll_timer)
add_subdirectory(memory)
add_subdirectory(precision)
add_subdirectory(assert)
add_subdirectory(constants)
add_subdirectory(utilities)
add_subdirectory(low_level_file_io)
add_subdirectory(timer)

#mid-level
set(MID_LEVEL_LIST_DIR numerical_utilities non_linear_solvers splines integration fftpack fft collective remap WENO poisson_solvers)
add_subdirectory(numerical_utilities)
add_subdirectory(non_linear_solvers)
add_subdirectory(splines)
add_subdirectory(WENO)
add_subdirectory(integration)
IF(FFTPACK_ENABLED)
	add_subdirectory(fftpack)
ELSE()
	LIST(REMOVE_ITEM MID_LEVEL_LIST_DIR fftpack)
ENDIF()
add_subdirectory(fft)
IF(DISABLED_MPI_MODULE)
# nothing to do
ELSE()
  add_subdirectory(collective)
  add_subdirectory(remap)
ENDIF()
IF(NOT STDF95)
add_subdirectory(poisson_solvers)
ENDIF()

#high-level
add_library(geometry_functions STATIC ${CMAKE_SOURCE_DIR}/fields/geometry_functions.F90)
add_dependencies(geometry_functions sll_utilities sll_assertion sll_constants sll_working_precision)
add_subdirectory(interpolators)
add_subdirectory(mapped_meshes)
add_subdirectory(fields)
add_subdirectory(ode_solvers)
add_subdirectory(distribution_function)
add_subdirectory(advection_field)
add_subdirectory(simulation)
IF(FISHPACK_ENABLED)
	add_subdirectory(fishpack)
ENDIF()



IF(MUDPACK_ENABLED)
	add_subdirectory(mudpack)
	add_subdirectory(multigrid)
ENDIF()
add_subdirectory(LinRood)


#add_subdirectory(maxwell_solvers) #use old mesh_types 
#add_subdirectory(BSL) top
add_subdirectory(CSL)
#add_subdirectory(quasi_neutral_solvers) don't compile problem with fftw3_thread! top



#add_subdirectory(selalib_level EXCLUDE_FROM_ALL)
#add_custom_target(install make selalib_low_level)



######################################################################
#                              TESTS                                 #
######################################################################

enable_testing()

add_test(NAME memory COMMAND test_memory)
set_tests_properties(memory PROPERTIES TIMEOUT 20)

#add_test(precision precision/test_precision) #too longer

set(passRegex "Assertion error triggered in file")
add_test(NAME assert COMMAND test_assert)
set_tests_properties(assert PROPERTIES PASS_REGULAR_EXPRESSION "${passRegex}")

add_test(NAME constants COMMAND test_constants)
add_test(NAME utilities COMMAND test_utilities)
add_test(NAME low_level_file_io COMMAND test_io)
set_tests_properties(low_level_file_io PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
add_test(NAME timer COMMAND test_timer)
add_test(NAME tridiagonal COMMAND test_tridiagonal)
add_test(NAME newton_raphson COMMAND test_newton_raphson)
add_test(NAME splines COMMAND test_splines #test_non_unif_splines 
  test_arbitrary_degree_splines)
set_tests_properties(splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
add_test(NAME integration COMMAND test_integration)
add_test(NAME fft COMMAND test_fft)
IF(FISHPACK_ENABLED)
	add_test(fishpack fishpack/thwscyl fishpack/thwsplr fishpack/thwsssp fishpack/tpois3d)
ENDIF()
IF(MUDPACK_ENABLED)
	add_test(NAME mudpack COMMAND tmud34sp tmud24sp )
        set_tests_properties(mudpack PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
        set(PROCS 4)
        set(ARGS "")
	add_test(NAME multigrid COMMAND 
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_multigrid
	${MPIEXEC_POSTFLAGS} ${ARGS})
        set_tests_properties(multigrid PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
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

add_test(NAME WENO COMMAND test_WENO_interp test_WENO_recon)
add_test(NAME interpolators COMMAND test_interpolators_1d
				    test_interpolators_2d)
add_test(NAME mapped_meshes COMMAND test_mapped_meshes_1d
				    test_mapped_meshes_2d)
add_test(NAME fields COMMAND test_scalar_field)
add_test(NAME ode_solvers COMMAND test_implicit_ode_nonuniform)
add_test(NAME distribution_function COMMAND test_distribution_function)
add_test(NAME advection_field COMMAND test_advection_field)

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


set(FORTRANCL_ENABLED CACHE BOOL OFF)
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

project (sll)
enable_language (Fortran)
cmake_minimum_required(VERSION 2.8)
#INCLUDE(CMakeForceCompiler)

# Determine how-to install the modules
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
# Add the modules directory to the list of include directories
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

set(CMAKE_BINARY_DIR "${CMAKE_BINARY_DIR}/bin")

# Add new directory to find the FIND<PACKAGE>.cmake files
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules ${CMAKE_MODULE_PATH} )

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

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
set(HDF5_ENABLED YES CACHE BOOL " ")

# recherche de HDF5
IF(HDF5_ENABLED)
	IF( DEFINED ENV{HDF5_ROOT} )
		set(HDF5_ROOT $ENV{HDF5_ROOT})
	ENDIF()
	find_package(HDF5 QUIET)
	if (HDF5_FOUND)
		message(STATUS "HDF5 FOUND")
		include_directories(${HDF5_INCLUDE_DIRS})
	else (HDF5_FOUND)
		add_definitions(-DNOHDF5)
		set(HDF5_ENABLED NO)
		set(HDF5_LIBRARIES "")
	endif(HDF5_FOUND)
ELSE(HDF5_ENABLED)
	add_definitions(-DNOHDF5)
ENDIF(HDF5_ENABLED)
##########################################################

##########################################################
# add the cache entry FFTPACK_ENABLED for enable/disable fftpack
set(FFTPACK_ENABLED NO CACHE BOOL " ")

# recherche le package FFTW
set(FFTW_ENABLED YES CACHE BOOL " ")
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
		set(FFTW_ENABLED NO)
	ENDIF(FFTW_FOUND)
ENDIF(FFTW_ENABLED)
##########################################################

##########################################################
# add the cache entry PASTIX_ENABLED for enable/disable pastix
set(PASTIX_ENABLED NO CACHE BOOL " ")

# recherche le package PASTIX
set(PASTIX ENABLED YES CACHE BOOL " ")
IF(PASTIX_ENABLED)
	IF( DEFINED ENV{PASTIX_ROOT} )
		set(PASTIX_ROOT $ENV{PASTIX_ROOT})
	ENDIF()
	find_package(PASTIX QUIET)
	IF(PASTIX_FOUND)
		message(STATUS "PASTIX FOUND")
	ELSE(PASTIX_FOUND)
		message(STATUS "PASTIX NOT FOUND")
		set(PASTIX_ENABLED NO)
	ENDIF(PASTIX_FOUND)
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

#add flag to disabled MPI stuff for debug
set(DISABLED_MPI_MODULE NO CACHE BOOL " ")


# recherche de MPI
IF( DEFINED ENV{MPI_ROOT} )
	set(MPI_ROOT $ENV{MPI_ROOT})
	set(_MPI_BASE_DIR ${MPI_ROOT})
ENDIF()
find_package(MPI)
IF(MPI_FOUND)
	message(STATUS "MPI FOUND")
	include_directories(${MPI_Fortran_INCLUDE_PATH})
	#CMAKE_FORCE_Fortran_COMPILER(mpif90 ${Fortran_COMPILER_NAME})
	#CMAKE_FORCE_Fortran_COMPILER(${MPI_Fortran_COMPILER} ${Fortran_COMPILER_NAME})
ELSE(MPI_FOUND)
	message(STATUS "MPI NOT FOUND")
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

set(FISHPACK_ENABLED CACHE BOOL OFF)
set(MUDPACK_ENABLED CACHE BOOL OFF)


#low-level
set(LOW_LEVEL_LIST_DIR memory precision assert constants utilities low_level_file_io timer)
add_subdirectory(memory)
add_subdirectory(precision)
add_subdirectory(assert)
add_subdirectory(constants)
add_subdirectory(utilities)
add_subdirectory(low_level_file_io)
add_subdirectory(timer)

#mid-level
set(MID_LEVEL_LIST_DIR numerical_utilities non_linear_solvers splines integration fftpack fft collective remap WENO )#poisson_solvers)
add_subdirectory(numerical_utilities)
add_subdirectory(non_linear_solvers)
add_subdirectory(splines)
add_subdirectory(integration)
IF(FFTPACK_ENABLED)
	add_subdirectory(fftpack)
ELSE()
	LIST(REMOVE_ITEM MID_LEVEL_LIST_DIR fftpack)
ENDIF()
add_subdirectory(fft)
IF(NOT DISABLED_MPI_MODULE)
  add_subdirectory(collective)
  add_subdirectory(remap)
ENDIF()
add_subdirectory(WENO)
add_subdirectory(poisson_solvers)

#hight-level
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
	add_subdirectory(mgd3)
ENDIF()


#add_subdirectory(maxwell_solvers) #use old mesh_types 
#add_subdirectory(BSL) top
#add_subdirectory(CSL) don't compile error in SConstruct! top
#add_subdirectory(quasi_neutral_solvers) don't compile problem with fftw3_thread! top



add_subdirectory(selalib_level EXCLUDE_FROM_ALL)
#add_custom_target(install make selalib_low_level)



######################################################################
#                              TESTS                                 #
######################################################################

enable_testing()

IF(NOT Fortran_COMPILER_NAME STREQUAL "ifort")
add_test(memory memory/test_memory)
ENDIF()

#add_test(precision precision/test_precision) #too longer

set(passRegex "Assertion error triggered in file")
add_test(assert assert/test_assert)
set_tests_properties(assert PROPERTIES PASS_REGULAR_EXPRESSION "${passRegex}")

add_test(constants constants/test_constants)
add_test(utilities utilities/test_utilities)
add_test(low_level_file_io low_level_file_io/test_io)
add_test(timer timer/test_timer)
add_test(tridiagonal numerical_utilities/test_tridiagonal)
add_test(newton_raphson non_linear_solvers/test_newton_raphson)
add_test(splines splines/test_splines)
add_test(splines splines/test_non_unif_splines)
add_test(integration integration/test_integration)
add_test(fft fft/test_fft)
IF(FISHPACK_ENABLED)
	add_test(fishpack fishpack/thwscyl fishpack/thwsplr fishpack/thwsssp fishpack/tpois3d)
ENDIF()
IF(MUDPACK_ENABLED)
	add_test(mudpack mudpack/tmud34sp mudpack/tmud24sp )
	add_test(mgd3 mgd3/tmgd3)
ENDIF()

# we can call ${MPIEXEC_MAX_NUMPROCS}
set(PROCS 2)
set(ARGS "")
add_test(collective ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS} collective/test_collective ${MPIEXEC_POSTFLAGS} ${ARGS})
#set_tests_properties(collective PROPERTIES FAIL_REGULAR_EXPRESSION "${failRegex}")

set(PROCS 8)
set(ARGS "")
add_test(remap ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS} remap/test_remap ${MPIEXEC_POSTFLAGS} ${ARGS})

add_test(WENO WENO/test_WENO)
add_test(interpolators interpolators/test_interpolators_1d)
add_test(interpolators interpolators/test_interpolators_2d)
add_test(mapped_meshes mapped_meshes/test_mapped_meshes)
add_test(mapped_meshes_1d mapped_meshes/test_mapped_meshes_1d)
add_test(fields fields/test_scalar_field)
add_test(ode_solvers ode_solvers/test_implicit_ode_nonuniform)
add_test(distribution_function distribution_function/test_distribution_function)
add_test(advection_field advection_field/test_advection_field)

#add_test(maxwell_solvers maxwell_solvers/test_maxwell_2d)
   #add_test(quasi_neutral_solvers quasi_neutral_solvers/test_quasi_neutral_solver)
IF(FFTPACK_ENABLED)
	add_test(poisson_solvers poisson_solvers/test_poisson_1d)
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

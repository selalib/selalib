project (sll)
enable_language (Fortran)
cmake_minimum_required(VERSION 2.8)
INCLUDE(CMakeForceCompiler)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

add_definitions(-DGFORTRAN -DDEBUG -DMPIF90)

# use MPIF90 as compiler
# set(CMAKE_Fortran_COMPILER mpif90)

get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if(Fortran_COMPILER_NAME STREQUAL "gfortran")
	set(CMAKE_Fortran_FLAGS "-g -Wall -cpp -pedantic -ffree-line-length-none -std=f2003 -fall-intrinsics")
	message(STATUS "COMPILER GFORTRAN")
else()
	message(STATUS "NO FORTRAN COMPILER FOUND")
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

#ajoute les modules locaux
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules ${CMAKE_MODULE_PATH} )	
message(STATUS "CMAKE_SOURCE_DIR = " ${CMAKE_SOURCE_DIR})

#set(PATH_WITH_SPACE $ENV{PATH})
#string(REPLACE ":" ";" PATH_WITH_SPACE ${PATH_WITH_SPACE})
#message(STATUS "PATH_WITH_SPACE = " ${PATH_WITH_SPACE})
#set( CMAKE_PREFIX_PATH ${PATH_WITH_SPACE})
#set( CMAKE_PREFIX_PATH /opt/local )

#set(CMAKE_LIBRARY_PATH "/opt/local;${CMAKE_LIBRARY_PATH}")

# recherche de MPI
IF( DEFINED ENV{MPI_ROOT} )
  set(MPI_ROOT $ENV{MPI_ROOT})
  set(_MPI_BASE_DIR ${MPI_ROOT})
ENDIF()
find_package(MPI)
if (MPI_FOUND)
  message(STATUS "MPI FOUND")
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  message(STATUS "PARALLEL FORTRAN IS ${MPI_Fortran_COMPILER}")
  #CMAKE_FORCE_Fortran_COMPILER(mpif90 ${Fortran_COMPILER_NAME})
  CMAKE_FORCE_Fortran_COMPILER(${MPI_Fortran_COMPILER} ${Fortran_COMPILER_NAME})
  message(STATUS "Compiler : ${CMAKE_Fortran_COMPILER}")
  #set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
else (MPI_FOUND)
  message(STATUS "MPI NOT FOUND")
endif(MPI_FOUND)

# recherche de HDF5
IF( DEFINED ENV{HDF5_ROOT} )
  set(HDF5_ROOT $ENV{HDF5_ROOT})
ENDIF()
find_package(HDF5 QUIET)
if (HDF5_FOUND)
    message(STATUS "HDF5 FOUND")
    message(STATUS "libraries = " ${HDF5_LIBRARIES}) 
    message(STATUS "include = " ${HDF5_INCLUDE_DIRS}) 
    include_directories(${HDF5_INCLUDE_DIRS})
else (HDF5_FOUND)
    add_definitions(-DNOHDF5)
    set(HDF5_LIBRARIES "")
endif(HDF5_FOUND)

# recherche le package FFTW
IF( DEFINED ENV{FFTW_ROOT} )
  set(FFTW_ROOT $ENV{FFTW_ROOT})
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
include_directories(remap poisson_solvers utilities memory assert precision fields) 


# pour la doc 
find_package(Doxygen)
if(DOXYGEN_FOUND)
   configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in 
   ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
   add_custom_target(doc ${DOXYGEN_EXECUTABLE} 
   ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY 
   ${CMAKE_CURRENT_BINARY_DIR} 
   COMMENT "Generating API documentation with Doxygen" VERBATIM)
else(DOXYGEN_FOUND)
   MESSAGE(STATUS "DOXYGEN NOT FOUND")
endif(DOXYGEN_FOUND)



add_subdirectory(memory)
add_subdirectory(precision)
add_subdirectory(assert)
add_subdirectory(constants)
add_subdirectory(utilities)
add_subdirectory(low_level_file_io)
add_subdirectory(timer)
add_subdirectory(numerical_utilities)
add_subdirectory(non_linear_solvers)
add_subdirectory(splines)
add_subdirectory(integration)
add_subdirectory(fftpack)
add_subdirectory(fft)
add_subdirectory(collective)
add_subdirectory(remap)
add_subdirectory(interpolators)
add_subdirectory(WENO)
add_subdirectory(mapped_meshes)
add_subdirectory(fields)
add_subdirectory(ode_solvers)
add_subdirectory(poisson_solvers)
add_subdirectory(distribution_function)
add_subdirectory(advection_field)
add_subdirectory(simulation)
add_subdirectory(fishpack)


#add_subdirectory(maxwell_solvers)
#add_subdirectory(poisson_solvers)
#add_subdirectory(BSL)

#add_subdirectory(CSL) don't compile error in SConstruct!
#add_subdirectory(quasi_neutral_solvers) don't compile problem with fftw3_thread!

######################################################################
#                              TESTS                                 #
######################################################################

enable_testing()

add_test(memory memory/test_memory)

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
add_test(fishpack fishpack/thwscyl fishpack/thwsplr fishpack/thwsssp fishpack/tpois3d)


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
add_test(poisson_solvers poisson_solvers/test_poisson_1d)
#add_test(poisson_solvers poisson_solvers/test_poisson_2d)
#add_test(poisson_solvers poisson_solvers/test_poisson_3d_seq)
#add_test(poisson_solvers mpirun -np 2 poisson_solvers/test_poisson_3d_par)
   #add_test(poisson_solvers poisson_solvers/unit_test)
#add_test(BSL BSL/test_bsl)
#add_test(BSL BSL/test_landau_4d)
#add_test(BSL BSL/test_landau_2d)
   #add_test(CSL CSL/test_csl)
add_test(simulation_1d simulation/test_simulation_1d)

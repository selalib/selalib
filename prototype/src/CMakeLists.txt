project (sll)
enable_language (Fortran)
cmake_minimum_required(VERSION 2.8)

add_definitions(-DGFORTRAN -DDEBUG -DMPIF90)

# use MPIF90 as compiler
set(CMAKE_Fortran_COMPILER mpif90)
# some useful flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -Wall -cpp -ffree-line-length-none -std=f2003 -fall-intrinsics")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

#ajoute les modules locaux
set (CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR}/CMakeModules ${CMAKE_MODULE_PATH} )	
message(STATUS "CMAKE_SOURCE_DIR = " ${CMAKE_SOURCE_DIR})

set(PATH_WITH_SPACE $ENV{PATH})
string(REPLACE ":" ";" PATH_WITH_SPACE ${PATH_WITH_SPACE})
#message(STATUS "PATH_WITH_SPACE = " ${PATH_WITH_SPACE})
set( CMAKE_PREFIX_PATH ${PATH_WITH_SPACE})
#set( CMAKE_PREFIX_PATH /opt/local )

# recherche de HDF5
find_package(HDF5 COMPONENTS Fortran)

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
find_package(FFTW)

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

# FFTW_INCLUDE_DIR = fftw3.f03
# FFTW_LIBRARIES = libfftw3.a
# FFTW_FOUND = true if FFTW3 is found

#set(FFTW_FOUND TRUE)
#IF( DEFINED ENV{FFTWINCPATH} )
#  set(FFTW_INCLUDE_DIR $ENV{FFTWINCPATH})
#ELSE( DEFINED ENV{FFTWINCPATH} )
#  set(FFTW_INCLUDE_DIR "")
#  set(FFTW_FOUND FALSE)
#ENDIF( DEFINED ENV{FFTWINCPATH} )

#IF( FFTW_FOUND AND DEFINED ENV{FFTWLIBPATH} )
#  FIND_LIBRARY(FFTW_LIBRARIES fftw3 $ENV{FFTWLIBPATH})
#ELSE()
#  set(FFTW_LIBRARIES "")
#  set(FFTW_FOUND FALSE)
#ENDIF()

#IF(FFTW_FOUND)
#  MESSAGE(STATUS "FFTW FOUND")
#  MESSAGE(STATUS "FFTW_INCLUDE_DIR=${FFTW_INCLUDE_DIR}")
#  MESSAGE(STATUS "FFTW_LIBRARIES=${FFTW_LIBRARIES}")
#  set(FFT_ADD_MODULE fftpack_module fftw_module)
#  include_directories(${FFTW_INCLUDE_DIR})
#ELSE()
#  MESSAGE(STATUS "FFTW NOT FOUND")
#  set(FFT_ADD_MODULE fftpack_module)
#  add_definitions(-D_NOFFTW)
#ENDIF()

# include's
include_directories(remap poisson_solvers utilities memory mesh_types assert precision fft fftw) 

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

#add_subdirectory(interpolators)
#add_subdirectory(WENO)
#add_subdirectory(mapped_meshes)
#add_subdirectory(fields)
#add_subdirectory(ode_solvers)
#add_subdirectory(distribution_function)
#add_subdirectory(maxwell_solvers)
#add_subdirectory(poisson_solvers)
#add_subdirectory(BSL)

#add_subdirectory(CSL) don't compile error in SConstruct!
#add_subdirectory(quasi_neutral_solvers) don't compile problem with fftw3_thread!


# add tests
enable_testing()
add_test(memory memory/test_memory)
add_test(precision precision/test_precision)
add_test(assert assert/test_assert)
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
add_test(collective mpirun -np 4 collective/test_collective)
add_test(remap mpirun -np 8 remap/test_remap)

#add_test(interpolators interpolators/test_interpolators)
#add_test(WENO WENO/test_WENO)
#add_test(mapped_meshes mapped_meshes/test_mapped_meshes)
#add_test(fields fields/test_mesh_types)
#add_test(distribution_function distribution_function/test_distribution_function)
#add_test(maxwell_solvers maxwell_solvers/test_maxwell_2d)
#add_test(ode_solvers ode_solvers/test_ode_solvers)
   #add_test(quasi_neutral_solvers quasi_neutral_solvers/test_quasi_neutral_solver)
#add_test(poisson_solvers poisson_solvers/test_poisson_1d)
#add_test(poisson_solvers poisson_solvers/test_poisson_2d)
#add_test(poisson_solvers poisson_solvers/test_poisson_3d_seq)
#add_test(poisson_solvers mpirun -np 2 poisson_solvers/test_poisson_3d_par)
   #add_test(poisson_solvers poisson_solvers/unit_test)
#add_test(BSL BSL/test_bsl)
#add_test(BSL BSL/test_landau_4d)
#add_test(BSL BSL/test_landau_2d)
   #add_test(CSL CSL/test_csl)

project (sll)
enable_language (Fortran)
cmake_minimum_required(VERSION 2.8)

add_definitions(-DGFORTRAN -DDEBUG -DMPIF90)

# use MPIF90 as compiler
set(CMAKE_Fortran_COMPILER mpif90)
# some useful flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -Wall -cpp -ffree-line-length-none -std=f2003 -fall-intrinsics")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# recherche de HDF5
find_package(HDF5 COMPONENTS Fortran Fortran_HL)

if (HDF5_FOUND)
    message(STATUS "HDF5 trouve")
    message(STATUS "libraries = " ${HDF5_LIBRARIES}) 
    message(STATUS "include = " ${HDF5_INCLUDE_DIRS}) 
    include_directories(${HDF5_INCLUDE_DIRS})
else (HDF5_FOUND)
    add_definitions(-DNOHDF5)
endif(HDF5_FOUND)

#ajoute les modules locaux
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMakeModules)	

# recherche le package GMP
find_package(FFTW)
if (FFTW_FOUND)
	include_directories(${FFTW_INCLUDE_DIR})
endif()


# include's
include_directories(remap poisson_solvers utilities memory mesh_types assert precision) 

add_subdirectory(precision)
add_subdirectory(memory)
add_subdirectory(collective)
add_subdirectory(utilities)
add_subdirectory(remap)
add_subdirectory(constants)
add_subdirectory(numerical_utilities)
add_subdirectory(diagnostics)
add_subdirectory(splines)
add_subdirectory(mesh_types)
add_subdirectory(timer)
add_subdirectory(fft)
add_subdirectory(assert)

# add tests
enable_testing()
add_test(timer timer/test_timer)
add_test(memory memory/test_memory)
add_test(collective collective/test_collective)
add_test(utilities utilities/test_utilities)
add_test(remap mpirun -np 8 remap/test_remap)
add_test(constants constants/test_constants)
add_test(tridiagonal numerical_utilities/test_tridiagonal)
add_test(diagnostics diagnostics/test_diagnostics)
add_test(splines-1 splines/test_splines)
add_test(splines-2 splines/test_non_unif_splines)
add_test(mesh_types mesh_types/test_mesh_types)
add_test(fft fft/test_fft)
add_test(assert assert/test_assert)

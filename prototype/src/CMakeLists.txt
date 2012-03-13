project (sll)
enable_language (Fortran)
cmake_minimum_required(VERSION 2.8)

add_definitions(-DGFORTRAN -DDEBUG -DMPIF90 -DNOHDF5)

# use MPIF90 as compiler
set(CMAKE_Fortran_COMPILER mpif90)
# some useful flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O5 -Wall -cpp -ffree-line-length-none -std=f2003 -fall-intrinsics")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# include's
include_directories(remap poisson_solvers utilities memory mesh_types assert precision) 

add_subdirectory(memory)
add_subdirectory(precision)
add_subdirectory(assert)
add_subdirectory(collective)
add_subdirectory(utilities)
add_subdirectory(remap)
add_subdirectory(constants)
add_subdirectory(numerical_utilities)

# add tests
enable_testing()
add_test(memory memory/test_memory)
add_test(collective collective/test_collective)
add_test(utilities utilities/test_utilities)
add_test(remap mpirun -np 2 remap/test_remap)
add_test(constants constants/test_constants)
add_test(tridiagonal numerical_utilities/test_tridiagonal)
#add_test(assert assert/test_assert)

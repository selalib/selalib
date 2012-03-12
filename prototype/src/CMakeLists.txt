project (sll)
enable_language (Fortran)
cmake_minimum_required(VERSION 2.8)

add_definitions(-DGFORTRAN -DDEBUG -DMPIF90 -DNOHDF5)

# FLAGS
message (STATUS ${CMAKE_Fortran_COMPILER})
set(CMAKE_Fortran_COMPILER mpif90)
#set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS .h .F90)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -cpp -ffree-line-length-none -std=f2003")

# pour les includes
include_directories(remap poisson_solvers utilities memory mesh_types assert precision) 

# pour les sources de test
add_subdirectory(memory)

# pour les tests
enable_testing()
add_test(memory memory/unit_test)

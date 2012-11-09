
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

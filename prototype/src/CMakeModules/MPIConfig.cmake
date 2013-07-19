SET(MPI_MODULE_ENABLED ON CACHE BOOL " ")

IF($ENV{HOSTNAME} MATCHES "hydra")
   SET(MPI_Fortran_COMPILER "mpiifort")
ENDIF()

IF(MPI_MODULE_ENABLED)
   FIND_PACKAGE(MPI REQUIRED Fortran)
ENDIF(MPI_MODULE_ENABLED)

IF(MPI_FOUND)
   MESSAGE(STATUS "MPI FOUND")
   INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
ELSE(MPI_FOUND)
   MESSAGE(STATUS "MPI NOT FOUND")
   SET(MPI_MODULE_ENABLED OFF CACHE BOOL " " FORCE)
ENDIF(MPI_FOUND)


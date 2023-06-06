if($ENV{HOSTNAME} MATCHES "hydra*")
  set(MPI_Fortran_COMPILER "mpiifort")
endif()

if(MPI_ENABLED)
  find_package(MPI REQUIRED Fortran CXX)
endif(MPI_ENABLED)

if(MPI_FOUND)
  message(STATUS "MPI FOUND")
  find_path(
    MPI_Fortran_MOD_DIR
    NAMES mpi.mod
    PATHS $ENV{MPI_FORTRAN_MOD_DIR} ${MPI_Fortran_INCLUDE_PATH})
  if(MPI_Fortran_MOD_DIR)
    set(MPI_Fortran_INCLUDE_PATH ${MPI_Fortran_MOD_DIR}
                                 ${MPI_Fortran_INCLUDE_PATH})
  endif(MPI_Fortran_MOD_DIR)
else(MPI_FOUND)
  message(STATUS "MPI NOT FOUND")
  set(MPI_ENABLED
      OFF
      CACHE BOOL " " FORCE)
endif(MPI_FOUND)

mark_as_advanced(MPI_EXTRA_LIBRARY MPI_LIBRARY MPI_Fortran_MOD_DIR)

mark_as_advanced(CLEAR MPI_Fortran_COMPILER)

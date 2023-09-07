if(MPI_ENABLED)

  find_package(MPI REQUIRED Fortran CXX)

  if(MPI_FOUND)
    message(STATUS "MPI FOUND")
    find_path(
      MPI_Fortran_MOD_DIR
      NAMES mpi.mod
      PATHS $ENV{MPI_Fortran_MOD_DIR})
    if(MPI_Fortran_MOD_DIR)
      set(MPI_Fortran_ADDITIONAL_INCLUDE_DIRS ${MPI_Fortran_MOD_DIR})
    endif(MPI_Fortran_MOD_DIR)
  else(MPI_FOUND)
    message(STATUS "MPI NOT FOUND")
    set(MPI_ENABLED
        OFF
        CACHE BOOL " " FORCE)
  endif(MPI_FOUND)

endif(MPI_ENABLED)

ADD_LIBRARY(dagmg_seq STATIC dagmg.f90 dagmg_mumps.f90)

IF(MPI_MODULE_ENABLED)
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  ADD_LIBRARY(dagmg_par STATIC dagmg.f90 dagmg_par.f90)
ENDIF(MPI_MODULE_ENABLED)

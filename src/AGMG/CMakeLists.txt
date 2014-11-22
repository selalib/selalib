IF(Fortran_COMPILER_NAME MATCHES "gfortran")
   set(CMAKE_Fortran_FLAGS_DEBUG "-g -w ")
   set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -w ")
ELSEIF(Fortran_COMPILER_NAME MATCHES "ifort")
   SET(CMAKE_Fortran_FLAGS_RELEASE "-r8 -nowarn -O3 -xHost -ip -openmp")
   SET(CMAKE_Fortran_FLAGS_DEBUG "-g -r8 -nowarn")
ENDIF()

ADD_LIBRARY(blas_agmg STATIC blas_agmg.f)
ADD_LIBRARY(lapack_agmg STATIC lapack_agmg.f)
ADD_LIBRARY(mumps_agmg STATIC dagmg_mumps.f90)

ADD_LIBRARY(dagmg_seq STATIC dagmg.f90 )
TARGET_LINK_LIBRARIES(dagmg_seq mumps_agmg lapack_agmg blas_agmg) 

ADD_EXECUTABLE(test_dagmg_seq Example_seq.F90)
TARGET_LINK_LIBRARIES(test_dagmg_seq dagmg_seq sll_working_precision)

IF(MPI_MODULE_ENABLED)

  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

  FIND_LIBRARY(MUMPS_COMMON_LIBRARY
    NAMES
    mumps_common
    PATHS
    /opt/local
    $ENV{PETSC_DIR}
    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}
    $ENV{MUMPS_DIR}
    SUFFIXES lib
  )

  FIND_LIBRARY(DMUMPS_LIBRARY
    NAMES
    dmumps
    PATHS
    /opt/local
    $ENV{PETSC_DIR}
    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}
    $ENV{MUMPS_DIR}
    SUFFIXES lib
  )

  IF (MUMPS_COMMON_LIBRARY AND DMUMPS_LIBRARY)
    SET(MUMPS_LIBRARIES ${DMUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY})
  ENDIF()

  INCLUDE(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(MUMPS "Could not find MUMPS " MUMPS_LIBRARIES)

  IF(MUMPS_FOUND)
    ADD_LIBRARY(dagmg_par STATIC dagmg_par.f90)
    ADD_EXECUTABLE(test_dagmg_par Example_par.f90)
    TARGET_LINK_LIBRARIES(test_dagmg_par dagmg_par ${MUMPS_LIBRARIES} lapack_agmg blas_agmg)
    MARK_AS_ADVANCED(MUMPS_LIBRARIES)
  ENDIF(MUMPS_FOUND)

ENDIF()

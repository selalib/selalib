IF(CMAKE_Fortran_COMPILER MATCHES "gfortran")

   ADD_DEFINITIONS(-DGFORTRAN)
   SET(CMAKE_Fortran_FLAGS_RELEASE "-w -ffree-line-length-none -fall-intrinsics -O3 -fopenmp" CACHE STRING "opt flags" FORCE)
   SET(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -ffree-line-length-none -fall-intrinsics -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow -O0" CACHE STRING "debug flags" FORCE)

ELSEIF(CMAKE_Fortran_COMPILER MATCHES "ifort")

      SET(CMAKE_Fortran_FLAGS_RELEASE "-nowarn -O3 -xHost -ip -openmp" CACHE STRING "opt flags" FORCE)
      SET(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check all,noarg_temp_created -fpe0 -traceback -ftrapuv " CACHE STRING "debug flags" FORCE)

ELSEIF(CMAKE_Fortran_COMPILER MATCHES "xlf")

   ADD_DEFINITIONS(-DIBM)

   SET(CMAKE_Fortran_FLAGS_DEBUG "-qextname=flush -qthreaded -qhalt=e -qxlf2003=polymorphic" CACHE STRING "opt flags" FORCE)
   SET(CMAKE_Fortran_FLAGS_RELEASE "-qextname=flush -qthreaded -qhalt=e -qxlf2003=polymorphic" CACHE STRING "debug flags" FORCE)

ELSE()

   MESSAGE(SEND_ERROR "NO KNOWN FORTRAN COMPILER FOUND")

ENDIF()

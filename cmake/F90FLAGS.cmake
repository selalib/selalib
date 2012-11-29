IF(CMAKE_Fortran_COMPILER MATCHES "gfortran")

   ADD_DEFINITIONS(-DGFORTRAN)

   SET(CMAKE_Fortran_FLAGS_RELEASE "-fomit-frame-pointer -Wall -cpp -ffree-line-length-none -fall-intrinsics -O3")
   SET(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -cpp -ffree-line-length-none -fall-intrinsics -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow -O0")

ELSEIF(CMAKE_Fortran_COMPILER MATCHES "ifort")

   ADD_DEFINITIONS(-DINTEL)

      SET(CMAKE_Fortran_FLAGS_DEBUG "-O3 -xHost -ip")
      SET(CMAKE_Fortran_FLAGS_RELEASE "-O0 -check all -fpe0 -traceback -ftrapuv ")

ELSEIF(CMAKE_Fortran_COMPILER MATCHES "xlf")

   ADD_DEFINITIONS(-DIBM)

   SET(CMAKE_Fortran_FLAGS_DEBUG "-qextname=flush -qthreaded -qhalt=e -qxlf2003=polymorphic")
   SET(CMAKE_Fortran_FLAGS_RELEASE "-qextname=flush -qthreaded -qhalt=e -qxlf2003=polymorphic")

ELSE()

   MESSAGE(STATUS "NO KNOWN FORTRAN COMPILER FOUND")

ENDIF()

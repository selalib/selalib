IF(CMAKE_Fortran_COMPILER MATCHES "gfortran")

   ADD_DEFINITIONS(-DGFORTRAN)

   IF(CMAKE_BUILD_TYPE MATCHES Release)
      SET(CMAKE_Fortran_FLAGS "-fomit-frame-pointer -Wall -cpp -ffree-line-length-none -fall-intrinsics -O3")
   ELSE()
      SET(CMAKE_Fortran_FLAGS "-g -Wall -cpp -ffree-line-length-none -fall-intrinsics -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow -O0")
   ENDIF()

ELSEIF(CMAKE_Fortran_COMPILER MATCHES "ifort")

   ADD_DEFINITIONS(-DINTEL)

   IF(CMAKE_BUILD_TYPE MATCHES Release)
      SET(CMAKE_Fortran_FLAGS "-O3 -xHost -ip")
   ELSE()
      SET(CMAKE_Fortran_FLAGS "-O0 -check all -fpe0 -traceback -ftrapuv ")
   ENDIF()

ELSEIF(CMAKE_Fortran_COMPILER MATCHES "xlf")

   ADD_DEFINITIONS(-DIBM)

   SET(CMAKE_Fortran_FLAGS "-qextname=flush -qthreaded -qhalt=e -qxlf2003=polymorphic")

ELSE()

   MESSAGE(STATUS "NO KNOWN FORTRAN COMPILER FOUND")

ENDIF()

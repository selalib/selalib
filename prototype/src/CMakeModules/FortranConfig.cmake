# Determine how-to install the modules. CMAKE_BINARY_DIR is the directory
# in which the make command is invoked.
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# Add the modules directory to the list of include directories
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

GET_FILENAME_COMPONENT(Fortran_COMPILER_NAME "${CMAKE_Fortran_COMPILER}" NAME)

IF (CMAKE_Fortran_COMPILER MATCHES "ifort")
   EXEC_PROGRAM(${CMAKE_Fortran_COMPILER} ARGS "-v" OUTPUT_VARIABLE source_path)
   STRING(REGEX MATCH "1[0-9]\\.[0-9]\\.[0-9]" Fortran_COMPILER_VERSION ${source_path})
ELSE()
   EXEC_PROGRAM(${CMAKE_Fortran_COMPILER} ARGS "--version" OUTPUT_VARIABLE source_path)
   STRING(REGEX MATCH "4\\.[0-9]\\.[0-9]" Fortran_COMPILER_VERSION ${source_path})
ENDIF()
MESSAGE(STATUS "Fortran compiler : ${Fortran_COMPILER_NAME}-${Fortran_COMPILER_VERSION}")


IF(CMAKE_Fortran_COMPILER MATCHES "gfortran*")
   ADD_DEFINITIONS(-DGFORTRAN)
   SET(CMAKE_Fortran_FLAGS_RELEASE "-w -ffree-line-length-none -fall-intrinsics -O3 -fopenmp")
   SET(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -cpp -pedantic -ffree-line-length-none -std=f2003 -fall-intrinsics -fbounds-check -fbacktrace -ffpe-trap=zero,overflow -O0")

ELSEIF(CMAKE_Fortran_COMPILER MATCHES "ifort")

   SET(CMAKE_Fortran_FLAGS_RELEASE "-nowarn -O3 -xHost -ip -openmp")
   SET(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check all,noarg_temp_created -fpe0 -traceback -ftrapuv -fpic")

ELSEIF(Fortran_COMPILER MATCHES "IBM")

   SET(CMAKE_Fortran_FLAGS_RELEASE "-WF,-qnotrigraph -qextname=flush -qthreaded -qhalt=e -qxlf2003=polymorphic")
   SET(CMAKE_Fortran_FLAGS_DEBUG "-WF,-qnotrigraph -qextname=flush -qthreaded -qhalt=e -qxlf2003=polymorphic")

ELSE()

   MESSAGE(STATUS "NO KNOWN FORTRAN COMPILER FOUND")

ENDIF()

MARK_AS_ADVANCED(CLEAR CMAKE_Fortran_COMPILER)

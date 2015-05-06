# Determine how-to install the modules. CMAKE_BINARY_DIR is the directory
# in which the make command is invoked.
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# Add the modules directory to the list of include directories
#INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

GET_FILENAME_COMPONENT(Fortran_COMPILER_NAME "${CMAKE_Fortran_COMPILER}" NAME)

IF (Fortran_COMPILER_NAME MATCHES ifort)
  EXEC_PROGRAM(${CMAKE_Fortran_COMPILER} ARGS "-v" OUTPUT_VARIABLE source_path)
  MESSAGE(STATUS "${source_path}")
  STRING(REGEX MATCH "1[0-9]\\.[0-9]\\.[0-9]" Fortran_COMPILER_VERSION ${source_path})
ELSE()
  EXEC_PROGRAM(${CMAKE_Fortran_COMPILER} ARGS "--version" OUTPUT_VARIABLE source_path)
  STRING(REGEX MATCH "4\\.[0-9]\\.[0-9]" Fortran_COMPILER_VERSION ${source_path})
ENDIF()
MESSAGE(STATUS "Fortran ${Fortran_COMPILER_NAME}-${Fortran_COMPILER_VERSION}")


IF(Fortran_COMPILER_NAME MATCHES gfortran)
  ADD_DEFINITIONS(-DGFORTRAN)
  SET(CMAKE_Fortran_FLAGS_RELEASE "-w -ffree-line-length-none -fall-intrinsics -O3 -fPIC")
  SET(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -cpp -pedantic -ffree-line-length-none -std=f2008 -fall-intrinsics -fbounds-check -fbacktrace -ffpe-trap=zero,overflow -O0 -fcheck-array-temporaries")

ELSEIF(Fortran_COMPILER_NAME MATCHES ifort)

  SET(CMAKE_Fortran_FLAGS_RELEASE "-nowarn -O3 -xHost -ip -fpic")
  SET(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check all,noarg_temp_created -fpe0 -traceback -ftrapuv -fpic")

ELSEIF(Fortran_COMPILER MATCHES "IBM")

  SET(CMAKE_Fortran_FLAGS_DEBUG "-qextname=flush -qxlf2003=polymorphic")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-qnosave -qextname=flush -qxlf2003=polymorphic")

ELSE()

  MESSAGE(SEND_ERROR "NO KNOWN FORTRAN COMPILER FOUND")

ENDIF()


IF(OPENMP_ENABLED)
  FIND_PACKAGE(OpenMP_Fortran)
  SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenMP_Fortran_FLAGS}")
  SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_Fortran_FLAGS}")
ENDIF()

SET(F2008 FALSE)

if(Fortran_COMPILER_NAME STREQUAL "gfortran")
  if(Fortran_COMPILER_VERSION VERSION_LESS "4.8.0")
    message(STATUS "Insufficient gfortran version for F2008")
  else()
    message(STATUS "GNU fortran version OK for F2008")
    SET(F2008 TRUE)
  endif()
elseif(Fortran_COMPILER_NAME STREQUAL "ifort")
  if(Fortran_COMPILER_VERSION VERSION_LESS "14.0.0")
    message(STATUS "Insufficient ifort version for F2008")
  else()
    message(STATUS "Intel fortran version OK for F2008")
    SET(F2008 TRUE)
  endif()
endif()

MARK_AS_ADVANCED(CLEAR CMAKE_Fortran_COMPILER)
MARK_AS_ADVANCED(CLEAR CMAKE_C_COMPILER)
MARK_AS_ADVANCED(CLEAR CMAKE_CXX_COMPILER)



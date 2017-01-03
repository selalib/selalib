# Determine how-to install the modules. CMAKE_BINARY_DIR is the directory
# in which the make command is invoked.
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# Add the modules directory to the list of include directories
#INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

GET_FILENAME_COMPONENT(Fortran_COMPILER_NAME "${CMAKE_Fortran_COMPILER}" NAME)
MESSAGE(STATUS "CMAKE_Fortran_COMPILER_ID:${CMAKE_Fortran_COMPILER_ID}")
SET(FULL_FORTRAN2003 FALSE)

IF (CMAKE_Fortran_COMPILER_ID MATCHES Intel)

  EXEC_PROGRAM(${CMAKE_Fortran_COMPILER} ARGS "-v" OUTPUT_VARIABLE source_path)
  MESSAGE(STATUS "${source_path}")
  STRING(REGEX MATCH "1[0-9]\\.[0-9]\\.[0-9]" Fortran_COMPILER_VERSION ${source_path})
  SET(CMAKE_Fortran_FLAGS_RELEASE "-nowarn -O3 -xHost -ip -fpic")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -check all,noarg_temp_created -fpe0 -traceback -ftrapuv -fpic")
  SET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "-shared-intel")

  if(Fortran_COMPILER_VERSION VERSION_LESS "14.0.0")
    message(STATUS "Insufficient ifort version for the F2003 standard")
  else()
    message(STATUS "Intel fortran version OK for the F2003 standard")
    SET(FULL_FORTRAN2003 TRUE)
  endif()

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES PGI)

  EXEC_PROGRAM(${CMAKE_Fortran_COMPILER} ARGS "--version" OUTPUT_VARIABLE source_path)
  STRING(REGEX MATCH "1[0-9]\\.[0-9][0-9]\\-[0-9]" Fortran_COMPILER_VERSION ${source_path})
  SET(CMAKE_Fortran_FLAGS_DEBUG "-Mbounds -O0 -g")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-acc -Minfo=accel -fast ")
  SET(FULL_FORTRAN2003 TRUE)

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES IBM)

  SET(CMAKE_Fortran_FLAGS_DEBUG   "-qextname=flush -qxlf2003=polymorphic")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-qnosave -qextname=flush -qxlf2003=polymorphic")
  SET(FULL_FORTRAN2003 TRUE)

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES GNU)

  EXEC_PROGRAM(${CMAKE_Fortran_COMPILER} ARGS "--version" OUTPUT_VARIABLE source_path)
  STRING(REGEX MATCH "[4-7]\\.[0-9]\\.[0-9]" Fortran_COMPILER_VERSION ${source_path})

  ADD_DEFINITIONS(-DGFORTRAN)
  SET(CMAKE_Fortran_FLAGS_RELEASE "-w -ffree-line-length-none -fall-intrinsics -O3 -fPIC -march=native ")
  IF(APPLE)
    SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -mno-avx")
  ENDIF(APPLE)
  SET(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -Wall -cpp -ffree-line-length-none -std=f2008 -pedantic -Wconversion -Wconversion-extra -Wintrinsics-std -fcheck=array-temps,bounds,do,pointer,recursion -fall-intrinsics -finit-real=snan -finit-integer=-9999 -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck-array-temporaries")

  SET(UNUSED_FUNCTION_WARNING_ENABLED OFF CACHE BOOL "Add -Wunused-function flag to gfortran")
  IF(NOT UNUSED_FUNCTION_WARNING_ENABLED)
    SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wno-unused-function")
  ENDIF()

  SET(UNUSED_DUMMY_WARNING_ENABLED OFF CACHE BOOL   "Add -Wunused-dummy-argument flag to gfortran")
  IF(NOT UNUSED_DUMMY_WARNING_ENABLED)
    SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wno-unused-dummy-argument")
  ENDIF()

  if(Fortran_COMPILER_VERSION VERSION_LESS "4.8.0")
    message(STATUS "Insufficient gfortran version for the Fortran 2003 Standard")
  else()
    message(STATUS "GNU fortran version OK for the Fortran 2003 standard")
    SET(FULL_FORTRAN2003 TRUE)
  endif()

ELSE()

  MESSAGE(SEND_ERROR "NO KNOWN FORTRAN COMPILER FOUND")

ENDIF()

MESSAGE(STATUS "Fortran ${Fortran_COMPILER_NAME}-${Fortran_COMPILER_VERSION}")

SET(ADDITIONAL_COMPILER_FLAGS "" CACHE STRING "The user can define additional compiler flags here")
SET(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} ${ADDITIONAL_COMPILER_FLAGS}")
SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${ADDITIONAL_COMPILER_FLAGS}")

IF(OPENMP_ENABLED)
  FIND_PACKAGE(OpenMP_Fortran)
  SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenMP_Fortran_FLAGS}")
  SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_Fortran_FLAGS}")
ENDIF()

IF(FULL_FORTRAN2003)
  ADD_DEFINITIONS(-DFULL_FORTRAN2003)
ENDIF(FULL_FORTRAN2003)

MARK_AS_ADVANCED(CLEAR CMAKE_Fortran_COMPILER)
MARK_AS_ADVANCED(CLEAR CMAKE_C_COMPILER)
MARK_AS_ADVANCED(CLEAR CMAKE_CXX_COMPILER)



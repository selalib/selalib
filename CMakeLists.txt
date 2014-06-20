CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING 
   "Choose the type of build, options are: Debug Release.")
ELSE()
   SET(CMAKE_BUILD_TYPE Debug CACHE STRING 
   "Choose the type of build, options are: Debug Release.")
ENDIF()

IF(OPTIONS_FILE)
INCLUDE(${OPTIONS_FILE})
ENDIF()

PROJECT(selalib)
ENABLE_LANGUAGE(Fortran)

# local modules path
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake
                      ${CMAKE_SOURCE_DIR}/prototype/src/CMakeModules)

INCLUDE(F90FLAGS)
# Tests+Dashboard configuration
INCLUDE(CTest)

SET(CTEST_DROP_SITE_CDASH TRUE)
ENABLE_TESTING()
INCLUDE(Pipol)

MESSAGE(STATUS "CMAKE_Fortran_COMPILER:${CMAKE_Fortran_COMPILER}")
IF (${CMAKE_Fortran_COMPILER} MATCHES "ifort")
  EXEC_PROGRAM(${CMAKE_Fortran_COMPILER} ARGS "-v" OUTPUT_VARIABLE source_path)
  STRING(REGEX MATCH "1[0-9]\\.[0-9]\\.[0-9]" Fortran_COMPILER_VERSION ${source_path})
ELSE()
  EXEC_PROGRAM(${CMAKE_Fortran_COMPILER} ARGS "--version" OUTPUT_VARIABLE source_path)
  STRING(REGEX MATCH "4\\.[0-9]\\.[0-9]" Fortran_COMPILER_VERSION ${source_path})
ENDIF()
MESSAGE(STATUS "Fortran_COMPILER_VERSION:${Fortran_COMPILER_VERSION}")
GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

IF(PIPOL_IMAGE)
  SET(BUILDNAME "${Fortran_COMPILER_NAME}-${Fortran_COMPILER_VERSION}-${PIPOL_IMAGE_NAME}")
  SET(SITE "${PIPOL_SITE}")
ELSE(PIPOL_IMAGE)
  SET(BUILDNAME "${Fortran_COMPILER_NAME}-${Fortran_COMPILER_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
ENDIF(PIPOL_IMAGE)

MESSAGE(STATUS "BUILDNAME:${BUILDNAME}")

SET(HDF5_PARALLEL_ENABLED ON CACHE BOOL "Use Parallel HDF5" FORCE)
SET(FFTW_ENABLED          ON CACHE BOOL "Use fftw library " FORCE)
ADD_DEFINITIONS(-D_FFTW)
SET(FFT_LIB "FFTW" CACHE STRING "fftw library used" FORCE)

SET(SLL_BUILD_PACKAGE OFF CACHE BOOL "Build package only")

IF(SLL_BUILD_PACKAGE)

   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
   INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)
   FIND_PACKAGE(MPI REQUIRED Fortran)
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   INCLUDE(HDF5Config)
   IF(NOT ${HDF5_IS_PARALLEL})
      MESSAGE(SEND_ERROR "HDF5 must be parallel")
   ENDIF()
   FIND_PACKAGE(FFTW)
   ADD_SUBDIRECTORY (prototype/src/fftpack)
   ADD_SUBDIRECTORY (prototype/src/package)
   # Packages generation
   INCLUDE(CPackConfig) 

ELSE()

   ADD_SUBDIRECTORY (prototype/src)

   IF(NOT STDF95)
      #ADD_SUBDIRECTORY (slv2d/src)
      #ADD_SUBDIRECTORY (slv2d/vlasov_4d)
      #ADD_SUBDIRECTORY (pic2d/src)
      #ADD_SUBDIRECTORY (gysela_public/src)
   ENDIF(NOT STDF95)
   
   # Pipol rc and nightly files are installed under share/pipol
   FILE(GLOB rc_files  $ENV{HOME}/.pipol/rc $ENV{HOME}/rc.*)
   FOREACH(_f ${rc_files})
     INSTALL(FILES ${_f} DESTINATION share/pipol/selalib)
   ENDFOREACH(_f ${rc_files})
   
   FILE(GLOB nightly_files $ENV{HOME}/.pipol/nightly/*)
   FOREACH(_f ${nightly_files})
     INSTALL(FILES ${_f} DESTINATION share/pipol/selalib/nightly)
   ENDFOREACH(_f ${nightly_files})
   
   SET(CMAKE_INCLUDE_CURRENT_DIR ON)

   INSTALL(FILES AUTHORS ChangeLog COPYING.txt INSTALL README.txt
           DESTINATION share/doc/pipol/)
   
   INCLUDE(DartConfig)
   INCLUDE(Dart)

ENDIF(SLL_BUILD_PACKAGE)

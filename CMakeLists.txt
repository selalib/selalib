CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: Debug Release.")
ELSE()
   SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build,
options are: Debug Release.")
ENDIF()


PROJECT(SELALIB)
ENABLE_LANGUAGE(Fortran)

# local modules path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake
                      ${CMAKE_SOURCE_DIR}/prototype/src/CMakeModules)

include(F90FLAGS)
# Tests+Dashboard configuration
include(CTest)



SET(CTEST_DROP_SITE_CDASH TRUE)
enable_testing()
include(Pipol)

MESSAGE(STATUS "CMAKE_Fortran_COMPILER:${CMAKE_Fortran_COMPILER}")
if (${CMAKE_Fortran_COMPILER} MATCHES "ifort")
  EXEC_PROGRAM(${CMAKE_Fortran_COMPILER} ARGS "-v" OUTPUT_VARIABLE source_path)
  string(REGEX MATCH "1[0-9]\\.[0-9]\\.[0-9]" Fortran_COMPILER_VERSION ${source_path})
else()
  EXEC_PROGRAM(${CMAKE_Fortran_COMPILER} ARGS "--version" OUTPUT_VARIABLE source_path)
  string(REGEX MATCH "4\\.[0-9]\\.[0-9]" Fortran_COMPILER_VERSION ${source_path})
endif()
MESSAGE(STATUS "Fortran_COMPILER_VERSION:${Fortran_COMPILER_VERSION}")
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if(PIPOL_IMAGE)
  set(BUILDNAME "${Fortran_COMPILER_NAME}-${Fortran_COMPILER_VERSION}-${PIPOL_IMAGE_NAME}-${CMAKE_BUILD_TYPE}")
  set(SITE "${PIPOL_SITE}")
else(PIPOL_IMAGE)
  set(BUILDNAME "${Fortran_COMPILER_NAME}-${Fortran_COMPILER_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_BUILD_TYPE}")
endif(PIPOL_IMAGE)

MESSAGE(STATUS "BUILDNAME:${BUILDNAME}")

SET(HDF5_PARALLEL_ENABLED ON CACHE BOOL "Use Parallel HDF5" FORCE)
SET(FFTW_ENABLED          ON CACHE BOOL "Use fftw library " FORCE)
ADD_DEFINITIONS(-D_FFTW)

SET(SLL_BUILD_PACKAGE OFF CACHE BOOL "Build package only")

IF(SLL_BUILD_PACKAGE)
   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
   INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)
   FIND_PACKAGE(MPI REQUIRED Fortran)
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   INCLUDE(HDF5Config)
   FIND_PACKAGE(FFTW)
   add_subdirectory (prototype/src/package)
   # Packages generation
   include(CPackConfig) 
ELSE()
   add_subdirectory (prototype/src)

   IF(NOT STDF95)
      add_subdirectory (slv2d/src)
      add_subdirectory (pic2d/src)
      add_subdirectory (gysela_public/src)
   ENDIF(NOT STDF95)
   
   # Pipol rc and nightly files are installed under share/pipol
   file(GLOB rc_files  $ENV{HOME}/.pipol/rc $ENV{HOME}/rc.*)
   foreach(_f ${rc_files})
     install(FILES ${_f} DESTINATION share/pipol/selalib)
   endforeach(_f ${rc_files})
   
   file(GLOB nightly_files $ENV{HOME}/.pipol/nightly/*)
   foreach(_f ${nightly_files})
     install(FILES ${_f} DESTINATION share/pipol/selalib/nightly)
   endforeach(_f ${nightly_files})
   
   set(CMAKE_INCLUDE_CURRENT_DIR ON)
   install(FILES AUTHORS ChangeLog COPYING.txt INSTALL README.txt
           DESTINATION share/doc/pipol/)
   
   include(DartConfig)
   include(Dart)

ENDIF(SLL_BUILD_PACKAGE)

cmake_minimum_required(VERSION 2.8)

PROJECT(SELALIB)
ENABLE_LANGUAGE(Fortran)

# local modules path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake
                      ${CMAKE_SOURCE_DIR}/prototype/src/CMakeModules)

include(CTest)
set(CTEST_DROP_SITE_CDASH TRUE)
# Tests+Dashboard configuration
enable_testing()
include(Pipol)

EXEC_PROGRAM(${CMAKE_Fortran_COMPILER} ARGS --version OUTPUT_VARIABLE source_path)
string(REGEX MATCH "1[0-9]\\.[0-9]\\.[0-9]" Fortran_COMPILER_VERSION ${source_path})
MESSAGE(STATUS "Fortran_COMPILER_VERSION:${Fortran_COMPILER_VERSION}")

MESSAGE(STATUS "CMAKE_Fortran_COMPILER:${CMAKE_Fortran_COMPILER}")
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if(PIPOL_IMAGE)
  set(BUILDNAME "${Fortran_COMPILER_NAME}-${Fortran_COMPILER_VERSION}-${PIPOL_IMAGE_NAME}")
  set(SITE "${PIPOL_SITE}")
else(PIPOL_IMAGE)
  set(BUILDNAME "${Fortran_COMPILER_NAME}-${Fortran_COMPILER_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
endif(PIPOL_IMAGE)

MESSAGE(STATUS "BUILDNAME:${BUILDNAME}")

SET(HDF5_PARALLEL_ENABLED ON CACHE BOOL "Use Parallel HDF5" FORCE)
SET(FFTW_ENABLED          ON CACHE BOOL "Use fftw library " FORCE)
add_subdirectory (prototype/src)
add_subdirectory (slv2d/src)
add_subdirectory (pic2d/src)

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

# Packages generation
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "SeLaLib prototype")
set(CPACK_PACKAGE_DESCRIPTION 
  "This Pipol package steps are: 
     1. architecture reservation and OS deployement
     2. git update from gforge 
     3. compilation
     4. test and visualisation on a dashboard
     5. packages generation
     6. gforge packages upload 
  ")

set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING.txt")

set(CPACK_PACKAGE_NAME "selalib-prototype")

set(CPACK_PACKAGE_VERSION_MAJOR 1)
set(CPACK_PACKAGE_VERSION_MINOR 0)
set(CPACK_PACKAGE_VERSION_PATCH 0)

IF(PIPOL_IMAGE)
  SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-1.2-${PIPOL_IMAGE_NAME}")
ENDIF(PIPOL_IMAGE)

set(CPACK_PACKAGE_CONTACT "Pierre Navaro 0368850173")

include(CPack)

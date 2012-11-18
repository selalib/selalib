cmake_minimum_required(VERSION 2.8)

PROJECT(SELALIB)

# local modules path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake
                      ${CMAKE_SOURCE_DIR}/prototype/src/CMakeModules)

include(CTest)
# Tests+Dashboard configuration
enable_testing()
include(Pipol)


if(PIPOL_IMAGE)
  set(BUILDNAME "${CMAKE_GENERATOR}-${CMAKE_Fortran_COMPILER_ID}-${PIPOL_IMAGE_NAME}")
  set(SITE "${PIPOL_SITE}")
else(PIPOL_IMAGE)
  set(BUILDNAME "${CMAKE_Fortran_COMPILER}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
  message(STATUS "${CMAKE_Fortran_COMPILER}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
endif(PIPOL_IMAGE)

SET(HDF5_PARALLEL_ENABLED ON CACHE BOOL "Use Parallel HDF5" FORCE)
SET(FFTW_ENABLED          ON CACHE BOOL "Use fftw library " FORCE)
add_subdirectory (prototype/src)
add_subdirectory (slv2d/src)

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

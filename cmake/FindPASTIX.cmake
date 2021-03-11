# Pastix lib requires linking to a blas library. It is up to the user of this
# module to find a BLAS and link to it. Pastix requires SCOTCH or METIS
# (partitioning and reordering tools) as well

if(DEFINED ENV{PASTIX_ROOT})
  set(PASTIX_ROOT
      $ENV{PASTIX_ROOT}
      CACHE PATH "pastix location")
else()
  set(PASTIX_ROOT
      /usr/local
      CACHE PATH "pastix location")
endif()

find_path(
  PASTIX_INCLUDE_DIRS
  NAMES pastix_fortran.h
  HINTS ${PASTIX_ROOT}
  PATH_SUFFIXES include Include INCLUDE
  DOC "PATH TO pastix_fortran.h")

find_library(
  PASTIX_LIBRARY
  NAMES pastix
  HINTS ${PASTIX_ROOT}
  PATH_SUFFIXES lib Lib LIB
  DOC "PATH TO libpastix.a")

find_library(
  PASTIX_MATRIX_DRIVER_LIBRARY
  NAMES matrix_driver
  HINTS ${PASTIX_ROOT}
  PATH_SUFFIXES lib Lib LIB
  DOC "PATH TO libmatrix_driver.a")

find_library(
  HWLOC_LIBRARY
  NAMES hwloc
  HINTS $ENV{HWLOC_HOME} ${PASTIX_ROOT} /usr/local /opt/local
  PATH_SUFFIXES lib Lib LIB lib64
  DOC "PATH TO hwloc library")

set(PASTIX_LIBRARIES
    ${PASTIX_LIBRARY};${PASTIX_MATRIX_DRIVER_LIBRARY};${HWLOC_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PASTIX DEFAULT_MSG PASTIX_INCLUDE_DIRS
                                  PASTIX_LIBRARIES)

if(PASTIX_FOUND)
  message(STATUS "PASTIX_INCLUDE_DIRS:${PASTIX_INCLUDE_DIRS}")
  message(STATUS "PASTIX_LIBRARIES:${PASTIX_LIBRARIES}")
  add_definitions(-DPASTIX)
endif(PASTIX_FOUND)

find_path(
  MURGE_INCLUDE_DIRS
  NAMES murge.inc
  HINTS ${PASTIX_ROOT}
  PATH_SUFFIXES include Include INCLUDE
  DOC "PATH TO murge.inc")

find_library(
  MURGE_LIBRARIES
  NAMES pastix_murge
  HINTS ${PASTIX_ROOT}
  PATH_SUFFIXES lib Lib LIB
  DOC "PATH TO libpastix_murge.a")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MURGE DEFAULT_MSG MURGE_INCLUDE_DIRS
                                  MURGE_LIBRARIES)

mark_as_advanced(
  PASTIX_INCLUDE_DIRS
  PASTIX_LIBRARY
  PASTIX_MATRIX_DRIVER_LIBRARY
  PASTIX_LIBRARIES
  HWLOC_LIBRARY
  MURGE_INCLUDE_DIRS
  MURGE_LIBRARIES)

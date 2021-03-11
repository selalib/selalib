# MUMPS lib requires linking to a blas library. It is up to the user of this
# module to find a BLAS and link to it. MUMPS requires SCOTCH or METIS
# (partitioning and reordering tools) as well

if(DEFINED ENV{MUMPS_ROOT})
  set(MUMPS_ROOT
      $ENV{MUMPS_ROOT}
      CACHE PATH "mumps location")
else()
  set(MUMPS_ROOT
      /usr/local
      CACHE PATH "mumps location")
endif()

find_path(
  MUMPS_INCLUDE_DIRS
  NAMES dmumps_struc.h
  HINTS ${MUMPS_ROOT}
  PATH_SUFFIXES include Include INCLUDE
  DOC "PATH TO dmumps_struc.h")

find_library(
  MUMPS_LIBRARY
  NAMES dmumps
  HINTS ${MUMPS_ROOT}
  PATH_SUFFIXES lib Lib LIB
  DOC "PATH TO libdmumps.a")

find_library(
  MUMPS_COMMMON_LIBRARY
  NAMES mumps_common
  HINTS ${MUMPS_ROOT}
  PATH_SUFFIXES lib Lib LIB
  DOC "PATH TO libmumps_common.a")

find_library(
  MUMPS_SIMPLE_LIBRARY
  NAMES mumps_simple
  HINTS ${MUMPS_ROOT}
  PATH_SUFFIXES lib Lib LIB
  DOC "PATH TO libmumps_simple.a")

set(MUMPS_LIBRARIES
    ${MUMPS_LIBRARY};${MUMPS_COMMON_LIBRARY};${MUMPS_SIMPLE_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS DEFAULT_MSG MUMPS_INCLUDE_DIRS
                                  MUMPS_LIBRARIES)

if(MUMPS_FOUND)
  message(STATUS "MUMPS_INCLUDE_DIRS:${MUMPS_INCLUDE_DIRS}")
  message(STATUS "MUMPS_LIBRARIES:${MUMPS_LIBRARIES}")
  add_definitions(-DMUMPS)
endif(MUMPS_FOUND)

mark_as_advanced(MUMPS_INCLUDE_DIRS MUMPS_LIBRARY MUMPS_COMMON_LIBRARY
                 MUMPS_SIMPLE_LIBRARY MUMPS_LIBRARIES)

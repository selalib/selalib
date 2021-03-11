# * Try to find SUITESPARSE Once done this will define
#
# SUITESPARSE_FOUND            - system has SUITESPARSE SUITESPARSE_INCLUDE_DIRS
# - the SUITESPARSE include directory SUITESPARSE_LIBRARIES        - Link these
# to use SUITESPARSE SUITESPARSE_LIBRARY_DIR      - Library main directory
# containing suitesparse libs SUITESPARSE_LIBRARY_DIRS     - all Library
# directories containing suitesparse libs

set(SUITESPARSE_ROOT
    "/usr"
    CACHE PATH "Root directory for SuiteSParse library")

if(SUITESPARSE_INCLUDE_DIRS)
  # Already in cache, be silent
  set(SUITESPARSE_FIND_QUIETLY TRUE)
endif(SUITESPARSE_INCLUDE_DIRS)

find_path(
  CHOLMOD_INCLUDE_DIR cholmod.h
  PATHS ${SUITESPARSE_ROOT} /opt/local /usr/local
  PATH_SUFFIXES include include/suitesparse)

macro(FIND_SUITESPARSE_LIBRARY LIB)

  set(_LIB ${LIB}-NOTFOUND)
  find_library(
    _LIB
    NAMES ${LIB}
    PATHS ${SUITESPARSE_ROOT} /opt/local /usr/local
    PATH_SUFFIXES lib lib64)

  if(_LIB)
    list(APPEND SUITESPARSE_LIBRARIES ${_LIB})
  endif()

endmacro(FIND_SUITESPARSE_LIBRARY)

# Add cholmod include directory to collection include directories
if(CHOLMOD_INCLUDE_DIR)
  list(APPEND SUITESPARSE_INCLUDE_DIRS ${CHOLMOD_INCLUDE_DIR})
endif(CHOLMOD_INCLUDE_DIR)

find_suitesparse_library(umfpack)
find_suitesparse_library(cholmod)
find_suitesparse_library(amd)
find_suitesparse_library(colamd)
find_suitesparse_library(suitesparseconfig)
find_suitesparse_library(rt)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUITESPARSE DEFAULT_MSG SUITESPARSE_LIBRARIES
                                  SUITESPARSE_INCLUDE_DIRS)
mark_as_advanced(SUITESPARSE_LIBRARIES SUITESPARSE_INCLUDE_DIRS)

if(SUITESPARSE_FOUND)
  message(STATUS "SUITESPARSE_INCLUDE_DIRS:${SUITESPARSE_INCLUDE_DIRS}")
  message(STATUS "SUITESPARSE_LIBRARIES:${SUITESPARSE_LIBRARIES}")
  add_definitions(-DUMFPACK)
endif(SUITESPARSE_FOUND)

unset(_LIB CACHE)

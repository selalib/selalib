# Pastix lib requires linking to a blas library. It is up to the user of this
# module to find a BLAS and link to it. Pastix requires SCOTCH

if(DEFINED ENV{SCOTCH_ROOT})
  set(SCOTCH_ROOT
      $ENV{SCOTCH_ROOT}
      CACHE PATH "scotch library location")
else()
  set(SCOTCH_ROOT
      /usr/local
      CACHE PATH "scotch library location")
endif()

find_path(
  SCOTCH_INCLUDE_DIRS
  NAMES scotch.h scotchf.h
  HINTS ${SCOTCH_ROOT}
  PATH_SUFFIXES include Include INCLUDE
  DOC "PATH TO scotch.h and scotchf.h")

find_library(
  SCOTCH_LIBRARY
  NAMES scotch
  HINTS ${SCOTCH_ROOT}
  PATH_SUFFIXES lib Lib LIB
  DOC "PATH TO libscotch.a")

find_library(
  SCOTCHERR_LIBRARY
  NAMES scotcherr
  HINTS ${SCOTCH_ROOT}
  PATH_SUFFIXES lib Lib LIB
  DOC "PATH TO libscotcherr.a")

find_library(
  SCOTCHERREXIT_LIBRARY
  NAMES scotcherrexit
  HINTS ${SCOTCH_ROOT}
  PATH_SUFFIXES lib Lib LIB
  DOC "PATH TO libscotcherrexit.a")

set(SCOTCH_LIBRARIES ${SCOTCH_LIBRARY} ${SCOTCHERR_LIBRARY}
                     ${SCOTCHERREXIT_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCOTCH DEFAULT_MSG SCOTCH_INCLUDE_DIRS
                                  SCOTCH_LIBRARIES)

if(SCOTCH_FOUND)
  message(STATUS "SCOTCH_INCLUDE_DIRS:${SCOTCH_INCLUDE_DIRS}")
  message(STATUS "SCOTCH_LIBRARIES:${SCOTCH_LIBRARIES}")
endif(SCOTCH_FOUND)

find_path(
  PTSCOTCH_INCLUDE_DIRS
  NAMES ptscotch.h ptscotchf.h
  HINTS ${SCOTCH_ROOT}
  PATH_SUFFIXES include Include INCLUDE
  DOC "PATH TO ptscotch.h and ptscotchf.h")
find_library(
  PTSCOTCH_LIBRARY
  NAMES ptscotch
  HINTS ${SCOTCH_ROOT}
  PATH_SUFFIXES lib Lib LIB
  DOC "PATH TO libptscotch.a")
find_library(
  PTSCOTCHERR_LIBRARY
  NAMES ptscotcherr
  HINTS ${SCOTCH_ROOT}
  PATH_SUFFIXES lib Lib LIB
  DOC "PATH TO libptscotcherr.a")
find_library(
  PTSCOTCHERREXIT_LIBRARY
  NAMES ptscotcherrexit
  HINTS ${SCOTCH_ROOT}
  PATH_SUFFIXES lib Lib LIB
  DOC "PATH TO libptscotcherrexit.a")

set(PTSCOTCH_LIBRARIES ${PTSCOTCH_LIBRARY} ${PTSCOTCHERR_LIBRARY}
                       ${PTSCOTCHERREXIT_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PTSCOTCH DEFAULT_MSG PTSCOTCH_INCLUDE_DIRS
                                  PTSCOTCH_LIBRARIES)
if(PTSCOTCH_FOUND)
  message(STATUS "PTSCOTCH_INCLUDE_DIR:${PTSCOTCH_INCLUDE_DIRS}")
  message(STATUS "PTSCOTCH_LIBRARIES:${PTSCOTCH_LIBRARIES}")
endif(PTSCOTCH_FOUND)

mark_as_advanced(
  PTSCOTCHERREXIT_LIBRARY
  PTSCOTCHERR_LIBRARY
  PTSCOTCH_LIBRARY
  PTSCOTCH_INCLUDE_DIRS
  SCOTCHERREXIT_LIBRARY
  SCOTCHERR_LIBRARY
  SCOTCH_INCLUDE_DIRS
  SCOTCH_LIBRARY)

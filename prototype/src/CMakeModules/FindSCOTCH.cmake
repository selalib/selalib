# Pastix lib requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.
# Pastix requires SCOTCH 


find_path(SCOTCH_INCLUDE_DIRS
		NAMES scotch.h scotchf.h 
		HINTS ${SCOTCH_ROOT}
		PATH_SUFFIXES include Include INCLUDE
		DOC "PATH TO scotch.h and scotchf.h")
find_library(SCOTCH_LIBRARY
		NAMES scotch
		HINTS ${SCOTCH_ROOT}
		PATH_SUFFIXES lib Lib LIB
		DOC "PATH TO libscotch.a")
find_library(SCOTCHERR_LIBRARY
		NAMES scotcherr
		HINTS ${SCOTCH_ROOT}
		PATH_SUFFIXES lib Lib LIB
		DOC "PATH TO libscotcherr.a")
find_library(SCOTCHERREXIT_LIBRARY
		NAMES scotcherrexit
		HINTS ${SCOTCH_ROOT}
		PATH_SUFFIXES lib Lib LIB
		DOC "PATH TO libscotcherrexit.a")

set (SCOTCH_LIBRARIES ${SCOTCH_LIBRARY};${SCOTCHERR_LIBRARY};${SCOTCHERREXIT_LIBRARY})
if (SCOTCH_INCLUDE_DIRS AND SCOTCH_LIBRARIES)
  set(SCOTCH_FOUND YES)
endif(SCOTCH_INCLUDE_DIRS AND SCOTCH_LIBRARIES)

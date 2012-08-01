# Pastix lib requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.
# Pastix requires PTSCOTCH 


find_path(PTSCOTCH_INCLUDE_DIRS
		NAMES scotch.h scotchf.h 
		HINTS ${PTSCOTCH_ROOT}
		PATH_SUFFIXES include Include INCLUDE
		DOC "PATH TO scotch.h and scotchf.h")
find_library(PTSCOTCH_LIBRARY
		NAMES scotch
		HINTS ${PTSCOTCH_ROOT}
		PATH_SUFFIXES lib Lib LIB
		DOC "PATH TO libscotch.a")
find_library(PTSCOTCHERR_LIBRARY
		NAMES scotcherr
		HINTS ${PTSCOTCH_ROOT}
		PATH_SUFFIXES lib Lib LIB
		DOC "PATH TO libscotcherr.a")
find_library(PTSCOTCHERREXIT_LIBRARY
		NAMES scotcherrexit
		HINTS ${PTSCOTCH_ROOT}
		PATH_SUFFIXES lib Lib LIB
		DOC "PATH TO libscotcherrexit.a")

set (PTSCOTCH_LIBRARIES ${PTSCOTCH_LIBRARY};${PTSCOTCHERR_LIBRARY};${PTSCOTCHERREXIT_LIBRARY})
if (PTSCOTCH_INCLUDE_DIRS AND PTSCOTCH_LIBRARIES)
  set(PTSCOTCH_FOUND YES)
endif(PTSCOTCH_INCLUDE_DIRS AND PTSCOTCH_LIBRARIES)

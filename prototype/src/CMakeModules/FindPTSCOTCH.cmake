# Pastix lib requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.
# Pastix requires PTSCOTCH 


find_path(PTSCOTCH_INCLUDE_DIRS
		NAMES ptscotch.h ptscotchf.h 
		HINTS ${PTSCOTCH_ROOT}
		PATH_SUFFIXES include Include INCLUDE
		DOC "PATH TO ptscotch.h and ptscotchf.h")
find_library(PTSCOTCH_LIBRARY
		NAMES ptscotch
		HINTS ${PTSCOTCH_ROOT}
		PATH_SUFFIXES lib Lib LIB
		DOC "PATH TO libptscotch.a")
find_library(PTSCOTCHERR_LIBRARY
		NAMES ptscotcherr
		HINTS ${PTSCOTCH_ROOT}
		PATH_SUFFIXES lib Lib LIB
		DOC "PATH TO libptscotcherr.a")
find_library(PTSCOTCHERREXIT_LIBRARY
		NAMES ptscotcherrexit
		HINTS ${PTSCOTCH_ROOT}
		PATH_SUFFIXES lib Lib LIB
		DOC "PATH TO libptscotcherrexit.a")

set (PTSCOTCH_LIBRARIES @PTSCOTCH_LIBRARY@;@PTSCOTCHERR_LIBRARY@;@PTSCOTCHERREXIT_LIBRARY@)
if (PTSCOTCH_INCLUDE_DIRS AND PTSCOTCH_LIBRARIES)
  set(PTSCOTCH_FOUND YES)
endif(PTSCOTCH_INCLUDE_DIRS AND PTSCOTCH_LIBRARIES)

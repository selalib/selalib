# Pastix lib requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.
# Pastix requires MURGE 

find_path(MURGE_INCLUDE_DIRS
		NAMES murge.inc
		HINTS ${MURGE_ROOT}
		PATH_SUFFIXES include Include INCLUDE
		DOC "PATH TO murge.inc")
find_library(MURGE_LIBRARIES
		NAMES pastix_murge 
		HINTS ${MURGE_ROOT}
		PATH_SUFFIXES lib Lib LIB
		DOC "PATH TO libpastix_murge.a")
if (MURGE_INCLUDE_DIRS AND MURGE_LIBRARIES)
  set(MURGE_FOUND YES)
endif(MURGE_INCLUDE_DIRS AND MURGE_LIBRARIES)

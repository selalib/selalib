# Pastix lib requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.
# Pastix requires MURGE 

FIND_PATH(MURGE_INCLUDE_DIRS
	    NAMES murge.inc
	    HINTS ${PASTIX_ROOT}
	    PATH_SUFFIXES include Include INCLUDE
	    DOC "PATH TO murge.inc")

FIND_LIBRARY(MURGE_LIBRARIES
		 NAMES pastix_murge 
		 HINTS ${PASTIX_ROOT}
		 PATH_SUFFIXES lib Lib LIB
		 DOC "PATH TO libpastix_murge.a")

IF (MURGE_INCLUDE_DIRS AND MURGE_LIBRARIES)

  SET(MURGE_FOUND YES)

ENDIF(MURGE_INCLUDE_DIRS AND MURGE_LIBRARIES)

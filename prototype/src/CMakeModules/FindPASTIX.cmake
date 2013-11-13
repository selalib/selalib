# Pastix lib requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.
# Pastix requires SCOTCH or METIS (partitioning and reordering tools) as well


FIND_PATH(PASTIX_INCLUDE_DIRS
	    NAMES pastix_fortran.h
	    HINTS ${PASTIX_ROOT}
	    PATH_SUFFIXES include Include INCLUDE
	    DOC "PATH TO pastix_fortran.h")

FIND_LIBRARY(PASTIX_LIBRARY NAMES pastix
		 HINTS ${PASTIX_ROOT}
		 PATH_SUFFIXES lib Lib LIB
		 DOC "PATH TO libpastix.a")

FIND_LIBRARY(PASTIX_MATRIX_DRIVER_LIBRARY NAMES matrix_driver
		 HINTS ${PASTIX_ROOT}
		 PATH_SUFFIXES lib Lib LIB
		 DOC "PATH TO libmatrix_driver.a")

SET(PASTIX_LIBRARIES ${PASTIX_LIBRARY};${PASTIX_MATRIX_DRIVER_LIBRARY})

IF (PASTIX_INCLUDE_DIRS AND PASTIX_LIBRARIES)

  SET(PASTIX_FOUND YES)

ENDIF(PASTIX_INCLUDE_DIRS AND PASTIX_LIBRARIES)

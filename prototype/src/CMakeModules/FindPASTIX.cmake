# Pastix lib requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.
# Pastix requires SCOTCH or METIS (partitioning and reordering tools) as well


find_path(PASTIX_INCLUDE_DIRS
		NAMES pastix.h pastix_fortran.h
		HINTS ${PASTIX_ROOT}
		PATH_SUFFIXES include Include INCLUDE
		DOC "PATH TO pastix.h and pastix_fortran.h")
find_library(PASTIX_LIBRARIES
		NAMES pastix
		HINTS ${PASTIX_ROOT}
		PATH_SUFFIXES lib Lib LIB
		DOC "PATH TO libpastix.a")
if (PASTIX_INCLUDE_DIRS AND PASTIX_LIBRARIES)
  set(PASTIX_FOUND YES)
endif(PASTIX_INCLUDE_DIRS AND PASTIX_LIBRARIES)

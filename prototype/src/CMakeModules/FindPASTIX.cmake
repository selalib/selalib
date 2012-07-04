find_path(PASTIX_INCLUDE_DIRS NAMES pastix.h pastix_fortran.h pastix_nompi.h
	HINTS ${PASTIX_ROOT}
	PATH_SUFFIXES include
	DOC "PATH TO pastix.h pastix_fortran.h pastix_nompi.h")

find_library(PASTIX_PASTIX_LIBRARY NAMES libpastix.a and libpastix_murge.a
	HINTS ${PASTIX_ROOT}
	PATH_SUFFIXES lib 
	DOC "PATH TO libpastix.a libpastix_murge.a")

set(PASTIX_LIBRARIES @PASTIX_PASTIX_LIBRARY@)

if (PASTIX_INCLUDE_DIRS AND
	PASTIX_PASTIX_LIBRARY)
  set(PASTIX_FOUND YES)
  include_directories(${PASTIX_INCLUDE_DIRS})
endif()

# - Try to find SUITESPARSE
# Once done this will define
#  
#  SUITESPARSE_FOUND            - system has SUITESPARSE
#  SUITESPARSE_INCLUDE_DIRS     - the SUITESPARSE include directory
#  SUITESPARSE_LIBRARIES        - Link these to use SUITESPARSE
#  SUITESPARSE_LIBRARY_DIR      - Library main directory containing suitesparse libs
#  SUITESPARSE_LIBRARY_DIRS     - all Library directories containing suitesparse libs

IF (SUITESPARSE_INCLUDE_DIRS)
  # Already in cache, be silent
  SET(SUITESPARSE_FIND_QUIETLY TRUE)
ENDIF (SUITESPARSE_INCLUDE_DIRS)

FIND_PATH( CHOLMOD_INCLUDE_DIR cholmod.h
           PATHS /usr/local/include 
                 /usr/include 
        	 /usr/include/suitesparse/ 
                 /opt/local/include)

MESSAGE(STATUS "CHOLMOD_INCLUDE_DIR:${CHOLMOD_INCLUDE_DIR}")

MACRO(FIND_SUITESPARSE_LIBRARY LIB )

SET(_LIB ${LIB}-NOTFOUND)
FIND_LIBRARY( _LIB
              NAMES ${LIB}
              PATHS /usr/lib /usr/lib64 /usr/local/lib /opt/local/lib )

IF(_LIB)
   LIST ( APPEND SUITESPARSE_LIBRARIES ${_LIB})
ENDIF()

ENDMACRO(FIND_SUITESPARSE_LIBRARY)


# Add cholmod include directory to collection include directories
IF ( CHOLMOD_INCLUDE_DIR )
   LIST ( APPEND SUITESPARSE_INCLUDE_DIRS ${CHOLMOD_INCLUDE_DIR} )
ENDIF( CHOLMOD_INCLUDE_DIR )

FIND_SUITESPARSE_LIBRARY(umfpack)
FIND_SUITESPARSE_LIBRARY(cholmod)
FIND_SUITESPARSE_LIBRARY(amd)
FIND_SUITESPARSE_LIBRARY(colamd)
FIND_SUITESPARSE_LIBRARY(suitesparseconfig)
FIND_SUITESPARSE_LIBRARY(rt)
   
IF (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)
   SET(SUITESPARSE_FOUND TRUE)
   MESSAGE(STATUS "SUITESPARSE_INCLUDE_DIRS:${SUITESPARSE_INCLUDE_DIRS}")
   MESSAGE(STATUS "SUITESPARSE_LIBRARIES:${SUITESPARSE_LIBRARIES}")
ELSE (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)
   SET( SUITESPARSE_FOUND FALSE )
ENDIF (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)


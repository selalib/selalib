# - Try to find SUITESPARSE
# Once done this will define
#  
#  SUITESPARSE_FOUND            - system has SUITESPARSE
#  SUITESPARSE_INCLUDE_DIRS     - the SUITESPARSE include directory
#  SUITESPARSE_LIBRARIES        - Link these to use SUITESPARSE
#  SUITESPARSE_SPQR_LIBRARY     - name of spqr library (necessary due to error in debian package)
#  SUITESPARSE_SPQR_LIBRARY_DIR - name of spqr library (necessary due to error in debian package)
#  SUITESPARSE_LIBRARY_DIR      - Library main directory containing suitesparse libs
#  SUITESPARSE_LIBRARY_DIRS     - all Library directories containing suitesparse libs
#  SUITESPARSE_SPQR_VALID       - automatic identification whether or not spqr package is installed correctly

IF (SUITESPARSE_INCLUDE_DIRS)
  # Already in cache, be silent
  SET(SUITESPARSE_FIND_QUIETLY TRUE)
ENDIF (SUITESPARSE_INCLUDE_DIRS)

FIND_PATH( CHOLMOD_INCLUDE_DIR cholmod.h
  	        PATHS /usr/local/include 
        	        /usr/include 
        	        /usr/include/suitesparse/ 
         )

MESSAGE(STATUS "CHOLMOD_INCLUDE_DIR:${CHOLMOD_INCLUDE_DIR}")

FIND_LIBRARY( SUITESPARSE_LIBRARY_DIR
           NAMES cholmod
           PATHS /usr/lib 
                 /usr/lib64 
                 /usr/local/lib )
MESSAGE(STATUS "SUITESPARSE_LIBRARY_DIR:${SUITESPARSE_LIBRARY_DIR}")

# Add cholmod include directory to collection include directories
IF ( CHOLMOD_INCLUDE_DIR )
	list ( APPEND SUITESPARSE_INCLUDE_DIRS ${CHOLMOD_INCLUDE_DIR} )
ENDIF( CHOLMOD_INCLUDE_DIR )

# if we found the library, add it to the defined libraries
IF ( SUITESPARSE_LIBRARY_DIR )

    # Skipped, as this is set for apple in the block above
    list ( APPEND SUITESPARSE_LIBRARIES amd)
    list ( APPEND SUITESPARSE_LIBRARIES btf)
    list ( APPEND SUITESPARSE_LIBRARIES camd)
    list ( APPEND SUITESPARSE_LIBRARIES ccolamd)
    list ( APPEND SUITESPARSE_LIBRARIES cholmod)
    list ( APPEND SUITESPARSE_LIBRARIES colamd)
#   list ( APPEND SUITESPARSE_LIBRARIES csparse)
    list ( APPEND SUITESPARSE_LIBRARIES cxsparse)
    list ( APPEND SUITESPARSE_LIBRARIES klu)
#   list ( APPEND SUITESPARSE_LIBRARIES spqr)
    list ( APPEND SUITESPARSE_LIBRARIES umfpack)
    list ( APPEND SUITESPARSE_LIBRARIES suitesparseconfig)
   
    # Metis and spqr are optional
    FIND_LIBRARY( SUITESPARSE_METIS_LIBRARY
                  NAMES metis
                  PATHS ${SUITESPARSE_LIBRARY_DIR} )
    IF (SUITESPARSE_METIS_LIBRARY)			
	    list ( APPEND SUITESPARSE_LIBRARIES metis)
    ENDIF(SUITESPARSE_METIS_LIBRARY)

   if(EXISTS  "${CHOLMOD_INCLUDE_DIR}/SuiteSparseQR.hpp")
	  SET(SUITESPARSE_SPQR_VALID TRUE CACHE BOOL "SuiteSparseSPQR valid")
   else()
	  SET(SUITESPARSE_SPQR_VALID false CACHE BOOL "SuiteSparseSPQR valid")
   endif()

   if(SUITESPARSE_SPQR_VALID)
	  FIND_LIBRARY( SUITESPARSE_SPQR_LIBRARY
	  NAMES spqr
	  PATHS ${SUITESPARSE_LIBRARY_DIR} )
	IF (SUITESPARSE_SPQR_LIBRARY)			
	    list ( APPEND SUITESPARSE_LIBRARIES spqr)
	ENDIF (SUITESPARSE_SPQR_LIBRARY)
    endif()
       
ENDIF( SUITESPARSE_LIBRARY_DIR )  
   
IF (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)
   SET(SUITESPARSE_FOUND TRUE)
ELSE (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)
   SET( SUITESPARSE_FOUND FALSE )
ENDIF (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)

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
        	     /usr/include/suitesparse/ )

MESSAGE(STATUS "CHOLMOD_INCLUDE_DIR:${CHOLMOD_INCLUDE_DIR}")

FIND_LIBRARY( AMD_LIBRARY
              NAMES amd
              PATHS /usr/lib /usr/lib64 /usr/local/lib )

FIND_LIBRARY( CHOLMOD_LIBRARY
              NAMES cholmod
              PATHS /usr/lib /usr/lib64 /usr/local/lib )

FIND_LIBRARY( SUITESPARSECONFIG_LIBRARY
              NAMES suitesparseconfig
              PATHS /usr/lib /usr/lib64 /usr/local/lib )

MESSAGE(STATUS "CHOLMOD_LIBRARY:${CHOLMOD_LIBRARY}")

# Add cholmod include directory to collection include directories
IF ( CHOLMOD_INCLUDE_DIR )
   LIST ( APPEND SUITESPARSE_INCLUDE_DIRS ${CHOLMOD_INCLUDE_DIR} )
ENDIF( CHOLMOD_INCLUDE_DIR )

# if we found the library, add it to the defined libraries
IF ( CHOLMOD_LIBRARY AND SUITESPARSECONFIG_LIBRARY AND AMD_LIBRARY)

   # Skipped, as this is set for apple in the block above
   LIST ( APPEND SUITESPARSE_LIBRARIES amd)
#   LIST ( APPEND SUITESPARSE_LIBRARIES btf)
#   LIST ( APPEND SUITESPARSE_LIBRARIES camd)
#   LIST ( APPEND SUITESPARSE_LIBRARIES ccolamd)
   LIST ( APPEND SUITESPARSE_LIBRARIES cholmod)
   LIST ( APPEND SUITESPARSE_LIBRARIES colamd)
#  LIST ( APPEND SUITESPARSE_LIBRARIES csparse)
#  LIST ( APPEND SUITESPARSE_LIBRARIES cxsparse)
#   LIST ( APPEND SUITESPARSE_LIBRARIES klu)
#  LIST ( APPEND SUITESPARSE_LIBRARIES spqr)
   LIST ( APPEND SUITESPARSE_LIBRARIES umfpack)
   LIST ( APPEND SUITESPARSE_LIBRARIES suitesparseconfig)
   
   # Metis and spqr are optional
   FIND_LIBRARY( SUITESPARSE_METIS_LIBRARY
                 NAMES metis
                 PATHS ${SUITESPARSE_LIBRARY_DIR} )
   IF (SUITESPARSE_METIS_LIBRARY)			
      LIST ( APPEND SUITESPARSE_LIBRARIES metis)
   ENDIF(SUITESPARSE_METIS_LIBRARY)

   IF(EXISTS  "${CHOLMOD_INCLUDE_DIR}/SuiteSparseQR.hpp")
      SET(SUITESPARSE_SPQR_VALID TRUE CACHE BOOL "SuiteSparseSPQR valid")
   ELSE()
      SET(SUITESPARSE_SPQR_VALID false CACHE BOOL "SuiteSparseSPQR valid")
   ENDIF()

   IF(SUITESPARSE_SPQR_VALID)

      FIND_LIBRARY( SUITESPARSE_SPQR_LIBRARY
                    NAMES spqr
                    PATHS /usr/lib /usr/lib64 /usr/local/lib )

      IF (SUITESPARSE_SPQR_LIBRARY)			
         LIST ( APPEND SUITESPARSE_LIBRARIES spqr)
      ENDIF (SUITESPARSE_SPQR_LIBRARY)

   ENDIF()
       
ENDIF( CHOLMOD_LIBRARY AND SUITESPARSECONFIG_LIBRARY AND AMD_LIBRARY)
   
IF (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)
   SET(SUITESPARSE_FOUND TRUE)
ELSE (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)
   SET( SUITESPARSE_FOUND FALSE )
ENDIF (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)

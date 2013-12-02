SET(WSMP_ENABLED OFF CACHE BOOL "Use Watson Sparse Matrix Library ")
SET(WSMP_ROOT    $ENV{WSMP_ROOT} CACHE PATH "Watson Sparse Matrix location")

IF(WSMP_ENABLED)

   FIND_LIBRARY(WSMP_LIBRARY NAMES wsmp64
                             HINTS ${WSMP_ROOT} 
                             PATH_SUFFIXES lib
                             DOC "Watson Sparse Matrix Library"
               )

   MESSAGE(STATUS "WSMP_LIBRARY:${WSMP_LIBRARY}")

ENDIF()

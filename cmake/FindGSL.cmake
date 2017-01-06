IF (BUILD_GSL)

  INCLUDE(ExternalProject)

  EXTERNALPROJECT_ADD( gsl_project
     URL  http://mirror0.babylon.network/gnu/gsl/gsl-1.16.tar.gz
     SOURCE_DIR ${CMAKE_BINARY_DIR}/gsl
     BINARY_DIR ${CMAKE_BINARY_DIR}/gsl
     CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/gsl/configure --prefix=${CMAKE_BINARY_DIR}
     BUILD_COMMAND ${MAKE}
  )
  
  EXTERNALPROJECT_ADD( fgsl_project
     URL  http://www.lrz.de/services/software/mathematik/gsl/fortran/download/fgsl-1.0.0.tar.gz
     SOURCE_DIR ${CMAKE_BINARY_DIR}/fgsl
     BINARY_DIR ${CMAKE_BINARY_DIR}/fgsl
     CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/fgsl/configure --prefix=${CMAKE_BINARY_DIR} PKG_CONFIG_PATH=${CMAKE_BINARY_DIR}/lib/pkgconfig
     BUILD_COMMAND ${MAKE}
  )

  ADD_DEPENDENCIES(fgsl_project gsl_project)

  SET (FGSL_INCLUDES "${CMAKE_BINARY_DIR}/include/fgsl")
  SET (FGSL_LIBRARIES fgsl gsl gslcblas)

  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
  SET(FGSL_FOUND TRUE)

ELSE()

  # - Find GSL
  # Find the native GSL INCLUDEs and library
  #
  #  GSL_INCLUDES    - where to find gsl/gsl_*.h, etc.
  #  GSL_LIBRARIES   - List of libraries when using GSL.
  #  GSL_FOUND       - True if GSL found.
  
  IF (GSL_INCLUDES)
    SET (GSL_FIND_QUIETLY TRUE)
  ENDIF (GSL_INCLUDES)
  
  FIND_PATH (GSL_INCLUDES NAMES gsl/gsl_math.h HINTS ${CMAKE_BINARY_DIR} PATH_SUFFIXES include)
  
  FIND_LIBRARY (GSL_LIB NAMES gsl HINTS ${CMAKE_BINARY_DIR} PATH_SUFFIXES lib )
  
  SET (GSL_CBLAS_LIB "" CACHE FILEPATH "If your program fails to link
  (usually because GSL is not automatically linking a CBLAS and no other
  component of your project provides a CBLAS) then you may need to point
  this variable to a valid CBLAS.  Usually GSL is distributed with
  libgslcblas.{a,so} (next to GSL_LIB) which you may use if an optimized
  CBLAS is unavailable.")
  
  SET (GSL_LIBRARIES "${GSL_LIB}" "${GSL_CBLAS_LIB}")
  
  # handle the QUIETLY and REQUIRED arguments and set GSL_FOUND to TRUE if
  # all listed variables are TRUE
  INCLUDE (FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS (GSL DEFAULT_MSG GSL_LIBRARIES GSL_INCLUDES)
  
  MARK_AS_ADVANCED (GSL_LIB GSL_CBLAS_LIB GSL_INCLUDES)
  
  IF (GSL_FOUND)
  
    FIND_PATH (FGSL_INCLUDES NAMES fgsl.mod 
               HINTS /usr/local/include
               PATH_SUFFIXES fgsl 
               DOC "Path to fgsl.mod")
    
    FIND_LIBRARY (FGSL_LIB NAMES fgsl DOC "Path to libfgsl.a" )
    
    IF (FGSL_LIB) 
      SET (FGSL_LIBRARIES "${FGSL_LIB}" "${GSL_LIBRARIES}")
    ENDIF (FGSL_LIB) 
    
    INCLUDE (FindPackageHandleStandardArgs)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS (FGSL DEFAULT_MSG FGSL_LIBRARIES FGSL_INCLUDES)
    
    MARK_AS_ADVANCED (FGSL_LIB FGSL_INCLUDES)

  
  ENDIF (GSL_FOUND)

ENDIF (BUILD_GSL)


IF(FGSL_FOUND)
  MESSAGE(STATUS "FGSL_INCLUDES:${FGSL_INCLUDES}")
  MESSAGE(STATUS "FGSL_LIB:${FGSL_LIB}")
ENDIF(FGSL_FOUND)

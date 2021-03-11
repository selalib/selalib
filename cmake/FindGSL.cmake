if(BUILD_GSL)

  include(ExternalProject)

  ExternalProject_Add(
    gsl_project
    URL http://mirror0.babylon.network/gnu/gsl/gsl-1.16.tar.gz
    SOURCE_DIR ${CMAKE_BINARY_DIR}/gsl
    BINARY_DIR ${CMAKE_BINARY_DIR}/gsl
    CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/gsl/configure
                      --prefix=${CMAKE_BINARY_DIR}
    BUILD_COMMAND ${MAKE})

  ExternalProject_Add(
    fgsl_project
    URL http://www.lrz.de/services/software/mathematik/gsl/fortran/download/fgsl-1.0.0.tar.gz
    SOURCE_DIR ${CMAKE_BINARY_DIR}/fgsl
    BINARY_DIR ${CMAKE_BINARY_DIR}/fgsl
    CONFIGURE_COMMAND
      ${CMAKE_BINARY_DIR}/fgsl/configure --prefix=${CMAKE_BINARY_DIR}
      PKG_CONFIG_PATH=${CMAKE_BINARY_DIR}/lib/pkgconfig
    BUILD_COMMAND ${MAKE})

  add_dependencies(fgsl_project gsl_project)

  set(FGSL_INCLUDES "${CMAKE_BINARY_DIR}/include/fgsl")
  set(FGSL_LIBRARIES fgsl gsl gslcblas)

  link_directories(${CMAKE_BINARY_DIR}/lib)
  set(FGSL_FOUND TRUE)

else()

  # * Find GSL Find the native GSL INCLUDEs and library
  #
  # GSL_INCLUDES    - where to find gsl/gsl_*.h, etc. GSL_LIBRARIES   - List of
  # libraries when using GSL. GSL_FOUND       - True if GSL found.

  if(GSL_INCLUDES)
    set(GSL_FIND_QUIETLY TRUE)
  endif(GSL_INCLUDES)

  find_path(
    GSL_INCLUDES
    NAMES gsl/gsl_math.h
    HINTS ${CMAKE_BINARY_DIR}
    PATH_SUFFIXES include)

  find_library(
    GSL_LIB
    NAMES gsl
    HINTS ${CMAKE_BINARY_DIR}
    PATH_SUFFIXES lib)

  set(GSL_CBLAS_LIB
      ""
      CACHE
        FILEPATH
        "If your program fails to link
  (usually because GSL is not automatically linking a CBLAS and no other
  component of your project provides a CBLAS) then you may need to point
  this variable to a valid CBLAS.  Usually GSL is distributed with
  libgslcblas.{a,so} (next to GSL_LIB) which you may use if an optimized
  CBLAS is unavailable.")

  set(GSL_LIBRARIES "${GSL_LIB}" "${GSL_CBLAS_LIB}")

  # handle the QUIETLY and REQUIRED arguments and set GSL_FOUND to TRUE if all
  # listed variables are TRUE
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(GSL DEFAULT_MSG GSL_LIBRARIES GSL_INCLUDES)

  mark_as_advanced(GSL_LIB GSL_CBLAS_LIB GSL_INCLUDES)

  if(GSL_FOUND)

    find_path(
      FGSL_INCLUDES
      NAMES fgsl.mod
      HINTS /usr/local/include
      PATH_SUFFIXES fgsl
      DOC "Path to fgsl.mod")

    find_library(
      FGSL_LIB
      NAMES fgsl
      DOC "Path to libfgsl.a")

    if(FGSL_LIB)
      set(FGSL_LIBRARIES "${FGSL_LIB}" "${GSL_LIBRARIES}")
    endif(FGSL_LIB)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(FGSL DEFAULT_MSG FGSL_LIBRARIES
                                      FGSL_INCLUDES)

    mark_as_advanced(FGSL_LIB FGSL_INCLUDES)

  endif(GSL_FOUND)

endif(BUILD_GSL)

if(FGSL_FOUND)
  message(STATUS "FGSL_INCLUDES:${FGSL_INCLUDES}")
  message(STATUS "FGSL_LIB:${FGSL_LIB}")
endif(FGSL_FOUND)

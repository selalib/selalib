# FFTW_INCLUDE_DIR = fftw3.f03 FFTW_LIBRARIES = libfftw3.a FFTW_FOUND = true if
# FFTW3 is found

if(DEFINED ENV{FFTW_ROOT})
  set(FFTW_ROOT
      $ENV{FFTW_ROOT}
      CACHE PATH "FFTW location")
else()
  set(FFTW_ROOT
      "/usr/local"
      CACHE PATH "FFTW location")
endif()

set(TRIAL_PATHS
    ${FFTW_ROOT}
    $ENV{FFTW_HOME}
    $ENV{FFTW_DIR}
    $ENV{FFTW_BASE}
    $ENV{FFTW_ROOT_DIR}
    /usr
    /usr/local
    /usr/lib64/mpich2
    /usr/lib64/openmpi
    /opt/local)

set(FFTW_F2003
    ON
    CACHE BOOL "Use FFTW Fortran 2003 interface")

if(FFTW_F2003)
  find_path(
    FFTW_INCLUDE_DIRS
    NAMES fftw3.f03
    HINTS ${TRIAL_PATHS} $ENV{FFTW_INCLUDE}
    PATH_SUFFIXES include
    DOC "path to fftw3.f03")
  if(FFTW_INCLUDE_DIRS)
    add_definitions(-DFFTW_F2003)
  else()
    message(
      WARNING
        "Could not find FFTW F2003 header file, falling back to F77 interface..."
    )
    find_path(
      FFTW_INCLUDE_DIRS
      NAMES fftw3.f
      HINTS ${TRIAL_PATHS} $ENV{FFTW_INCLUDE}
      PATH_SUFFIXES include
      DOC "path to fftw3.f")
    set(FFTW_F2003
        OFF
        CACHE BOOL "Use FFTW Fortran 2003 interface" FORCE)
    remove_definitions(-DFFTW_F2003)
  endif(FFTW_INCLUDE_DIRS)
else()
  remove_definitions(-DFFTW_F2003)
endif(FFTW_F2003)

find_library(
  FFTW_LIBRARY
  NAMES fftw3
  HINTS ${TRIAL_PATHS} $ENV{FFTW_LIB}
  PATH_SUFFIXES lib lib64)

set(FFTW_LIBRARIES ${FFTW_LIBRARY})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_INCLUDE_DIRS
                                  FFTW_LIBRARIES)
if(USE_MKL AND NOT FFTW_FOUND)

  if(FFTW_F2003)
    message(
      WARNING
        "Intel MKL wrappers to FFTW in use. F2003 interface not available, falling back to FFTW F77 interface..."
    )
  endif()
  find_path(
    FFTW_INCLUDE_DIRS
    NAMES fftw3.f
    HINTS $ENV{MKLROOT}/include
    PATH_SUFFIXES fftw)
  set(FFTW_LIBRARIES ${LAPACK_LIBRARIES})
  set(FFTW_F2003
      OFF
      CACHE BOOL "Use FFTW Fortran 2003 interface" FORCE)
  remove_definitions(-DFFTW_F2003)

  find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_INCLUDE_DIRS
                                    FFTW_LIBRARIES)

endif()

if(FFTW_FOUND)

  message(STATUS "FFTW_INCLUDE_DIRS:${FFTW_INCLUDE_DIRS}")
  message(STATUS "FFTW_LIBRARIES:${FFTW_LIBRARIES}")
  mark_as_advanced(FFTW_INCLUDE_DIRS FFTW_LIBRARIES)

endif(FFTW_FOUND)

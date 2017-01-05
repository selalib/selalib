# FFTW_INCLUDE_DIR = fftw3.f03
# FFTW_LIBRARIES = libfftw3.a
# FFTW_FOUND = true if FFTW3 is found

IF(DEFINED ENV{FFTW_ROOT})
   SET(FFTW_ROOT $ENV{FFTW_ROOT} CACHE PATH "FFTW location")
ELSE()
   SET(FFTW_ROOT "/usr/local" CACHE PATH "FFTW location")
ENDIF()

SET(TRIAL_PATHS ${FFTW_ROOT}
                $ENV{FFTW_HOME}
                $ENV{FFTW_DIR}
                $ENV{FFTW_BASE}
                $ENV{FFTW_ROOT_DIR}
                /usr
                /usr/local
                /usr/lib64/mpich2
                /usr/lib64/openmpi
                /opt/local
 )


SET(FFTW_F2003 ON CACHE BOOL "Use FFTW Fortran 2003 interface")

IF(FFTW_F2003)
  FIND_PATH(FFTW_INCLUDE_DIRS NAMES fftw3.f03 
                              HINTS ${TRIAL_PATHS} $ENV{FFTW_INCLUDE}
                              PATH_SUFFIXES include DOC "path to fftw3.f03")
  IF(FFTW_INCLUDE_DIRS)
    ADD_DEFINITIONS(-DFFTW_F2003)
  ELSE()
    MESSAGE(WARNING "Could not find FFTW F2003 header file, falling back to F77 interface...")
    FIND_PATH(FFTW_INCLUDE_DIRS NAMES fftw3.f
                               HINTS ${TRIAL_PATHS} $ENV{FFTW_INCLUDE}
                               PATH_SUFFIXES include DOC "path to fftw3.f")
    SET(FFTW_F2003 OFF CACHE BOOL "Use FFTW Fortran 2003 interface" FORCE)
    REMOVE_DEFINITIONS(-DFFTW_F2003)
  ENDIF(FFTW_INCLUDE_DIRS)
ELSE()
   REMOVE_DEFINITIONS(-DFFTW_F2003)
ENDIF(FFTW_F2003)

FIND_LIBRARY(FFTW_LIBRARY NAMES fftw3 
                          HINTS ${TRIAL_PATHS} $ENV{FFTW_LIB}
                          PATH_SUFFIXES lib lib64)

SET(FFTW_LIBRARIES ${FFTW_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW 
                                  DEFAULT_MSG 
                                  FFTW_INCLUDE_DIRS FFTW_LIBRARIES)
IF(USE_MKL AND NOT FFTW_FOUND)

   IF(FFTW_F2003)
      MESSAGE(WARNING "Intel MKL wrappers to FFTW in use. F2003 interface not available, falling back to FFTW F77 interface...")
   ENDIF()
   FIND_PATH(FFTW_INCLUDE_DIRS NAMES fftw3.f 
                              HINTS $ENV{MKLROOT}/include
                              PATH_SUFFIXES fftw)
   SET(FFTW_LIBRARIES ${LAPACK_LIBRARIES})
   SET(FFTW_F2003 OFF CACHE BOOL "Use FFTW Fortran 2003 interface" FORCE)
   REMOVE_DEFINITIONS(-DFFTW_F2003)

   FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW 
                                     DEFAULT_MSG 
                                     FFTW_INCLUDE_DIRS FFTW_LIBRARIES)

ENDIF()


IF(FFTW_FOUND)

   # this include_directories is commented because we use 
   # target_include_directories in fft interface (PN)
   #INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIRS})
   MESSAGE(STATUS "FFTW_INCLUDE_DIRS:${FFTW_INCLUDE_DIRS}")
   MESSAGE(STATUS "FFTW_LIBRARIES:${FFTW_LIBRARIES}")
   MARK_AS_ADVANCED( FFTW_INCLUDE_DIRS FFTW_LIBRARIES)

ENDIF(FFTW_FOUND)

# FFTW_INCLUDE_DIR = fftw3.f03
# FFTW_LIBRARIES = libfftw3.a
# FFTW_FOUND = true if FFTW3 is found

SET(TRIAL_PATHS 
                $ENV{FFTW_HOME}
                $ENV{FFTW_ROOT}
                $ENV{FFTW_DIR}
                /usr
                /usr/local
                /usr/lib64/mpich2
                /usr/lib64/openmpi
                /opt/local
 )

FIND_PATH(FFTW_INCLUDE_DIRS NAMES fftw3.f03 
                            HINTS ${TRIAL_PATHS} $ENV{FFTW_INCLUDE}
                            PATH_SUFFIXES include DOC "path to fftw3.f03")

IF(FFTW_INCLUDE_DIRS)
   SET(FFTW_F2003 ON CACHE BOOL "FFTW Fortran 2003 header file available")
   ADD_DEFINITIONS(-DFFTW_F2003)
ELSE()
   FIND_PATH(FFTW_INCLUDE_DIRS NAMES fftw3.f
                               HINTS ${TRIAL_PATHS} $ENV{FFTW_INCLUDE}
                               PATH_SUFFIXES include DOC "path to fftw3.f")
ENDIF(FFTW_INCLUDE_DIRS)

IF( NOT FFTW_F2003) 
   REMOVE_DEFINITIONS(-DFFTW_F2003)
ENDIF()

#IF(FFTW_MPI_INCLUDE_DIR)
#   SET(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIRS} ${FFTW_MPI_INCLUDE_DIR})
#ENDIF(FFTW_MPI_INCLUDE_DIR)
#
FIND_LIBRARY(FFTW_LIBRARY NAMES fftw3 
                          HINTS ${TRIAL_PATHS} $ENV{FFTW_LIB}
                          PATH_SUFFIXES lib lib64)

#FIND_LIBRARY(FFTW_THREADS_LIBRARY NAMES fftw3_threads HINTS ${TRIAL_PATHS} PATH_SUFFIXES lib lib64)

IF(FFTW_LIBRARY)
   SET(FFTW_LIBRARIES ${FFTW_LIBRARY})
ENDIF()

#IF(FFTW_THREADS_LIBRARY)
#   SET(FFTW_LIBRARIES ${FFTW_THREADS_LIBRARY} ${FFTW_LIBRARY})
#ELSE()
#   MESSAGE(STATUS "No threaded fftw3 installation")
#ENDIF()

#FIND_PATH(FFTW_MPI_INCLUDE_DIR NAMES fftw3-mpi.f03 HINTS ${TRIAL_PATHS} PATH_SUFFIXES include DOC "path to fftw3-mpi.f03")
#FIND_LIBRARY(FFTW_MPI_LIBRARY NAMES fftw3_mpi HINTS ${TRIAL_PATHS} PATH_SUFFIXES lib lib64)
#IF(FFTW_MPI_LIBRARY)
#   SET(FFTW_LIBRARIES ${FFTW_MPI_LIBRARY} ${FFTW_LIBRARY})
#ELSE()
#   MESSAGE(STATUS "No mpi fftw3 installation")
#ENDIF()

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW DEFAULT_MSG
                                  FFTW_INCLUDE_DIRS FFTW_LIBRARIES)

MARK_AS_ADVANCED( FFTW_INCLUDE_DIRS
                  FFTW_LIBRARIES)

IF(FFTW_INCLUDE_DIRS)
   INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIRS})
ENDIF(FFTW_INCLUDE_DIRS)
MESSAGE(STATUS "FFTW_INCLUDE_DIRS:${FFTW_INCLUDE_DIRS}")
MESSAGE(STATUS "FFTW_LIBRARIES:${FFTW_LIBRARIES}")

IF($ENV{MKLROOT})

   FIND_PATH(FFTW_MKL_INCLUDE_DIR NAMES fftw3.f 
                                  HINTS $ENV{MKLROOT}/include
                                  PATH_SUFFIXES fftw)

   MESSAGE(STATUS "MKL FFTW FOUND (NOT LINKED)")
   MESSAGE(STATUS "FFTW_MKL_INCLUDE_DIR:${FFTW_MKL_INCLUDE_DIR}")

ENDIF($ENV{MKLROOT})

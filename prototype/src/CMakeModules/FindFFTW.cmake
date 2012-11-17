# FFTW_INCLUDE_DIR = fftw3.f03
# FFTW_LIBRARIES = libfftw3.a
# FFTW_FOUND = true if FFTW3 is found

SET(TRIAL_PATHS 
                $ENV{FFTW_ROOT}
                /usr
                /usr/local
                /usr/lib64/mpich2
                /usr/lib64/openmpi
                /opt/local
 )

FIND_PATH(FFTW_INCLUDE_DIRS NAMES fftw3.f03 HINTS ${TRIAL_PATHS} PATH_SUFFIXES include DOC "path tp fftw3.f03")
FIND_PATH(FFTW_MPI_INCLUDE_DIR NAMES fftw3-mpi.f03 HINTS ${TRIAL_PATHS} PATH_SUFFIXES include DOC "path to fftw3-mpi.f03")
IF(FFTW_MPI_INCLUDE_DIR)
   SET(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIRS} ${FFTW_MPI_INCLUDE_DIR})
ENDIF(FFTW_MPI_INCLUDE_DIR)

FIND_LIBRARY(FFTW_LIBRARY NAMES fftw3 HINTS ${TRIAL_PATHS} PATH_SUFFIXES lib lib64)
FIND_LIBRARY(FFTW_THREADS_LIBRARY NAMES fftw3_threads HINTS ${TRIAL_PATHS} PATH_SUFFIXES lib lib64)
FIND_LIBRARY(FFTW_MPI_LIBRARY NAMES fftw3_mpi HINTS ${TRIAL_PATHS} PATH_SUFFIXES lib lib64)

IF(FFTW_LIBRARY)
   SET(FFTW_LIBRARIES ${FFTW_LIBRARY})
ELSE()
   MESSAGE(SEND_ERROR "No fftw3 installation")
ENDIF()

IF(FFTW_THREADS_LIBRARY)
   SET(FFTW_LIBRARIES ${FFTW_THREADS_LIBRARY} ${FFTW_LIBRARY})
ELSE()
   MESSAGE(STATUS "No threaded fftw3 installation")
ENDIF()

IF(FFTW_MPI_LIBRARY)
   SET(FFTW_LIBRARIES ${FFTW_MPI_LIBRARY} ${FFTW_THREADS_LIBRARY} ${FFTW_LIBRARY})
ELSE()
   MESSAGE(STATUS "No mpi fftw3 installation")
ENDIF()

SET(FFTW_FOUND FALSE)

IF($ENV{HOSTNAME} MATCHES "hpc-f0*")
   SET(FFTW_INCLUDE_DIRS "/home/math/navaro/local/include")
   SET(FFTW_LIBRARIES "/home/math/navaro/local/lib/libfftw3.a" 
                      "/home/math/navaro/local/lib/libfftw3_threads.a" 
                      "/home/math/navaro/local/lib/libfftw3_mpi.a")
ENDIF()

IF(FFTW_INCLUDE_DIRS AND FFTW_LIBRARIES)
   MESSAGE(STATUS "FFTW_INCLUDE_DIRS=${FFTW_INCLUDE_DIRS}")
   MESSAGE(STATUS "FFTW_LIBRARIES=${FFTW_LIBRARIES}")
   SET(FFTW_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED( FFTW_INCLUDE_DIRS
                  FFTW_MPI_INCLUDE_DIR
                  FFTW_LIBRARIES
                  FFTW_FOUND              )

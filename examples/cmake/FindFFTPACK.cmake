SET(TRIAL_PATHS $ENV{SELALIB_ROOT}/usr
                /usr/local)

FIND_LIBRARY(FFTPACK_LIBRARIES NAMES dfftpack
                               HINTS ${TRIAL_PATHS}   
		                         PATH_SUFFIXES lib 
		                         DOC "PATH TO libdfftpack")

IF (FFTPACK_LIBRARIES)

   MESSAGE(STATUS "FFTPACK FOUND")
   SET(FFTPACK_FOUND TRUE)

ELSE()

   MESSAGE(STATUS "FFTPACK NOT FOUND")
   SET(FFTPACK_FOUND FALSE)

ENDIF(FFTPACK_LIBRARIES)



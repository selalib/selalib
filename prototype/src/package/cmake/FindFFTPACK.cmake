SET(TRIAL_PATHS /usr /usr/local /opt/local)

FIND_LIBRARY(FFTPACK_LIBRARIES NAMES dfftpack
                               HINTS ${TRIAL_PATHS}   
                               PATH_SUFFIXES lib 
                               DOC "PATH TO libdfftpack")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTPACK DEFAULT_MSG
                                  FFTPACK_INCLUDE_DIRS FFTPACK_LIBRARIES)


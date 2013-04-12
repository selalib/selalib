
##########################################################
# add the cache entry FFT_DEFAULT_LIBRARY to define the default library use by sll_fft
SET(FFT_DEFAULT_LIBRARY SLLFFT CACHE STRING "specify the library use in sll_fft, options are : FFTPACK, FFTW or SLLFFT")

IF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTPACK")
  add_definitions(-DFFTPACK)
ELSEIF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTW")
  add_definitions(-DFFTW)
ELSE()
  add_definitions(-DSLLFFT)
ENDIF()

# add the cache entry FFTPACK_ENABLED for enable/disable fftpack
SET(FFTPACK_ENABLED OFF CACHE BOOL " ")

# recherche le package FFTW
SET(FFTW_ENABLED ON CACHE BOOL " ")

IF(FFTW_ENABLED)
   IF( DEFINED ENV{FFTW_ROOT} )
      SET(FFTW_ROOT $ENV{FFTW_ROOT})
   ENDIF()
   FIND_PACKAGE(FFTW QUIET)
ENDIF(FFTW_ENABLED)

IF(FFTW_FOUND)
   MESSAGE(STATUS "FFTW FOUND")
   include_directories(${FFTW_INCLUDE_DIRS})
   SET(FFT_ADD_MODULE dfftpack fftw_module)
ELSE(FFTW_FOUND)
   MESSAGE(STATUS "FFTW NOT FOUND")
   SET(FFTW_ENABLED NO CACHE BOOL " " FORCE)
   add_definitions(-D_NOFFTW)
   SET(FFT_ADD_MODULE dfftpack)
ENDIF(FFTW_FOUND)

IF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTPACK" AND NOT FFTPACK_ENABLED)
   MESSAGE("Please put on FFTPACK_ENABLED to use fftpack")
ENDIF()

IF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTW" AND NOT FFTW_FOUND)
   MESSAGE("You can't use fftw library because it's not installed")
ELSEIF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTW" AND NOT FFTW_ENABLED)
   MESSAGE("Please put on FFTW_ENABLED to use fftw")
ENDIF()

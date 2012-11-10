
##########################################################
# add the cache entry FFT_DEFAULT_LIBRARY to define the default library use by sll_fft
set(FFT_DEFAULT_LIBRARY SLLFFT CACHE STRING "specify the library use in sll_fft")


# add the cache entry FFTPACK_ENABLED for enable/disable fftpack
set(FFTPACK_ENABLED OFF CACHE BOOL " ")

# recherche le package FFTW
set(FFTW_ENABLED ON CACHE BOOL " ")

IF(FFTW_ENABLED)
   IF( DEFINED ENV{FFTW_ROOT} )
      set(FFTW_ROOT $ENV{FFTW_ROOT})
   ENDIF()
   FIND_PACKAGE(FFTW QUIET)
ENDIF(FFTW_ENABLED)

IF(FFTW_FOUND)
   message(STATUS "FFTW FOUND")
   include_directories(${FFTW_INCLUDE_DIRS})
   set(FFT_ADD_MODULE fftpack_module fftw_module)
ELSE(FFTW_FOUND)
   message(STATUS "FFTW NOT FOUND")
   set(FFTW_ENABLED NO CACHE BOOL " " FORCE)
   add_definitions(-D_NOFFTW)
   set(FFT_ADD_MODULE fftpack_module)
ENDIF(FFTW_FOUND)

IF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTPACK" AND NOT FFTPACK_ENABLED)
	message("Please put on FFTPACK_ENABLED to use fftpack")

ENDIF()
IF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTW" AND NOT FFTW_FOUND)
	message("You can't use fftw library because it's not installed")
ELSEIF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTW" AND NOT FFTW_ENABLED)
	message("Please put on FFTW_ENABLED to use fftw")
ENDIF()

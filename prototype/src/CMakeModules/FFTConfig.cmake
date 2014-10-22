
SET(FFT_LIB SLLFFT CACHE STRING "fft library, options are : FFTPACK, FFTW or SLLFFT")

IF( DEFINED ENV{FFTW_ROOT} )

   SET(FFTW_ROOT $ENV{FFTW_ROOT})

ENDIF()

FIND_PACKAGE(FFTW)

IF(${FFT_LIB} STREQUAL "FFTPACK")

  ADD_DEFINITIONS(-DFFTPACK)

ELSEIF(${FFT_LIB} STREQUAL "FFTW")

  IF(FFTW_FOUND)
     ADD_DEFINITIONS(-DFFTW)
  ELSE(FFTW_FOUND)
     MESSAGE(SEND_ERROR "FFTW NOT FOUND, try set FFTW_ROOT or change FFT_LIB ")
  ENDIF(FFTW_FOUND)

  SET(FFTW_ENABLED ON)

ELSE()

  ADD_DEFINITIONS(-DSLLFFT)

ENDIF()

MESSAGE(STATUS "FFT_LIB=${FFT_LIB}")
IF(FFTW_FOUND)
   MESSAGE(STATUS "Set FFT_LIB=FFTW to use fftw library")
ENDIF()

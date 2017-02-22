# Extract environmental variable FFTW_ROOT
IF( DEFINED ENV{FFTW_ROOT} )
   SET(FFTW_ROOT $ENV{FFTW_ROOT})
ENDIF()

# Set default library to FFTW if available, SLLFFT otherwise
IF (FFTW_ENABLED)
  FIND_PACKAGE(FFTW)
  IF(FFTW_FOUND)
    SET(FFT_DEFAULT FFTW)
  ELSE()
    SET(FFT_DEFAULT SLLFFT)
  ENDIF()
ENDIF()

# Create cache variable FFT_LIB and give it default value
SET(FFT_LIB ${FFT_DEFAULT} CACHE STRING "fft library, options are : FFTW, SLLFFT or FFTPACK")

# Perform checks and add definitions whenever FFT_LIB is changed by user
IF(${FFT_LIB} STREQUAL "FFTPACK")

  MESSAGE(FATAL_ERROR "The fftpack interface is deprecated, sorry")
  ADD_DEFINITIONS(-DFFTPACK)

ELSEIF(${FFT_LIB} STREQUAL "FFTW")

  IF(FFTW_FOUND)
     ADD_DEFINITIONS(-DFFTW)
  ELSE(FFTW_FOUND)
     MESSAGE(SEND_ERROR "FFTW NOT FOUND, try set FFTW_ROOT or change FFT_LIB ")
  ENDIF(FFTW_FOUND)
  SET(FFTW_ENABLED ON)

ELSEIF(${FFT_LIB} STREQUAL "SLLFFT")

  MESSAGE(WARNING "Careful: SLLFFT only works on arrays with power of 2 number of elements")
  ADD_DEFINITIONS(-DSLLFFT)

ELSE()

  MESSAGE(FATAL_ERROR "Unrecognized option for FFT_LIB" )

ENDIF()

MESSAGE(STATUS "FFT_LIB=${FFT_LIB}")

MARK_AS_ADVANCED(FFTW_LIBRARY)

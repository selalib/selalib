# Set default library to FFTW if available, SLLFFT otherwise
if(FFTW_ENABLED)
  find_package(FFTW)
  if(FFTW_FOUND)
    set(FFT_DEFAULT FFTW)
  else()
    set(FFT_DEFAULT SLLFFT)
  endif()
endif()

# Create cache variable FFT_LIB and give it default value
set(FFT_LIB
    ${FFT_DEFAULT}
    CACHE STRING "fft library, options are : FFTW, SLLFFT or FFTPACK")

# Perform checks and add definitions whenever FFT_LIB is changed by user
if(${FFT_LIB} STREQUAL "FFTPACK")

  message(FATAL_ERROR "The fftpack interface is deprecated, sorry")
  add_definitions(-DFFTPACK)

elseif(${FFT_LIB} STREQUAL "FFTW")

  if(FFTW_FOUND)
    add_definitions(-DFFTW)
  else(FFTW_FOUND)
    message(SEND_ERROR "FFTW NOT FOUND, try set FFTW_ROOT or change FFT_LIB ")
  endif(FFTW_FOUND)
  set(FFTW_ENABLED ON)

elseif(${FFT_LIB} STREQUAL "SLLFFT")

  message(
    WARNING
      "Careful: SLLFFT only works on arrays with power of 2 number of elements")
  add_definitions(-DSLLFFT)

else()

  message(FATAL_ERROR "Unrecognized option for FFT_LIB")

endif()

message(STATUS "FFT_LIB=${FFT_LIB}")

mark_as_advanced(FFTW_LIBRARY)

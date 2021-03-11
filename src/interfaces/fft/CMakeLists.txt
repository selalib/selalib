if(FFT_LIB STREQUAL "FFTPACK")

  set(FFT_SRC sll_m_fft_fftpack.F90)
  set(DEPENDENCIES dfftpack)
  add_definitions(-D_DEFAULTFFTLIB=FFTPACK_MOD)

elseif(FFT_LIB STREQUAL "FFTW")

  if(FFTW_FOUND AND FFTW_ENABLED)
    set(FFT_SRC sll_m_fft_fftw3.F90)
    add_definitions(-D_DEFAULTFFTLIB=FFTW_MOD)
  endif()

else()

  if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    set(FFT_SRC sll_m_fft_sllfft_pgi.F90)
  else(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    set(FFT_SRC sll_m_fft_sllfft.F90)
  endif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)

  add_definitions(-D_DEFAULTFFTLIB=SLLFFT_MOD)

endif()

add_library(sll_fft STATIC ${FFT_SRC})

if(FFTW_FOUND AND FFTW_ENABLED)
  add_library(sll_fftw3 STATIC sll_m_fftw3.F90)
  target_include_directories(sll_fftw3 PRIVATE ${FFTW_INCLUDE_DIRS})
  list(APPEND DEPENDENCIES sll_fftw3 ${FFTW_LIBRARIES})
endif(FFTW_FOUND AND FFTW_ENABLED)

target_link_libraries(sll_fft sll_utilities sll_memory sll_assert
                      ${DEPENDENCIES})
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)

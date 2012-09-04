set(DEPENDENCIES
	sll_utilities
	sll_timer
	sll_constants
	sll_memory
	sll_assertion
	sll_working_precision
)

IF(FFTW_FOUND AND FFTW_ENABLED)
  add_library(fftw_module STATIC fftw_module.F90)
  LIST(APPEND DEPENDENCIES fftw_module)
  target_link_libraries(fftw_module ${FFTW_LIBRARIES})
  include_directories(${FFTW_INCLUDE_DIRS})
ENDIF()

IF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTPACK")
  set(FFT_SRC sll_fft_mod_through_fftpack.F90)
  set(DEPENDENCIES fftpack_module ${DEPENDENCIES})
  add_definitions(-D_DEFAULTFFTLIB=FFTPACK_MOD)
ELSEIF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTW")
  set(FFT_SRC sll_fft_mod_through_fftw.F90)
  add_definitions(-D_DEFAULTFFTLIB=FFTW_MOD)
ELSE()
  set(FFT_SRC sll_fft_selalib_implementation.F90)
  add_definitions(-D_DEFAULTFFTLIB=SLLFFT_MOD)
ENDIF()

add_library(sll_fft STATIC
	${FFT_SRC}
	fft_utils.F90

)
add_dependencies(sll_fft ${DEPENDENCIES})

IF(FFTW_FOUND AND FFTW_ENABLED)
	target_link_libraries(sll_fft fftw_module)
ENDIF()
IF(FFTPACK_ENABLED)
	target_link_libraries(sll_fft fftpack_module)
ENDIF()

add_executable(test_fft unit_test.F90)
target_link_libraries(test_fft
	sll_fft
	${DEPENDENCIES}
)


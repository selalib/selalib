add_definitions(-D_DEFAULTFFTLIB=FFTPACK_MOD)

IF(FFTW_FOUND AND FFTW_ENABLED)
  add_library(fftw_module STATIC fftw_module.F90)
  target_link_libraries(fftw_module ${FFTW_LIBRARIES})
  include_directories(${FFTW_INCLUDE_DIRS})
ELSE()
	add_definitions(-D_NOFFTW)
ENDIF()

IF(NOT FFTPACK_ENABLED)
	add_definitions(-D_NOFFTPACK)
ENDIF()

add_library(sll_fft STATIC sll_fft.F90)
IF(FFTW_FOUND AND FFTW_ENABLED)
	target_link_libraries(sll_fft fftw_module)
ENDIF()
IF(FFTPACK_ENABLED)
	target_link_libraries(sll_fft fftpack_module)
ENDIF()

add_executable(test_fft unit_test.F90)
target_link_libraries(test_fft
	sll_fft
	sll_utilities
	sll_timer
	sll_constants
	sll_memory
	sll_assertion
)


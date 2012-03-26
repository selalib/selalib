
IF(FFTW_FOUND)
  add_library(fftw_module STATIC fftw_module.F90)
  target_link_libraries(fftw_module ${FFTW_LIBRARIES})
  set(SRC_FILE sll_fft.F90 fftpack_lib.F90 fftw_lib.F90 my_lib.F90)
ELSE()
  set(SRC_FILE sll_fft.F90 fftpack_lib.F90 my_lib.F90)
ENDIF()

set(SRC_FILE sll_fft.F90)
add_library(sll_fft STATIC ${SRC_FILE})
target_link_libraries(sll_fft ${FFT_ADD_MODULE})
add_executable(test_fft unit_test.F90)
target_link_libraries(test_fft sll_fft sll_utilities sll_timer sll_constants sll_memory sll_assertion)


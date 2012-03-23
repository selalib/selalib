
IF(FFTW_FOUND)
  add_library(fftw_module STATIC fftw_module.F90)
  target_link_libraries(fftw_module ${FFTW_LIBRARIES})
ENDIF()

add_library(sll_fft STATIC sll_fft.F90)
target_link_libraries(sll_fft ${FFT_ADD_MODULE})
add_executable(test_fft unit_test.F90)
target_link_libraries(test_fft sll_fft sll_constants sll_memory
sll_working_precision sll_utilities sll_assertion sll_timer)


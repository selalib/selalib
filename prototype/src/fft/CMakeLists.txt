set(FFT_SRC sll_fft.F90)
add_library(sll_fft STATIC sll_fft.F90)
target_link_libraries( sll_fft )
add_executable(test_fft unit_test.F90)
target_link_libraries(test_fft sll_fft sll_constants sll_memory
sll_working_precision sll_utilities sll_assert sll_timer)

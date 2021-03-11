add_executable(test_fft test_fft.F90)
target_link_libraries(test_fft sll_fft)
add_test(NAME fft COMMAND test_fft)
set_tests_properties(fft PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_executable(test_fft_timings test_fft_timings.F90)
target_link_libraries(test_fft_timings sll_fft sll_timer)

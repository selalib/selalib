ADD_EXECUTABLE( test_fft test_fft.F90 )
TARGET_LINK_LIBRARIES( test_fft sll_fft )
ADD_TEST( NAME fft COMMAND test_fft )
SET_TESTS_PROPERTIES( fft PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_fft_timings test_fft_timings.F90 )
TARGET_LINK_LIBRARIES( test_fft_timings sll_fft sll_timer)

ADD_EXECUTABLE( test_fft test_fft.F90 )
TARGET_LINK_LIBRARIES( test_fft sll_fft )
ADD_TEST( NAME fft COMMAND test_fft )

ADD_EXECUTABLE( test_fft_timings test_fft_timings.F90 )
TARGET_LINK_LIBRARIES( test_fft_timings sll_fft sll_timer)

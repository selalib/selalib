ADD_EXECUTABLE( test_fft test_fft.F90 )
TARGET_LINK_LIBRARIES( test_fft sll_fft )
ADD_TEST( NAME fft COMMAND test_fft )

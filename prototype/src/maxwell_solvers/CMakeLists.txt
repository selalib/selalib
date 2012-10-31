add_library(sll_maxwell_2d STATIC sll_maxwell_2d_fdtd.F90)
add_dependencies(sll_maxwell_2d sll_assertion sll_constants sll_memory 
                                sll_collective sll_utilities sll_fft)

ADD_EXECUTABLE(test_maxwell_2d_fdtd unit_test.F90)

TARGET_LINK_LIBRARIES(test_maxwell_2d_fdtd sll_assertion 
                                           sll_constants 
                                           sll_memory 
                                           sll_collective 
                                           sll_utilities 
                                           sll_fft 
                                           sll_maxwell_2d)

ADD_EXECUTABLE(test_maxwell_2d_pstd sll_maxwell_2d_pstd.F90)

TARGET_LINK_LIBRARIES(test_maxwell_2d_pstd sll_assertion 
                                           sll_constants 
                                           sll_memory 
                                           sll_collective 
                                           sll_utilities 
                                           ${FFTW_LIBRARIES} -lfftw3_threads )

ADD_LIBRARY(sll_maxwell_2d STATIC sll_maxwell_2d.F90
                                  sll_maxwell_2d_fdtd.F90
                                  sll_maxwell_2d_pstd.F90)

ADD_DEPENDENCIES(sll_maxwell_2d sll_assertion sll_constants sll_memory 
                                sll_collective sll_utilities sll_fft)

ADD_EXECUTABLE(test_maxwell_2d_fdtd unit_test_2d_fdtd.F90)

TARGET_LINK_LIBRARIES(test_maxwell_2d_fdtd sll_maxwell_2d
                                           sll_assertion 
                                           sll_constants 
                                           sll_memory 
                                           sll_collective 
                                           sll_utilities 
                                           sll_fft )

ADD_EXECUTABLE(test_maxwell_2d_pstd unit_test_2d_pstd.F90 )

TARGET_LINK_LIBRARIES(test_maxwell_2d_pstd sll_maxwell_2d
                                           sll_assertion 
                                           sll_constants 
                                           sll_memory 
                                           sll_collective 
                                           sll_utilities 
                                           sll_fft)

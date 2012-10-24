add_library(sll_maxwell_2d STATIC sll_maxwell_2d_fdtd.F90)
add_executable(test_maxwell_2d_fdtd unit_test.F90)
target_link_libraries(test_maxwell_2d_fdtd sll_assertion 
                                           sll_constants 
                                           sll_memory 
                                           sll_collective 
                                           sll_utilities 
                                           sll_fft 
                                           sll_maxwell_2d)

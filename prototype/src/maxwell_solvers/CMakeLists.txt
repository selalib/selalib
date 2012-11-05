ADD_LIBRARY(sll_maxwell STATIC sll_maxwell.F90
                               sll_maxwell_2d_fdtd.F90
                               sll_maxwell_2d_pstd.F90)

ADD_DEPENDENCIES(sll_maxwell sll_assertion sll_constants sll_memory 
                              sll_collective sll_utilities sll_fft)

ADD_EXECUTABLE(test_maxwell_2d_fdtd unit_test_2d_fdtd.F90)

TARGET_LINK_LIBRARIES(test_maxwell_2d_fdtd sll_maxwell
                                           sll_assertion 
                                           sll_constants 
                                           sll_memory 
                                           sll_collective 
                                           sll_utilities 
                                           sll_fft )

IF(FFTW_FOUND)

   ADD_EXECUTABLE(test_maxwell_2d_pstd unit_test_2d_pstd.F90 )
   
   TARGET_LINK_LIBRARIES(test_maxwell_2d_pstd sll_maxwell
                                              sll_assertion 
                                              sll_constants 
                                              sll_memory 
                                              sll_collective 
                                              sll_utilities 
                                              ${FFTW_LIBRARIES})
ENDIF(FFTW_FOUND)

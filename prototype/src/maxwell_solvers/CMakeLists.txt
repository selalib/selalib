IF(FFTW_FOUND AND FFTW_ENABLED)

   ADD_LIBRARY(sll_maxwell_solvers STATIC sll_maxwell_solvers.F90
                                          sll_maxwell_2d_fdtd.F90
                                          sll_maxwell_2d_pstd.F90
                                          sll_maxwell_3d_pstd.F90)

   ADD_EXECUTABLE(test_maxwell_2d_pstd unit_test_2d_pstd.F90 )
   
   TARGET_LINK_LIBRARIES(test_maxwell_2d_pstd sll_maxwell_solvers
                                              sll_assertion 
                                              sll_constants 
                                              sll_memory 
                                              sll_utilities 
                                              ${FFTW_LIBRARIES})
   # ****************************************************************************
   #                           PARALLEL MAXWELL SOLVERS
   # ****************************************************************************

   IF(HDF5_ENABLE_PARALLEL AND HDF5_PARALLEL_ENABLED)
   
     ADD_LIBRARY(sll_maxwell_parallel STATIC 
                 sll_maxwell_2d_pstd_parallel.F90 )
   
     ADD_DEPENDENCIES(sll_maxwell_parallel sll_remap sll_collective)
     
     ADD_EXECUTABLE(test_maxwell_2d_pstd_parallel 
                    unit_test_2d_pstd_parallel.F90 )
   
     TARGET_LINK_LIBRARIES(test_maxwell_2d_pstd_parallel 
                           sll_maxwell_parallel
                           sll_remap
                           sll_file_io_parallel
                           sll_collective
                           sll_fft
                           sll_memory
                           sll_timer
                           sll_utilities
                           sll_assertion )
   ENDIF()
ELSE()

   ADD_LIBRARY(sll_maxwell_solvers STATIC sll_maxwell_solvers.F90 
                                          sll_maxwell_2d_fdtd.F90)

ENDIF(FFTW_FOUND AND FFTW_ENABLED)

ADD_EXECUTABLE(test_maxwell_2d_fdtd unit_test_2d_fdtd.F90)

ADD_DEPENDENCIES(sll_maxwell_solvers sll_assertion 
                                     sll_constants 
                                     sll_memory 
                                     sll_utilities )

TARGET_LINK_LIBRARIES(test_maxwell_2d_fdtd sll_maxwell_solvers
                                           sll_assertion 
                                           sll_constants 
                                           sll_memory 
                                           sll_utilities )

ADD_EXECUTABLE(test_allocate test.F90)
TARGET_LINK_LIBRARIES(test_allocate sll_memory)


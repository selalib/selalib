IF(FFTW_FOUND AND FFTW_ENABLED)
  # this is very strange. Why compile fftw3.F90 here? this compilation should
  # depend on the configuration flags only!
   ADD_LIBRARY(sll_maxwell_solvers STATIC ../fft/fftw3.F90
                                          sll_maxwell_solvers_base.F90
                                          sll_maxwell_2d_fdtd.F90
                                          sll_maxwell_2d_pstd.F90
                                          sll_maxwell_2d_diga.F90
                                          sll_maxwell_3d_pstd.F90)

   TARGET_LINK_LIBRARIES(sll_maxwell_solvers sll_utilities ${FFTW_LIBRARIES})

   ADD_EXECUTABLE(test_maxwell_2d_pstd unit_test_2d_pstd.F90 )
   TARGET_LINK_LIBRARIES(test_maxwell_2d_pstd sll_maxwell_solvers)

   # ****************************************************************
   #                           PARALLEL MAXWELL SOLVERS
   # ***************************************************************

   IF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)
   
     ADD_LIBRARY(sll_maxwell_parallel STATIC 
                 sll_maxwell_2d_pstd_parallel.F90 )
   
     TARGET_LINK_LIBRARIES(sll_maxwell_parallel sll_remap ${FFTW_LIBRARIES})
     
     ADD_EXECUTABLE(test_maxwell_2d_pstd_parallel 
                    unit_test_2d_pstd_parallel.F90 )
   
     TARGET_LINK_LIBRARIES(test_maxwell_2d_pstd_parallel 
                           sll_maxwell_parallel
                           sll_file_io_parallel)
   ENDIF()

ELSE()

   ADD_LIBRARY(sll_maxwell_solvers STATIC sll_maxwell_solvers_base.F90 
                                          sll_maxwell_2d_fdtd.F90
                                          sll_maxwell_2d_diga.F90
   )

ENDIF(FFTW_FOUND AND FFTW_ENABLED)

TARGET_LINK_LIBRARIES(sll_maxwell_solvers sll_dg_fields sll_utilities )

ADD_EXECUTABLE(test_maxwell_2d_fdtd unit_test_2d_fdtd.F90)
TARGET_LINK_LIBRARIES(test_maxwell_2d_fdtd sll_maxwell_solvers)

ADD_EXECUTABLE(test_maxwell_2d_diga unit_test_2d_diga.F90)
TARGET_LINK_LIBRARIES(test_maxwell_2d_diga sll_maxwell_solvers)


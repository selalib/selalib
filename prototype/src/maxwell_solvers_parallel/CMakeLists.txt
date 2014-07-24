IF(FFTW_FOUND AND FFTW_ENABLED)

   
   ADD_LIBRARY(sll_maxwell_parallel STATIC 
               sll_maxwell_2d_diga.F90
               sll_maxwell_2d_diga_parallel.F90
               sll_maxwell_2d_pstd_parallel.F90 )
   
   TARGET_LINK_LIBRARIES(sll_maxwell_parallel sll_remap ${FFTW_LIBRARIES})
     
   IF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)

     ADD_EXECUTABLE(test_maxwell_2d_pstd_parallel 
                    unit_test_2d_pstd_parallel.F90 )
   
     TARGET_LINK_LIBRARIES(test_maxwell_2d_pstd_parallel 
                           sll_maxwell_parallel
                           sll_file_io_parallel)


     ADD_EXECUTABLE(test_maxwell_2d_diga_parallel unit_test_2d_diga_parallel.F90 )
   
     TARGET_LINK_LIBRARIES(test_maxwell_2d_diga_parallel sll_maxwell_parallel)

   ENDIF()

ENDIF(FFTW_FOUND AND FFTW_ENABLED)

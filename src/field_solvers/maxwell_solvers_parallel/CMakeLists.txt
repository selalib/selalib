IF(FFTW_FOUND AND FFTW_ENABLED)

  ADD_LIBRARY(sll_maxwell_parallel STATIC 
              #sll_maxwell_2d_diga_parallel.F90 
              sll_maxwell_2d_pstd_parallel.F90 )
  
  TARGET_LINK_LIBRARIES(sll_maxwell_parallel sll_remap ${FFTW_LIBRARIES})
    
  IF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND BUILD_TESTING)
    ADD_SUBDIRECTORY(testing)
  ENDIF()

ENDIF(FFTW_FOUND AND FFTW_ENABLED)

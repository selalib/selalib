ADD_LIBRARY(sll_maxwell_parallel STATIC 
            sll_m_maxwell_2d_periodic_cartesian_par.F90 )

TARGET_LINK_LIBRARIES(sll_maxwell_parallel sll_remap sll_fft )
  
IF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF()

IF (NOT (FFT_LIB STREQUAL "SLLFFT"))

  SET( SPARSEGRID_SRC_FILES    
    sll_m_sparse_grid_interpolator.F90
    sll_m_sparse_grid_4d.F90
    sll_m_sparse_grid_3d.F90
    sll_m_sparse_grid_2d.F90
    sll_m_poisson_2d_sparse_grid_fft.F90
    sll_m_poisson_3d_sparse_grid_fft.F90)

  ADD_LIBRARY(sll_sparse_grid STATIC ${SPARSEGRID_SRC_FILES})
  TARGET_LINK_LIBRARIES(sll_sparse_grid sll_interpolators sll_periodic_interpolation sll_utilities ${LAPACK_LIBRARIES})

  IF (BUILD_TESTING)
    ADD_SUBDIRECTORY(testing)
  ENDIF()


ENDIF()

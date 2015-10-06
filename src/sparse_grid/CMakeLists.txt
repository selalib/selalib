IF (FFTW_FOUND AND FFTW_ENABLED)

  SET( SPARSEGRID_SRC_FILES    
    sll_sparsegrid_interpolator.F90
    sparse_grid_4d.F90
    sparse_grid_3d.F90
    sparse_grid_2d.F90
    sll_poisson_2d_sparsegrid_fft.F90
    sll_poisson_3d_sparsegrid_fft.F90)

  ADD_LIBRARY(sll_sparsegrid ${SPARSEGRID_SRC_FILES})
  TARGET_LINK_LIBRARIES(sll_sparsegrid sll_interpolators sll_periodic_interpolation sll_utilities ${LAPACK_LIBRARIES})

ENDIF()
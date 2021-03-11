if(NOT (FFT_LIB STREQUAL "SLLFFT"))

  set(SPARSEGRID_SRC_FILES
      sll_m_sparse_grid_interpolator.F90 sll_m_sparse_grid_4d.F90
      sll_m_sparse_grid_3d.F90 sll_m_sparse_grid_2d.F90
      sll_m_poisson_2d_sparse_grid_fft.F90 sll_m_poisson_3d_sparse_grid_fft.F90)

  add_library(sll_sparse_grid STATIC ${SPARSEGRID_SRC_FILES})
  target_link_libraries(
    sll_sparse_grid sll_interpolators sll_periodic_interpolation sll_utilities
    ${LAPACK_LIBRARIES})

  if(BUILD_TESTING)
    add_subdirectory(testing)
  endif()

endif()

if(HDF5_ENABLED AND HDF5_IS_PARALLEL)

  add_library(sll_fcisl STATIC sll_m_fcisl.F90 sll_m_derivative_2d_oblic.F90
                               sll_m_fcisl_toroidal.F90)

  target_link_libraries(sll_fcisl sll_advection sll_interpolators sll_utilities
                        sll_file_io sll_parallel_array_utilities)

  # Tests
  if(BUILD_TESTING)

    add_subdirectory(testing)

  endif(BUILD_TESTING)

endif(HDF5_ENABLED AND HDF5_IS_PARALLEL)

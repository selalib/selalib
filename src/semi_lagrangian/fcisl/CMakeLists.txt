IF(HDF5_ENABLED AND HDF5_IS_PARALLEL)

  ADD_LIBRARY( sll_fcisl STATIC sll_m_fcisl.F90
                                sll_m_derivative_2d_oblic.F90
                                sll_m_fcisl_toroidal.F90 )

  TARGET_LINK_LIBRARIES( sll_fcisl 
    sll_advection
    sll_interpolators
    sll_utilities
    sll_file_io
    sll_parallel_array_utilities )


  # Tests
  IF(BUILD_TESTING)

    ADD_SUBDIRECTORY(testing)

  ENDIF(BUILD_TESTING)

ENDIF(HDF5_ENABLED AND HDF5_IS_PARALLEL)

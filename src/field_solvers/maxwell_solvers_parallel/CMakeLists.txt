add_library(sll_maxwell_parallel STATIC
            sll_m_maxwell_2d_periodic_cartesian_par.F90)

target_link_libraries(sll_maxwell_parallel sll_remap sll_fft)

if(HDF5_IS_PARALLEL
   AND HDF5_PARALLEL_ENABLED
   AND BUILD_TESTING)
  add_subdirectory(testing)
endif()

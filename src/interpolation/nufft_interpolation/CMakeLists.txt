add_library(sll_nufft_interpolation STATIC sll_m_nufft_interpolation.F90)
target_link_libraries(sll_nufft_interpolation sll_fft nufft)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)

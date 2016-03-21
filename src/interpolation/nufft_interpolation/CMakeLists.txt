ADD_LIBRARY( sll_nufft_interpolation STATIC sll_m_nufft_interpolation.F90 )
TARGET_LINK_LIBRARIES( sll_nufft_interpolation sll_fft nufft )

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)

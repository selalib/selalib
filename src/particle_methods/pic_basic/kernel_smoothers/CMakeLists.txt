ADD_LIBRARY(sll_kernel_smoothers STATIC 
  sll_m_kernel_smoother_spline_1d.F90
  sll_m_kernel_smoother_spline_2d.F90)

TARGET_LINK_LIBRARIES( sll_kernel_smoothers
  sll_memory
  sll_assert
  sll_pic_interface
  sll_particle_types
  sll_splines
  )


# Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
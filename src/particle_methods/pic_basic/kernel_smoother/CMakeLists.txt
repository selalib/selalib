ADD_LIBRARY(sll_kernel_smoothers STATIC kernel_smoother_base.F90 
  sll_kernel_smoother_spline_2d.F90
  sll_kernel_smoother_spline_1d.F90)

TARGET_LINK_LIBRARIES( sll_kernel_smoothers
  sll_memory
  sll_assert
  sll_pic_interface
  sll_particle_types
  sll_splines
  )


# Tests
IF(BUILD_TESTING)
  ADD_EXECUTABLE( test_kernel_smoother_spline_1d unit_test_kernel_smoother_spline_1d.F90)
  TARGET_LINK_LIBRARIES( test_kernel_smoother_spline_1d sll_kernel_smoothers )
  ADD_TEST( NAME kernel_smoother_spline_1d COMMAND test_kernel_smoother_spline_1d)

ENDIF(BUILD_TESTING)
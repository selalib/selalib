ADD_LIBRARY(sll_kernel_smoothers STATIC kernel_smoother_base.F90 
  sll_kernel_smoother_linear_2d.F90)

TARGET_LINK_LIBRARIES( sll_kernel_smoothers
  sll_memory
  sll_assert
  sll_pic_interface
  sll_splines
  )
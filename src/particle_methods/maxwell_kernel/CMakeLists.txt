ADD_LIBRARY(sll_pic_maxwell STATIC sll_pic_maxwell_base.F90 
  sll_pic_maxwell_1d_fourier_modes.F90)

TARGET_LINK_LIBRARIES( sll_pic_maxwell
  sll_memory
  sll_assert
  sll_pic_interface
  sll_particle_types
  )
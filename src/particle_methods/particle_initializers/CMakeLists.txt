ADD_LIBRARY(sll_particle_initializers STATIC sll_particle_initializer_2x2v.F90 )

TARGET_LINK_LIBRARIES( sll_particle_initializers
  sll_memory
  sll_assert
  sll_pic_interface
  )
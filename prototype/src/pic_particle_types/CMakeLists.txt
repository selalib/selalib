ADD_LIBRARY( sll_pic_particles STATIC 
  particle_representation.F90
  particle_group_2d.F90 )

TARGET_LINK_LIBRARIES(sll_pic_particles
  sll_memory )


ADD_EXECUTABLE(test_pic_particles unit_test.F90  )

TARGET_LINK_LIBRARIES(test_pic_particles sll_pic_particles )


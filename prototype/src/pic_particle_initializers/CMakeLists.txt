ADD_LIBRARY( sll_pic_initializers STATIC 
  sll_particle_init4D.F90 )

TARGET_LINK_LIBRARIES( sll_pic_initializers
  sll_pic_particles
  sll_visu_pic
  sll_logical_meshes
  sll_working_precision
  sll_assertion
  sll_memory )

ADD_EXECUTABLE(test_pic_initializers unit_test.F90  )

TARGET_LINK_LIBRARIES(test_pic_initializers
  sll_pic_initializers
#  sll_visu_pic
 sll_collective
 # sll_file_io )
)
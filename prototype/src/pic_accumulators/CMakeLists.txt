ADD_LIBRARY( sll_accumulators STATIC 
  sll_accumulators.F90 )

TARGET_LINK_LIBRARIES(sll_accumulators
  sll_memory
  sll_assertion
  sll_working_precision )


ADD_EXECUTABLE(test_pic_accumulator unit_test.F90  )

TARGET_LINK_LIBRARIES(test_pic_accumulator
  sll_accumulators
  sll_pic_particles
  sll_pic_initializers
  sll_logical_meshes )

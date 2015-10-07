ADD_LIBRARY( sll_pic_particles STATIC 
  particle_representation.F90
  particle_group_2d.F90 
  particle_group_4d.F90
  )

TARGET_LINK_LIBRARIES( sll_pic_particles
  sll_memory
  sll_working_precision
  sll_assert
  sll_meshes )

# Ctest
IF(BUILD_TESTING)
  ADD_EXECUTABLE( test_pic_particles unit_test.F90 )
  TARGET_LINK_LIBRARIES( test_pic_particles sll_pic_particles )
  ADD_TEST( NAME pic_particles COMMAND test_pic_particles )
  SET_TESTS_PROPERTIES( pic_particles PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(BUILD_TESTING)

ADD_LIBRARY( sll_simple_pic_particles STATIC
  simple_pic_4d_particle.F90
  simple_pic_4d_group.F90
  )

TARGET_LINK_LIBRARIES( sll_simple_pic_particles
  sll_memory
  sll_working_precision
  sll_assert
  sll_meshes
  sll_pic_interface )

# Ctest
IF(BUILD_TESTING)
  ADD_EXECUTABLE( test_simple_pic_particles unit_test.F90 )
  TARGET_LINK_LIBRARIES( test_simple_pic_particles sll_simple_pic_particles )
  ADD_TEST( NAME simple_pic_particles COMMAND test_simple_pic_particles )
  SET_TESTS_PROPERTIES( simple_pic_particles PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(BUILD_TESTING)

ADD_LIBRARY( sll_simple_pic STATIC
  simple_pic_4d_particle.F90
  simple_pic_4d_group.F90
  )

TARGET_LINK_LIBRARIES( sll_simple_pic
  sll_memory
  sll_working_precision
  sll_assert
  sll_meshes
  sll_remapped_pic
  # sll_pic_interface
  sll_pic_accumulators
  sll_remapped_pic_utilities
  )

# Ctest
IF(BUILD_TESTING)
  ADD_EXECUTABLE( test_simple_pic unit_test.F90 )
  TARGET_LINK_LIBRARIES( test_simple_pic sll_simple_pic )
  ADD_TEST( NAME simple_pic COMMAND test_simple_pic )
  SET_TESTS_PROPERTIES( simple_pic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(BUILD_TESTING)

# CTest
ADD_EXECUTABLE( test_pic_accumulators test_pic_accumulators.F90 )
TARGET_LINK_LIBRARIES( test_pic_accumulators
  sll_pic_utilities
  sll_pic_initializers
  sll_pic_accumulators
  sll_pic_particles
  sll_meshes
  sll_assert)
ADD_TEST( NAME pic_accumulators COMMAND test_pic_accumulators )
SET_TESTS_PROPERTIES( pic_accumulators PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

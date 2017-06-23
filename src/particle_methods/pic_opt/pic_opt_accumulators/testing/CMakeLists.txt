# CTest
ADD_EXECUTABLE( test_pic_accumulators test_pic_accumulators.F90 )
TARGET_LINK_LIBRARIES( test_pic_accumulators
  sll_pic_opt_utilities
  sll_pic_opt_initializers
  sll_pic_opt_accumulators
  sll_pic_opt_particles
  sll_meshes
  sll_assert)
ADD_TEST( NAME pic_accumulators COMMAND test_pic_accumulators )
SET_TESTS_PROPERTIES( pic_accumulators PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

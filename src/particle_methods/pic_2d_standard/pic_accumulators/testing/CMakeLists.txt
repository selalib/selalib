# CTest
ADD_EXECUTABLE( test_pic_accumulator unit_test.F90 )
TARGET_LINK_LIBRARIES( test_pic_accumulator
  sll_pic_utilities
  sll_pic_initializers
  sll_pic_accumulators
  sll_pic_particles
  sll_meshes
  sll_assert
  sll_constants )
ADD_TEST( NAME pic_accumulator COMMAND test_pic_accumulator )
SET_TESTS_PROPERTIES( pic_accumulator PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

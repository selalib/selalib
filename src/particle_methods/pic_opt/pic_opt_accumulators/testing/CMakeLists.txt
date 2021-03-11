# CTest
add_executable(test_pic_accumulators test_pic_accumulators.F90)
target_link_libraries(
  test_pic_accumulators
  sll_pic_opt_utilities
  sll_pic_opt_initializers
  sll_pic_opt_accumulators
  sll_pic_opt_particles
  sll_meshes
  sll_assert)
add_test(NAME pic_accumulators COMMAND test_pic_accumulators)
set_tests_properties(pic_accumulators PROPERTIES PASS_REGULAR_EXPRESSION
                                                 "PASSED")

add_executable(test_pic_initializers_2d test_pic_initializers_2d.F90)
target_link_libraries(test_pic_initializers_2d sll_pic_opt_initializers
                      sll_collective)
add_test(NAME pic_initializers_2d COMMAND test_pic_initializers_2d)
set_tests_properties(pic_initializers_2d PROPERTIES PASS_REGULAR_EXPRESSION
                                                    "PASSED")

add_executable(test_pic_initializers_4d test_pic_initializers_4d.F90)
target_link_libraries(test_pic_initializers_4d sll_pic_opt_initializers
                      sll_collective)
add_test(NAME pic_initializers_4d COMMAND test_pic_initializers_4d)
set_tests_properties(pic_initializers_4d PROPERTIES PASS_REGULAR_EXPRESSION
                                                    "PASSED")

add_executable(test_pic_initializers_6d test_pic_initializers_6d.F90)
target_link_libraries(test_pic_initializers_6d sll_pic_opt_initializers
                      sll_collective)
add_test(NAME pic_initializers_6d COMMAND test_pic_initializers_6d)
set_tests_properties(pic_initializers_6d PROPERTIES PASS_REGULAR_EXPRESSION
                                                    "PASSED")

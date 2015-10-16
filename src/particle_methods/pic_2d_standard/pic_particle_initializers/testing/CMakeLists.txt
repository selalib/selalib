# CTest
add_executable( test_pic_initializers unit_test.F90 )
target_link_libraries( test_pic_initializers sll_pic_initializers sll_collective )
add_test( NAME pic_initializers COMMAND test_pic_initializers )
set_tests_properties( pic_initializers PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

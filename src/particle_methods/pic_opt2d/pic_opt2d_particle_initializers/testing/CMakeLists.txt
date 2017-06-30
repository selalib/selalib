#These tests are disabled for PGI compiler because of macros used inside. (PN)
if(NOT PGI)

  add_executable( test_pic_initializers_2d test_pic_initializers_2d.F90 )
  target_link_libraries( test_pic_initializers_2d sll_pic_opt2d_initializers sll_collective )
  add_test( NAME pic_initializers_2d COMMAND test_pic_initializers_2d )
  set_tests_properties( pic_initializers_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  add_executable( test_pic_initializers_4d test_pic_initializers_4d.F90 )
  target_link_libraries( test_pic_initializers_4d sll_pic_opt2d_initializers sll_collective )
  add_test( NAME pic_initializers_4d COMMAND test_pic_initializers_4d )
  set_tests_properties( pic_initializers_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

endif(NOT PGI)

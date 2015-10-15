ADD_EXECUTABLE( test_simple_pic test_simple_pic.F90 )
TARGET_LINK_LIBRARIES( test_simple_pic sll_simple_pic )

ADD_TEST( NAME simple_pic COMMAND test_simple_pic )
SET_TESTS_PROPERTIES( simple_pic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_bsl_lt_pic_4d test_bsl_lt_pic_4d.F90 )
TARGET_LINK_LIBRARIES( test_bsl_lt_pic_4d
  sll_timer
  sll_assert
  sll_remapped_pic
  sll_bsl_lt_pic
  )

ADD_TEST( NAME bsl_lt_pic_4d COMMAND test_bsl_lt_pic_4d )
SET_TESTS_PROPERTIES( bsl_lt_pic_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

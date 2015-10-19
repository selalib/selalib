  ADD_EXECUTABLE( test_reduction test_reduction.F90 )
  TARGET_LINK_LIBRARIES( test_reduction sll_reduction )
  ADD_TEST( NAME reduction COMMAND test_reduction )
  SET_TESTS_PROPERTIES( reduction PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

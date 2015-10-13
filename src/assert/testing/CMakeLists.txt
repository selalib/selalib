# Ctest
IF(CMAKE_BUILD_TYPE STREQUAL DEBUG)
   ADD_EXECUTABLE( test_assert unit_test.F90 )
   TARGET_LINK_LIBRARIES( test_assert sll_assert )
   ADD_TEST( NAME assert COMMAND test_assert )
   SET( passRegex "Assertion error triggered in file" )
   SET_TESTS_PROPERTIES( assert PROPERTIES PASS_REGULAR_EXPRESSION "${passRegex}" )
ENDIF(CMAKE_BUILD_TYPE STREQUAL DEBUG)

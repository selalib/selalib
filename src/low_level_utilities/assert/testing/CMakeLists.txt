# [YG] As assert causes abnormal termination with abort(), there is no way to
# tell CTest that the unit test should in fact fail catastrophically.
#
# # Ctest
# IF(CMAKE_BUILD_TYPE STREQUAL Debug )
#    ADD_EXECUTABLE( test_assert test_assert.F90 )
#    TARGET_LINK_LIBRARIES( test_assert sll_assert )
#    ADD_TEST( NAME assert COMMAND test_assert )
#    SET( passRegex "Assertion error triggered in file" )
#    SET_TESTS_PROPERTIES( assert PROPERTIES PASS_REGULAR_EXPRESSION "${passRegex}" )
# ENDIF(CMAKE_BUILD_TYPE STREQUAL Debug )

# Unit test to be run manually: it must fail
if( CMAKE_BUILD_TYPE STREQUAL Debug )
  add_executable( test_assert test_assert.F90 )
  target_link_libraries( test_assert sll_assert )
endif()

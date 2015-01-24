ADD_LIBRARY( sll_constants STATIC sll_constants.F90 )
TARGET_LINK_LIBRARIES( sll_constants sll_working_precision )

# Ctest
IF(BUILD_TESTING)
  ADD_EXECUTABLE( test_constants unit_test.F90 )
  TARGET_LINK_LIBRARIES(test_constants sll_constants)
  ADD_TEST( NAME constants COMMAND test_constants )
ENDIF(BUILD_TESTING)

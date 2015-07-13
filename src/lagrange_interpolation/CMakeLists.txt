ADD_LIBRARY( sll_lagrange_interpolation STATIC
             sll_lagrange_interpolation.F90)

TARGET_LINK_LIBRARIES( sll_lagrange_interpolation sll_utilities
                                                  sll_boundary_condition_descriptors)

# Ctest
IF(BUILD_TESTING)

  ADD_EXECUTABLE( test_lagrange_interpolation unit_test.F90 )
  TARGET_LINK_LIBRARIES( test_lagrange_interpolation 
    sll_lagrange_interpolation sll_constants )
  
  ADD_TEST( NAME lagrange_interpolation COMMAND test_lagrange_interpolation )
  SET_TESTS_PROPERTIES( lagrange_interpolation PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ENDIF(BUILD_TESTING)

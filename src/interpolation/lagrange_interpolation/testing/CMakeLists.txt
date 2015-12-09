  ADD_EXECUTABLE( test_lagrange_interpolation test_lagrange_interpolation.F90 )
  TARGET_LINK_LIBRARIES( test_lagrange_interpolation 
    sll_lagrange_interpolation  )
  
  ADD_TEST( NAME lagrange_interpolation COMMAND test_lagrange_interpolation )
  SET_TESTS_PROPERTIES( lagrange_interpolation PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

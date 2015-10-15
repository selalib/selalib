  ADD_EXECUTABLE( test_periodic_interpolation test_periodic_interpolation.F90 )
  TARGET_LINK_LIBRARIES( test_periodic_interpolation sll_periodic_interpolation )
  ADD_TEST( NAME periodic_interpolation COMMAND test_periodic_interpolation )
  SET_TESTS_PROPERTIES( periodic_interpolation PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

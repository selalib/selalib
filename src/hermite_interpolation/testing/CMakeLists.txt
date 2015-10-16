  ADD_EXECUTABLE( test_hermite_interpolation_1d test_hermite_interpolation_1d.F90 )
  TARGET_LINK_LIBRARIES( test_hermite_interpolation_1d
    sll_hermite_interpolation
    sll_constants )

  ADD_EXECUTABLE( test_hermite_interpolation_2d test_hermite_interpolation_2d.F90 )
  TARGET_LINK_LIBRARIES( test_hermite_interpolation_2d
    sll_hermite_interpolation
    sll_constants )

  ADD_EXECUTABLE( test_hermite_aligned_interpolation_2d 
    test_hermite_aligned_interpolation_2d.F90 )
  TARGET_LINK_LIBRARIES( test_hermite_aligned_interpolation_2d
    sll_hermite_aligned_interpolation
    sll_constants )

  # Ctest
  ADD_TEST( NAME hermite_interpolation_1d COMMAND test_hermite_interpolation_1d )
  SET_TESTS_PROPERTIES( hermite_interpolation_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  ADD_TEST( NAME hermite_interpolation_2d COMMAND test_hermite_interpolation_2d )
  SET_TESTS_PROPERTIES( hermite_interpolation_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  ADD_TEST( NAME hermite_aligned_interpolation_2d COMMAND test_hermite_aligned_interpolation_2d )
  SET_TESTS_PROPERTIES( hermite_aligned_interpolation_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

IF(BUILD_TESTING)
  ADD_EXECUTABLE( test_integration test_integration.F90
                                   test_func_module.F90 )
  TARGET_LINK_LIBRARIES( test_integration sll_quadrature )
  # Ctest
  ADD_TEST( NAME integration COMMAND test_integration )
  SET_TESTS_PROPERTIES( integration PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF(BUILD_TESTING)

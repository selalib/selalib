IF(BUILD_TESTING)

  ADD_LIBRARY( quadrature_helper STATIC test_func_module.F90 )

  ADD_EXECUTABLE( test_fekete test_fekete.F90 )
  TARGET_LINK_LIBRARIES( test_fekete sll_quadrature )
  ADD_TEST( NAME quadrature_fekete COMMAND test_fekete )
  SET_TESTS_PROPERTIES( quadrature_fekete PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

  ADD_EXECUTABLE( test_gauss_legendre test_gauss_legendre.F90 )
  TARGET_LINK_LIBRARIES( test_gauss_legendre quadrature_helper sll_quadrature )
  ADD_TEST( NAME quadrature_gauss_legendre COMMAND test_gauss_legendre )
  SET_TESTS_PROPERTIES( quadrature_gauss_legendre PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

  ADD_EXECUTABLE( test_gauss_lobatto test_gauss_lobatto.F90 )
  TARGET_LINK_LIBRARIES( test_gauss_lobatto quadrature_helper sll_quadrature )
  ADD_TEST( NAME quadrature_gauss_lobatto COMMAND test_gauss_lobatto )
  SET_TESTS_PROPERTIES( quadrature_gauss_lobatto PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

  ADD_EXECUTABLE( test_rectangle test_rectangle.F90 )
  TARGET_LINK_LIBRARIES( test_rectangle quadrature_helper sll_quadrature )
  ADD_TEST( NAME quadrature_rectangle COMMAND test_rectangle )
  SET_TESTS_PROPERTIES( quadrature_rectangle PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

  ADD_EXECUTABLE( test_trapz test_trapz.F90 )
  TARGET_LINK_LIBRARIES( test_trapz quadrature_helper sll_quadrature )
  ADD_TEST( NAME quadrature_trapz COMMAND test_trapz )
  SET_TESTS_PROPERTIES( quadrature_trapz PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ENDIF(BUILD_TESTING)

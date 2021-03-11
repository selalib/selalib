if(BUILD_TESTING)

  add_library(quadrature_helper STATIC test_func_module.F90)
  target_link_libraries(quadrature_helper sll_working_precision)

  add_executable(test_fekete test_fekete.F90)
  target_link_libraries(test_fekete sll_quadrature)
  add_test(NAME quadrature_fekete COMMAND test_fekete)
  set_tests_properties(quadrature_fekete PROPERTIES PASS_REGULAR_EXPRESSION
                                                    "PASSED")

  add_executable(test_gauss_legendre test_gauss_legendre.F90)
  target_link_libraries(test_gauss_legendre quadrature_helper sll_quadrature)
  add_test(NAME quadrature_gauss_legendre COMMAND test_gauss_legendre)
  set_tests_properties(quadrature_gauss_legendre
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

  add_executable(test_gauss_lobatto test_gauss_lobatto.F90)
  target_link_libraries(test_gauss_lobatto quadrature_helper sll_quadrature)
  add_test(NAME quadrature_gauss_lobatto COMMAND test_gauss_lobatto)
  set_tests_properties(quadrature_gauss_lobatto
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

  add_executable(test_rectangle test_rectangle.F90)
  target_link_libraries(test_rectangle quadrature_helper sll_quadrature)
  add_test(NAME quadrature_rectangle COMMAND test_rectangle)
  set_tests_properties(quadrature_rectangle PROPERTIES PASS_REGULAR_EXPRESSION
                                                       "PASSED")

  add_executable(test_trapz test_trapz.F90)
  target_link_libraries(test_trapz quadrature_helper sll_quadrature)
  add_test(NAME quadrature_trapz COMMAND test_trapz)
  set_tests_properties(quadrature_trapz PROPERTIES PASS_REGULAR_EXPRESSION
                                                   "PASSED")

  add_executable(test_primitives test_primitives.F90)
  target_link_libraries(test_primitives sll_quadrature)
  add_test(NAME quadrature_primitives COMMAND test_primitives)
  set_tests_properties(quadrature_primitives PROPERTIES PASS_REGULAR_EXPRESSION
                                                        "PASSED")

endif(BUILD_TESTING)

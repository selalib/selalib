ADD_LIBRARY( sll_integration STATIC gauss_legendre.F90
                                    gauss_lobatto.F90
                                    rectangle.F90
                                    trapz.F90
                                    gauss_triangle.F90
                                    fekete.F90)

TARGET_LINK_LIBRARIES( sll_integration sll_splines
                                       sll_meshes
                                       sll_memory
                                       sll_utilities)

# Ctest
IF(BUILD_TESTING)
  ADD_EXECUTABLE( test_integration unit_test.F90 test_func_module.F90 )
  TARGET_LINK_LIBRARIES( test_integration sll_integration )
  ADD_TEST( NAME integration COMMAND test_integration )
  SET_TESTS_PROPERTIES( integration PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF(BUILD_TESTING)

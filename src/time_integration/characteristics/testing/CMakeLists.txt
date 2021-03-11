add_library(characteristics_1d_test_helper STATIC
            m_characteristics_1d_test_helper.F90)
target_link_libraries(characteristics_1d_test_helper sll_characteristics)
macro(CHARACTERISTICS_TEST TEST_NAME)
  set(EXECUTABLE test_${TEST_NAME})
  add_executable(${EXECUTABLE} test_${TEST_NAME}.F90)
  target_link_libraries(${EXECUTABLE} characteristics_1d_test_helper
                        sll_characteristics sll_constants sll_working_precision)
  add_test(NAME ${TEST_NAME} COMMAND ${EXECUTABLE})
  set_tests_properties(${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
endmacro(CHARACTERISTICS_TEST)

characteristics_test(characteristics_1d_explicit_euler)
characteristics_test(characteristics_1d_explicit_euler_conservative)
characteristics_test(characteristics_1d_trapezoid)
characteristics_test(characteristics_1d_trapezoid_conservative)
characteristics_test(characteristics_2d_explicit_euler)
characteristics_test(characteristics_2d_verlet)
characteristics_test(characteristics_2d_explicit_euler_conservative)

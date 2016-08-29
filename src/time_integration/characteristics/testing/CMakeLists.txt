add_library( characteristics_1d_test_helper OBJECT
  m_characteristics_1d_test_helper.F90 )

  MACRO( CHARACTERISTICS_TEST TEST_NAME )
    SET( EXECUTABLE test_${TEST_NAME} )
    ADD_EXECUTABLE( ${EXECUTABLE} test_${TEST_NAME}.F90 
      $<TARGET_OBJECTS:characteristics_1d_test_helper> )
    TARGET_LINK_LIBRARIES( ${EXECUTABLE}
      sll_characteristics
      sll_constants
      sll_working_precision )
    ADD_TEST( NAME ${TEST_NAME} COMMAND ${EXECUTABLE} )
    SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES
      PASS_REGULAR_EXPRESSION "PASSED" )
  ENDMACRO( CHARACTERISTICS_TEST )

  CHARACTERISTICS_TEST( characteristics_1d_explicit_euler              )
  CHARACTERISTICS_TEST( characteristics_1d_explicit_euler_conservative )
  CHARACTERISTICS_TEST( characteristics_1d_trapezoid                   ) 
  CHARACTERISTICS_TEST( characteristics_1d_trapezoid_conservative      )
  CHARACTERISTICS_TEST( characteristics_2d_explicit_euler              )
  CHARACTERISTICS_TEST( characteristics_2d_verlet                      )
  CHARACTERISTICS_TEST( characteristics_2d_explicit_euler_conservative )


#  ADD_EXECUTABLE(
#    test_characteristics_1d_explicit_euler
#    test_characteristics_1d_explicit_euler.F90)
#  TARGET_LINK_LIBRARIES(
#    test_characteristics_1d_explicit_euler
#    sll_characteristics)
#
#  ADD_EXECUTABLE(
#    test_characteristics_1d_trapezoid
#    test_characteristics_1d_trapezoid.F90)
#  TARGET_LINK_LIBRARIES(
#    test_characteristics_1d_trapezoid
#    sll_characteristics)
#
#  ADD_EXECUTABLE(
#    test_characteristics_1d_explicit_euler_conservative
#    test_characteristics_1d_explicit_euler_conservative.F90)
#  TARGET_LINK_LIBRARIES(
#    test_characteristics_1d_explicit_euler_conservative
#    sll_characteristics)
#
#  ADD_EXECUTABLE(
#    test_characteristics_1d_trapezoid_conservative
#    test_characteristics_1d_trapezoid_conservative.F90)
#  TARGET_LINK_LIBRARIES(
#    test_characteristics_1d_trapezoid_conservative
#    sll_characteristics)
#
#  ADD_EXECUTABLE(
#    test_characteristics_2d_explicit_euler
#    test_characteristics_2d_explicit_euler.F90)
#  TARGET_LINK_LIBRARIES(
#    test_characteristics_2d_explicit_euler
#    sll_characteristics)
#
#  ADD_EXECUTABLE(
#    test_characteristics_2d_verlet
#    test_characteristics_2d_verlet.F90)
#  TARGET_LINK_LIBRARIES(
#    test_characteristics_2d_verlet
#    sll_characteristics)
#
#  ADD_EXECUTABLE(
#    test_characteristics_2d_explicit_euler_conservative
#    test_characteristics_2d_explicit_euler_conservative.F90)
#  TARGET_LINK_LIBRARIES(
#    test_characteristics_2d_explicit_euler_conservative
#    sll_characteristics)
#
## Ctests
#ADD_TEST( NAME characteristics_1d_explicit_euler              COMMAND test_characteristics_1d_explicit_euler              )
#ADD_TEST( NAME characteristics_1d_explicit_euler_conservative COMMAND test_characteristics_1d_explicit_euler_conservative )
#ADD_TEST( NAME characteristics_1d_trapezoid                   COMMAND test_characteristics_1d_trapezoid                   )
#ADD_TEST( NAME characteristics_1d_trapezoid_conservative      COMMAND test_characteristics_1d_trapezoid_conservative      )
#ADD_TEST( NAME characteristics_2d_explicit_euler              COMMAND test_characteristics_2d_explicit_euler              )
#ADD_TEST( NAME characteristics_2d_explicit_euler_conservative COMMAND test_characteristics_2d_explicit_euler_conservative )
#ADD_TEST( NAME characteristics_2d_verlet                      COMMAND test_characteristics_2d_verlet                      )
#
#SET_TESTS_PROPERTIES( characteristics_1d_explicit_euler              PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
#SET_TESTS_PROPERTIES( characteristics_1d_explicit_euler_conservative PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
#SET_TESTS_PROPERTIES( characteristics_1d_trapezoid                   PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
#SET_TESTS_PROPERTIES( characteristics_1d_trapezoid_conservative      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
#SET_TESTS_PROPERTIES( characteristics_2d_explicit_euler              PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
#SET_TESTS_PROPERTIES( characteristics_2d_explicit_euler_conservative PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
#SET_TESTS_PROPERTIES( characteristics_2d_verlet                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

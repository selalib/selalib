  SET( CHARACTERISTICS_SRC_FILES 
	sll_characteristics_1d_explicit_euler.F90
	sll_characteristics_1d_explicit_euler_conservative.F90
	sll_characteristics_1d_trapezoid.F90
	sll_characteristics_1d_trapezoid_conservative.F90
	sll_characteristics_1d_base.F90
	sll_characteristics_2d_verlet.F90
	sll_characteristics_2d_explicit_euler.F90
	sll_characteristics_2d_explicit_euler_conservative.F90
	sll_characteristics_2d_base.F90)

  ADD_LIBRARY(
    sll_characteristics
    STATIC
    ${CHARACTERISTICS_SRC_FILES})
  TARGET_LINK_LIBRARIES( 
    sll_characteristics
    sll_utilities
    sll_interpolators)

  ADD_EXECUTABLE(
    test_characteristics_1d_explicit_euler
    unit_test_characteristics_1d_explicit_euler.F90)
  TARGET_LINK_LIBRARIES(
    test_characteristics_1d_explicit_euler
    sll_characteristics)

  ADD_EXECUTABLE(
    test_characteristics_1d_trapezoid
    unit_test_characteristics_1d_trapezoid.F90)
  TARGET_LINK_LIBRARIES(
    test_characteristics_1d_trapezoid
    sll_characteristics)

  ADD_EXECUTABLE(
    test_characteristics_1d_explicit_euler_conservative
    unit_test_characteristics_1d_explicit_euler_conservative.F90)
  TARGET_LINK_LIBRARIES(
    test_characteristics_1d_explicit_euler_conservative
    sll_characteristics)

  ADD_EXECUTABLE(
    test_characteristics_1d_trapezoid_conservative
    unit_test_characteristics_1d_trapezoid_conservative.F90)
  TARGET_LINK_LIBRARIES(
    test_characteristics_1d_trapezoid_conservative
    sll_characteristics)

  ADD_EXECUTABLE(
    test_characteristics_2d_explicit_euler
    unit_test_characteristics_2d_explicit_euler.F90)
  TARGET_LINK_LIBRARIES(
    test_characteristics_2d_explicit_euler
    sll_characteristics)

  ADD_EXECUTABLE(
    test_characteristics_2d_verlet
    unit_test_characteristics_2d_verlet.F90)
  TARGET_LINK_LIBRARIES(
    test_characteristics_2d_verlet
    sll_characteristics)

  ADD_EXECUTABLE(
    test_characteristics_2d_explicit_euler_conservative
    unit_test_characteristics_2d_explicit_euler_conservative.F90)
  TARGET_LINK_LIBRARIES(
    test_characteristics_2d_explicit_euler_conservative
    sll_characteristics)

# Ctests
ADD_TEST( NAME characteristics_1d_explicit_euler              COMMAND test_characteristics_1d_explicit_euler              )
ADD_TEST( NAME characteristics_1d_explicit_euler_conservative COMMAND test_characteristics_1d_explicit_euler_conservative )
ADD_TEST( NAME characteristics_1d_trapezoid                   COMMAND test_characteristics_1d_trapezoid                   )
ADD_TEST( NAME characteristics_1d_trapezoid_conservative      COMMAND test_characteristics_1d_trapezoid_conservative      )
ADD_TEST( NAME characteristics_2d_explicit_euler              COMMAND test_characteristics_2d_explicit_euler              )
ADD_TEST( NAME characteristics_2d_explicit_euler_conservative COMMAND test_characteristics_2d_explicit_euler_conservative )
ADD_TEST( NAME characteristics_2d_verlet                      COMMAND test_characteristics_2d_verlet                      )

SET_TESTS_PROPERTIES( characteristics_1d_explicit_euler              PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
SET_TESTS_PROPERTIES( characteristics_1d_explicit_euler_conservative PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
SET_TESTS_PROPERTIES( characteristics_1d_trapezoid                   PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
SET_TESTS_PROPERTIES( characteristics_1d_trapezoid_conservative      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
SET_TESTS_PROPERTIES( characteristics_2d_explicit_euler              PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
SET_TESTS_PROPERTIES( characteristics_2d_explicit_euler_conservative PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
SET_TESTS_PROPERTIES( characteristics_2d_verlet                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )



ADD_LIBRARY(characteristics_1d_test_helper STATIC  m_characteristics_1d_test_helper.F90 )
MACRO( CHARACTERISTICS_TEST TEST_NAME )
  SET( EXECUTABLE test_${TEST_NAME} )
  ADD_EXECUTABLE( ${EXECUTABLE} test_${TEST_NAME}.F90 )
  TARGET_LINK_LIBRARIES( ${EXECUTABLE}
    characteristics_1d_test_helper
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

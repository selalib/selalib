ADD_EXECUTABLE( test_ode_integrators test_ode_integrators.F90 m_ode_collection.F90 )
TARGET_LINK_LIBRARIES( test_ode_integrators sll_ode_integrators sll_assert )

IF( PYTHON3_FOUND )
  SET( PYTHON python3 )
ELSE()
  SET( PYTHON python )
ENDIF()

FILE( COPY test_ode_integrators.py
  DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )
ADD_TEST( NAME ode_integrators 
  COMMAND ${PYTHON} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_ode_integrators.py )
SET_TESTS_PROPERTIES( ode_integrators
  PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )


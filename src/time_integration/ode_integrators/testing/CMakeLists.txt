add_executable(test_ode_integrators test_ode_integrators.F90
                                    m_ode_collection.F90)
target_link_libraries(test_ode_integrators sll_ode_integrators sll_assert)

file(COPY test_ode_integrators.py DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
add_test(NAME ode_integrators
	COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_ode_integrators.py)
set_tests_properties(ode_integrators PROPERTIES PASS_REGULAR_EXPRESSION
                                                "PASSED")

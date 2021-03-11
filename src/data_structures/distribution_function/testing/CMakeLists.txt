add_executable(test_distribution_function test_distribution_function.F90)
target_link_libraries(test_distribution_function sll_distribution_function)
add_test(NAME distribution_function COMMAND test_distribution_function)

add_executable(test_initial_distribution test_initial_distribution.F90)
target_link_libraries(test_initial_distribution sll_distribution_function)

set(NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/initial_distribution.nml)
add_test(NAME initial_distribution COMMAND test_initial_distribution
                                           ${NML_FILE})
set_tests_properties(initial_distribution PROPERTIES PASS_REGULAR_EXPRESSION
                                                     "PASSED")
